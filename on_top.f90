! ontop_module.f90 - Calculates on-top pair density
!
! This file contains routines for calculating the on-top pair density
! from a 2-electron reduced density matrix (2-RDM) file.
module ontop_module
   use wfxinfo     ! Contains cartes (atomic positions), TMN, Alpha arrays
   use geninfo     ! Contains general info like pi
   use inputdat    ! Contains n_atoms_scan, n_directions, scan_directions, atom_indices, dm2_file
   use scan_utils  ! Contains generic scanning utility
   use cube_module ! For cube file output
   implicit none

   private
   public :: ontop, ontop_calculation

contains

   subroutine ontop(dm2_file_name)
      implicit none

      ! Input parameters
      character(len=*), intent(in) :: dm2_file_name

      ! Local variables for normalization factor computation
      integer :: i, max_primitive_index, orbital
      double precision :: norm_factor_1
      double precision, external :: dfact_dp
      double precision, allocatable :: omg(:)

      ! Variables for DM2 caching
      integer :: dm2_count, idx, j, k, l
      double precision :: DMval
      integer, allocatable, dimension(:,:) :: dm2_indices
      double precision, allocatable, dimension(:) :: dm2_values

      ! Determine max primitive index
      max_primitive_index = size(Alpha)
      allocate(omg(noccmo))

      ! Allocate and compute normalized primitives
      if (allocated(normalized_primitives)) deallocate(normalized_primitives)
      allocate(normalized_primitives(max_primitive_index))

      ! Precompute normalization factors for primitives
      do i = 1, max_primitive_index
         ! Compute double factorial product
         norm_factor_1 = dfact_dp(2*TMN(i,1)-1) * &
                        dfact_dp(2*TMN(i,2)-1) * &
                        dfact_dp(2*TMN(i,3)-1)

         ! Complete normalization factor
         if (ontop_hf) then
            normalized_primitives(i) = 1.0d0
         else
            normalized_primitives(i) = &
               (2.d0*Alpha(i)/pi)**(0.75d0) * &
               sqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))) / &
               sqrt(norm_factor_1)
         end if
      end do

      ! Determine orbital configuration and compute occupation factors
      if (uhf .or. opsh) then
         ! Open shell case
         ! Determine number of alpha and beta orbitals
         nalfaorb = 0
         nbetaorb = 0

         ! Initialize omg array (this should ideally come from module data)
         do i = 1, noccmo
            omg(i) = 1  ! Default to alpha (this is a placeholder - real data should come from wfxinfo)
         end do

         do orbital = 1, noccmo
            if (omg(orbital) == 1) then
               nalfaorb = nalfaorb + 1
            else
               nbetaorb = nbetaorb + 1
            end if
         end do

         ! Allocate and compute occupation factors
         if (allocated(occupation_factors)) deallocate(occupation_factors)
         allocate(occupation_factors(max(nalfaorb, nbetaorb)))
         if (allocated(orbital_type)) deallocate(orbital_type)
         allocate(orbital_type(max(nalfaorb, nbetaorb)))

         ! Reset counters
         nalfaorb = 0
         nbetaorb = 0

         ! Compute occupation factors for open shell
         do orbital = 1, noccmo
            if (omg(orbital) == 1) then
               ! Alpha orbital
               nalfaorb = nalfaorb + 1
               occupation_factors(nalfaorb) = &
                   (0.5d0*Occ(orbital)*(1.0d0 - 0.5d0*Occ(orbital)))**(0.5d0) - &
                   2.0d0*0.5d0*Occ(orbital)*(1.0d0 - 0.5d0*Occ(orbital))
               orbital_type(nalfaorb) = 1  ! Alpha
            else
               ! Beta orbital
               nbetaorb = nbetaorb + 1
               occupation_factors(nbetaorb) = &
                   (0.5d0*Occ(orbital)*(1.0d0 - 0.5d0*Occ(orbital)))**(0.5d0) - &
                   2.0d0*0.5d0*Occ(orbital)*(1.0d0 - 0.5d0*Occ(orbital))
               orbital_type(nbetaorb) = 2  ! Beta
            end if
         end do
      else
         ! Closed shell case
         if (allocated(occupation_factors)) deallocate(occupation_factors)
         allocate(occupation_factors(noccmo))

         ! Compute occupation factors for closed shell
         do orbital = 1, noccmo
            occupation_factors(orbital) = &
                (0.5d0*Occ(orbital)*(1.0d0 - 0.5d0*Occ(orbital)))**(0.5d0) - &
                2.0d0*0.5d0*Occ(orbital)*(1.0d0 - 0.5d0*Occ(orbital))
         end do
      endif

      ! Cache DM2 elements for performance
      write(*,*) "Reading and caching DM2 elements..."
      ! First count the number of significant DM2 elements
      dm2_count = 0
      open(unit=5, file=dm2_file_name, form='unformatted', status='old')
      do while (.true.)
          read(5, end=99) i, j, k, l, DMval
          if (i.eq.0) exit
          if (abs(DMval) > 1.0D-12) then  ! Only cache significant elements
              dm2_count = dm2_count + 1
          endif
      end do
99    close(5)

      write(*,*) "Found", dm2_count, "significant DM2 elements"

      ! Allocate arrays to store DM2 elements
      allocate(dm2_indices(dm2_count, 4))
      allocate(dm2_values(dm2_count))

      ! Read and store DM2 elements
      open(unit=5, file=dm2_file_name, form='unformatted', status='old')
      idx = 0
      do while (.true.)
          read(5, end=100) i, j, k, l, DMval
          if (i.eq.0) exit
          if (abs(DMval) > 1.0D-12) then  ! Only cache significant elements
              idx = idx + 1
              dm2_indices(idx, 1) = i
              dm2_indices(idx, 2) = j
              dm2_indices(idx, 3) = k
              dm2_indices(idx, 4) = l

              ! Apply normalization if needed
              if (ontop_hf) then
                  dm2_values(idx) = DMval
              else
                  dm2_values(idx) = normalized_primitives(i) * normalized_primitives(j) * &
                                   normalized_primitives(k) * normalized_primitives(l) * DMval
              endif
          endif
      end do
100   close(5)

      ! Save cached DM2 data for use in property calculations
      call set_cached_dm2(dm2_indices, dm2_values, idx)

      ! Call generic scanning utility for 1D scan output
      call scan_and_print(dm2_file_name)

      ! Debug output
      write(*,*) "Primitive normalization factors:" !JD
      do i=1,min(5,max_primitive_index)
        write(*,*) "Primitive", i, "norm:", normalized_primitives(i) !JD
      end do

      ! Cleanup
      if (allocated(omg)) deallocate(omg)
      if (allocated(dm2_indices)) deallocate(dm2_indices)
      if (allocated(dm2_values)) deallocate(dm2_values)

      ! Note: Don't deallocate these as they need to persist for the cube generation
      ! if (allocated(normalized_primitives)) deallocate(normalized_primitives)
      ! if (allocated(occupation_factors)) deallocate(occupation_factors)
      ! if (allocated(orbital_type)) deallocate(orbital_type)

   end subroutine ontop

   ! Separate entry point for on-top calculations with cube output
   subroutine ontop_calculation(dm2_file_name, output_type)
      implicit none

      character(len=*), intent(in) :: dm2_file_name
      character(len=*), intent(in) :: output_type  ! "density" or "all"

      ! Process the DM2 file and set up calculations
      call ontop(dm2_file_name)

      ! Generate cube file for on-top density
      if (output_type == "all" .or. output_type == "on_top") then
         write(*,*) "Generating on-top density cube file..."
         call write_cube_file("ontop_density.cube", "on_top", &
                             get_ontop_value)
      end if

   contains
      ! Adapter function for cube file generation
      function get_ontop_value(x, y, z) result(value)
         double precision, intent(in) :: x, y, z
         double precision :: value

         value = get_cached_ontop_value(x, y, z)
      end function get_ontop_value
   end subroutine ontop_calculation

end module ontop_module

