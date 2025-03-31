! pair_density.f90
! Module for calculating and analyzing pair density (two-electron density) properties
! This module handles full pair density, and its components (C1 and C2)

module pair_density_module
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use wfxinfo
   use geninfo
   use wfxinfo_hf
   use properties, only: register_pair_property, pair_property_func
   implicit none

   ! Module variables for normalization
   double precision, allocatable, private :: norm_prim(:)
   logical, private :: initialized = .false.

   ! DM2 file names and caching variables
   character(len=40), private :: dm2_file, dm2hf_file, dm2hfl_file

   ! Caching variables for DM2 elements
   integer, private :: cached_count = 0
   integer, allocatable, private, dimension(:,:) :: cached_indices
   double precision, allocatable, private, dimension(:) :: cached_values

   ! External function declarations
   double precision, external :: Prim, dfact_dp

contains

   ! Initialize the module
   subroutine initialize_pair_density(dm2_file, dm2hfname, dm2hflname)
      character(len=40), intent(in) :: dm2_file, dm2hfname, dm2hflname
      integer :: i
      procedure(pair_property_func), pointer :: full_pd_func => null()
      procedure(pair_property_func), pointer :: c1_pd_func => null()
      procedure(pair_property_func), pointer :: c2_pd_func => null()

      if (initialized) return

      ! Store file names
      dm2hf_file = dm2hfname
      dm2hfl_file = dm2hflname

      ! Calculate normalization of primitives
      allocate(norm_prim(nprim))
      do i = 1, nprim
         norm_prim(i) = (2.d0 * Alpha(i) / pi)**0.75d0 * &
                        sqrt((4.d0 * Alpha(i))**(dble(TMN(i, 1) + TMN(i, 2) + TMN(i, 3))) / &
                        (dble(dfact_dp(2 * TMN(i, 1) - 1) * dfact_dp(2 * TMN(i, 2) - 1) * dfact_dp(2 * TMN(i, 3) - 1))))
      end do

      ! Register pair density properties
      full_pd_func => full_pair_density
      c1_pd_func => c1_pair_density
      c2_pd_func => c2_pair_density

      call register_pair_property("Full Pair Density", &
                                 "Full two-electron density ρ₂(r₁,r₂)", full_pd_func)
      call register_pair_property("C1 Pair Density", &
                                 "C1 component of pair density (static correlation)", c1_pd_func)
      call register_pair_property("C2 Pair Density", &
                                 "C2 component of pair density (dynamic correlation)", c2_pd_func)

      initialized = .true.
      write(*,*) "Pair density module initialized."
   end subroutine initialize_pair_density

   ! Cache DM2 elements for a specific file
   subroutine cache_dm2_file(filename)
      character(len=40), intent(in) :: filename
      integer :: i, j, k, l, count_elements, idx
      double precision :: DMval

      write(*,*) "Caching DM2 elements from ", trim(filename)

      ! Count significant elements
      count_elements = 0
      open(unit=5, file=filename, form='unformatted', status='old')
      do while (.true.)
         read(5, end=99) i, j, k, l, DMval
         if (i == 0) exit
         if (abs(DMval) > 1.0D-12) then
            count_elements = count_elements + 1
         endif
      end do
99    close(5)

      write(*,*) "Found ", count_elements, " significant elements."

      ! Allocate arrays
      if (allocated(cached_indices)) deallocate(cached_indices)
      if (allocated(cached_values)) deallocate(cached_values)
      allocate(cached_indices(count_elements, 4))
      allocate(cached_values(count_elements))

      ! Read and store elements
      idx = 0
      open(unit=5, file=filename, form='unformatted', status='old')
      do while (.true.)
         read(5, end=100) i, j, k, l, DMval
         if (i == 0) exit
         if (abs(DMval) > 1.0D-12) then
            idx = idx + 1
            cached_indices(idx, 1:4) = [i, j, k, l]
            cached_values(idx) = DMval
         endif
      end do
100   close(5)

      cached_count = idx
      write(*,*) "Successfully cached ", cached_count, " DM2 elements."
   end subroutine cache_dm2_file

   ! Implementation of full pair density for property system
   function full_pair_density(r1, r2) result(value)
      real(dp), dimension(3), intent(in) :: r1, r2
      real(dp) :: value

      ! Cache DM2 elements if not already done
      call cache_dm2_file(dm2_file)

      ! Calculate pair density
      value = compute_pair_density_cached(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3))
   end function full_pair_density

   ! Implementation of C1 pair density (static correlation)
   function c1_pair_density(r1, r2) result(value)
      real(dp), dimension(3), intent(in) :: r1, r2
      real(dp) :: value
      real(dp) :: pd_hfl, pd_hf

      ! First get HFL pair density
      call cache_dm2_file(dm2hfl_file)
      pd_hfl = compute_pair_density_cached(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3))

      ! Then get HF pair density
      call cache_dm2_file(dm2hf_file)
      pd_hf = compute_pair_density_cached(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3))

      ! C1 = HFL - HF
      value = pd_hfl - pd_hf
   end function c1_pair_density

   ! Implementation of C2 pair density (dynamic correlation)
   function c2_pair_density(r1, r2) result(value)
      real(dp), dimension(3), intent(in) :: r1, r2
      real(dp) :: value
      real(dp) :: pd_full, pd_hfl

      ! First get full pair density
      call cache_dm2_file(dm2_file)
      pd_full = compute_pair_density_cached(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3))

      ! Then get HFL pair density
      call cache_dm2_file(dm2hfl_file)
      pd_hfl = compute_pair_density_cached(r1(1), r1(2), r1(3), r2(1), r2(2), r2(3))

      ! C2 = Full - HFL
      value = pd_full - pd_hfl
   end function c2_pair_density

   ! Main calculation function that uses cached DM2 elements
   function compute_pair_density_cached(x1, y1, z1, x2, y2, z2) result(pd_value)
      double precision, intent(in) :: x1, y1, z1, x2, y2, z2
      double precision :: pd_value
      integer :: idx, i, j, k, l
      double precision :: prim_i, prim_j, prim_k, prim_l, n_prim_t

      pd_value = 0.d0

      do idx = 1, cached_count
         i = cached_indices(idx, 1)
         j = cached_indices(idx, 2)
         k = cached_indices(idx, 3)
         l = cached_indices(idx, 4)

         n_prim_t = norm_prim(i) * norm_prim(j) * norm_prim(k) * norm_prim(l)

         prim_i = Prim(x1, y1, z1, i)
         prim_j = Prim(x2, y2, z2, j)
         prim_k = Prim(x1, y1, z1, k)
         prim_l = Prim(x2, y2, z2, l)

         pd_value = pd_value + cached_values(idx) * n_prim_t * &
                   prim_i * prim_j * prim_k * prim_l
      end do

      ! Apply threshold for very small values
      if (abs(pd_value) < 1.0D-11) then
         pd_value = 0.d0
      endif
   end function compute_pair_density_cached

   ! Function for direct use that takes separate coordinates
   function compute_pair_density(x1, y1, z1, x2, y2, z2) result(pd_value)
      double precision, intent(in) :: x1, y1, z1, x2, y2, z2
      double precision :: pd_value

      ! Make sure we're initialized
      if (.not. initialized) then
         write(*,*) "Error: Pair density module not initialized"
         pd_value = 0.d0
         return
      end if

      ! Use cached version
      pd_value = compute_pair_density_cached(x1, y1, z1, x2, y2, z2)
   end function compute_pair_density

   ! Adapter function for cube2 file output
   function full_pd_adapter(x1, y1, z1, x2, y2, z2) result(value)
      double precision, intent(in) :: x1, y1, z1, x2, y2, z2
      double precision :: value

      value = compute_pair_density(x1, y1, z1, x2, y2, z2)
   end function full_pd_adapter

   ! Adapter function for cube2 file output
   function c1_pd_adapter(x1, y1, z1, x2, y2, z2) result(value)
      double precision, intent(in) :: x1, y1, z1, x2, y2, z2
      double precision :: value
      real(dp) :: r1(3), r2(3)

      r1 = [real(dp) :: x1, y1, z1]
      r2 = [real(dp) :: x2, y2, z2]
      value = c1_pair_density(r1, r2)
   end function c1_pd_adapter

   ! Adapter function for cube2 file output
   function c2_pd_adapter(x1, y1, z1, x2, y2, z2) result(value)
      double precision, intent(in) :: x1, y1, z1, x2, y2, z2
      double precision :: value
      real(dp) :: r1(3), r2(3)

      r1 = [real(dp) :: x1, y1, z1]
      r2 = [real(dp) :: x2, y2, z2]
      value = c2_pair_density(r1, r2)
   end function c2_pd_adapter


   ! Clean up the module
   subroutine cleanup_pair_density()
      if (allocated(norm_prim)) deallocate(norm_prim)
      if (allocated(cached_indices)) deallocate(cached_indices)
      if (allocated(cached_values)) deallocate(cached_values)
      initialized = .false.
      write(*,*) "Pair density module cleaned up."
   end subroutine cleanup_pair_density

   ! Entry point for pair density calculations
   subroutine pair_density_calculation(dm2_file, dm2hfname, dm2hflname, output_type, grid, grid_2)
      use cube_module        ! For cube file output
      use properties  ! For property access

      character(len=40), intent(in) :: dm2_file, dm2hfname, dm2hflname
      character(len=10), intent(in) :: output_type  ! "full", "c1", "c2", "all"
      type(grid), intent(in) :: grid, grid_2

      ! Initialize pair density module
      call initialize_pair_density(dm2_file, dm2hfname, dm2hflname)

      ! Enable properties based on output_type
      if (output_type == "full" .or. output_type == "all") then
         call enable_property("Full Pair Density")
      end if

      if (output_type == "c1" .or. output_type == "all") then
         call enable_property("C1 Pair Density")
      end if

      if (output_type == "c2" .or. output_type == "all") then
         call enable_property("C2 Pair Density")
      end if

      ! Generate cube files based on enabled properties
      if (is_property_enabled("Full Pair Density")) then
         write(*,*) "Generating full pair density cube2 file..."
         call write_cube_file_pair("full_pair_density.cube2", "Full Pair Density", &
                              full_pd_adapter, grid, grid_2)
      end if

      if (is_property_enabled("C1 Pair Density")) then
         write(*,*) "Generating C1 pair density cube2 file..."
         call write_cube_file_pair("c1_pair_density.cube2", "C1 Pair Density", &
                              c1_pd_adapter, grid, grid_2)
      end if

      if (is_property_enabled("C2 Pair Density")) then
         write(*,*) "Generating C2 pair density cube2 file..."
         call write_cube_file_pair("c2_pair_density.cube2", "C2 Pair Density", &
                              c2_pd_adapter, grid, grid_2)
      end if

      ! Cleanup
      call cleanup_pair_density()
   end subroutine pair_density_calculation

   subroutine pair_density_nucleus_calculation(dm2_file, dm2hfname, dm2hflname, output_type, grid)
      use cube_module        ! For cube file output
      use properties  ! For property access

      character(len=40), intent(in) :: dm2_file, dm2hfname, dm2hflname
      character(len=10), intent(in) :: output_type  ! "full", "c1", "c2", "all"
      type(grid), intent(in) :: grid
      type(grid) :: grid_2
      integer ::i

      do i = 1,3  !
         grid_2%max_vals(i) = 0.0_dp
         grid_2%step_sizes(i) = 0.0_dp
      enddo

      ! Initialize pair density module
      call initialize_pair_density(dm2_file, dm2hfname, dm2hflname)

      ! Enable properties based on output_type
      if (output_type == "full" .or. output_type == "all") then
         call enable_property("Full Pair Nucleus Density")
      end if

      if (output_type == "c1" .or. output_type == "all") then
         call enable_property("C1 Pair Nucleus Density")
      end if

      if (output_type == "c2" .or. output_type == "all") then
         call enable_property("C2 Pair Nucleus Density")
      end if

      ! Generate cube files based on enabled properties
      if (is_property_enabled("Full Pair Nucleus Density")) then
         write(*,*) "Generating full pair Nucleus Density cube file..."
         call write_cube_file_pair("full_pair_density_nucleus.cube", "Full Pair Nucleus Density", &
                              full_pd_adapter, grid, grid_2)
      end if

      if (is_property_enabled("C1 Pair Nucleus Density")) then
         write(*,*) "Generating C1 pair Nucleus Density cube file..."
         call write_cube_file_pair("c1_pair_density_nucleus.cube", "C1 Pair Nucleus Density", &
                              c1_pd_adapter, grid, grid_2)
      end if

      if (is_property_enabled("C2 Pair Nucleus Density")) then
         write(*,*) "Generating C2 pair Nucleus Density cube file..."
         call write_cube_file_pair("c2_pair_density_nucleus.cube", "C2 Pair Nucleus Density", &
                              c2_pd_adapter, grid, grid_2)
      end if

      ! Cleanup
      call cleanup_pair_density()
   end subroutine pair_density_nucleus_calculation


end module pair_density_module
