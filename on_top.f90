subroutine ontop(dm2_file)
   use wfxinfo     ! Contains cartes (atomic positions), TMN, Alpha arrays
   use geninfo     ! Contains general info like pi
   use inputdat    ! Contains n_atoms_scan, n_directions, scan_directions, atom_indices
   use scan_utils  ! Will contain generic scanning utility
   use on_top_cube  ! Add this import
   implicit none
   
   ! Input parameters
   character*40, intent(in) :: dm2_file
   
   ! Local variables for normalization factor computation
   integer :: i, max_primitive_index, orbital
   double precision :: norm_factor_1

   double precision :: dfact_dp
   double precision, allocatable :: omg(:) 
   
   ! Determine max primitive index
   max_primitive_index = size(Alpha)
   allocate(omg(noccmo))
   
   ! Allocate and compute normalized primitives
   if (allocated(normalized_primitives)) deallocate(normalized_primitives)
   allocate(normalized_primitives(max_primitive_index))
   
   ! Precompute normalization factors for primitives
   do i = 1, max_primitive_index
      ! Compute double factorial product
      norm_factor_1 = dble(dfact_dp(2*TMN(i,1)-1) * &
                           dfact_dp(2*TMN(i,2)-1) * &
                           dfact_dp(2*TMN(i,3)-1))
      
      ! Complete normalization factor
      if (ontop_hf) then
          normalized_primitives(i) = 1.0d0
     else
      normalized_primitives(i) = &
         (2.d0*Alpha(i)/pi)**(0.75d0) * &
         sqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))) / &
         norm_factor_1
     end if
   end do
   
   ! Determine orbital configuration and compute occupation factors
   if (uhf .or. opsh) then
      ! Open shell case
      ! Determine number of alpha and beta orbitals
      nalfaorb = 0
      nbetaorb = 0
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
   if (ontop_cube_output) then
      call ontop_cube_write(dm2_file)
   endif

   ! Call generic scanning utility
   call scan_and_print(dm2_file)
        write(*,*) "Primitive normalization factors:" !JD
        do i=1,5
        write(*,*) "Primitive", i, "norm:", normalized_primitives(i) !JD
        end do   
   ! Cleanup
   if (allocated(normalized_primitives)) deallocate(normalized_primitives)
   if (allocated(occupation_factors)) deallocate(occupation_factors)
   if (allocated(orbital_type)) deallocate(orbital_type)

end subroutine ontop

