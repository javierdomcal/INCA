module scan_utils
    use wfxinfo
    use geninfo
    use inputdat
    implicit none
    
    double precision, external :: Density, Prim, MO_a, MO_b, MoOr, Dens_a,Dens_b
    ! Module-level variables to be used in calculations
    double precision, allocatable, public :: normalized_primitives(:)
    double precision, allocatable, public :: occupation_factors(:)
    integer, allocatable, public :: orbital_type(:)
    
! Add to scan_utils.f90 module variables
logical, private :: cached_dm2 = .false.
integer, private :: cached_count = 0
integer, allocatable, private, dimension(:,:) :: cached_indices
double precision, allocatable, private, dimension(:) :: cached_values
contains
    subroutine scan_and_print(dm2_file)
        ! Input parameters
        character(len=*), intent(in) :: dm2_file
        
        ! Local variables
        integer :: iatom, idir, i1, n_points, scan_dim
        double precision :: scan_pos(3), atom_pos(3)
        double precision :: ontop_val, density_val, indicator_val
        
        ! Calculate number of points to scan
        n_points = nint(end_distance/step_size)
        
        ! Open output file and write header
        open(unit=6, file="ontop.dat")
        write(6,'(A)') "atom,direction,r,ontop,density,indicator_dynamic"
        
        ! Loop over atoms to scan
        do iatom = 1, n_atoms_scan
            ! Get atom position
            atom_pos = cartes(atom_indices(iatom),:)
            
            ! Loop through each specified direction
            do idir = 1, n_directions
                ! Determine which dimension to scan
                select case(scan_directions(idir))
                    case('x')
                        scan_dim = 1
                    case('y')
                        scan_dim = 2
                    case('z')
                        scan_dim = 3
                end select
                
                ! Scan in the specified direction
                do i1 = -n_points, n_points
                    ! Initialize scan position at atom position
                    scan_pos = atom_pos
                    
                    ! Modify the coordinate for the current scanning direction
                    scan_pos(scan_dim) = atom_pos(scan_dim) + i1*step_size
                    
                    ! Calculate ontop value
                    ontop_val = generate_ontop_values(scan_pos, dm2_file)
                    
                    ! Calculate density
                    if (uhf .or. opsh) then
                         density_val = Dens_a(scan_pos(1), scan_pos(2),scan_pos(3)) + Dens_b(scan_pos(1), scan_pos(2),scan_pos(3))
                    else
                         density_val = Density(scan_pos(1), scan_pos(2), scan_pos(3))
                    end if
                    ! Calculate indicator value
                    indicator_val = generate_indicator_values(scan_pos)
                    
                    ! Write results for this point
                    ! Calculate relative coordinate value
                    write(6,'(I0,A,A,A,F10.5,A,E16.8,A,E16.8,A,E16.8)') &
                        atom_indices(iatom), ',', scan_directions(idir), ',', &
                        (scan_pos(scan_dim) - atom_pos(scan_dim))*0.529177249d0, ',', &
                        ontop_val, ',', density_val, ',', indicator_val
                end do
            end do
        end do
        
        ! Close output file
        close(6)

       write(*,*) "=== Reference Point Values ===" !JD
scan_pos = (/0.0d0, 0.0d0, 0.0d0/)  ! Nuclear position !JD
ontop_val = generate_ontop_values(scan_pos, dm2_file)
write(*,*) "At nuclear position:", scan_pos !JD
write(*,*) "On-top value: ", ontop_val !JD
write(*,*) "Density value: ", Density(scan_pos(1), scan_pos(2), scan_pos(3)) !JD 
    end subroutine scan_and_print
    
    ! Function to generate on-top values
    function generate_ontop_values(scan_pos, dm2_file) result(on_top_val)
        ! Input parameters
        double precision, intent(in) :: scan_pos(3)
        character*40, intent(in) :: dm2_file
        double precision :: on_top_val
        
        ! Local variables
        integer :: i, j, k, l
        double precision :: DMval
! Keep track of how many elements we've printed
integer :: debug_count = 0 !JD
! Max elements to print
integer, parameter :: MAX_DEBUG = 10 !JD
        
        ! Initialize on-top value
        on_top_val = 0.d0
        
        ! Open DM2 file
        open(unit=5, file=dm2_file, form='unformatted', status='old')
      ! In generate_ontop_values function
write(*,*) "=== DM2 File Reading Debug ===" !JD
write(*,*) "Reading from file:", trim(dm2_file) !JD
write(*,*) "First few DM2 elements:" !JD



        ! Read and process DM2 file elements
        do while (.true.)
            read(5, end=99) i, j, k, l, DMval
            if (i.eq.0) exit
            
    if (ontop_hf) then
       DMval = DMval
    else
            ! Precomputed normalized primitives
            DMval = normalized_primitives(i) * normalized_primitives(j) * &
                    normalized_primitives(k) * normalized_primitives(l) * DMval
            
    end if
    ! Add debug prints for first few elements
    if (debug_count < MAX_DEBUG) then !JD
        write(*,*) "DM2(",i,j,k,l,") = ", DMval !JD
        debug_count = debug_count + 1 !JD
    endif

            ! Accumulate on-top value
            on_top_val = on_top_val + DMval * &
                        Prim(scan_pos(1),scan_pos(2),scan_pos(3),i) * &
                        Prim(scan_pos(1),scan_pos(2),scan_pos(3),j) * &
                        Prim(scan_pos(1),scan_pos(2),scan_pos(3),k) * &
                        Prim(scan_pos(1),scan_pos(2),scan_pos(3),l)
        end do
        
99      continue
        close(5)
        
    end function generate_ontop_values
    
    ! Function to generate indicator values
    function generate_indicator_values(scan_pos) result(indicator)
        ! Input parameters
        double precision, intent(in) :: scan_pos(3)
        double precision :: indicator
        
        ! Local variables
        integer :: orbital
        double precision :: orbital_value
        
        ! Initialize indicator
        indicator = 0.0d0
        
        ! Handle different electron configuration cases
        if (uhf .or. opsh) then
            ! Open shell case with separate alpha and beta orbitals
            do orbital = 1, size(occupation_factors)
                if (orbital_type(orbital) == 1) then
                    ! Alpha orbital
                    orbital_value = MO_a(scan_pos(1), scan_pos(2), scan_pos(3), orbital)
                else
                    ! Beta orbital
                    orbital_value = MO_b(scan_pos(1), scan_pos(2), scan_pos(3), orbital)
                endif
                
                ! Accumulate indicator
                indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
            end do
        else
            ! Closed shell case
            do orbital = 1, noccmo
                orbital_value = MoOr(scan_pos(1), scan_pos(2), scan_pos(3), orbital)
                indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
            end do
        endif
        
        ! Threshold check
        if (abs(indicator) .lt. 1.0D-11) then
            indicator = 0.d0
        endif


    end function generate_indicator_values


! Add to scan_utils.f90 module subprograms
subroutine set_cached_dm2(indices, values, count)
   implicit none
   integer, intent(in) :: count
   integer, dimension(count,4), intent(in) :: indices
   double precision, dimension(count), intent(in) :: values
   
   if (allocated(cached_indices)) deallocate(cached_indices)
   if (allocated(cached_values)) deallocate(cached_values)
   
   allocate(cached_indices(count,4))
   allocate(cached_values(count))
   
   cached_indices = indices
   cached_values = values
   cached_count = count
   cached_dm2 = .true.
end subroutine set_cached_dm2

subroutine clear_cached_dm2()
   implicit none
   if (allocated(cached_indices)) deallocate(cached_indices)
   if (allocated(cached_values)) deallocate(cached_values)
   cached_count = 0
   cached_dm2 = .false.
end subroutine clear_cached_dm2

function get_cached_ontop_value(x, y, z) result(ontop_val)
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: ontop_val
   integer :: i
   double precision :: prim_val_i, prim_val_j, prim_val_k, prim_val_l
   
   ontop_val = 0.0d0
   
   do i = 1, cached_count
      ! Calculate primitive values once for each point
      prim_val_i = Prim(x, y, z, cached_indices(i,1))
      prim_val_j = Prim(x, y, z, cached_indices(i,2))
      prim_val_k = Prim(x, y, z, cached_indices(i,3))
      prim_val_l = Prim(x, y, z, cached_indices(i,4))
      
      ! Calculate contribution to on-top value
      ontop_val = ontop_val + cached_values(i) * &
                 prim_val_i * prim_val_j * prim_val_k * prim_val_l
   end do
end function get_cached_ontop_value

function calculate_indicator_at_point(x, y, z) result(indicator)
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: indicator
   integer :: orbital
   double precision :: orbital_value
   double precision, external :: MO_a, MO_b, MoOr
   
   indicator = 0.0d0
   
   if (uhf .or. opsh) then
      do orbital = 1, size(occupation_factors)
         if (orbital_type(orbital) == 1) then
            orbital_value = MO_a(x, y, z, orbital)
         else
            orbital_value = MO_b(x, y, z, orbital)
         endif
         indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
      end do
   else
      do orbital = 1, noccmo
         orbital_value = MoOr(x, y, z, orbital)
         indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
      end do
   endif
   
   if (abs(indicator) < 1.0D-11) then
      indicator = 0.d0
   endif
end function calculate_indicator_at_point
end module scan_utils

