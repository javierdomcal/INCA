! ontop_cube.f90
! Module for generating cube files for on-top pair density, density, and indicator dynamic
! This module provides functions to generate standard cube files containing 3D volumetric data
! for visualization and analysis in molecular modeling software

module on_top_cube
    use wfxinfo
    use geninfo
    use cubeinfo
    use inputdat
    use scan_utils, only: normalized_primitives, occupation_factors, orbital_type
    implicit none

    ! External function declarations
    double precision, external :: Density, Prim, MO_a, MO_b, MoOr, Dens_a, Dens_b

    ! Module variables for caching DM2 elements
    logical, private :: cached_dm2 = .false.
    integer, private :: cached_count = 0
    integer, allocatable, private, dimension(:,:) :: cached_indices
    double precision, allocatable, private, dimension(:) :: cached_values

contains

    !----------------------------------------------------------------------
    ! Main subroutine to generate all cube files
    !----------------------------------------------------------------------
    subroutine ontop_cube_write(dm2_file)
        implicit none
        
        character*40, intent(in) :: dm2_file
        
        ! Local variables
        double precision, dimension(3) :: cube_center, cube_step
        integer, dimension(3) :: cube_npoints
        
        ! For caching DM2 elements
        integer :: i, j, k, l, dm2_count, idx
        double precision :: DMval
        integer, allocatable, dimension(:,:) :: dm2_indices
        double precision, allocatable, dimension(:) :: dm2_values
        
        ! Get grid parameters (fixed medium resolution)
        call get_cube_parameters(cube_center, cube_step, cube_npoints)
        
        ! Set up cube file parameters
        center = cube_center
        step = cube_step
        np = cube_npoints
        
        write(*,*) "Preparing cube files with grid dimensions:", np(1), "x", np(2), "x", np(3)
        write(*,*) "Grid spacing:", step(1), "bohr"
        
        ! Always cache DM2 elements for performance
        ! First count the number of significant DM2 elements
        write(*,*) "Reading and caching DM2 elements..."
        dm2_count = 0
        open(unit=5, file=dm2_file, form='unformatted', status='old')
        do while (.true.)
            read(5, end=99) i, j, k, l, DMval
            if (i.eq.0) exit
            if (abs(DMval) > 1.0D-12) then  ! Only cache significant elements
                dm2_count = dm2_count + 1
            endif
        end do
99      close(5)
        
        write(*,*) "Found", dm2_count, "significant DM2 elements"
        
        ! Allocate arrays to store DM2 elements
        allocate(dm2_indices(dm2_count, 4))
        allocate(dm2_values(dm2_count))
        
        ! Read and store DM2 elements
        open(unit=5, file=dm2_file, form='unformatted', status='old')
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
100     close(5)
        
        ! Set cached DM2 data for use by calculation functions
        call set_cached_dm2(dm2_indices, dm2_values, idx)
        
        ! Generate cube files for each property
        write(*,*) "Generating on-top pair density cube file..."
        call write_ontop_cube("ontop.cube", dm2_file)
        
        write(*,*) "Generating electron density cube file..."
        call write_density_cube("density.cube")
        
        write(*,*) "Generating indicator dynamic cube file..."
        call write_indicator_cube("ID.cube")
        
        ! Cleanup cached data
        call clear_cached_dm2()
        deallocate(dm2_indices)
        deallocate(dm2_values)
        
        write(*,*) "Cube files generated successfully."
    end subroutine ontop_cube_write

    !----------------------------------------------------------------------
    ! Get appropriate cube file grid parameters
    !----------------------------------------------------------------------
    subroutine get_cube_parameters(center, step, npoints)
        implicit none
        
        double precision, dimension(3), intent(out) :: center, step
        integer, dimension(3), intent(out) :: npoints
        
        double precision :: molecule_size, max_extent
        integer :: i
        
        ! Calculate molecule center and size
        center = 0.0d0
        max_extent = 0.0d0
        
        do i = 1, natoms
            center = center + cartes(i,:)
        end do
        center = center / dble(natoms)
        
        ! Find maximum extent of molecule from center
        do i = 1, natoms
            max_extent = max(max_extent, sqrt(sum((cartes(i,:) - center)**2)))
        end do
        
        ! Add buffer around molecule
        max_extent = max_extent + 5.0d0  ! Add 5 bohr buffer
        
        ! Set parameters for medium resolution (fixed)
        step = (/0.25d0, 0.25d0, 0.25d0/)
        
        ! Calculate number of points needed to cover the molecule
        do i = 1, 3
            npoints(i) = ceiling(2.0d0 * max_extent / step(i))
            ! Make sure it's odd for centered grid
            if (mod(npoints(i), 2) == 0) npoints(i) = npoints(i) + 1
        end do
    end subroutine get_cube_parameters

    !----------------------------------------------------------------------
    ! Write on-top pair density cube file
    !----------------------------------------------------------------------
    subroutine write_ontop_cube(filename, dm2_file)
        implicit none
        
        character*(*), intent(in) :: filename
        character*40, intent(in) :: dm2_file
        
        ! Local variables
        integer :: i, j, k, ia, icount, imod
        double precision :: x, y, z, ontop_value
        double precision, dimension(3) :: xm
        double precision, parameter :: zero = 0.0d0
        
        ! Open cube file
        open(unit=2, file=filename)
        
        ! Write header
        write(2,*) "CUBE FILE: On-top pair density"
        write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
        write(2,*) natoms, center(1), center(2), center(3)
        write(2,*) np(1), step(1), zero, zero
        write(2,*) np(2), zero, step(2), zero
        write(2,*) np(3), zero, zero, step(3)
        
        ! Write atom positions
        do ia = 1, natoms
            write(2,*) an(ia), chrg(ia), cartes(ia,1), cartes(ia,2), cartes(ia,3)
        end do
        
        ! Calculate grid starting point
        do i = 1, 3
            if (mod(np(i), 2) .eq. 0) then ! even number
                xm(i) = center(i) - (step(i)/2.0d0) - ((np(i)-2)/2) * step(i)
            else ! odd number
                xm(i) = center(i) - ((np(i)-1)/2) * step(i)
            end if
        end do
        
        ! Write volumetric data
        icount = 0
        imod = max(1, np(1)*np(2)*np(3)/20) ! For progress reporting
        
        do i = 1, np(1)
            x = xm(1) + step(1) * (i-1)
            do j = 1, np(2)
                y = xm(2) + step(2) * (j-1)
                do k = 1, np(3)
                    z = xm(3) + step(3) * (k-1)
                    
                    ! Progress reporting
                    icount = icount + 1
                    if (mod(icount, imod) == 0) then
                        write(*,'(A,F6.2,A)') "  Progress: ", &
                                100.0d0*dble(icount)/(np(1)*np(2)*np(3)), "%"
                    endif
                    
                    ! Calculate on-top value at this point
                    ontop_value = get_cached_ontop_value(x, y, z)
                    
                    ! Write value to cube file with proper format
                    write(2, '(E16.6E3)') ontop_value
                end do
            end do
        end do
        
        close(2)
    end subroutine write_ontop_cube

    !----------------------------------------------------------------------
    ! Write electron density cube file
    !----------------------------------------------------------------------
    subroutine write_density_cube(filename)
        implicit none
        
        character*(*), intent(in) :: filename
        
        ! Local variables
        integer :: i, j, k, ia, icount, imod
        double precision :: x, y, z, density_value
        double precision, dimension(3) :: xm
        double precision, parameter :: zero = 0.0d0
        
        ! Open cube file
        open(unit=2, file=filename)
        
        ! Write header
        write(2,*) "CUBE FILE: Electron density"
        write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
        write(2,*) natoms, center(1), center(2), center(3)
        write(2,*) np(1), step(1), zero, zero
        write(2,*) np(2), zero, step(2), zero
        write(2,*) np(3), zero, zero, step(3)
        
        ! Write atom positions
        do ia = 1, natoms
            write(2,*) an(ia), chrg(ia), cartes(ia,1), cartes(ia,2), cartes(ia,3)
        end do
        
        ! Calculate grid starting point
        do i = 1, 3
            if (mod(np(i), 2) .eq. 0) then ! even number
                xm(i) = center(i) - (step(i)/2.0d0) - ((np(i)-2)/2) * step(i)
            else ! odd number
                xm(i) = center(i) - ((np(i)-1)/2) * step(i)
            end if
        end do
        
        ! Write volumetric data
        icount = 0
        imod = max(1, np(1)*np(2)*np(3)/20) ! For progress reporting
        
        do i = 1, np(1)
            x = xm(1) + step(1) * (i-1)
            do j = 1, np(2)
                y = xm(2) + step(2) * (j-1)
                do k = 1, np(3)
                    z = xm(3) + step(3) * (k-1)
                    
                    ! Progress reporting
                    icount = icount + 1
                    if (mod(icount, imod) == 0) then
                        write(*,'(A,F6.2,A)') "  Progress: ", &
                                100.0d0*dble(icount)/(np(1)*np(2)*np(3)), "%"
                    endif
                    
                    ! Calculate density value at this point
                    if (uhf .or. opsh) then
                        density_value = Dens_a(x, y, z) + Dens_b(x, y, z)
                    else
                        density_value = Density(x, y, z)
                    endif
                    
                    ! Write value to cube file
                    write(2, '(E16.6E3)') density_value
                end do
            end do
        end do
        
        close(2)
    end subroutine write_density_cube

    !----------------------------------------------------------------------
    ! Write indicator dynamic cube file
    !----------------------------------------------------------------------
    subroutine write_indicator_cube(filename)
        implicit none
        
        character*(*), intent(in) :: filename
        
        ! Local variables
        integer :: i, j, k, ia, icount, imod
        double precision :: x, y, z, indicator_value
        double precision, dimension(3) :: xm
        double precision, parameter :: zero = 0.0d0
        
        ! Open cube file
        open(unit=2, file=filename)
        
        ! Write header
        write(2,*) "CUBE FILE: Indicator dynamic"
        write(2,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
        write(2,*) natoms, center(1), center(2), center(3)
        write(2,*) np(1), step(1), zero, zero
        write(2,*) np(2), zero, step(2), zero
        write(2,*) np(3), zero, zero, step(3)
        
        ! Write atom positions
        do ia = 1, natoms
            write(2,*) an(ia), chrg(ia), cartes(ia,1), cartes(ia,2), cartes(ia,3)
        end do
        
        ! Calculate grid starting point
        do i = 1, 3
            if (mod(np(i), 2) .eq. 0) then ! even number
                xm(i) = center(i) - (step(i)/2.0d0) - ((np(i)-2)/2) * step(i)
            else ! odd number
                xm(i) = center(i) - ((np(i)-1)/2) * step(i)
            end if
        end do
        
        ! Write volumetric data
        icount = 0
        imod = max(1, np(1)*np(2)*np(3)/20) ! For progress reporting
        
        do i = 1, np(1)
            x = xm(1) + step(1) * (i-1)
            do j = 1, np(2)
                y = xm(2) + step(2) * (j-1)
                do k = 1, np(3)
                    z = xm(3) + step(3) * (k-1)
                    
                    ! Progress reporting
                    icount = icount + 1
                    if (mod(icount, imod) == 0) then
                        write(*,'(A,F6.2,A)') "  Progress: ", &
                                100.0d0*dble(icount)/(np(1)*np(2)*np(3)), "%"
                    endif
                    
                    ! Calculate indicator value at this point
                    indicator_value = calculate_indicator_at_point(x, y, z)
                    
                    ! Write value to cube file
                    write(2, '(E16.6E3)') indicator_value
                end do
            end do
        end do
        
        close(2)
    end subroutine write_indicator_cube

    !----------------------------------------------------------------------
    ! Caching functions for performance optimization
    !----------------------------------------------------------------------
    
    ! Set cached DM2 data
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

    ! Clear cached DM2 data
    subroutine clear_cached_dm2()
        implicit none
        if (allocated(cached_indices)) deallocate(cached_indices)
        if (allocated(cached_values)) deallocate(cached_values)
        cached_count = 0
        cached_dm2 = .false.
    end subroutine clear_cached_dm2

    !----------------------------------------------------------------------
    ! Property calculation functions
    !----------------------------------------------------------------------
    
    ! Calculate on-top pair density value at a given point using cached DM2 data
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
        
        ! Apply threshold for very small values
        if (abs(ontop_val) < 1.0D-11) then
            ontop_val = 0.d0
        endif
    end function get_cached_ontop_value

    ! Calculate indicator dynamic value at a given point
    function calculate_indicator_at_point(x, y, z) result(indicator)
        implicit none
        double precision, intent(in) :: x, y, z
        double precision :: indicator
        integer :: orbital
        double precision :: orbital_value
        
        indicator = 0.0d0
        
        if (uhf .or. opsh) then
            ! Open shell case with separate alpha and beta orbitals
            do orbital = 1, size(occupation_factors)
                if (orbital_type(orbital) == 1) then
                    ! Alpha orbital
                    orbital_value = MO_a(x, y, z, orbital)
                else
                    ! Beta orbital
                    orbital_value = MO_b(x, y, z, orbital)
                endif
                
                ! Accumulate indicator
                indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
            end do
        else
            ! Closed shell case
            do orbital = 1, noccmo
                orbital_value = MoOr(x, y, z, orbital)
                indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
            end do
        endif
        
        ! Apply threshold for very small values
        if (abs(indicator) < 1.0D-11) then
            indicator = 0.d0
        endif
    end function calculate_indicator_at_point

end module on_top_cube
