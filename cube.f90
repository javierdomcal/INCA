! This module provides utilities for generating cube files, which are commonly used
! in computational chemistry and physics to represent volumetric data on a grid.
! It includes subroutines for writing cube files based on scalar and pairwise property functions.

module cube_module
    use inputdat
    use geninfo
    implicit none


contains

    subroutine write_cube_file(filename, property_name, property_function, grid_info)
        character(*), intent(in) :: filename, property_name
        interface
            function property_function(x, y, z) result(value)
                double precision, intent(in) :: x, y, z
                double precision :: value
            end function property_function
        end interface
        type(grid_parameters), intent(in) :: grid_info

        integer :: i, j, k, icount, iatom
        double precision :: x, y, z
        double precision, allocatable, dimension(:) :: origin
        integer, allocatable, dimension(:) :: dims
        integer :: unit_number

        ! Check if grid_info has the correct dimensionality
        if (grid_info%ndims /= 3) then
            write(*,*) "Error: write_cube_file requires 3D grid parameters"
            return
        end if

        ! Convert grid info to dimensions and origin
        allocate(dims(3), origin(3))
        dims = 2*int(grid_info%max_vals / grid_info%step_sizes)+1
        origin = -grid_info%max_vals

        open(newunit=unit_number, file=filename)

        ! Write header
        write(unit_number, '(A)') trim(property_name)
        write(unit_number, '(A)') "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
        write(unit_number, '(I5,3F12.6)') natoms, origin(1), origin(2), origin(3)

        do i = 1, 3
            if (grid_info%active_dims(i)) then
                if (i == 1) then
                    write(unit_number, '(I5,3F12.6)') dims(i), grid_info%step_sizes(i), 0.0d0, 0.0d0
                else if (i == 2) then
                    write(unit_number, '(I5,3F12.6)') dims(i), 0.0d0, grid_info%step_sizes(i), 0.0d0
                else
                    write(unit_number, '(I5,3F12.6)') dims(i), 0.0d0, 0.0d0, grid_info%step_sizes(i)
                end if
            else
                write(unit_number, '(I5,3F12.6)') 1, 0.0d0, 0.0d0, 0.0d0
            end if
        end do

        ! Write atom information
        do iatom = 1, natoms
            write(unit_number, '(I5,4F12.6)') an(iatom), chrg(iatom), &
                                             cartes(iatom,1), cartes(iatom,2), cartes(iatom,3)
        end do

        ! Write volumetric data
        icount = 0
        do i = 1, dims(1)
            x = origin(1) + (i-1) * grid_info%step_sizes(1)
            do j = 1, dims(2)
                y = origin(2) + (j-1) * grid_info%step_sizes(2)
                do k = 1, dims(3)
                    z = origin(3) + (k-1) * grid_info%step_sizes(3)

                    ! Calculate value at this point
                    write(unit_number, '(E16.8)', advance="no") property_function(x, y, z)
                    icount = icount + 1
                    if (mod(icount, 6) == 0) then
                        write(unit_number, *)  ! New line after every 6 values
                    end if
                end do
            end do
        end do
        if (mod(icount, 6) /= 0) write(unit_number, *)  ! Ensure final newline

        close(unit_number)
        deallocate(dims, origin)
    end subroutine write_cube_file

    subroutine write_cube2_file(filename, property_name, pair_property_function, grid_info)
        character(*), intent(in) :: filename, property_name
        interface
            function pair_property_function(x1, y1, z1, x2, y2, z2) result(value)
                double precision, intent(in) :: x1, y1, z1, x2, y2, z2
                double precision :: value
            end function pair_property_function
        end interface
        type(grid_parameters), intent(in) :: grid_info

        integer :: i1, j1, k1, i2, j2, k2, icount, iatom
        double precision :: x1, y1, z1, x2, y2, z2
        double precision, allocatable, dimension(:) :: origin
        integer, allocatable, dimension(:) :: dims
        integer :: unit_number


        ! Check if grid_info has the correct dimensionality
        if (grid_info%ndims /= 6) then
            write(*,*) "Error: write_cube2_file requires 6D grid parameters"
            return
        end if

        ! Convert grid info to dimensions and origin
        allocate(dims(6), origin(6))
        dims = 2*int(grid_info%max_vals / grid_info%step_sizes)+1
        origin = -grid_info%max_vals

        open(newunit=unit_number, file=filename)

        ! Write header
        write(unit_number, '(A)') trim(property_name)
        write(unit_number, '(A)') "OUTER LOOP: X1, MIDDLE LOOP: Y1, INNER LOOP: Z1, X2, Y2, Z2"
        write(unit_number, '(I5,3F12.6)') natoms, origin(1), origin(2), origin(3), origin(4), origin(5), origin(6)

        ! Write grid specifications for the first 3 dimensions (x1,y1,z1)
        do i1 = 1, 3
            if (grid_info%active_dims(i1)) then
                if (i1 == 1) then
                    write(unit_number, '(I5,3F12.6)') dims(i1), grid_info%step_sizes(i1), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
                else if (i1 == 2) then
                    write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, grid_info%step_sizes(i1), 0.0d0, 0.0d0, 0.0d0, 0.0d0
                else if (i1 == 3) then
                    write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, 0.0d0, grid_info%step_sizes(i1), 0.0d0, 0.0d0, 0.0d0
                else if (i1 == 4) then
                    write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, 0.0d0, 0.0d0, grid_info%step_sizes(i1), 0.0d0, 0.0d0
                else if (i1 == 5) then
                    write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, 0.0d0, 0.0d0, 0.0d0, grid_info%step_sizes(i1), 0.0d0
                else
                    write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, grid_info%step_sizes(i1)
                end if
            else
                write(unit_number, '(I5,3F12.6)') 1, 0.0d0, 0.0d0, 0.0d0
            end if
        end do

        ! Write atom information
        do iatom = 1, natoms
            write(unit_number, '(I5,4F12.6)') an(iatom), chrg(iatom), &
                                             cartes(iatom,1), cartes(iatom,2), cartes(iatom,3)
        end do

        ! Write additional header for second set of coordinates (optional)
        write(unit_number, '(A)') "SECOND ELECTRON COORDINATES"
        write(unit_number, '(3F12.6)') origin(4), origin(5), origin(6)

        ! Write grid specifications for the second 3 dimensions (x2,y2,z2)
        do i2 = 4, 6
            if (grid_info%active_dims(i2)) then
                write(unit_number, '(I5,3F12.6)') dims(i2), &
                                     grid_info%step_sizes(i2), 0.0d0, 0.0d0
            else
                write(unit_number, '(I5,3F12.6)') 1, 0.0d0, 0.0d0, 0.0d0
            end if
        end do

        ! Write volumetric data
        icount = 0
        do i1 = 1, dims(1)
            x1 = origin(1) + (i1-1) * grid_info%step_sizes(1)
            do j1 = 1, dims(2)
                y1 = origin(2) + (j1-1) * grid_info%step_sizes(2)
                do k1 = 1, dims(3)
                    z1 = origin(3) + (k1-1) * grid_info%step_sizes(3)
                    do i2 = 1, dims(4)
                        x2 = origin(4) + (i2-1) * grid_info%step_sizes(4)
                        do j2 = 1, dims(5)
                            y2 = origin(5) + (j2-1) * grid_info%step_sizes(5)
                            do k2 = 1, dims(6)
                                z2 = origin(6) + (k2-1) * grid_info%step_sizes(6)

                                ! Calculate value at this point
                                write(unit_number, '(E16.8)', advance="no") pair_property_function(x1, y1, z1, x2, y2, z2)
                                icount = icount + 1
                                if (mod(icount, 6) == 0) then
                                    write(unit_number, *)  ! New line after every 6 values
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
        if (mod(icount, 6) /= 0) write(unit_number, *)  ! Ensure final newline

        close(unit_number)
        deallocate(dims, origin)
    end subroutine write_cube2_file

end module cube_module
