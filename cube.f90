! This module provides utilities for generating cube files, which are commonly used
! in computational chemistry and physics to represent volumetric data on a grid.
! It includes subroutines for writing cube files based on scalar and pairwise property functions.

module cube_module
    use inputdat
    use geninfo
    implicit none


contains

    subroutine write_cube_file(filename, property_name, property_function)
        character(*), intent(in) :: filename, property_name
        interface
            function property_function(x, y, z) result(value)
                double precision, intent(in) :: x, y, z
                double precision :: value
            end function property_function
        end interface


        integer :: i, j, k, icount, iatom
        double precision :: x, y, z
        double precision, dimension(3) :: origin
        integer, dimension(3) :: dims
        integer :: unit_number



        ! Convert grid info to dimensions and origin
        do i = 1, 3
            origin(i) = -grid%max_vals(i)
            if (grid%step_sizes(i) .le. 0.00001d0) then
                dims(i) = 1
            else
                dims(i) = 2*int(grid%max_vals(i) / grid%step_sizes(i))+1
            end if
        end do


        open(newunit=unit_number, file=filename)

        ! Write header
        write(unit_number, '(A)') trim(property_name)
        write(unit_number, '(A)') "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
        write(unit_number, '(I5,3F12.6)') natoms, origin(1), origin(2), origin(3)

        do i = 1, 3
            if (i == 1) then
                write(unit_number, '(I5,3F12.6)') dims(i), grid%step_sizes(i), 0.0d0, 0.0d0
            else if (i == 2) then
                write(unit_number, '(I5,3F12.6)') dims(i), 0.0d0, grid%step_sizes(i), 0.0d0
            else
                write(unit_number, '(I5,3F12.6)') dims(i), 0.0d0, 0.0d0, grid%step_sizes(i)
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
            x = origin(1) + (i-1) * grid%step_sizes(1)
            do j = 1, dims(2)
                y = origin(2) + (j-1) * grid%step_sizes(2)
                do k = 1, dims(3)
                    z = origin(3) + (k-1) * grid%step_sizes(3)

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
    end subroutine write_cube_file

    subroutine write_cube_file_pair(filename, property_name, pair_property_function)
        character(*), intent(in) :: filename, property_name
        interface
            function pair_property_function(x1, y1, z1, x2, y2, z2) result(value)
                double precision, intent(in) :: x1, y1, z1, x2, y2, z2
                double precision :: value
            end function pair_property_function
        end interface

        integer :: i1, j1, k1, i2, j2, k2, icount, iatom
        double precision :: x1, y1, z1, x2, y2, z2
        double precision,  dimension(3) :: origin, origin_2
        integer,  dimension(3) :: dims, dims_2
        integer :: unit_number




        ! Convert grid info to dimensions and origin
        do i1 = 1, 3
            origin(i1) = -grid%max_vals(i1)
            origin_2(i1) = -grid_2%max_vals(i1)
            if (grid%step_sizes(i1) .le. 0.00001d0) then
                dims(i1) = 1
            else
                dims(i1) = 2*int(grid%max_vals(i1) / grid%step_sizes(i1))+1
            end if
            if (grid_2%step_sizes(i1) .le. 0.00001d0) then
                dims_2(i1) = 1
            else
                dims_2(i1) = 2*int(grid_2%max_vals(i1) / grid_2%step_sizes(i1))+1
            end if
        end do

        open(newunit=unit_number, file=filename)

        ! Write header
        write(unit_number, '(A)') trim(property_name)
        write(unit_number, '(A)') "OUTER LOOP: X1, MIDDLE LOOP: Y1, INNER LOOP: Z1, X2, Y2, Z2"
        write(unit_number, '(I5,3F12.6)') natoms, origin(1), origin(2), origin(3)

        ! Write grid specifications for the first 3 dimensions (x1,y1,z1)
        do i1 = 1, 3
            if (i1 == 1) then
                write(unit_number, '(I5,3F12.6)') dims(i1), grid%step_sizes(i1), 0.0d0, 0.0d0
            else if (i1 == 2) then
                write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, grid%step_sizes(i1), 0.0d0
            else if (i1 == 3) then
                write(unit_number, '(I5,3F12.6)') dims(i1), 0.0d0, 0.0d0, grid%step_sizes(i1)
            end if
        end do

        ! Write atom information
        do iatom = 1, natoms
            write(unit_number, '(I5,4F12.6)') an(iatom), chrg(iatom), &
                                             cartes(iatom,1), cartes(iatom,2), cartes(iatom,3)
        end do



        ! Write volumetric data
        icount = 0
        do i1 = 1, dims(1)
            x1 = origin(1) + (i1-1) * grid%step_sizes(1)
            do j1 = 1, dims(2)
                y1 = origin(2) + (j1-1) * grid%step_sizes(2)
                do k1 = 1, dims(3)
                    z1 = origin(3) + (k1-1) * grid%step_sizes(3)
                    do i2 = 1, dims_2(1)
                        x2 = origin_2(1) + (i2-1) * grid_2%step_sizes(1)
                        do j2 = 1, dims_2(2)
                            y2 = origin_2(2) + (j2-1) * grid_2%step_sizes(2)
                            do k2 = 1, dims_2(3)
                                z2 = origin_2(3) + (k2-1) * grid_2%step_sizes(3)

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
    end subroutine write_cube_file_pair

end module cube_module
