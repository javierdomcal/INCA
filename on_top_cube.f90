module on_top_cube
    use cube_utils
    use wfxinfo, only: uhf, opsh, noccmo ! Ensure these variables are defined in the module
    use geninfo
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

    subroutine ontop_cube_write(dm2_file)
        character(len=*), intent(in) :: dm2_file

        ! Generate cube files for each property
        call write_cube_file("ontop.cube", "On-top pair density", get_ontop_value, &
                                                 cube_origin, cube_dims, cube_spacing, cube_units, cube_format, include_atoms)
        call write_cube_file("density.cube", "Electron density", get_density_value, &
                                                 cube_origin, cube_dims, cube_spacing, cube_units, cube_format, include_atoms)
        call write_cube_file("ID.cube", "Indicator dynamic", get_indicator_value, &
                                                 cube_origin, cube_dims, cube_spacing, cube_units, cube_format, include_atoms)

        ! Cleanup cached data
        call clear_cached_dm2()

        write(*,*) "Cube files generated successfully."
    end subroutine ontop_cube_write

    function get_ontop_value(x, y, z) result(ontop_val)
        double precision, intent(in) :: x, y, z
        double precision :: ontop_val
        ontop_val = get_cached_ontop_value(x, y, z)
    end function get_ontop_value

    function get_density_value(x, y, z) result(density_val)
        double precision, intent(in) :: x, y, z
        double precision :: density_val
        if (associated(uhf) .or. associated(opsh)) then
            density_val = Dens_a(x, y, z) + Dens_b(x, y, z)
        else
            density_val = Density(x, y, z)
        endif
    end function get_density_value

    function get_indicator_value(x, y, z) result(indicator_val)
        double precision, intent(in) :: x, y, z
        double precision :: indicator_val
        indicator_val = calculate_indicator_at_point(x, y, z)
    end function get_indicator_value

    !----------------------------------------------------------------------
    ! Caching functions for performance optimization
    !----------------------------------------------------------------------

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

    !----------------------------------------------------------------------
    ! Property calculation functions
    !----------------------------------------------------------------------

    function get_cached_ontop_value(x, y, z) result(ontop_val)
        implicit none
        double precision, intent(in) :: x, y, z
        double precision :: ontop_val
        integer :: i
        double precision :: prim_val_i, prim_val_j, prim_val_k, prim_val_l

        ontop_val = 0.0d0

        do i = 1, cached_count
                prim_val_i = Prim(x, y, z, cached_indices(i,1))
                prim_val_j = Prim(x, y, z, cached_indices(i,2))
                prim_val_k = Prim(x, y, z, cached_indices(i,3))
                prim_val_l = Prim(x, y, z, cached_indices(i,4))

                ontop_val = ontop_val + cached_values(i) * &
                                        prim_val_i * prim_val_j * prim_val_k * prim_val_l
        end do

        if (abs(ontop_val) < 1.0D-11) ontop_val = 0.d0
    end function get_cached_ontop_value

    function calculate_indicator_at_point(x, y, z) result(indicator)
        implicit none
        double precision, intent(in) :: x, y, z
        double precision :: indicator
        integer :: orbital
        double precision :: orbital_value

        indicator = 0.0d0

        if (uhf .or. opsh) then
                do orbital = 1, size(occupation_factors)
                        orbital_value = merge(MO_a(x, y, z, orbital), MO_b(x, y, z, orbital), orbital_type(orbital) == 1)
                        indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
                end do
        else
                do orbital = 1, size(occupation_factors)
                        orbital_value = MoOr(x, y, z, orbital)
                        indicator = indicator + occupation_factors(orbital) * orbital_value**2.0d0
                end do
        endif

        if (abs(indicator) < 1.0D-11) indicator = 0.d0
    end function calculate_indicator_at_point

end module on_top_cube
