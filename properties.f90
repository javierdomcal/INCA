! properties.f90
! Unified properties management module for INCA program
! Simplifies the handling of properties by using a single registry

module properties
     implicit none

    ! Maximum number of properties that can be registered

    ! Property type enum
    integer, parameter :: PROP_SINGLE = 1
    integer, parameter :: PROP_PAIR = 2

    ! Property descriptor type

    ! Function pointer interfaces
    abstract interface
        function single_prop_func(r) result(value)
            import :: dp
            double precision, dimension(3), intent(in) :: r
            double precision :: value
        end function

        function pair_prop_func(r1, r2) result(value)
            import :: dp
            double precision, dimension(3), intent(in) :: r1, r2
            double precision :: value
        end function
    end interface

    ! Property registry
    abstract interface
    procedure(single_prop_func), pointer, nopass :: single_func => null()
    procedure(pair_prop_func), pointer, nopass :: pair_func => null()
    end interface



contains



    ! Calculate a single property value
    function calculate_single_property(single_func, r) result(value)
        double precision, dimension(3), intent(in) :: r
        procedure(single_prop_func), pointer, intent(in) :: single_func
        double precision :: value

        value = 0.0_dp
        if (associated(single_func)) then
            value = single_func(r)
        else
            write(*,*) "Error: Function pointer not associated for property:"
        end if
    end function calculate_single_property

    ! Calculate a pair property value
    function calculate_pair_property(pair_func, r1, r2) result(value)
         double precision, dimension(3), intent(in) :: r1, r2
        procedure(pair_prop_func), pointer, intent(in) :: pair_func
        double precision :: value

        value = 0.0_dp
        if (associated(pair_func)) then
            value = pair_func(r1, r2)
        else
            write(*,*) "Error: Function pointer not associated for property:"
        end if
    end function calculate_pair_property



end module properties