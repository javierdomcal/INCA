module properties
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    ! Maximum number of properties that can be registered (renamed to avoid conflict)
    integer, parameter :: PROP_MAX_SINGLE = 20
    integer, parameter :: PROP_MAX_PAIR = 10

    ! Arrays to store property data (renamed to avoid conflict)
    character(len=80), allocatable :: prop_args(:)
    character(len=40), allocatable :: prop_names(:)
    integer :: prop_count = 0

    ! Property descriptor type
    type :: property_descriptor
        character(len=40) :: name
        character(len=80) :: description
        integer :: dimensionality  ! 1=single, 2=pair
        logical :: enabled = .false.
    end type

    ! Single property type
    type :: single_property
        type(property_descriptor) :: descriptor
        procedure(single_property_func), pointer, nopass :: calculate => null()
    end type

    ! Pair property type
    type :: pair_property
        type(property_descriptor) :: descriptor
        procedure(pair_property_func), pointer, nopass :: calculate => null()
    end type

    ! Interfaces for property functions
    interface
        function single_property_func(r) result(value)
            import :: dp
            real(dp), dimension(3), intent(in) :: r
            real(dp) :: value
        end function

        function pair_property_func(r1, r2) result(value)
            import :: dp
            real(dp), dimension(3), intent(in) :: r1, r2
            real(dp) :: value
        end function
    end interface

    ! Registry for properties
    type(single_property), allocatable :: single_properties(:)
    type(pair_property), allocatable :: pair_properties(:)

    ! Tracking current number of registered properties
    integer :: num_single_properties = 0
    integer :: num_pair_properties = 0

contains

    ! Set property information from main program
    subroutine set_property_info(names, args, count)
        character(len=40), intent(in) :: names(:)
        character(len=80), intent(in) :: args(:)
        integer, intent(in) :: count

        ! Store the count
        prop_count = count

        ! Allocate and copy to module-level variables
        if (allocated(prop_names)) deallocate(prop_names)
        if (allocated(prop_args)) deallocate(prop_args)

        allocate(prop_names(count))
        allocate(prop_args(count))

        prop_names(1:count) = names(1:count)
        prop_args(1:count) = args(1:count)
    end subroutine set_property_info

    ! Initialize properties
    subroutine initialize_properties(max_single)
        integer, intent(in), optional :: max_single
        integer :: max_to_use

        ! Use provided value or default
        if (present(max_single)) then
            max_to_use = max_single
        else
            max_to_use = PROP_MAX_SINGLE
        endif

        ! Allocate property arrays based on input parameters
        if (allocated(single_properties)) deallocate(single_properties)
        if (allocated(pair_properties)) deallocate(pair_properties)

        allocate(single_properties(max_to_use))
        allocate(pair_properties(PROP_MAX_PAIR))

        ! Reset property counters
        num_single_properties = 0
        num_pair_properties = 0
    end subroutine

    ! Register a single property
    subroutine register_single_property(name, description, calculate)
        character(len=*), intent(in) :: name, description
        procedure(single_property_func), pointer, intent(in) :: calculate

        ! Increment counter and check bounds
        num_single_properties = num_single_properties + 1
        if (num_single_properties > size(single_properties)) then
            write(*,*) "Error: Too many single properties registered"
            return
        end if

        single_properties(num_single_properties)%descriptor%name = name
        single_properties(num_single_properties)%descriptor%description = description
        single_properties(num_single_properties)%descriptor%dimensionality = 1
        single_properties(num_single_properties)%calculate => calculate
    end subroutine

    ! Register a pair property
    subroutine register_pair_property(name, description, calculate)
        character(len=*), intent(in) :: name, description
        procedure(pair_property_func), pointer, intent(in) :: calculate

        ! Increment counter and check bounds
        num_pair_properties = num_pair_properties + 1
        if (num_pair_properties > size(pair_properties)) then
            write(*,*) "Error: Too many pair properties registered"
            return
        end if

        pair_properties(num_pair_properties)%descriptor%name = name
        pair_properties(num_pair_properties)%descriptor%description = description
        pair_properties(num_pair_properties)%descriptor%dimensionality = 2
        pair_properties(num_pair_properties)%calculate => calculate
    end subroutine

    ! Enable a property by name
    subroutine enable_property(name)
        character(len=*), intent(in) :: name
        integer :: i

        ! First check single properties
        do i = 1, num_single_properties
            if (trim(single_properties(i)%descriptor%name) == trim(name)) then
                single_properties(i)%descriptor%enabled = .true.
                return
            end if
        end do

        ! Then check pair properties
        do i = 1, num_pair_properties
            if (trim(pair_properties(i)%descriptor%name) == trim(name)) then
                pair_properties(i)%descriptor%enabled = .true.
                return
            end if
        end do

        write(*,*) "Warning: Property not found for enabling: ", trim(name)
    end subroutine

    ! Get property argument by name
    function get_property_argument(name) result(arg)
        character(len=*), intent(in) :: name
        character(len=80) :: arg
        integer :: i

        arg = ""

        ! Search for the property in the prop_names array
        if (allocated(prop_names)) then
            do i = 1, prop_count
                if (trim(prop_names(i)) == trim(name)) then
                    if (allocated(prop_args) .and. i <= size(prop_args)) then
                        arg = prop_args(i)
                    end if
                    return
                end if
            end do
        end if
    end function

    ! Get single property by name
    function get_single_property_by_name(name) result(prop_index)
        character(len=*), intent(in) :: name
        integer :: prop_index
        integer :: i

        prop_index = 0
        do i = 1, num_single_properties
            if (trim(single_properties(i)%descriptor%name) == trim(name)) then
                prop_index = i
                return
            end if
        end do
    end function

    ! Get pair property by name
    function get_pair_property_by_name(name) result(prop_index)
        character(len=*), intent(in) :: name
        integer :: prop_index
        integer :: i

        prop_index = 0
        do i = 1, num_pair_properties
            if (trim(pair_properties(i)%descriptor%name) == trim(name)) then
                prop_index = i
                return
            end if
        end do
    end function

    ! Check if a property is enabled
    function is_property_enabled(name) result(enabled)
        character(len=*), intent(in) :: name
        logical :: enabled
        integer :: single_index, pair_index

        single_index = get_single_property_by_name(name)
        if (single_index > 0) then
            enabled = single_properties(single_index)%descriptor%enabled
            return
        end if

        pair_index = get_pair_property_by_name(name)
        if (pair_index > 0) then
            enabled = pair_properties(pair_index)%descriptor%enabled
            return
        end if

        enabled = .false.
    end function

    ! Calculate a single property
    function calculate_single_property(name, coordinates) result(value)
        character(len=*), intent(in) :: name
        real(dp), dimension(3), intent(in) :: coordinates
        real(dp) :: value
        integer :: prop_index

        prop_index = get_single_property_by_name(name)
        if (prop_index > 0) then
            value = single_properties(prop_index)%calculate(coordinates)
        else
            value = 0.0_dp
            write(*,*) "Warning: Single property not found: ", trim(name)
        end if
    end function

    ! Calculate a pair property
    function calculate_pair_property(name, coordinates1, coordinates2) result(value)
        character(len=*), intent(in) :: name
        real(dp), dimension(3), intent(in) :: coordinates1, coordinates2
        real(dp) :: value
        integer :: prop_index

        prop_index = get_pair_property_by_name(name)
        if (prop_index > 0) then
            value = pair_properties(prop_index)%calculate(coordinates1, coordinates2)
        else
            value = 0.0_dp
            write(*,*) "Warning: Pair property not found: ", trim(name)
        end if
    end function

    ! Get list of all enabled single properties
    subroutine get_enabled_single_properties(names, count)
        character(len=40), allocatable, intent(out) :: names(:)
        integer, intent(out) :: count
        integer :: i, idx

        ! Count enabled properties
        count = 0
        do i = 1, num_single_properties
            if (single_properties(i)%descriptor%enabled) then
                count = count + 1
            end if
        end do

        ! Allocate and fill names array
        if (allocated(names)) deallocate(names)
        allocate(names(count))

        idx = 0
        do i = 1, num_single_properties
            if (single_properties(i)%descriptor%enabled) then
                idx = idx + 1
                names(idx) = single_properties(i)%descriptor%name
            end if
        end do
    end subroutine

    ! Get list of all enabled pair properties
    subroutine get_enabled_pair_properties(names, count)
        character(len=40), allocatable, intent(out) :: names(:)
        integer, intent(out) :: count
        integer :: i, idx

        ! Count enabled properties
        count = 0
        do i = 1, num_pair_properties
            if (pair_properties(i)%descriptor%enabled) then
                count = count + 1
            end if
        end do

        ! Allocate and fill names array
        if (allocated(names)) deallocate(names)
        allocate(names(count))

        idx = 0
        do i = 1, num_pair_properties
            if (pair_properties(i)%descriptor%enabled) then
                idx = idx + 1
                names(idx) = pair_properties(i)%descriptor%name
            end if
        end do
    end subroutine

end module properties