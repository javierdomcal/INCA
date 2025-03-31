module correlation_indicators
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use wfxinfo
   use geninfo
   use properties, only: register_single_property, single_property_func
   implicit none

   ! External function declarations
   double precision, external :: MO_a, MO_b, MoOr

   ! Module variables shared with other modules
   double precision, allocatable, public :: occupation_factors(:)
   integer, allocatable, public :: orbital_type(:)

contains

   ! Initialize the indicator module - call at program start
   subroutine initialize_indicators()
      procedure(single_property_func), pointer :: indicator_func => null()

      ! Register indicator dynamic as a property
      indicator_func => indicator_dynamic
      call register_single_property("Indicator Dynamic", &
                                   "Measure of electron correlation effects", indicator_func)
   end subroutine initialize_indicators

   ! Implementation of indicator dynamic property
   function indicator_dynamic(r) result(value)
      real(dp), dimension(3), intent(in) :: r
      real(dp) :: value

      value = calculate_indicator_at_point(r(1), r(2), r(3))
   end function indicator_dynamic

   ! Calculate indicator dynamic value at a given point
   function calculate_indicator_at_point(x, y, z) result(indicator)
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

   ! Function for use in API that takes separate coordinates
   function indicator_at_point(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value

      value = calculate_indicator_at_point(x, y, z)
   end function indicator_at_point

end module correlation_indicators