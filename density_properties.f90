! density_properties.f90
! Module for electron density and related properties
! This module handles registration and calculation of various density properties

module density_properties
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use wfxinfo
   use geninfo
   use wfxinfo_hf
   use properties, only: register_single_property, single_property_func
   implicit none

   ! External function declarations
   double precision, external :: Density, Dens_a, Dens_b, Dens_hf, MoOr, MO_a, MO_b

   ! Module state
   logical, private :: initialized = .false.

contains

   ! Initialize the module
   subroutine initialize_density_properties()
      procedure(single_property_func), pointer :: total_density_func => null()
      procedure(single_property_func), pointer :: alpha_density_func => null()
      procedure(single_property_func), pointer :: beta_density_func => null()
      procedure(single_property_func), pointer :: spin_density_func => null()
      procedure(single_property_func), pointer :: hf_density_func => null()

      if (initialized) return

      ! Register electron density properties
      total_density_func => total_electron_density
      alpha_density_func => alpha_electron_density
      beta_density_func => beta_electron_density
      spin_density_func => spin_density
      hf_density_func => hartree_fock_density

      call register_single_property("Electron Density", &
                                  "Total electron density ρ(r)", total_density_func)

      call register_single_property("Alpha Density", &
                                  "Alpha electron density ρα(r)", alpha_density_func)

      call register_single_property("Beta Density", &
                                  "Beta electron density ρβ(r)", beta_density_func)

      call register_single_property("Spin Density", &
                                  "Spin density ρα(r)-ρβ(r)", spin_density_func)

      call register_single_property("HF Density", &
                                  "Hartree-Fock electron density", hf_density_func)

      initialized = .true.
      write(*,*) "Density properties module initialized."
   end subroutine initialize_density_properties

   ! Implementation of total electron density for property system
   function total_electron_density(r) result(value)
      real(dp), dimension(3), intent(in) :: r
      real(dp) :: value

      value = Density(r(1), r(2), r(3))
   end function total_electron_density

   ! Implementation of alpha electron density for property system
   function alpha_electron_density(r) result(value)
      real(dp), dimension(3), intent(in) :: r
      real(dp) :: value

      value = Dens_a(r(1), r(2), r(3))
   end function alpha_electron_density

   ! Implementation of beta electron density for property system
   function beta_electron_density(r) result(value)
      real(dp), dimension(3), intent(in) :: r
      real(dp) :: value

      value = Dens_b(r(1), r(2), r(3))
   end function beta_electron_density

   ! Implementation of spin density for property system
   function spin_density(r) result(value)
      real(dp), dimension(3), intent(in) :: r
      real(dp) :: value

      value = Dens_a(r(1), r(2), r(3)) - Dens_b(r(1), r(2), r(3))
   end function spin_density

   ! Implementation of Hartree-Fock density for property system
   function hartree_fock_density(r) result(value)
      real(dp), dimension(3), intent(in) :: r
      real(dp) :: value

      value = Dens_hf(r(1), r(2), r(3))
   end function hartree_fock_density

   ! Adapter functions for cube file generation
   function density_adapter(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value

      value = Density(x, y, z)
   end function density_adapter

   function alpha_density_adapter(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value

      value = Dens_a(x, y, z)
   end function alpha_density_adapter

   function beta_density_adapter(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value

      value = Dens_b(x, y, z)
   end function beta_density_adapter

   function spin_density_adapter(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value

      value = Dens_a(x, y, z) - Dens_b(x, y, z)
   end function spin_density_adapter

   function hf_density_adapter(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value

      value = Dens_hf(x, y, z)
   end function hf_density_adapter

   ! Entry point for density calculations
   subroutine density_calculation(output_type, grid_params)
      use cube_module        ! For cube file output
      use properties  ! For property access

      character(len=10), intent(in) :: output_type  ! "total", "alpha", "beta", "spin", "hf", "all"
      type(grid_parameters), intent(in) :: grid_params

      ! Initialize density module
      call initialize_density_properties()

      ! Enable properties based on output_type
      if (output_type == "total" .or. output_type == "all") then
         call enable_property("Electron Density")
      end if

      if (output_type == "alpha" .or. output_type == "all") then
         call enable_property("Alpha Density")
      end if

      if (output_type == "beta" .or. output_type == "all") then
         call enable_property("Beta Density")
      end if

      if (output_type == "spin" .or. output_type == "all") then
         call enable_property("Spin Density")
      end if

      if (output_type == "hf" .or. output_type == "all") then
         call enable_property("HF Density")
      end if

      ! Generate cube files based on enabled properties
      if (is_property_enabled("Electron Density")) then
         write(*,*) "Generating electron density cube file..."
         call write_cube_file("density.cube", "Electron Density", &
                             density_adapter, grid_params)
      end if

      if (is_property_enabled("Alpha Density")) then
         write(*,*) "Generating alpha density cube file..."
         call write_cube_file("alpha_density.cube", "Alpha Electron Density", &
                             alpha_density_adapter, grid_params)
      end if

      if (is_property_enabled("Beta Density")) then
         write(*,*) "Generating beta density cube file..."
         call write_cube_file("beta_density.cube", "Beta Electron Density", &
                             beta_density_adapter, grid_params)
      end if

      if (is_property_enabled("Spin Density")) then
         write(*,*) "Generating spin density cube file..."
         call write_cube_file("spin_density.cube", "Spin Density", &
                             spin_density_adapter, grid_params)
      end if

      if (is_property_enabled("HF Density")) then
         write(*,*) "Generating Hartree-Fock density cube file..."
         call write_cube_file("hf_density.cube", "Hartree-Fock Density", &
                             hf_density_adapter, grid_params)
      end if
   end subroutine density_calculation

end module density_properties
