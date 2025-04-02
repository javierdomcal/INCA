module density_properties
   use wfxinfo
   use geninfo
   use wfxinfo_hf
   use cube_module, only: write_cube_file
   implicit none

   ! External function declarations
   double precision, external :: Density, Dens_a, Dens_b, Dens_hf

contains
   ! Density calculation functions compatible with single_prop_func interface
   function density_calc(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value
      value = Density(x, y, z)
   end function density_calc

   function alpha_density(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value
      value = Dens_a(x, y, z)
   end function alpha_density

   function beta_density(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value
      value = Dens_b(x, y, z)
   end function beta_density

   function spin_density(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value
      value = Dens_a(x, y, z) - Dens_b(x, y, z)
   end function spin_density

   function hf_density(x, y, z) result(value)
      double precision, intent(in) :: x, y, z
      double precision :: value
      value = Dens_hf(x, y, z)
   end function hf_density

   ! Entry point for density calculations
   subroutine density_calculation(output_type)
      character(len=10), intent(in) :: output_type

      select case(output_type)
         case('total')
            write(*,*) "Generating electron density cube file..."
            call write_cube_file("density.cube", "Electron Density", &
                                 density_calc)

         case('alpha')
            write(*,*) "Generating alpha density cube file..."
            call write_cube_file("alpha_density.cube", "Alpha Electron Density", &
                                 alpha_density)

         case('beta')
            write(*,*) "Generating beta density cube file..."
            call write_cube_file("beta_density.cube", "Beta Electron Density", &
                                 beta_density)

         case('spin')
            write(*,*) "Generating spin density cube file..."
            call write_cube_file("spin_density.cube", "Spin Density", &
                                 spin_density)

         case('hf')
            write(*,*) "Generating Hartree-Fock density cube file..."
            call write_cube_file("hf_density.cube", "Hartree-Fock Density", &
                                 hf_density)

         case('all')
            write(*,*) "Generating all density cube files..."

            call write_cube_file("density.cube", "Electron Density", &
                                 density_calc)

            call write_cube_file("alpha_density.cube", "Alpha Electron Density", &
                                 alpha_density)

            call write_cube_file("beta_density.cube", "Beta Electron Density", &
                                 beta_density)

            call write_cube_file("spin_density.cube", "Spin Density", &
                                 spin_density)

            call write_cube_file("hf_density.cube", "Hartree-Fock Density", &
                                 hf_density)

         case default
            write(*,*) "Invalid density output type: ", output_type
      end select
   end subroutine density_calculation

end module density_properties