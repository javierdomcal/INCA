module cube_utils
   implicit none
   logical :: generate_cube
   double precision, dimension(3) :: cube_origin
   double precision, dimension(3) :: cube_spacing
   integer, dimension(3) :: cube_points
   character*40 :: cube_filename
end module cube_utils

module on_top_cube
   implicit none
   logical :: compute_on_top
   double precision, allocatable, dimension(:,:) :: on_top_density
   character*40 :: on_top_cube_filename
end module on_top_cube

module pair_density
   implicit none
   logical :: compute_pair_density
   double precision, allocatable, dimension(:,:) :: pair_density_matrix
   character*40 :: pair_density_filename
end module pair_density

module properties_module
   type :: property_info
    character(len=40) :: name
    character(len=80) :: argument
    logical :: enabled
end type property_info
type(property_info), allocatable :: properties(:)
end module properties_module

module inputdat
implicit none
 logical :: readwfx, readlog, cube_flag, primcube, aocube, MOcube, denscube, gradient, laplacian, intracalc, id, ontop_hf
 logical :: dmn
 character*80 :: option
 !!!!!!!!!!!input filenames!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 character*80 :: nameinp
 character*40 :: wfxfilename
 character*40 :: logfilename
 character*40 :: wfxhfname !to compute coulomb hole
 !!!!!!!!!!!output_cubefile_names!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 character*40 :: nameprim
 character*40 :: nameao
 character*40 :: namemo   !name of the cubefile we will generate
 character*40 :: namedens
 character*40 :: namegrad
 character*40 :: namelap
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 character*40 :: nameint
 character*40 :: namec1
 character*40 :: dm2name

 character*40 :: dm2_file
 character*40 :: dm2hf_file
 character*40 :: dm2hfl_file

 ! Property system variables
 integer :: num_properties = 0
 integer :: max_properties = 20
 character(len=40), allocatable :: requested_properties(:)
 character(len=80), allocatable :: property_arguments(:)

 ! Scanning variables
 integer, allocatable :: atom_indices(:)
 character(len=1), allocatable :: scan_directions(:)
 integer :: n_atoms_scan = 0
 integer :: n_directions = 0
 real*8 :: end_distance = 0.0d0
 real*8 :: step_size = 0.0d0
 logical :: ontop_cube_output = .false.
 logical :: Grid_2_flag = .false.

 ! Density properties flags
 logical, allocatable :: density_properties_enabled(:)  ! Array for enabled density properties
 logical, allocatable :: pair_density_properties_enabled(:)  ! Array for enabled pair density properties
 logical :: indicator_dynamic_enabled = .false.  ! Flag for indicator dynamic

 ! Output type flags
 character(len=10) :: density_output_type = "total"  ! "total", "alpha", "beta", "spin", "hf", "all"
 character(len=10) :: pair_density_output_type = "full"  ! "full", "c1", "c2", "all"

 ! Grid definition
 type :: grid_type
     double precision, dimension(3) :: max_vals = 0.0d0
     double precision, dimension(3) :: step_sizes = 0.0d0
 end type grid_type

 ! Grid instances
 type(grid_type) :: grid
 type(grid_type) :: grid_2

 ! Maximum property settings
 integer :: max_pair_properties = 10

end module inputdat

!----------------------------------------------------------------------------------------------------

module cubeinfo
   integer :: mo !molecular orbital to represent in the cubefile
   integer :: npr  !primitive to represent in the cubefile
   integer :: cao  !atomic orbital to represent in the cubefile
   double precision, dimension(3) :: center !cube centered in (x,y,z)
   double precision, dimension(3) :: step !distance between points for each axis
   integer, dimension(3) :: np !number of points in the cube for each axis
end module cubeinfo

!-----------------------------------------------------------------------------------------------------

module quadinfo1
implicit none  ! Angular and radial integral
logical :: autoc                                       !calculate centres automatically (Becke)
logical :: radial_integral                    !=.true. integral of radial intra from 0 to infty
character*40 :: int_mom_name
integer :: nquad !number of quadratures
double precision, allocatable, dimension(:,:) :: cent  !integration centre
integer, allocatable, dimension(:) :: nelec_cent !number of electron per centre
double precision, allocatable, dimension(:) :: sfalpha !scaling factor for Gauss-Legendre

logical :: dif_nodes                      !different nodes per centre
integer :: gpt,npr !number of angular nodes, number of radial nodes
integer, allocatable, dimension(:) :: nradc
integer, allocatable, dimension(:) :: nangc  !angular grid points (diff nodes per centre)
double precision :: a,b                                !limits of the integral
double precision :: snglalpha                     !alpha parameter (same for all nodes)
end module quadinfo1

module quadinfo2 ! Angular integral at several distances
logical :: radial_plot                                !=.true if we want I vs R plot
character*40 :: r_plot_name                           !name of the radial plot
integer :: nblock                                     !number of blocks for the radius
integer, allocatable, dimension(:) :: n_an_per_part   !number of angluar points per block
double precision, allocatable, dimension(:,:) :: tart !space between radius
integer, allocatable, dimension(:) :: stp    !number of points for each block
end module quadinfo2

module intracube !vectorial plot with a cubefile!
logical :: cubeint
character*40 :: cubeintraname
double precision, dimension(3) :: center_i            !cube centered in (x,y,z)
double precision, dimension(3) :: step_i              !distance between points for each axis
integer, dimension(3) :: np3       !number of points in the cube for each axis
integer :: np3t
double precision, allocatable, dimension(:,:) :: r3
end module intracube

!-----------------------------------------------------------------------------------------------------

module geninfo           !general information, can be obtained from .wfx files.
  implicit none
  logical :: opsh  !open or closed shell
  integer :: natoms      !number of atoms
  integer :: nelec
  integer :: nprim       !number of primitives
  double precision :: netch !net charge
  integer :: mspin       !electronic spin multiplicity
  integer :: nalfae      !number of alpha electrons
  integer :: nbetae      !number of beta electrons
  integer, allocatable, dimension(:) :: Ra    !atomic centers of primitives
  integer, allocatable, dimension(:) :: Ptyp  !primitive type
  double precision, allocatable, dimension(:) :: Alpha !Primitive exponents
  integer, allocatable, dimension(:,:) :: TMN !matrix with t,m,n coeficients
  double precision, allocatable,  dimension(:,:) :: cartes !nuclear cartesian coordinates
  integer, allocatable, dimension(:) :: an !atomic number
  double precision, allocatable, dimension(:) :: chrg !nuclear charge
  double precision, parameter :: pi=4.d0*datan(1.d0)
end module

!--------------------------------------------------------------------------------------------------------

module wfxinfo !especific data of wfx files
  implicit none
  logical :: corr   !correlated method
  logical :: uhf   !unrestricted-restricted Hartree Fock???
  logical :: udens !relaxed or unrelaxed density
  integer :: noccmo !number of occupied MOs
  double precision, allocatable, dimension(:,:) :: T   !MOs in primitives coeficients (matriu T)
  double precision, allocatable, dimension(:) :: Occ   !ocupancies of orbitals
  !for open-shell
  double precision, allocatable, dimension(:,:) :: T_a !Alpha MOs coeficients
  double precision, allocatable, dimension(:,:) :: T_b !Beta MOs coeficients
  double precision, allocatable, dimension(:) :: Occ_a, Occ_b
  integer :: nalfaorb, nbetaorb
end module wfxinfo

!-------------------------------------------------------------------------------------------------------

module wfxinfo_hf
  integer :: noccmo2 !number of occupied MOs
  double precision, allocatable, dimension(:,:) :: T2   !MOs in primitives coeficients (matriu T)
  double precision, allocatable, dimension(:) :: Occ2   !ocupancies of orbitals
  !for open-shell
  double precision, allocatable, dimension(:,:) :: T_a2 !Alpha MOs coeficients
  double precision, allocatable, dimension(:,:) :: T_b2 !Beta MOs coeficients
  integer :: nalfaorb2, nbetaorb2
end module wfxinfo_hf

!---------------------------------------------------------------------------------------

module holes
   implicit none
   double precision,allocatable,dimension(:) :: hfintra,sdintra,corrintra
end module holes

!----------------------------------------------------------------------------------------

 module loginfo
   implicit none
   double precision, allocatable, dimension(:) :: Flg, N_prim !fixed expansion coeficient for the contraction of AOs into primitives, normalization constant for the primitives
   double precision, allocatable, dimension(:,:) :: Ckalk !expansion coeficient of MOs in AOs
   integer, allocatable, dimension(:) :: npao !number of primitives for each atomic orbital
   integer :: nao !number of atomic orbitals in the molecule
   integer, allocatable, dimension(:) :: aotyp !type of atomic orbital: s->0, px, py, pz ->1,  ...
end module loginfo

!------------------------------------------------------------------------------------------

module thresholds
implicit none
double precision :: thresh !threshold for tau
end module

!--------------------------------------------------------------------------------------------------------

module total_integral
implicit none
integer :: np1 !number of quadrature points
double precision, allocatable, dimension(:,:) :: r1 !quadrature points
double precision, allocatable, dimension(:) :: w1, w1_1, w1_2,w1_11,w1_22
double precision, allocatable, dimension(:) :: I_vec1
end module

module radial_scan
implicit none
integer, allocatable, dimension(:) :: smrd
integer :: np2 !number of grid points
double precision, allocatable, dimension(:,:) :: r2 !grid points
double precision, allocatable, dimension(:) :: radi !distances of the scan
integer :: nradi
double precision, allocatable, dimension(:) :: w2 !weight
double precision, allocatable, dimension(:) :: I_vec2
end module

!-------------------------------------------------------------------------------------------------------------

module qhermite
double precision, allocatable, dimension(:) :: rh, w_r !nodes and weights for gauss hermite(coef.)
end module

!-------------------------------------------------------------------------------------------------

module fractions
   double precision, parameter :: hr=1.d0*(3.d0**(-1.d0))
   double precision, parameter :: boh=5.d0*(3.d0**(-1.d0))
   double precision, parameter :: bs=1.d0*(6.d0**(-1.d0))
   double precision, parameter :: bih=2.d0*(3.d0**(-1.d0))
end module fractions

!-------------------------------------------------------------------------------------------------

module radis
   double precision, parameter, dimension(92) :: BL=(/0.327d0, 0.320d0, &
    + 1.219d0, 0.911d0, 0.793d0, 0.766d0, 0.699d0, 0.658d0,&
    + 0.900d0, 0.690d0, 1.545d0, 1.333d0, 1.199d0, 1.123d0, 1.110d0, 1.071d0, 1.039d0,&
    + 0.970d0, 1.978d0, 1.745d0, 1.337d0, 1.274d0, 1.236d0, 1.128d0, 1.180d0, 1.091d0,&
    + 1.089d0, 1.077d0, 1.146d0, 1.187d0, 1.199d0, 1.179d0, 1.209d0, 1.201d0, 1.201d0,&
    + 1.100d0, 2.217d0, 1.928d0, 1.482d0, 1.377d0, 1.353d0, 1.240d0, 1.287d0, 1.212d0,&
    + 1.229d0, 1.240d0, 1.362d0, 1.429d0, 1.385d0, 1.380d0, 1.421d0, 1.400d0, 1.397d0,&
    + 1.300d0, 2.442d0, 2.149d0, 1.653d0, 1.600d0, 1.600d0, 1.600d0, 1.600d0, 1.600d0,&
    + 1.600d0, 1.600d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0, 1.500d0,&
    + 1.364d0, 1.346d0, 1.256d0, 1.258d0, 1.222d0, 1.227d0, 1.227d0, 1.273d0, 1.465d0,&
    + 1.531d0, 1.434d0, 1.496d0, 1.500d0, 1.500d0, 1.450d0, 1.500d0, 1.500d0, 1.500d0,&
    + 1.500d0, 1.500d0,1.500d0/)
end module radis

!----------------------------------------------------------------------------------------

module located
        contains
                subroutine locate(iunit,string)
                integer :: ii
                integer iunit
                character string*(*)
                character*80 linia
                rewind(iunit)
                ii=0
                do while(ii.eq.0)
                  read(iunit,"(a80)", end=10)linia
                  if(index(linia,string).ne.0) then
                  ii=1
                  return
                 end if
                end do
                10 write(*,*) string,"", 'section not found'
                !return
                end subroutine
end module located

!---------------------------------------------------------------------------------------

module coeff
        double precision, allocatable, dimension(:) :: C_x,C_y,C_z !coefficients of V
        double precision, dimension(3) :: Lrtot
end module coeff

!---------------------------------------------------------------------------------------




