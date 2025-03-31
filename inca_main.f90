!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program wavefunction !main program
!!!!!reads .wfx and .log files and performs calculations about this info
use inputdat        !information about the calculations we want to do
use properties      !property management system (if enabled in compile)
use density_properties     !For density property calculations
use pair_density_module    !For pair density calculations
use correlation_indicators !For indicator dynamic calculations
use cube_module            !For cube file output
implicit none
integer :: a !defines to subroutine cubefile what function we want to represent: prim, ao, mo, dens
double precision, allocatable, dimension(:,:) :: r1,r2,r3 !nodes
character*80 :: dm2p,dm2phf,dm2phfl, dm2pc1 !name of the DM2 files

!provisional stuff
character*80 :: dm1name, namerad

readwfx=.false.
readlog=.false.

call readinput()

write(*,*) "=== Starting Main Program ==="
write(*,*) "Option selected: ", option
write(*,*) "ontop_hf flag: ", ontop_hf
write(*,*) "dm2name: ", trim(dm2name)

! Read wave function files if specified
if (readwfx) call filewfx(wfxfilename)  !reads info from a wfx file
if (readlog) call filelog(logfilename)  !reads info from a log file

! Initialize properties system if needed
if (option == "property" .or. any(density_properties_enabled) .or. &
    any(pair_density_properties_enabled) .or. indicator_dynamic_enabled) then
    write(*,*) "Initializing property system..."
    call initialize_properties(max_properties, max_pair_properties)
endif

! If needed, read intracule information from input file for legacy or C1/C2 pair density
if ((option.ne.'no').and.(option.ne.'ontop').and.(option.ne.'property')) then
    call readintra(nameinp, dm2p, dm2phf, dm2phfl)
endif

! Handle property-based calculations
if (option == "property" .or. any(density_properties_enabled) .or. &
    any(pair_density_properties_enabled) .or. indicator_dynamic_enabled) then

    write(*,*) "=== Running property-based calculations ==="

    ! Process density properties
    if (any(density_properties_enabled)) then
        write(*,*) "Calculating density properties..."

        ! Initialize density properties module
        call initialize_density_properties()

        ! Enable selected density properties
        if (density_properties_enabled(1)) then
            call enable_property("Electron Density")
            write(*,*) "  Electron density enabled"
        endif

        if (density_properties_enabled(2)) then
            call enable_property("Alpha Density")
            write(*,*) "  Alpha density enabled"
        endif

        if (density_properties_enabled(3)) then
            call enable_property("Beta Density")
            write(*,*) "  Beta density enabled"
        endif

        if (density_properties_enabled(4)) then
            call enable_property("Spin Density")
            write(*,*) "  Spin density enabled"
        endif

        if (density_properties_enabled(5)) then
            call enable_property("HF Density")
            write(*,*) "  HF density enabled"
        endif

        ! Generate cube files for enabled density properties
        if (is_property_enabled("Electron Density")) then
            write(*,*) "Generating electron density cube file..."
            call write_cube_file("density.cube", "Electron Density", &
                               density_adapter, grid)
        endif

        if (is_property_enabled("Alpha Density")) then
            write(*,*) "Generating alpha density cube file..."
            call write_cube_file("alpha_density.cube", "Alpha Electron Density", &
                               alpha_density_adapter, grid)
        endif

        if (is_property_enabled("Beta Density")) then
            write(*,*) "Generating beta density cube file..."
            call write_cube_file("beta_density.cube", "Beta Electron Density", &
                               beta_density_adapter, grid)
        endif

        if (is_property_enabled("Spin Density")) then
            write(*,*) "Generating spin density cube file..."
            call write_cube_file("spin_density.cube", "Spin Density", &
                               spin_density_adapter, grid)
        endif

        if (is_property_enabled("HF Density")) then
            write(*,*) "Generating Hartree-Fock density cube file..."
            call write_cube_file("hf_density.cube", "Hartree-Fock Density", &
                               hf_density_adapter, grid)
        endif
    endif

    ! Process pair density properties
    if (any(pair_density_properties_enabled)) then
        write(*,*) "Calculating pair density properties..."

        ! Get required DM2 file names if not already set
        if (dm2p == '') then
            call readintra(nameinp, dm2p, dm2phf, dm2phfl)
        endif

        ! Initialize pair density module
        call initialize_pair_density(dm2p, dm2phf, dm2phfl)

        ! Enable selected pair density properties
        if (pair_density_properties_enabled(1)) then
            call enable_property("Full Pair Density")
            write(*,*) "  Full pair density enabled"
        endif

        if (pair_density_properties_enabled(2)) then
            call enable_property("C1 Pair Density")
            write(*,*) "  C1 pair density enabled"
        endif

        if (pair_density_properties_enabled(3)) then
            call enable_property("C2 Pair Density")
            write(*,*) "  C2 pair density enabled"
        endif

        ! Generate cube files for enabled pair density properties
        if (is_property_enabled("Full Pair Density")) then
            write(*,*) "Generating full pair density cube file..."
            call write_cube2_file("full_pair_density.cube", "Full Pair Density", &
                                full_pd_adapter, grid, grid_2)
        endif

        if (is_property_enabled("C1 Pair Density")) then
            write(*,*) "Generating C1 pair density cube file..."
            call write_cube2_file("c1_pair_density.cube", "C1 Pair Density", &
                                c1_pd_adapter, grid, grid_2)
        endif

        if (is_property_enabled("C2 Pair Density")) then
            write(*,*) "Generating C2 pair density cube file..."
            call write_cube2_file("c2_pair_density.cube", "C2 Pair Density", &
                                c2_pd_adapter, grid, grid_2)
        endif

        ! Cleanup pair density module
        call cleanup_pair_density()
    endif

    ! Process indicator dynamic property
    if (indicator_dynamic_enabled) then
        write(*,*) "Calculating indicator dynamic..."

        ! Initialize indicator dynamic module
        call initialize_indicators()

        ! Enable indicator dynamic property
        call enable_property("Indicator Dynamic")

        ! Generate cube file for indicator dynamic
        write(*,*) "Generating indicator dynamic cube file..."
        call write_cube_file("indicator_dynamic.cube", "Indicator Dynamic", &
                           indicator_at_point, grid)
    endif

    write(*,*) "=== Property-based calculations completed ==="
endif

! Handle legacy cube file generation if requested
if (cube_flag) then
    write(*,*) "=== Generating legacy cube files ==="

    if (primcube) then
        a=1
        call cubefile(a,nameprim) !Generate cubefile with a Primitive
    end if

    if (aocube) then
        a=2
        call cubefile(a,nameao)  !Generate cubefile with AO
    end if

    if (mocube) then
        a=3
        call cubefile(a,namemo) !Generate cubefile with a MO
    end if

    if (denscube) then
        a=5
        call cubefile(a,namedens) !Generate cubefile with density from MO
    end if

    if (laplacian) then
        a=6
        call cubefile(a,namelap) !generate a cubefile with laplacian
    end if

    write(*,*) "=== Legacy cube files generated ==="
end if

! Handle on-top calculations (can be via property system or direct option)
if (option == "ontop" .or. ontop_cube_output) then
    write(*,*) "=== Entering On-top Calculation ==="

    if (ontop_hf) then
        write(*,*) "Using HF path"
        write(*,*) "WFX file: ", trim(wfxfilename)
        call wfxhf(wfxfilename)       ! Read HF wfx file
        call wrthfdm2()               ! Compute dm2 in primitives for HF
        call ontop('dm2hf.dat')       ! Calculate ontop with HF density
    else
        write(*,*) "Using CASSCF path"
        write(*,*) "DM2 file: ", trim(dm2name)
        call ontop(dm2name)           ! Original behavior for non-HF
    endif

    write(*,*) "=== On-top calculation completed ==="
endif

! Handle legacy Coulomb Hole calculations
if (option.eq."c1") then !calculate only c1 part of the Coulomb Hole
    write(*,*) "=== Calculating C1 part of Coulomb Hole ==="

    if (dmn) then
        !call dm2diff(dm2phfl,dm2phf,dm2pc1) WORKING ON IT
        call gridpoints() !for total integral
        call gridpoints2() !for radial scan
        call intracule(dm2pc1, .true.)
    else
        call wfxhf(wfxhfname) !read hf wfx file
        call wrthfdm2() !compute dm2 in primitives for HF, FCI (dens, |1rdm|^2, tot) and C1hole
        call gridpoints() !for total integral
        call gridpoints2() !for radial scan
        call intracule(dm2pc1, .false.)
    end if
    call angular_int() !radial scan of c1
    call radial_angular()  !integral of c1 and its moments

    write(*,*) "=== C1 calculation completed ==="
end if

if (option.eq."all") then !compute total Coulomb Hole
    write(*,*) "=== Calculating total Coulomb Hole ==="

    call wfxhf(wfxhfname)
    call wrthfdm2()
    call gridpoints()
    call gridpoints2()
    call intracule(dm2phf,.false.)
    call intracule(dm2p,.true.)
    call intracule(dm2phfl,.false.)
    ! call hole_ops(I_v1,I_v2) !subtract intracules
    call angular_int()
    call radial_angular()

    write(*,*) "=== Total Coulomb Hole calculation completed ==="
end if

if (option.eq."c1approx") then !compute approx c1hole with Becke-Roussel model
    write(*,*) "=== Calculating approximate C1 hole with Becke-Roussel model ==="

    call wfxhf(wfxhfname)    !read wfx of HF calculation (needed to compute dm2hf_prim)
    call wrthfdm2() !compute dm2 in primitives for HF, FCI (dens, |1rdm|^2, tot) and C1hole
    call gridpoints() !for total integral
    call gridpoints2() !for radial scan
    call intracule("dm1hf.dat",.false.) !compute intracule of |1rdm_hf|^2
    !compute c1hole
    call c1hole(90,1.d0)    !compute approximate c1 hole I(s,|1rdm_hf|^2)-\int_{r=0,inf}BRHOLE(r,s)*rho(r)dr
    call angular_int()
    call radial_angular()

    write(*,*) "=== Approximate C1 hole calculation completed ==="
end if

if (option.eq."intra") then !compute the intracule
    write(*,*) "=== Calculating intracule ==="

    call gridpoints()
    call gridpoints2()
    call intracule(dm2p,.true.)
    call angular_int()
    call radial_angular()

    write(*,*) "=== Intracule calculation completed successfully ==="
end if

write(*,*) "=== Program execution completed ==="

end program wavefunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   In this subroutine we read an input file where we specify what type of
!   calculations we want to do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinput()
use inputdat
use cubeinfo
use intracube
use located
use cube_module  ! For grid parameters
implicit none
character*80 :: name, line, keyword
integer :: i, j, ndims, io_stat
double precision :: max_val
logical :: active_dim
character*40 :: property_name, property_type

! Get input file name from command line
call getarg(1, name)
name = trim(name)
open(unit=3, file=name, status='OLD')

nameinp = name

! Set default values
readwfx = .false.
readlog = .false.
cube_flag = .false.
primcube = .false.
aocube = .false.
MOcube = .false.
denscube = .false.
gradient = .false.
laplacian = .false.
dmn = .false.
id = .false.
ontop_hf = .false.

! Basic file inputs (keeping at the beginning for compatibility)
call locate(3, "$wfxfile")
read(3, *) wfxfilename
if (wfxfilename.ne.'no') then
    readwfx = .true.
endif

call locate(3, "$logfile")
read(3, *) logfilename
if (logfilename.ne.'no') then
    readlog = .true.
endif

! Process either legacy options or new property-based format



    ! Check for any inactive dimensions
    rewind 3
    call locate(3, "$Grid")
    do i = 1, 3
        read(3, *, iostat=io_stat) max_val, step_size
            grid%max_val(i) = max_val
            grid%step_size(i) = step_size
    enddo

    rewind 3
    call locate(3, "$Grid_2")
    read(3, *) Grid_2

    rewind 3
    call locate(3, "$Grid_2")
    if (Grid_2) then
        do i = 1,3  !
            read(3, *, iostat=io_stat) max_val, step_size
            grid_2%max_val(i) = max_val
            grid_2%step_size(i) = step_size
        enddo
    else
        do i = 1, 3
            grid_2%max_val(i) = grid%max_val(i)
            grid_2%step_size(i) = grid%step_size(i)
        enddo
    endif


    ! Read property specifications
    if (scan_to_location(3, "$Properties")) then
        ! Allocate property flags with reasonable defaults
        max_properties = 100  ! Default maximum

        allocate(density_properties_enabled(5))  ! total, alpha, beta, spin, hf
        allocate(pair_density_properties_enabled(3))  ! full, c1, c2
        density_properties_enabled = .false.
        pair_density_properties_enabled = .false.
        indicator_dynamic_enabled = .false.

        ! Process each property line
        do
            read(3, '(A)', iostat=io_stat) line
            if (io_stat /= 0 .or. line(1:1) == '$') exit

            ! Remove leading/trailing spaces
            line = adjustl(trim(line))
            if (len_trim(line) == 0) cycle  ! Skip empty lines

            ! Convert to lowercase for easier parsing
            call to_lowercase(line)

            ! Check for single properties
            if (index(line, "density") > 0) then
                ! Enable total density by default
                density_properties_enabled(1) = .true.
                density_output_type = "total"

                ! Check for spin keyword
                if (index(line, "spin") > 0) then
                    ! Enable all three density types for spin
                    density_properties_enabled(1:3) = .true.  ! total, alpha, beta
                    density_properties_enabled(4) = .true.    ! spin density
                    density_output_type = "all"
                else if (index(line, "alpha") > 0) then
                    density_properties_enabled(2) = .true.
                else if (index(line, "beta") > 0) then
                    density_properties_enabled(3) = .true.
                else if (index(line, "hf") > 0) then
                    density_properties_enabled(5) = .true.
                end if
            end if

            ! Check for on-top pair density
            if (index(line, "on top") > 0 .or. index(line, "ontop") > 0) then
                ontop_cube_output = .true.
                option = "ontop"  ! This helps maintain legacy functionality
            end if

            ! Check for indicator dynamic
            if (index(line, "indicator") > 0) then
                indicator_dynamic_enabled = .true.
            end if

            ! Check for pair density
            if (index(line, "pair density") > 0 .or. index(line, "pair_density") > 0) then
                ! Enable full pair density by default
                pair_density_properties_enabled(1) = .true.
                pair_density_output_type = "full"

                ! Check for components keyword
                if (index(line, "component") > 0) then
                    ! Enable all pair density components
                    pair_density_properties_enabled = .true.  ! full, c1, c2
                    pair_density_output_type = "all"
                else if (index(line, "c1") > 0) then
                    pair_density_properties_enabled(2) = .true.
                else if (index(line, "c2") > 0) then
                    pair_density_properties_enabled(3) = .true.
                end if
            end if
        end do
    end if

    ! Handle remaining input sections for other functionality
    rewind 3
    call locate(3, "$On top")
    read(3, *) dm2name
    if (dm2name.eq."no") ontop_hf = .true.

    ! Only read these if needed (scanning properties)
    if (scan_to_location(3, "$scanning_props")) then
        read(3, *) n_atoms_scan
        allocate(atom_indices(n_atoms_scan))
        read(3, *) (atom_indices(i), i=1, n_atoms_scan)

        ! Read number of directions
        read(3, *) n_directions
        allocate(scan_directions(n_directions))
        read(3, *) (scan_directions(i), i=1, n_directions)

        ! Read distance parameters and convert to atomic units
        read(3, *) end_distance
        end_distance = end_distance/0.529177249d0  ! Convert Å to au

        read(3, *) step_size
        step_size = step_size/0.529177249d0  ! Convert Å to au
    endif
else
    ! Legacy input format
    ! Process cube file parameters
    call locate(3, "$cubefile")
    read(3, *) cube_flag
    if (cube_flag) then
        read(3, *) (center(i), i=1,3)
        read(3, *) (step(i), i=1,3)
        read(3, *) (np(i), i=1,3)
        rewind 3

        call locate(3, "$primitive")
        read(3, *) nameprim
        if (nameprim.ne.'no') primcube = .true.
        rewind 3

        if (readlog) then
            call locate(3, "$AO")
            read(3, *) nameao
            if (nameao.ne.'no') aocube = .true.
            read(3, *) cao
        endif

        rewind 3
        call locate(3, "$MO")
        read(3, *) namemo
        if (namemo.ne.'no') MOcube = .true.
        read(3, *) mo

        rewind 3
        call locate(3, "$density")
        read(3, *) namedens
        if (namedens.ne.'no') denscube = .true.

        rewind 3
        call locate(3, "$gradient")
        read(3, *) namegrad
        if (namegrad.ne.'no') gradient = .true.

        rewind 3
        call locate(3, "$laplacian")
        read(3, *) namelap
        if (namelap.ne.'no') laplacian = .true.

        rewind 3
    endif

    ! Legacy Coulomb Hole parameters
    rewind 3
    call locate(3, "$Coulomb Hole")
    read(3, *) option
    read(3, *) wfxhfname
    if (wfxhfname.ne."no") dmn = .true.

    ! Other legacy parameters
    call locate(3, "$ID")
    read(3, *) id
endif

close(3)

contains
    ! Helper function to scan to a specific section
    function scan_to_location(unit, section) result(found)
        integer, intent(in) :: unit
        character(*), intent(in) :: section
        logical :: found

        character(len=80) :: line

        rewind(unit)
        found = .false.
        do
            read(unit, '(A)', end=110, err=110) line
            if (index(line, trim(section)) > 0) then
                found = .true.
                exit
            endif
        enddo

        110 continue
        if (.not. found) rewind(unit)
    end function scan_to_location

    ! Helper function to check if one string contains another (case-insensitive)
    function scan_words(line, word) result(found)
        character(*), intent(in) :: line, word
        logical :: found
        character(len=80) :: line_lower, word_lower
        integer :: i

        ! Convert to lowercase for case-insensitive comparison
        line_lower = line
        word_lower = word
        do i = 1, len_trim(line_lower)
            if (line_lower(i:i) >= 'A' .and. line_lower(i:i) <= 'Z') then
                line_lower(i:i) = achar(iachar(line_lower(i:i)) + 32)
            endif
        enddo

        do i = 1, len_trim(word_lower)
            if (word_lower(i:i) >= 'A' .and. word_lower(i:i) <= 'Z') then
                word_lower(i:i) = achar(iachar(word_lower(i:i)) + 32)
            endif
        enddo

        found = index(trim(line_lower), trim(word_lower)) > 0
    end function scan_words



    ! Helper function to convert string to lowercase
    subroutine to_lowercase(str)
        character(*), intent(inout) :: str
        integer :: i

        do i = 1, len_trim(str)
            if (str(i:i) >= 'A' .and. str(i:i) <= 'Z') then
                str(i:i) = achar(iachar(str(i:i)) + 32)
            end if
        end do
    end subroutine to_lowercase
end subroutine readinput




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readintra(nameinp,dm2p,dm2phf,dm2phfl)
use thresholds
use quadinfo1
use quadinfo2
use intracube
use located
implicit none
character*80, intent(in) :: nameinp
character*80, intent(out) :: dm2p, dm2phf, dm2phfl
integer :: i, j
character*20 :: ccent
logical :: hf,hfl
logical :: beckw
   open(unit=3,file=nameinp,status='OLD')
   !read input for intracule
   call locate(3,'$DM2files')
   read(3,*) dm2p, dm2phf, dm2phfl
   if (dm2phf.eq.'no') dm2phf="dm2hf.dat"
   if (dm2phfl.eq.'no') dm2phfl="dm2phf.dat"
   call locate(3,'$Threshold') !threshold for integral screenings
   read(3,*) thresh
   rewind(3)
   call locate(3,'$radial_integral')
   read(3,*) radial_integral
   if (radial_integral) then   !read info to perform radial integral
      call centercalc() !calculate integration centers
      allocate(nradc(nquad))
      allocate(nangc(nquad))
      allocate(sfalpha(nquad))
      call locate(3,'Gauss-Legendre')
      !calculate nrad for each centre
      read(3,*) nradc(1)!reference number of radial nodes for each center
      read(3,*) a,b       !Choose limits of the radial integral
                          !if 0 the integral is computed from 0 to infty
      read(3,*) sfalpha !HEre is a fixed number but I will calculate it
      !calculate sfalpha for each centre
      rewind(3)
      call locate(3,'$Gauss-Lebedev')
      read(3,*) nangc(1) !fix it to 590
      rewind(3)
      !call auto()

      !call pdint() !compute Becke centres weights automatically (TO IMPLEMENT)
   end if !end radial integral

   call locate(3,'$radial_plot')
   read(3,*) radial_plot
   if (radial_plot) then
      read(3,*) r_plot_name
      read(3,*) nblock
      allocate(n_an_per_part(nblock)) !angular nodes per part
      allocate(tart(2,nblock))        !how many points per part
      allocate(stp(nblock))           !stepsize
      do i=1,nblock
         read(3,*) tart(:,i), stp(i), n_an_per_part(i)
      end do
   end if
   rewind(3) !end radial plot

   call locate (3,'$Vectorial_plot')
   read(3,*) cubeint
   if (cubeint) then
      read(3,*) cubeintraname
      read(3,*) (center_i(i), i=1,3) !cube centered in (x,y,z)
      read(3,*) (step_i(i), i=1,3) !distance between points in the axis
      read(3,*) (np3(i), i=1,3) !number of points for each axis
   end if

   close(3)
end subroutine readintra


subroutine centercalc()
use geninfo !contains cartes and natoms
use quadinfo1 !contains nquad and cent
implicit none
integer :: i,j,k, nqmax, sm, sm1, smp
double precision, allocatable, dimension(:,:) :: c1
integer, allocatable, dimension(:) :: nelec_c1
double precision, parameter :: trsh=1.d-1, zero=1.d-15, trs=5.d-1
double precision :: z
logical :: equal
!integer, allocatable, dimension(:) :: eqc
nqmax=(natoms*(natoms-1))+1
allocate(c1(3,nqmax))
allocate(nelec_c1(nqmax))
c1(:,1)=0.d0 !first center always at zero
sm=1
do j=1,natoms-1
   do k=j+1,natoms
      if (j.ne.k) then
         !if ((an(j).gt.1).and.(an(k).gt.1)) then !only for non-hydrogen atoms
         sm=sm+1
         c1(:,sm)=cartes(j,:)-cartes(k,:)
         nelec_c1(sm)=(an(i)-chrg(i))*(an(j)-chrg(j))
         z=abs(c1(3,sm))
         if (z.lt.trs) c1(3,sm)=0.d0 !to exploit symmetry
         smp=sm
         sm=sm+1 !count centres
         c1(3,sm)=-c1(3,smp)
         c1(1,sm)=-c1(1,smp) !I(x,y,z)=I(-x,-y,-z)
         c1(2,sm)=-c1(2,smp)
         nelec_c1(sm)=nelec_c1(smp)
      end if
   end do
end do
!check if there are equal centers
write(*,*) "I found",sm, "centres.", "Now I will discard equal centres"
nqmax=sm
sm=0
do i=2,nqmax-1
   do j=i+1,nqmax
      sm1=sm1+1
      z=abs(c1(1,i))+abs(c1(2,i))+abs(c1(3,i)) !sum x,y,z axis (if 0->equivalent centres exist)
      if (z.ne.0.d0) then                      !this center was not repeated
         equal=all(abs(c1(:,i)-c1(:,j)).lt.trs)
         if (equal) then             !center j is equal to center i
            sm=sm+1
            c1(:,j)=0.d0                     !set center j to zero
         end if
      end if
   end do
end do

nquad=nqmax-sm
allocate(cent(3,nquad))
allocate(nelec_cent(nquad))
sm=1
sm1=1
do i=2,nqmax
   z=abs(c1(1,i))+abs(c1(2,i))+abs(c1(3,i))
   if (z.ne.0.d0) then  !store all non-equivalent centres in a lower-dimensional array
      sm1=sm1+1
      cent(:,sm1)=c1(:,i)
      nelec_cent(sm1)=nelec_c1(i)
   end if
end do

write(*,*) "INTEGRATION CENTRES:", nquad
end subroutine centercalc

subroutine autocalc()
use geninfo
use wfxinfo
use quadinfo1
implicit none
!character*20 :: param
integer :: i
double precision :: nrad
double precision :: factor
do i=1,nquad
   factor=nelec_cent(i)/maxval(nelec_cent)
   nradc(i)=nrad*factor
end do





end subroutine autocalc






