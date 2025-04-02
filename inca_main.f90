program wavefunction
use inputdat
use properties
use cube_module
use density_properties
use ontop_module
use correlation_indicators
use pair_density_module
use geninfo
use wfxinfo
use wfxinfo_hf
implicit none
integer :: i

call readinput()

! Initialize wavefunction
if (readwfx) call filewfx(wfxfilename)
if (readlog) call filelog(logfilename)

! Initialize properties module with max_properties from inputdat
call initialize_properties(max_properties)

! Make property information available to the properties module
call set_property_info(requested_properties, property_arguments, num_properties)

write(*,*) num_properties
write(*,*) requested_properties
write(*,*) property_arguments

! Enable requested properties
do i = 1, num_properties
    if (requested_properties(i) == "density") then
        call density_calculation(property_arguments(i))
    else if (requested_properties(i) == "on_top") then
        call ontop_calculation(dm2_file, property_arguments(i))
    else if (requested_properties(i) == "indicator") then
        call indicator_calculation(property_arguments(i))
    else if (requested_properties(i) == "pair_density") then
        call pair_density_calculation(dm2_file, dm2hf_file, dm2hfl_file, &
                                      property_arguments(i))
    else if (requested_properties(i) == "pair_density_nucleus") then
        call pair_density_nucleus_calculation(dm2_file, dm2hf_file, dm2hfl_file, &
                                             property_arguments(i))
end do



! Clean up memory
if (allocated(requested_properties)) deallocate(requested_properties)
if (allocated(property_arguments)) deallocate(property_arguments)

end program wavefunction

subroutine readinput()
use inputdat
use cubeinfo
use intracube
use located
implicit none
character*80 :: name, line
integer :: i, io

! Get input file name from command line
call getarg(1, name)
name = trim(name)
open(unit=3, file=name, status='OLD')
nameinp = name

! Read wavefunction file names
call locate(3, "$wfxfile")
read(3, *) wfxfilename
readwfx = trim(wfxfilename) /= 'no'

call locate(3, "$logfile")
read(3, *) logfilename
readlog = trim(logfilename) /= 'no'

! Read grid parameters
rewind(3)
call locate(3, "$Grid")
do i = 1, 3
    read(3, *) grid%max_vals(i), grid%step_sizes(i)
end do

rewind(3)
call locate(3, "$Grid_2")
read(3, *) Grid_2_flag
if (Grid_2_flag) then
    do i = 1, 3
        read(3, *) grid_2%max_vals(i), grid_2%step_sizes(i)
    end do
else
    grid_2 = grid
end if

! Read properties
rewind(3)
call locate(3, "$Properties")
read(3, *) num_properties

allocate(requested_properties(num_properties))
allocate(property_arguments(num_properties))

do i = 1, num_properties
    read(3, '(A)', iostat=io) line
    if (io /= 0) exit
    line = adjustl(trim(line))

    ! Split line into property name and arguments
    read(line, *) requested_properties(i), property_arguments(i)
end do

! Read DM2 file names
rewind(3)
call locate(3, "$DM2files")
read(3,*) dm2_file, dm2hf_file, dm2hfl_file
if (dm2hf_file == 'no') dm2hf_file = "dm2hf.dat"
if (dm2hfl_file == 'no') dm2hfl_file = "dm2hfl.dat"

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
sm1=0 ! Initialize sm1 to avoid undefined variable
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