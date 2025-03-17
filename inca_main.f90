!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program wavefunction !main program
!!!!!reads .wfx and .log files and performs calculations about this info
use inputdat  !information about the calculations we want to do
implicit none
integer :: a !defines to subroutine cubefile what function we want to represent: prim, ao, mo, dens
double precision, allocatable, dimension(:,:) :: r1,r2,r3!nodes
character*80 :: dm2p,dm2phf,dm2phfl, dm2pc1 !name of the DM2 files

!provisional stuff
character*80 :: dm1name, namerad

readwfx=.false. 
readlog=.false. 

call readinput()

write(*,*) "=== Starting Main Program ===" !JD
write(*,*) "Option selected: ", option !JD
write(*,*) "ontop_hf flag: ", ontop_hf !JD
write(*,*) "dm2name: ", trim(dm2name) !JD

! Add just before the ontop section
write(*,*) "=== Entering On-top Calculation ===" !JD
if (ontop_hf) then
    write(*,*) "Using HF path" !JD
    write(*,*) "WFX file: ", trim(wfxfilename) !JD
else
    write(*,*) "Using CASSCF path" !JD
    write(*,*) "DM2 file: ", trim(dm2name) !JD
endif

if (readwfx) call filewfx(wfxfilename)  !reads info from a wfx file 
if (readlog) call filelog(logfilename)  !reads info from a log file
if ((option.ne.'no').and.(option.ne.'ontop')) call readintra(nameinp,dm2p, dm2phf,dm2phfl)    !read intracule info from input file

if (cube) then
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
end if

if (option.eq."c1") then !calculate only c1 part of the Coulomb Hole
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
end if        


if (option.eq."all") then !compute total Coulomb Hole
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
end if        

if (option.eq."c1approx") then !compute approx c1hole with Becke-Roussel model
   call wfxhf(wfxhfname)    !read wfx of HF calculation (needed to compute dm2hf_prim)
   call wrthfdm2() !compute dm2 in primitives for HF, FCI (dens, |1rdm|^2, tot) and C1hole
   call gridpoints() !for total integral
   call gridpoints2() !for radial scan
   call intracule("dm1hf.dat",.false.) !compute intracule of |1rdm_hf|^2
   !compute c1hole
   call c1hole(90,1.d0)    !compute approximate c1 hole I(s,|1rdm_hf|^2)-\int_{r=0,inf}BRHOLE(r,s)*rho(r)dr
   call angular_int()
   call radial_angular()
end if


if (option.eq."intra") then !compute the intracule
   call gridpoints()
   call gridpoints2()  
   call intracule(dm2p,.true.)  
   call angular_int()
   call radial_angular()
   write(*,*) "intracule computed succesfully"
end if

if (option.eq."ontop") then !compute on top pair density
    if (ontop_hf) then  ! This is a HF calculation
        call wfxhf(wfxfilename)       ! Read HF wfx file
        call wrthfdm2()               ! Compute dm2 in primitives for HF
        call ontop('dm2hf.dat')           ! Calculate ontop with HF density
    else                  ! Non-HF calculation
        call ontop(dm2name)           ! Original behavior for non-HF
    end if
end if   

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
implicit none
character*80 :: name
integer :: i

call getarg(1,name)  !gets the name of the input file (writen in the prompt)
name=trim(name)      !remove the blank spaces of the string 'name'
open(unit=3,file=name,status='OLD') 

nameinp=name
!set up logical variables as false
readwfx=.false. 
readlog=.false.
primcube=.false.
aocube=.false.
MOcube=.false.
denscube=.false.
gradient=.false.
laplacian=.false.
dmn=.false.
id=.false.
ontop_hf=.false.

call locate(3,"$wfxfile") 
read(3,*) wfxfilename 
   if (wfxfilename.ne.'no') then
      readwfx=.true.
   end if

call locate(3,"$logfile")
read(3,*) logfilename
if (logfilename.ne.'no') then
      readlog=.true.
end if 

rewind 3

call locate(3,"$cubefile")
read(3,*) cube
if (cube) then 
   read(3,*) (center(i), i=1,3) !cube centered in (x,y,z)
   read(3,*) (step(i), i=1,3) !distance between points in the axis
   read(3,*) (np(i), i=1,3) !number of points for each axis  
   rewind 3
   call locate(3,"$primitive")
   read(3,*) nameprim
   if (nameprim.ne.'no') primcube=.true.
   !read(3,*) npr !number of primitive we want to print
   rewind 3
   if (readlog) then    
      call locate(3,"$AO")
      read(3,*) nameao
      if (nameao.ne.'no') aocube=.true. 
      read(3,*) cao  !number of ao we want to print
   end if 
   rewind 3
   call locate(3,"$MO")
   read(3,*) namemo
   if (namemo.ne.'no') MOcube=.true.    
   read(3,*) mo  !number of MO we want to print
   rewind 3
   call locate(3,"$density")
   read(3,*) namedens 
   if (namedens.ne.'no') denscube=.true.       
   rewind 3
   call locate(3,"$gradient")
   read(3,*) namegrad
   if (namegrad.ne.'no') gradient=.true.
   rewind 3
   call locate(3,"$laplacian")
   read(3,*) namelap
   if (namelap.ne.'no') laplacian=.true.
   rewind 3   
end if
rewind 3
call locate(3,"$Coulomb Hole") !compute c1hole with HF
read(3,*) option
read(3,*) wfxhfname
if (wfxhfname.ne."no") dmn=.true.
rewind 3
call locate(3,"$On top")
read(3,*) dm2name
if (dm2name.eq."no") ontop_hf=.true.

rewind 3
call locate(3,"$scanning_props")
read(3,*) n_atoms_scan
allocate(atom_indices(n_atoms_scan))
read(3,*) (atom_indices(i), i=1,n_atoms_scan)

! Read number of directions
read(3,*) n_directions
allocate(scan_directions(n_directions))
read(3,*) (scan_directions(i), i=1,n_directions)

! Read distance parameters and convert to atomic units
read(3,*) end_distance
end_distance = end_distance/0.529177249d0  ! Convert Å to au

read(3,*) step_size  
step_size = step_size/0.529177249d0  ! Convert Å to au

call locate(3,"$ID")
read(3,*) id


call locate(3,"$Ontop_Output_format")
read(3,*) ontop_cube_output
close(3)
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





