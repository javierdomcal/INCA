!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gridpoints() 
use quadinfo1   
use total_integral
! Compute grid points to perform a multicenter radial integral from 0 to infinity    
!or from a to b (specify in the input file).
! Radial quadrature:Gauss-Legendre ; Angular quadrature: Gauss-Levedeb.
! Using Becke's weights to compute multicenter integrals.
! Neglects grid points with negative z using intracule symmetry, I(\Vec{r})=I(\Vec{-r})
! Neglects grid points with 0 weight.
! Output: r1(grid points), w1 (weights), and moments of I(r)*r**(n), where n{-2,2}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!global variables in quadinfo and total integral
!local variables
!1-levedeb quadrature
integer :: nang
double precision, allocatable, dimension(:,:) :: r_lb
double precision, allocatable, dimension(:) :: Wlb
double precision, allocatable, dimension(:) :: theta, phi
!2-legendre quadrature
integer :: nrad
double precision, allocatable, dimension(:) :: xl_i, wl_i
double precision, allocatable, dimension(:) :: radius !variable change from 0 to infty
!Legendre+Lebedev combination
double precision, allocatable, dimension(:) :: weight, srweight !total weight for each point (before neglecting points)
double precision, allocatable, dimension(:) :: we1,we2,we11,we22 !weight with the moments (n=-2 to n=2
!3-becke vi
double precision, allocatable, dimension(:) :: w_beck
!-starting grid point
integer, allocatable, dimension(:) :: smn
integer :: rgrid, rrgrid !number of sym reduced grid points
double precision, allocatable, dimension(:,:) :: fgr, brrg, rrg 
integer :: i, ia, ir, sm, smp, sma, smnn, j, k, smr, smpr, ngrid, i1
integer :: np !number of points
double precision, parameter :: pi=4.d0*atan(1.d0)
double precision, parameter :: trsh=1.d-15, trsh2=1.d-16
double precision :: xs


!count maximum number of points
allocate(smn(nquad)) !count number of points per quadrature centre
np=0
do i=1,nquad
  np=np+nradc(i)*nAngc(i)
end do
allocate(weight(np))
allocate(we1(np)); allocate(we2(np)); allocate(we11(np)); allocate(we22(np))
allocate(fgr(3,np)); allocate(brrg(3,np))
sm=0 !count total grid points (1-->np)
smp=0 !count last grid point of previous centre
smr=0  !count sym reduced grid points
smpr=0  !count sym reduced gp of previous centre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i1=1,nquad       !loop over centres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nrad=nradc(i1)
    nAng=nAngc(i1)
    !!Obtain radial nodes and weights (Gauss-Legendre)!!!!!!!!!!!!!!!!!!!!!
    allocate(xl_i(nrad))
    allocate(wl_i(nrad))
    call sub_GauLeg(-1.d0,1.d0,xl_i,wl_i,nrad)
    !!!Obtain angular nodes and weights (Gauss-Lebedev)!!!!!!!!!!!!!!!!!!!!!!!  
    allocate(r_lb(3,nang))
    allocate(Wlb(nang))
    Wlb=0.d0
    r_lb=0.d0
    if (nAng.eq.6) call LD0006(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !6 grid points
    if (nAng.eq.14) call LD0014(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !14 grid points
    if (nAng.eq.26) call LD0026(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !26 grid points 
    if (nAng.eq.38) call LD0038(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !38 grid points
    if (nAng.eq.50) call LD0050(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !50 grid points
    if (nAng.eq.74) call LD0074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !74 grid points
    if (nAng.eq.86) call LD0086(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !86 grid points
    if (nAng.eq.110) call LD0110(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !110 grid points
    if (nAng.eq.146) call LD0146(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng)  !146 grid points
    if (nAng.eq.170) call LD0170(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !170 grid points
    if (nAng.eq.194) call LD0194(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !194 grid points 
    if (nAng.eq.230) call LD0230(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !230 grid points
    if (nAng.eq.266) call LD0266(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !266 grid points 
    if (nAng.eq.302) call LD0302(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !302 grid points
    if (nAng.eq.350) call LD0350(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !350 grid points
    if (nAng.eq.434) call LD0434(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !434
    if (nAng.eq.590) call LD0590(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !590
    if (nAng.eq.770) call LD0770(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !770
    if (nAng.eq.974) call LD0974(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !974
    if (nAng.eq.1202) call LD1202(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !1202
    if (nAng.eq.1454) call LD1454(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !1454
    if (nAng.eq.1730) call LD1730(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !1730
    if (nAng.eq.2030) call LD2030(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !2030
    if (nAng.eq.2354) call LD2354(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !2354
    if (nAng.eq.2702) call LD2702(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !2702
    if (nAng.eq.3074) call LD3074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !3074
    if (nAng.eq.3470) call LD3470(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !3470
    if (nAng.eq.3890) call LD3890(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !3890
    if (nAng.eq.4334) call LD4334(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !4334
    if (nAng.eq.4802) call LD4802(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !4802
    if (nAng.eq.5294) call LD5294(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !5294
    if (nAng.eq.5810) call LD5810(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,nAng) !5810
    allocate(theta(nang))
    allocate(phi(nang))
    do i=1,nAng   !compute angles     
        theta(i)= acos(r_lb(3,i)) 
        if (sin(theta(i)).ne.0.d0) then
             xs=r_lb(1,i)/sin(theta(i))
             if (xs.gt.1.d0) xs=1.d0
             if (xs.lt.-1.d0) xs=-1.d0
              phi(i)=acos(xs)
              if (r_lb(2,i).lt.0.d0) phi(i)=-phi(i)
        else
              phi(i)=0.d0
        end if  
        Wlb(i)=Wlb(i)*2.d0*pi !store 2pi factor on the weight (see solid angle integral)
    end do                     !should be 4pi but we have a factor of 1/2 in the intracule!!
    deallocate(r_lb)
    !!!!!!!!!!!!!Compute grid points of quadrature!!!!!!!!!!!
    nGrid=nAng*nrad
    allocate(radius(nrad))   
    !compute radial points for the quadrature from 0 to infty or a to b
    if ((abs(a)+abs(b)).ne.0.d0) then
            write(*,*) "Integral is done from",a,"to", b    
       do ir=1,nrad    
          radius(ir)=(b-a)*0.5d0*xl_i(ir)+(a+b)*0.5d0
          wl_i(ir)=(b-a)*0.5d0*wl_i(ir) !radial integration (see pdf) 
       end do   
    else
       do ir=1,nrad    
          radius(ir)=(1.d0+xl_i(ir))/(1.d0-xl_i(ir))*sfalpha(i1) 
          wl_i(ir)=(2.d0*sfalpha(i1)/((1.d0-xl_i(ir))**2.d0))*wl_i(ir)*radius(ir)**2.d0    
        end do   
    end if    
    deallocate(xl_i)
 
    !compute grid points for all the becke centers

    smn(i1)=0
    do ir=1,nrad
        do ia=1,nAng
               sm=sm+1   !count total grid points    
               fgr(1,sm)=radius(ir)*cos(phi(ia))*sin(theta(ia))+cent(1,i1)
               fgr(2,sm)=radius(ir)*sin(phi(ia))*sin(theta(ia))+cent(2,i1)
               fgr(3,sm)=radius(ir)*cos(theta(ia))+cent(3,i1)
               !store total grid points (all quadratures)              
          !     if (fgr(3,sm).ge.0.d0) then
                  smn(i1)=smn(i1)+1 !sum the number of sym reduced points of each quadrature                  
          !     end if                  
        end do   
    end do 
    deallocate(theta); deallocate(phi)
    !write(*,*) "Total number of gp after center",i1,"=", sm, ngrid
    !neglect points by symmetry!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !write(*,*) nGrid-(smn(i1)), "Points will be neglected in center", i1     
    !rGrid=smn(i1)  !reduced grid points of current center
    !store reduced points and their total weight
    sm=smp   !start again from 1st point of center
    smr=smpr !start from 1st reduced point of center
    smnn=0   
!    write(*,*) "*********Neglecting points by symmetry*******************"
    do j=1,nrad
       do k=1,nAng  
          sm=sm+1                     !sum total grid points
          !if (fgr(3,sm).ge.0.d0) then !reduce points by symmetry
              smr=smr+1               !sum sym reduced grid points   
              brrg(:,smr)=fgr(:,sm)   !store sym reduced points
         !     if (fgr(3,sm).lt.trsh) then    !z is 0, do not multiply by 2
                   weight(smr)=Wlb(k)*wl_i(j)
         !     else
         !          weight(smr)=2.d0*Wlb(k)*wl_i(j) !z is positive, use sym (I(z)=I(-z))
         !     end if
              !store moments at the weight
              we1(smr)=weight(smr)*radius(j)**(-1.d0)
              we2(smr)=weight(smr)*radius(j)**(-2.d0)
              we11(smr)=weight(smr)*radius(j)**(1.d0)
              we22(smr)=weight(smr)*radius(j)**(2.d0)
        !  else if (fgr(3,sm).lt.0.d0) then 
              !sym neglected point    
       !       smnn=smnn+1
       !   end if  
       end do             
    end do   
    smpr=smr !store last reduced point of the quadrature
    smp=sm   !store last point of the quadrature
!    write(*,*) smnn, "points have been neglected in center number", i1  
    deallocate(wl_i) !deallocate legendre weigths
    deallocate(radius) !deallocate legendre nodes
    deallocate(wlb)   !deallocate levedeb weights
!    write(*,*) "Total number of gp after center", i1, "=", sm
!    write(*,*) "Total number of reduced gp after center", i1, "=", smr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end do !end loop over quadratures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rgrid=sum(smn) !set total number of grid points
write(*,*) "There are", rgrid, "symmetry reduced grid points for total integration"
!open(unit=3,file="all_points.dat")
!do i=1,np
!    write(3,*) fgr(1,i), fgr(2,i), fgr(3,i)
!end do
!close(3)
!open(unit=3,file="sym_red_gp.dat")
!do i=1,rgrid
!    write(3,*) brrg(1,i), brrg(2,i), brrg(3,i)
!end do
!close(3)
allocate(rrg(3,rgrid)); allocate(srweight(rgrid))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sm=0
do i=1,nquad                   !store the sym reduced points in a smaller matrix 
 do j=1,smn(i)                 !this is useless
   sm=sm+1
   rrg(:,sm)=brrg(:,sm)
   srweight(sm)=weight(sm)
 end do  
end do 

deallocate(brrg); deallocate(fgr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nquad.gt.1) then         !we have more than one quadrature centre
        !write(*,*) "number of points per centre", smn(:)
        allocate(w_beck(rgrid)) !compute becke weights for each point    
        call becke(rrg,smn,rgrid,nquad,cent,w_beck,nelec_cent)
        do i=1,rgrid
            srweight(i)=srweight(i)*w_beck(i) !store becke weight into the total one
        end do
        deallocate(w_beck)
       !reduce points with 0 weight
        sm=0
        sma=0
        do i=1,rgrid
            if ((srweight(i).le.trsh2)) then !not 1st center and little weight
               sm=sm+1                       !sum number of neglected points  
            end if        
        end do
        rrgrid=rgrid-sm !number of sym and weight reduced points
else !only 1 quadrature centre--> No Becke    
        rrgrid=rgrid
end if
!write(*,*) sm, "neglected points of zero weight"
deallocate(smn)

allocate(w1(rrgrid)) !reduced weight
allocate(r1(3,rrgrid))  !doubly reduced points
allocate(w1_1(rrgrid));allocate(w1_2(rrgrid)) !moments
allocate(w1_11(rrgrid));allocate(w1_22(rrgrid))

sm=0
!write(*,*) "Final number of grid points=", rrgrid
do i=1,rgrid
     if ((srweight(i).gt.trsh2)) then !remove points with low (zero) weight
           sm=sm+1 
           w1(sm)=srweight(i) 
           r1(:,sm)=rrg(:,i)           
           w1_1(sm)=we1(i) ; w1_2(sm)=we2(i)
           w1_11(sm)=we11(i) ; w1_22(sm)=we22(i)
     end if
end do

np1=rrgrid

!open(unit=3,file="sym_weight_red_gp.dat")
!do i=1,rrgrid
!    write(3,*) rrrg(1,i), rrrg(2,i), rrrg(3,i)
!end do
!close(3)
deallocate(weight);deallocate(rrg); deallocate(srweight)
deallocate(we1);deallocate(we2);deallocate(we11); deallocate(we22)
close(4)
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine becke(rrg,sumq,rgrid,nquad,cent,w_beck,nelec_cent)                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Based on A. D. Becke's paper from 1987, A multicenter numerical integration scheme for   !
! polyatomic molecules.                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes Becke's weights for each (sym reduced) point of the grid, with a weight for each!
! center. Then normalizes this weight so that the sum of all the weights of a single point !
! is 1. Finally stores the weight of the corresponding grid point in the array w_beck,     !
! since we will only use the weight of the quadrature centre from where the curren pòints  !
! have been originated.                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! As input: grid points, number of points per quadrature center, number of quad. centers   !
! As output: becke weights for each grid point.                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 implicit none
 integer, intent(in) :: rgrid, nquad !total grid, number of quadratures
 double precision, intent(in), dimension(3,nquad) :: cent !center of the quadratures
 double precision, intent(in), dimension(3,rgrid) :: rrg  !(reduced) grid points
 integer, intent(in), dimension(nquad) :: nelec_cent !weight of each center
 integer, dimension(nquad) :: sumq            !number of grid points per quad.
 double precision, intent(out), dimension(rgrid) :: w_beck !becke weights for each point
 !local variables
 double precision, dimension(rgrid,nquad) :: w_becke !weights of different quad. centers at same point
 integer :: i1,i,j,sm
 double precision, dimension(nquad,nquad) :: Rij, Xi !distance between centers, equation A4 xi
 double precision, dimension(nquad) :: P             !equation 13
 double precision :: ri, rj, s_ij, Ptot
 double precision ::  mu_ij, v_ij, a_ij, u_ij !eq.11, eq.a2, eq.a5, eq.a6 

 !compute R_ij
 do i=1,nquad
  do j=1,nquad
     Rij(i,j)=sqrt(sum((cent(:,i)-cent(:,j))**2.d0))  !distance between centers
     Xi(i,j)=nelec_cent(i)/nelec_cent(j)        !relacio entre centres
  end do
 end do     
 !compute Becke weights from all centres for all the grid points
 w_becke=0.d0     
 sm=0     
 do i1=1,rgrid !for all the grid points
   sm=sm+1     !count grid points
   do i=1,nquad   
        P(i)=1.d0
        do j=1,nquad
           if (i.ne.j) then        
             ri=sqrt(sum((rrg(:,i1)-cent(:,i))**2.d0)) !distance to center i (from point r)
             rj=sqrt(sum((rrg(:,i1)-cent(:,j))**2.d0)) !distance to center j (from point r)
             mu_ij=(ri-rj)*(Rij(i,j)**(-1.d0))
             u_ij=(xi(i,j)-1.d0)/(xi(i,j)+1.d0)  
             a_ij=u_ij*((u_ij**2.d0)-1.d0)**(-1.d0)
             !set limit for a_ij!!!!!!!!!!!!!
             !you can change this if you want
             if (a_ij.lt.-0.5d0) a_ij=-0.5d0 
             if (a_ij.gt.0.5d0) a_ij=0.5d0
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             v_ij=mu_ij+a_ij*(1.d0-mu_ij**2.d0) 
             s_ij=0.5d0*(1.d0-f_k(v_ij))
             P(i)=P(i)*s_ij
           end if              
        end do
       w_becke(i1,i)=P(i) !store weight of quadrature i at point i1     
   end do 
   Ptot=sum(P)
   w_becke(i1,:)=w_becke(i1,:)*(Ptot**(-1.d0)) !normalize weight to fulfill equation 3   
 end do
 
 !store single becke weight for each grid point in a single array
 sm=0
 w_beck=0.d0
 do i=1,nquad
  do j=1,sumq(i)
    sm=sm+1   
    w_beck(sm)=w_becke(sm,i)
  end do       
 end do

 contains        
     function f_k(val)   !equation 20, with value k=3
     implicit none
     double precision :: f_k
     double precision, intent(in) :: val
     double precision :: vl
     integer :: i0
     integer, parameter :: k=3
     vl=val
         do i0=1,k         
            f_k=pf(vl)
            vl=f_k
         end do       
     end function
     function pf(vl)  !equation 19
     implicit none
     double precision :: pf
     double precision, intent(in) :: vl
         pf=1.5d0*vl-0.5d0*vl**(3.d0)     
     end function
end subroutine becke         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gridpoints2()                                                    !
    ! Compute grid points for with manually given radius (in the input file)!   
    use quadinfo2 
    use radial_scan
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    !local variables
    !1-levedeb quadrature
    double precision, allocatable, dimension(:,:) :: r_lb
    double precision, allocatable, dimension(:) :: wlb
    double precision, allocatable, dimension(:) :: theta, phi
    double precision, allocatable, dimension(:) :: pweight
    !grid points
    double precision, allocatable, dimension(:,:) :: gpt !grid points without reducing sym
    double precision :: z, step
    integer :: i, ia, ir, sm, smnn, j, n_an, smrad
    double precision, parameter :: pi=4.d0*atan(1.d0)
    double precision, parameter :: trsh=1.d-15, trsh2=1.d-16
    double precision :: xs
    
    nradi=sum(stp)
    allocate(radi(nradi))
    
    np2=0; sm=0
    do i=1,nblock
        step=abs(tart(2,i)-tart(1,i))/dble(stp(i))
        do j=1,stp(i)
            sm=sm+1
            radi(sm)=tart(1,i)+step*dble(j) !compute radial points
        end do 
        np2=np2+stp(i)*n_an_per_part(i) !compute number of grid points
    end do
    
    allocate(smrd(nradi)) !count number of grid points per radi
    allocate(gpt(3,np2)); allocate(pweight(np2)) 
    
    smrad=0 !sum over radius
    sm=0    !sum for grid points
    smnn=0  !sum symmetry reduced grid points
    do i=1,nblock !loop for each radius fragment
        n_an=n_an_per_part(i)
        allocate(r_lb(3,n_an))
        allocate(Wlb(n_an))
        !!!Obtain angular nodes and weights (Gauss-Lebedev)!!!!!!!!!!!!!!!!!!!!!!!   
        if (n_an.eq.6) call LD0006(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !6 grid points
        if (n_an.eq.14) call LD0014(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !14 grid points
        if (n_an.eq.26) call LD0026(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !26 grid points 
        if (n_an.eq.38) call LD0038(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !38 grid points
        if (n_an.eq.50) call LD0050(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !50 grid points
        if (n_an.eq.74) call LD0074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !74 grid points
        if (n_an.eq.86) call LD0086(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !86 grid points
        if (n_an.eq.110) call LD0110(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !110 grid points
        if (n_an.eq.146) call LD0146(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an)  !146 grid points
        if (n_an.eq.170) call LD0170(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !170 grid points
        if (n_an.eq.194) call LD0194(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !194 grid points 
        if (n_an.eq.230) call LD0230(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !230 grid points
        if (n_an.eq.266) call LD0266(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !266 grid points 
        if (n_an.eq.302) call LD0302(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !302 grid points
        if (n_an.eq.350) call LD0350(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !350 grid points
        if (n_an.eq.434) call LD0434(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !434
        if (n_an.eq.590) call LD0590(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !590
        if (n_an.eq.770) call LD0770(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !770
        if (n_an.eq.974) call LD0974(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !974
        if (n_an.eq.1202) call LD1202(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !1202
        if (n_an.eq.1454) call LD1454(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !1454
        if (n_an.eq.1730) call LD1730(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !1730
        if (n_an.eq.2030) call LD2030(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !2030
        if (n_an.eq.2354) call LD2354(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !2354
        if (n_an.eq.2702) call LD2702(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !2702
        if (n_an.eq.3074) call LD3074(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !3074
        if (n_an.eq.3470) call LD3470(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !3470
        if (n_an.eq.3890) call LD3890(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !3890
        if (n_an.eq.4334) call LD4334(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !4334
        if (n_an.eq.4802) call LD4802(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !4802
        if (n_an.eq.5294) call LD5294(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !5294
        if (n_an.eq.5810) call LD5810(r_lb(1,:),r_lb(2,:),r_lb(3,:),Wlb,n_an) !5810
     
        allocate(theta(n_an))
        allocate(phi(n_an))
        do j=1,n_an   !compute angles     
            theta(j)= acos(r_lb(3,j))
            if (sin(theta(j)).ne.0.d0) then
                xs=r_lb(1,j)/sin(theta(j))
                if (xs.gt. 1.d0) xs=1.d0
                if (xs.lt.-1.d0) xs=-1.d0
                phi(j)=acos(xs)
                if (r_lb(2,j).lt. 0.d0) phi(j)=-phi(j)
            else
                phi(j)=0.d0
            end if
            Wlb(j)=Wlb(j)*2.d0*pi !store 2pi factor on the weight (see solid angle integral)
        end do   
        deallocate(r_lb)
        !angles and weights are computed for each block
    
        !!!!!!!!Compute grid points!!!!!!!!!!!!!!!!!!!!!!!!!! 
        do ir=1,stp(i) !loop for each radius inside the block
            smrad=smrad+1 !sum over radius
            smrd(smrad)=0
            do ia=1,n_an  !sum over angular points
                sm=sm+1 !sum grid points
                z=radi(smrad)*cos(theta(ia)) !compute z
                if (z.ge.0.d0) then
                    smnn=smnn+1 !sum the number of reduced points of each quadrature
                    smrd(smrad)=smrd(smrad)+1 !sum number of points per radius
                    gpt(1,smnn)=radi(smrad)*cos(phi(ia))*sin(theta(ia))
                    gpt(2,smnn)=radi(smrad)*sin(phi(ia))*sin(theta(ia))
                    gpt(3,smnn)=z
                    if (z.lt.trsh) then !Don't use 0.d0!!!!!!!!!!
                        pweight(smnn)=Wlb(ia)
                    else
                        pweight(smnn)=2.d0*Wlb(ia)
                    end if   
                end if
            end do   
         end do  
         deallocate(wlb); deallocate(theta); deallocate(phi) 
    end do   !end loop for radius fragment
     
    write(*,*) "Original number of grid points", np2,"=",sm
    write(*,*) "Number of grid points after I(r)=I(-r)", smnn
    write(*,*) "Number of grid points per radi", smrd(:)
    np2=smnn
    allocate(r2(3,smnn))
    allocate(w2(smnn))
    do i=1,np2  !store the data in reduced size matrix
        r2(:,i)=gpt(:,i) 
        w2(i)=pweight(i)          
    end do
    deallocate(gpt)
end subroutine gridpoints2


subroutine gridpoints3()
    use intracube
    implicit none      
    !local variables
    integer :: i,j,k,sm
    double precision, dimension(3) :: xm  

    np3t=np3(1)*np3(2)*np3(3)
    allocate(r3(3,np3t))
    do i=1,3
        if (mod(np3(i),2).eq.0) then !odd number
            xm(i)=center_i(i)-(step_i(i)/2)-((np3(i)-2)/2)*step_i(i)
        else !even number  !Calculates the starting point acording to cubeinfo data
             xm(i)=center_i(i)-((np3(i)-1)/2)*step_i(i)
        end if
    end do
    sm=0
    do i=1,np3(1)         !Depending on a calculates the point with a different function
        do j=1,np3(2)
            do k=1,np3(3)
                sm=sm+1
                r3(1,sm)=xm(1)+step_i(1)*(i-1)
                r3(2,sm)=xm(2)+step_i(2)*(j-1)
                r3(3,sm)=xm(3)+step_i(3)*(k-1)
            end do   
        end do
    end do
    write(*,*) "xm", xm(:)
    write(*,*) "np_i", np3(:)
    do i=1,np3t
        write(7,*) r3(:,i)
    end do  
end subroutine gridpoints3    


