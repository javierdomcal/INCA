!!!!!!!!!!!!!!!!!!!!!!!!!!integrals of the intracule!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine angular_int()
use radial_scan
use quadinfo2
implicit none
   double precision, dimension(nradi) :: r_intra
   double precision :: in1,in2,in11,in22
   integer :: i,k,sm
   open(unit=3, file=r_plot_name)   
      sm=0
      r_intra=0.d0 !this is stored in quadratures module 
      do i=1,nradi        
            do k=1,smrd(i) !sum reduced angular points per radi  
                sm=sm+1            
                r_intra(i)=r_intra(i)+w2(sm)*I_vec2(sm)!Perform the angular quadrature 
            end do
            r_intra(i)=radi(i)**(2.d0)*r_intra(i)
            in1=r_intra(i)*radi(i)
            in2=r_intra(i)*radi(i)**2.d0
            in11=r_intra(i)/radi(i)
            in22=r_intra(i)/(radi(i)**2.d0)
            write(3,'(6(F14.8,1X))') radi(i), r_intra(i), in1,in2,in11,in22
      end do 
      close(3)
      write(*,*) "Angular integral done", "Radial scran file is at", r_plot_name
end subroutine
 !Takes I_vec from intracule subroutine, and performs its radial and/or angular integrals  
 
subroutine radial_angular()
use total_integral
use quadinfo1
implicit none
 double precision :: r_integral, m1,m2,m11,m22
 integer :: i
     WRITE(*, '(A, I4, A)') "Total integral using", nquad, " integration centres"
     r_integral=0.d0; m1=0.d0; m2=0.d0; m11=0.d0; m22=0.d0
     do i=1,np1
         r_integral=r_integral+w1(i)*I_vec1(i)
         m1=m1+w1_1(i)*I_vec1(i)
         m2=m2+w1_2(i)*I_vec1(i)
         m11=m11+w1_11(i)*I_vec1(i)
         m22=m22+w1_22(i)*I_vec1(i)
     end do 
     open(unit=3, file=int_mom_name)
     write(3,*) "n=0", "n=-1", "n=-2", "n=1", "n=2"
     write(3,*) r_integral, m1,m2,m11,m22
     close(3)
     write(*,*) "Angular and radial integral done", "Output file is", int_mom_name
end subroutine

subroutine cubeintra(ngrid,I_vec)
use geninfo
use intracube 
implicit none
integer, intent(in) :: ngrid
double precision, intent(in), dimension(ngrid) :: I_vec
integer :: i

     open(unit=3, file=cubeintraname)
     write(3,*) "CUBE FILE"
     write(3,*) "OUTER LOOP:X, MIDDLE LOOP:Y, INNER LOOP:Z"
     write(3,*) natoms, center_i(1), center_i(2), center_i(3)
     write(3,*) np3(1), step_i(1), 0.d0, 0.d0
     write(3,*) np3(2), 0.d0, step_i(2), 0.d0
     write(3,*) np3(3), 0.d0, 0.d0, step_i(3)
     do i=1,natoms
            write(3,*) an(i), chrg(i), cartes(i,1), cartes(i,2), cartes(i,3)
     end do
     do i=1,ngrid
           write(3,*) I_vec(i)
     end do  
     close(3)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
end subroutine
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SUBROUTINES TO OBTAIN QUADRATURE WEIGHTS AND NODES!!!!!!!!!!!!!!!!!!!!!!! 
subroutine gauherm(Lrtot,n) 
use qhermite !module contains gauss-hermite nodes and weights
!performs gauss-hermite quadrature                          
!we will obtain the nodes and weights depending on the degree of the polinomial (2n-1) --> (n) 
!(n is nx ny or nz in intracule.f90)
 implicit none
 !global variable
 integer, intent(in) :: Lrtot !degree of polynomial
 integer, intent(out) :: n    !number of gauss-hermite nodes
 !local variables
 integer :: factn, i2
 double precision :: Cons
 double precision, parameter :: pi=4.d0*atan(1.d0)
 
 if (Lrtot.gt.0) then
 !compute the number of nodes for exact Hermite quadrature
     if(MOD(Lrtot,2).eq.0) then !even 
           n=(Lrtot/2) +1 !nx->number of nodes
     else 
           n=(Lrtot+1)/2 
     end if
 else
     !only one node
           n=1
 end if  
 
 allocate(rh(n))  !nodes
 allocate(w_r(n)) !weights
 
 !store rh values depending on pol. degree
  
  if (n.eq.1) then !1 node-->degree of the pol is 0 or 1 
   rh(1)=0.0000000000000d0   
  else if (n.eq.2) then !2 nodes-->degree of the pol is 3 or 2                
   rh(2)=0.7071067811865475 
   rh(1)=-0.7071067811865475 
  else if (n.eq.3) then !3 nodes-->degree of pol is 5 or 4
   rh(2)=0.d0
   rh(3)= 1.224744871391589 
   rh(1)=-1.224744871391589  
  else if (n.eq.4) then   !7 or 6
   rh(3)=0.5246476232752903 
   rh(4)=1.650680123885785 
   rh(2)=-0.5246476232752903 
   rh(1)=-1.650680123885785 
  else if (n.eq.5) then    !9 or 8
   rh(1)=-2.020182870456086
   rh(2)=-0.9585724646138185
   rh(3)=0.000000000000000d0
   rh(4)=2.020182870456086
   rh(5)=0.9585724646138185
  else if (n.eq.6) then    !11 or 10
   rh(4)=0.4360774119276165 
   rh(5)=1.335849074013697 
   rh(6)=2.350604973674492   
   rh(3)=-0.4360774119276165 
   rh(2)=-1.335849074013697 
   rh(1)=-2.350604973674492
  else if (n.eq.7) then   !13 or 12
   rh(4)=0.d0
   rh(5)=0.8162878828589647 
   rh(6)=1.673551628767471 
   rh(7)=2.651961356835233
   rh(3)=-0.8162878828589647 
   rh(2)=-1.673551628767471 
   rh(1)=-2.651961356835233
  else if (n.eq.8) then    !15 or 14
   rh(5)=0.3811869902073221 
   rh(6)=1.157193712446780 
   rh(7)=1.981656756695843 
   rh(8)=2.930637420257244
   rh(4)=-0.3811869902073221 
   rh(3)=-1.157193712446780 
   rh(2)=-1.981656756695843 
   rh(1)=-2.930637420257244
  else if (n.eq.9) then     !17 or 16
   rh(6)=0.7235510187528376 
   rh(7)=1.468553289216668 
   rh(8)=2.266580584531843 
   rh(9)=3.190993201781528 
   rh(5)=0.d0
   rh(4)=-0.7235510187528376 
   rh(3)=-1.468553289216668 
   rh(2)=-2.266580584531843 
   rh(1)=-3.190993201781528 
  else if (n.eq.10) then     !19 or 18
    rh(6)=0.3429013272237046 
    rh(7)=1.036610829789514 
    rh(8)=1.756683649299882 
    rh(9)=2.532731674232790 
    rh(10)=3.436159118837738 
    rh(5)=-0.3429013272237046 
    rh(4)=-1.036610829789514 
    rh(3)=-1.756683649299882 
    rh(2)=-2.532731674232790 
    rh(1)=-3.436159118837738 
  else if (n.eq.11) then     !degree of pol is 21 or 20
    rh(6)=0.d0
    rh(7)=0.6568095668820998 
    rh(8)=1.326557084494933 
    rh(9)=2.025948015825755 
    rh(10)=2.783290099781652 
    rh(11)=3.668470846559583 
    rh(5)=-0.6568095668820998 
    rh(4)=-1.326557084494933 
    rh(3)=-2.025948015825755 
    rh(2)=-2.783290099781652 
    rh(1)=-3.668470846559583 
  else 
   write(*,*) "Error, total ang. momenta greater than 3"
  end if
  
  factn=1     !compute weights using the equation
  do i2=1,n   
   factn=factn*i2  
  end do
  Cons=((2.d0)**(dble(n-1)) * dble(factn)* sqrt(pi))/dble(n**2)
   
  do i2=1,n
   w_r(i2)=cons*(1.d0/(Hermite(rh(i2),n))**2.d0)  
  end do  

  contains
    function Hermite(x,n1) !
    double precision :: Hermite
    integer, intent(in) :: n1 !numer of nodes of pol
    double precision, intent(in) :: x !roots of Hermite pol
   
   if (n1.eq.1) then !0 deg
     Hermite=1.d0
   else if (n1.eq.2) then
     Hermite=2.d0*x 
   else if (n1.eq.3) then
     Hermite=4.d0*(x**2.d0) -2.d0
   else if (n1.eq.4) then
     Hermite=8.d0*(x**3.d0) -12.d0*x
   else if (n1.eq.5) then
     Hermite=16.d0*(x**4.d0)-48.d0*(x**2.d0)+12.d0  
   else if (n1.eq.6) then
     Hermite=32.d0*(x**5.d0)-160.d0*(x**3.d0)+120.d0*x  
   else if (n1.eq.7) then
     Hermite=64.d0*(x**6.d0)-480.d0*(x**4.d0)+720.d0*(x**2.d0)-120.d0
   else if (n1.eq.8) then
     Hermite=128.d0*(x**7.d0)-1344.d0*(x**5.d0)+3360.d0*(x**3.d0)-1680.d0*x    
   else if (n1.eq.9) then
     Hermite=256.d0*(x**8.d0)-3584.d0*(x**6.d0)+13440.d0*(x**4.d0)-13440.d0*(x**2.d0)+1680.d0  
   else if (n1.eq.10) then
     Hermite=512.d0*(x**9.d0)-9216.d0*(x**7.d0)+48384.d0*(x**5.d0)-80640.d0*(x**3.d0)+30240.d0*x
   else if (n1.eq.11) then
   Hermite=1024.d0*(x**10.d0)-23040.d0*(x**8.d0)+161280.d0*(x**6.d0)-403200.d0*(x**4.d0)+302400.d0*(x**2.d0)-30240.d0  
   end if       
   end function

end subroutine gauherm 

!******************************************************************************

!Gauss-Legendreren subroutine. t eta w balioak ematen diguna

!Modification of 
!https://github.com/NREL/OpenWARP/blob/master/source/NemohImproved/Nemoh/Solver/Core/Gaussm3.f90
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays t(1:n) and  w(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula. 
!* For detailed explanations finding weights & abscissas, see  "Numerical Recipes in Fortran */

 SUBROUTINE sub_GauLeg(X1,X2,t,w,n)
	!beste modu bat X1 eta X2 kendu, eta XL eta XM
	!bukaera t(i) = -z eta t(n+1-i)= +z izango dira

      IMPLICIT NONE
      
      INTEGER :: m,i,j
      INTEGER,intent(in) :: n   !Number of Gaussian points
      REAL*8, dimension(n),intent(out) :: W,t
      REAL*8 :: XM,XL,X1,X2,EPS,P1,P2,P3,pi,Z1,Z,PP

	!Relative precision
   	EPS = 1.D-14

  	!double precision arccosine. Pi value=3.14159
 	pi = DACOS(-1.D0)

	!N = number of Gauss Points
	!Roots are symmetric in the interval - so only need to find half of them  
   	m = (n + 1) / 2
	
	!The coats are going to be X1 = -1 and X2 = 1, Gauss-Legendre 
      	XM=0.5D0*(X1+X2)
      	XL=0.5D0*(X2-X1)


	!Loop over the desired roots
      	DO i = 1,m
         Z = DCOS (pi * (i - 0.25D0)/(n + 0.5D0))
	!Starting with the above approximation to the i-th root,
	!we enter the main loop of refinement by NEWTON'S method   
 10      P1 = 1.D0
         P2 = 0.D0

	!Loop up the recurrence relation to get the Legendre
	!polynomial evaluated at z                
         DO j = 1,n
            P3 = P2
            P2 = P1
            P1 = ((2.D0 * j - 1.D0) * Z * P2 - (j - 1.D0) * P3)/j
         END DO
!p1 is now the desired Legendre polynomial.
!We next compute pp, its derivative, by a standard relation involving also p2, 
!the polynomial of one lower order. 
         PP = n * (Z * P1 - P2)/(Z * Z - 1.D0)
         Z1 = Z
         Z = Z1 - P1/PP	      ! Newton's Method  */

         IF (DABS(Z-Z1) .GT. EPS) GO TO 10

	! Roots will be symmetric about the origin  
         t(i) = XM - XL * Z
         t(n + 1 - i) = XM + XL * Z
	!Compute the weight and its symmetric counterpart 
         W(i) = 2.D0 * XL/((1.D0 - Z * Z) * PP * PP)
         W(n + 1 - i) = W(i)
      END DO  

	RETURN   !not neccesary
END SUBROUTINE Sub_GauLeg 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine angular_nodes(gpt,nAng)                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Converts angular node option to number of angular points (nAng)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 
integer, intent(in) :: gpt
integer, intent(out) :: nAng
 if (gpt.eq.1) nAng=6
 if (gpt.eq.2) nAng=14
 if (gpt.eq.3) nAng=26
 if (gpt.eq.4) nAng=38
 if (gpt.eq.5) nAng=50
 if (gpt.eq.6) nAng=74
 if (gpt.eq.7) nAng=86
 if (gpt.eq.8) nAng=110
 if (gpt.eq.9) nAng=146
 if (gpt.eq.10) nAng=170
 if (gpt.eq.11) nAng=194
 if (gpt.eq.12) nAng=230
 if (gpt.eq.13) nAng=266
 if (gpt.eq.14) nAng=302
 if (gpt.eq.15) nAng=350
 if (gpt.eq.16) nAng=434
 if (gpt.eq.17) nAng=590
 if (gpt.eq.18) nAng=770
 if (gpt.eq.19) nAng=974
 if (gpt.eq.20) nAng=1202
 if (gpt.eq.21) nAng=1454
 if (gpt.eq.22) nAng=1730
 if (gpt.eq.23) nAng=2030
 if (gpt.eq.24) nAng=2354
 if (gpt.eq.25) nAng=2702
 if (gpt.eq.26) nAng=3074
 if (gpt.eq.27) nAng=3470
 if (gpt.eq.28) nAng=3890
 if (gpt.eq.29) nAng=4334
 if (gpt.eq.30) nAng=4802
 if (gpt.eq.31) nAng=5294
 if (gpt.eq.32) nAng=5810
end subroutine angular_nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


