!*******************************************************************************
function Prim(x,y,z,npr) 
   use geninfo
   implicit none
   double precision :: Prim
   double precision, intent(in) :: x, y, z
   integer, intent(in) :: npr
   double precision :: dx, dy, dz, r2, expterm
   
   
   dx = x-cartes(Ra(npr),1)
   dy = y-cartes(Ra(npr),2)
   dz = z-cartes(Ra(npr),3)
   r2 = dx**2.d0 + dy**2.d0 + dz**2.d0
   expterm = -Alpha(npr)*r2
   
   
   Prim = (dx)**(TMN(npr,1))* &
          (dy)**(TMN(npr,2))* &
          (dz)**(TMN(npr,3))* &
          exp(expterm)
   
end function

  function expr(x,y,z,xa,ya,za,npr) 
   use geninfo
   implicit none
   double precision :: expr
   double precision, intent(in) :: x, y, z,xa,ya,za
   integer, intent(in) :: npr
   expr=exp(-Alpha(npr)*((x-xa)**2.d0+(y-ya)**2.d0+(z-za)**2.d0))
   !write(*,*) "Expr=-------------", expr
 end function

   
!**************************************************************************
function sd(x,y,z,npr) !second derivative of primitives analitically
use geninfo
implicit none
double precision :: sd, Prim
double precision, intent(in) :: x, y, z
integer, intent(in) :: npr
double precision :: laplx, laply, laplz
double precision :: t,m,n
double precision :: xa,ya,za

t=0.d0 
m=0.d0 
n=0.d0
xa=0.d0
ya=0.d0 
za=0.d0
laplx=0.d0 
laply=0.d0
laplz=0.d0
sd=0.d0

t=TMN(npr,1)
m=TMN(npr,2)
n=TMN(npr,3)
xa=cartes(Ra(npr),1)
ya=cartes(Ra(npr),2)
za=cartes(Ra(npr),3)

laplx=(Prim(x,y,z,npr)/(x**t))*(t*(t-1)*x**(t-2) &
+4*t*(x**(t-1))*Alpha(npr)*(xa-x)+((x**t)*Alpha(npr) &
*(4*Alpha(npr)*(xa-x)-2)))

laply=(Prim(x,y,z,npr)/(y**m))*(m*(m-1)*y**(m-2) &
+4*m*(y**(m-1))*Alpha(npr)*(ya-y)+((y**m)*Alpha(npr) &
*(4*Alpha(npr)*(ya-y)-2)))

laplz=(Prim(x,y,z,npr)/(z**n))*(n*(n-1)*z**(n-2) &
+4*n*(z**(n-1))*Alpha(npr)*(za-z)+((z**n)*Alpha(npr) &
*(4*Alpha(npr)*(za-z)-2)))

sd=laplx+laply+laplz

end function
    
      
!**************************************************************************
   function AO(x,y,z,cao)
   use geninfo
   use loginfo      !cao from 1 to nao, 1->1S, 2->2S, 3->2Px, 4->2Py, ...
   implicit none
   double precision :: AO
   double precision :: Prim
   double precision, intent(in) :: x, y, z
   integer, intent(in) :: cao !choosen ao in the input file
   integer :: smm !sum for all the primitives
   integer :: i
   integer :: sm, kk, k, kp
   smm=0
   AO=0.d0
   sm=0
   k=1
   if (cao.eq.1) then 
     k=1
   else if (cao.eq.2) then
     kp=1
     k=k+npao(1)
   else  
     kp=cao-1    !kp=previous AO 
     do i=1,kp  !set k starting primitive, Flg(k)
      k=k+npao(i)  !sum the primitives until the current AO
     end do
   end if     
   kk=k+npao(cao)-1 !kk=last primitive of the current AO (cao)
    do i=k,kk        
      AO= AO + Flg(i) * Prim(x,y,z,i)    
    end do 
  end function AO
   
!*************************************************************************************   
   function MoOr(x,y,z,mo) !Closed shell molecular orbitals
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr
   double precision :: Prim
   integer, intent(in) :: mo
   integer :: i    
   MoOr=0.d0
   do i=1,nprim
      MoOr= MoOr+ T(mo,i)*Prim(x,y,z,i)
   end do   
   end function

   function MoOr_hf(x,y,z,mo)
   use wfxinfo_hf
   use geninfo
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr_hf
   double precision :: Prim
   integer, intent(in) :: mo
   integer :: i
   MoOr_hf=0.d0
   do i=1,nprim
      MoOr_hf=MoOr_hf+T2(mo,i)*Prim(x,y,z,i)
   end do   
   end function

 !**************************************************************************
  function MO_a(x,y,z,mo) !Apha molecular orbitals
   use wfxinfo
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_a
   double precision :: Prim
   integer :: i
   MO_a=0.d0
   do i=1,nprim
     MO_a= MO_a + T_a(mo,i) * Prim(x,y,z,i)
   end do
   end function
   
   function MO_b(x,y,z,mo) !Beta molecular orbitals
   use wfxinfo
   use geninfo
   implicit none
    
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_b
   double precision :: Prim
   integer :: i
   MO_b=0.d0
   do i=1,nprim
     MO_b=MO_b+ T_b(mo,i)*Prim(x,y,z,i)
   end do
   end function   
 !**************************************************************************
  function MO_aHF(x,y,z,mo) !Apha molecular orbitals
   use wfxinfo_hf
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_aHF
   double precision :: Prim
   integer :: i
   MO_aHF=0.d0
   do i=1,nprim
     MO_aHF= MO_aHF + T_a2(mo,i) * Prim(x,y,z,i)
   end do
   end function
   
   function MO_bHF(x,y,z,mo) !Beta molecular orbitals
   use wfxinfo_hf
   use geninfo
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: mo
   double precision :: MO_bHF
   double precision :: Prim
   integer :: i
   MO_bHF=0.d0
   do i=1,nprim
     MO_bHF=MO_bHF+ T_b2(mo,i)*Prim(x,y,z,i)
   end do
   end function   
 !**********************************************************************
   function Density(x,y,z)   !for closed shell
   use wfxinfo
   use wfxinfo_hf
   use geninfo 
   implicit none
   double precision, intent(in) :: x, y, z
   double precision :: MoOr, Mo_a, Mo_b, MoOr_hf
   double precision :: Density, dens_a, dens_b
   integer :: i
    Density=0.d0
    if (corr) then
      do i=1,noccmo
        Density= Density + Occ(i) * MoOr(x,y,z,i)**2.d0  
      end do
    else 
      do i =1,noccmo2
        Density = Density + Occ2(i) * MoOr_hf(x,y,z,i)**2.d0
      end do
   end if      
    end function
    
    function Dens_hf(x,y,z) !for restricted HF
    use wfxinfo_hf
    use geninfo 
    implicit none
    double precision :: Dens_hf, MoOr_hf
    double precision, intent(in) :: x,y,z
    integer :: i
    Dens_hf=0.d0
    do i=1,noccmo2
       Dens_hf=Dens_hf+Occ2(i)*MoOr_hf(x,y,z,i)**2.d0
    end do
    end function
!*******************************************************************     
function Dens_a(x,y,z) !alpha density (for open-shell)
use wfxinfo
use wfxinfo_hf
use geninfo
implicit none
double precision, intent(in) :: x, y, z
double precision :: MO_a, MoOr, MO_aHF  ! Added MO_aHF for HF
double precision :: Dens_a
integer :: i

Dens_a=0.d0
if (corr) then
  do i=1,nalfaorb  ! Using nalfaorb instead of nalfae
    Dens_a=Dens_a+MO_a(x,y,z,i)**2.d0
  end do
else
  ! Handle HF case
  do i=1,nalfaorb2
    Dens_a=Dens_a+MO_aHF(x,y,z,i)**2.d0
  end do
end if
end function  

!****************************************************************
function Dens_b(x,y,z) !beta density (for open-shell)
use wfxinfo
use wfxinfo_hf
use geninfo
implicit none
double precision, intent(in) :: x, y, z
double precision :: MO_b, MoOr, MO_bHF  ! Added MO_bHF for HF
double precision :: Dens_b
integer :: i

Dens_b=0.d0
if (corr) then
  do i=1,nbetaorb  ! Using nbetaorb instead of nalfae
    Dens_b=Dens_b+MO_b(x,y,z,i)**2.d0
  end do
else
  ! Handle HF case
  do i=1,nbetaorb2
    Dens_b=Dens_b+MO_bHF(x,y,z,i)**2.d0
  end do
end if
end function  
!*****************************************************************

    
function lapl(x,y,z) !computes Laplacian numerically
use wfxinfo
implicit none 
double precision :: lapl                           !this could be a function!!
double precision, intent(in) :: x,y,z
double precision :: spx, spy, spz !laplacian, square partial x, s.p.y, spz
double precision :: density,h
h=0.0000001d0
spx=(Density(x+h,y,z)+Density(x-h,y,z)-2.d0*Density(x,y,z))/(h**2)
spy=(Density(x,y+h,z)+Density(x,y-h,z)-2.d0*Density(x,y,z))/(h**2)
spz=(Density(x,y,z+h)+Density(x,y,z-h)-2.d0*Density(x,y,z))/(h**2)
lapl=spx+spy+spz
end function

!************************************************************************* 

!function Laplacian(x,y,z)  !laplacian analitically using fd and sd functions (first derivative and second derivative of primitives  !no funtziona
!use wfxinfo
!use geninfo
!implicit none
!double precision :: Laplacian, Laplx, Laply, Laplz, dp, d2p
!double precision :: prim
!double precision, intent(in) :: x, y, z  
!integer :: i,j,k
!Laplacian=0.d0                      !consider the UHF case!
!Laplx=0.d0
!Laply=0.d0
!Laplz=0.d0
!do i=1,noccmo    
!          !current molecular orbital value (also valid for UHF)
!   do j=1,nprim
!    do k=1,nprim    
!     Laplx=Laplx+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,1)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,1))**2)   
!     Laply=Laply+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,2)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,2))**2)   
!     Laplz=Laplz+2*(T(i,j)*T(i,k)*d2p(x,y,z,j,3)*Prim(x,y,z,k)+(T(i,j)*dp(x,y,z,j,3))**2)            
!     end do                
!    end do   
!end do
!Laplacian=laplx+laply+laplz
!end function              
!*******************************************************************
subroutine gradient(x,y,z,gradx,grady,gradz) !computes the gradient numerically
                                             !now we won't use that but could be useful
 implicit none
 double precision, intent(in) :: x,y,z
 double precision, intent(out) :: gradx, grady, gradz
 double precision :: density, h

 h=0.0000000000001d0 !change in the derivative

 gradx=(Density(x+h,y,z)-Density(x,y,z))/h
 grady=(Density(x,y+h,z)-Density(x,y,z))/h
 gradz=(Density(x,y,z+h)-Density(x,y,z))/h  
end subroutine
    
!*****************************************************

  function rDM1(x1,y1,z1,x2,y2,z2) !1-electron reduced density matrix (HF case)
                                   !for alpha and beta electrons
   use wfxinfo                     
   use geninfo
   implicit none
   double precision :: rDM1
   double precision :: rDM1_a, rDM1_b
   double precision :: MoOr, MO_a, MO_b
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i  
       if (uhf.or.opsh) then
         if (nelec.eq.1) then   !only one electron 
           rDM1=0.d0
           do i=1,noccmo
                rDM1=rDM1+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
           end do
         else        
           rDM1_a=0.d0
           rDM1_b=0.d0  
           if (dble(nalfae).gt.0.1d0) then  
             do i=1,nalfae
                rDM1_a=rDM1_a+MO_a(x1,y1,z1,i)*MO_a(x2,y2,z2,i)
             end do
           end if
           if (dble(nbetae).gt.0.1d0) then  
              do i=1,nbetae
                rDM1_b=rDM1_b+MO_b(x1,y1,z1,i)*MO_b(x2,y2,z2,i)
              end do
           end if  
           rDM1=rDM1_a+rDM1_b
         end if  
      else   !restricted Hartree Fock
           do i=1,noccmo
             if (occ(i).gt.0.d0) then
                rDM1=rDM1+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
             end if   
           end do
       end if    
   end function

   function rDM1_alf(x1,y1,z1,x2,y2,z2) !only for HF?
   use wfxinfo
   use geninfo
   implicit none
   double precision :: rdm1_alf
   double precision :: MoOr, MO_a, MO_b
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i
    rDM1_alf=0.d0
   if (.not.(opsh)) then !closed shell (You should consider the case of Closed Shell UHF)
          rDM1_alf=0.d0     
          do i=1,noccmo
            rDM1_alf=rDM1_alf+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
          end do
          rDM1_alf=rDM1_alf/2.d0  
         
   else !opsh         
        if (nelec.eq.1) then   !only one electron
           rDM1_alf=0.d0
           do i=1,noccmo
                rDM1_alf=rDM1_alf+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
           end do
        else
           rDM1_alf=0.d0
           if (dble(nalfae).gt.0.1d0) then
             do i=1,nalfae
                rDM1_alf=rDM1_alf+MO_a(x1,y1,z1,i)*MO_a(x2,y2,z2,i)
             end do
           else 
                write(*,*) "No alpha electrons, rdm1alf=0"   
                rDM1_alf=0.d0
           end if
        end if  
   end if  
   end function   
        
   function rDM1_bet(x1,y1,z1,x2,y2,z2)
   use wfxinfo
   use geninfo
   implicit none
   double precision :: rdm1_bet
   double precision :: MoOr, MO_a, MO_b
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i
   if (uhf) then
        if (nelec.eq.1) then   !only one electron
           rDM1_bet=0.d0
           do i=1,noccmo
                rDM1_bet=rDM1_bet+Occ(i)*MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
           end do
        else
           rDM1_bet=0.d0
           if (dble(nalfae).gt.0.1d0) then
             do i=1,nalfae
                rDM1_bet=rDM1_bet+MO_a(x1,y1,z1,i)*MO_a(x2,y2,z2,i)
             end do
           else
                write(*,*) "No alpha electrons, rdm1alf=0"
                rDM1_bet=0.d0
           end if
        end if
    else
        rdm1_bet=0.d0
        do i=1,noccmo
              if (occ(i).gt.0.d0) then
                 rDM1_bet=rDM1_bet+MoOr(x1,y1,z1,i)*MoOr(x2,y2,z2,i)
              end if
        end do
    end if
    end function

    function rDM1_alf_hf(x1,y1,z1,x2,y2,z2)
    use wfxinfo_hf
    use geninfo
   implicit none
   double precision :: rdm1_alf_hf
   double precision :: MoOr_hf
   double precision, intent(in) :: x1, y1, z1, x2, y2, z2
   integer :: i
   rdm1_alf_hf=0.d0
   do i=1,noccmo2
       rdm1_alf_hf=rdm1_alf_hf+Occ2(i)*MoOr_hf(x1,y1,z1,i)*MoOr_hf(x2,y2,z2,i)
   end do
   rdm1_alf_hf=rdm1_alf_hf/2.d0
   end function

   function rho2(x1,y1,z1,x2,y2,z2,dm2p)
       use wfxinfo
       use geninfo
       implicit none
       integer :: i,j,k,l
       double precision :: DMval
       double precision :: rho2
       double precision :: prim
       double precision, intent(in) :: x1,y1,z1,x2,y2,z2
       character, intent(in) :: dm2p
       open(unit=5,file=dm2p, form='unformatted') !open binary file
       rewind(5)
       rho2=0.d0
       do while (.true.)
           !call cpu_time(T2)
           read(5,end=100) i,j,k,l,DMval !read a line from binary file .dm2
           if (i.eq.nprim+1) then
              goto 100 !file is finished
           end if
           rho2=rho2+DMval*Prim(x1,y1,z1,i)*Prim(x2,y2,z2,j) &
                          *Prim(x1,y1,z1,k)*Prim(x2,y2,z2,l)
       end do   
       100 continue
       close(5)
   end function

                    
