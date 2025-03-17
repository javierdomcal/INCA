!This bunch of subroutines are used for Becke-Roussel stuff
subroutine dens_ops(x,y,z,Dens,E_kin,Grad,Lapl)
!Compute Density, kinetic energy density, gradient, and laplacian
!at x,y,z point
use geninfo
use wfxinfo
implicit none
!!!global variables!!!!!!!!!!!!!!!!!!
double precision, intent(in) :: x,y,z
double precision, intent(out) :: Dens,E_kin,Grad,Lapl
!!!local variables!!!!!!!!!!!!!!!!!!!!!!
!primitive derivatives (matrices)
double precision, dimension(nprim) :: pval !value of primitive
double precision, dimension(nprim) ::  dpxval, dpyval, dpzval !value of derivative of primitive
double precision, dimension(nprim) :: fxval,fyval,fzval !value of fragment of derivative of prim.
double precision, dimension(nprim) :: dpxx,dpyy,dpzz !value of double derivative of primitive
!primitive derivatives (on the fly values)
double precision :: dpx,dpy,dpz !fragment of derivative of primitives
double precision :: dx,dy,dz !derivative of primitives
double precision :: dxx,dyy,dzz !double derivative of primitive
double precision :: xa,ya,za    !center of primitives
double precision :: pr !primitive
!functions
double precision :: density, prim
double precision :: E_kn
double precision :: gradx,grady,gradz
double precision :: laplx,laply,laplz, laplx1, laplx2, laply1, laply2, laplz1, laplz2
double precision :: Ek_x,Ek_y,Ek_z   !x,y,z components of kinetic energy density
double precision :: Ekin_x,Ekin_y,Ekin_z
double precision :: Tmult, Tmult2
integer :: i,j,k,l,m,i2
double precision :: grad2, gx1,gx2,gy1,gy2,gz1,gz2
Dens=Density(x,y,z) !compute 1spin density
do i=1,nprim !compute and store value of primitives and its derivatives
    xa=cartes(Ra(i),1); ya=cartes(Ra(i),2); za=cartes(Ra(i),3)    
    pr=Prim(x,y,z,i)
    pval(i)=pr !value of primitive i   
    call dprim(i,x,y,z,xa,ya,za,pr,dpx,dpy,dpz,dx,dy,dz) 
    fxval(i)=dpx; fyval(i)=dpy; fzval(i)=dpz    
    dpxval(i)=dx; dpyval(i)=dy; dpzval(i)=dz   
    call ddprim(i,x,y,z,xa,ya,za,pval(i),dx,dy,dz,dpx,dpy,dpz,dxx,dyy,dzz)
    dpxx(i)=dxx; dpyy(i)=dyy; dpzz(i)=dzz 
end do
E_kin=0.d0; Ekin_x=0.d0; Ekin_y=0.d0; Ekin_z=0.d0
Gradx=0.d0; Grady=0.d0; Gradz=0.d0; Grad2=0.d0
Lapl=0.d0; Laplx=0.d0; Laply=0.d0; Laplz=0.d0
do i=1,noccmo
   if (Occ(i).gt.0.d0) then
    do j=1,nprim
      do k=1,nprim
         Tmult=T(i,j)*T(i,k)*Occ(i)
         Gradx=Gradx+Tmult*(dpxval(j)*pval(k)+pval(j)*dpxval(k))
         Grady=Grady+Tmult*(dpyval(j)*pval(k)+pval(j)*dpyval(k))
         Gradz=Gradz+Tmult*(dpzval(j)*pval(k)+pval(j)*dpzval(k))
         Laplx=Laplx+Tmult*((2.d0*dpxval(j)*dpxval(k))+pval(j)*dpxx(k) &
                    +dpxx(j)*pval(k))
         Laply=Laply+Tmult*((2.d0*dpyval(j)*dpyval(k))+pval(j)*dpyy(k) &
                    +dpyy(j)*pval(k))
         Laplz=Laplz+Tmult*((2.d0*dpzval(j)*dpzval(k))+pval(j)*dpzz(k) &
                    +dpzz(j)*pval(k))
         Ekin_x=Ekin_x+Tmult*dpxval(j)*dpxval(k)
         Ekin_y=Ekin_y+Tmult*dpyval(j)*dpyval(k)
         Ekin_z=Ekin_z+Tmult*dpzval(j)*dpzval(k)
     end do
   end do
  end if 
end do 
Grad=(gradx/2.d0)**(2.d0)+(grady/2.d0)**(2.d0)+(gradz/2.d0)**(2.d0) !scalar product of gradient (one spin)
Lapl=Laplx+Laply+Laplz
E_kin=(Ekin_x+Ekin_y+Ekin_z)
!compute all one spin
if (.not.(opsh)) then
 !       write(*,*) "Closed shell, Alpha=Beta"
        Dens=Dens/2.d0
       !Grad=Grad/2.d0
        Lapl=Lapl/2.d0
        E_kin=E_kin/2.d0
else
        write(*,*) "Caution: Open shell"
end if        
end subroutine

subroutine dprim(mu,x,y,z,xa,ya,za,pr,dpx,dpy,dpz,dx,dy,dz) 
use geninfo
implicit none
double precision, intent(in) :: x,y,z,xa,ya,za,pr
double precision, intent(out) :: dpx, dpy, dpz,dx,dy,dz
integer, intent(in) :: mu
double precision :: fx1,fy1,fz1
double precision :: expr
dpx=0.d0
dpy=0.d0
dpz=0.d0
if (abs(x-xa).gt.0d-5) then
        fx1=tmn(mu,1)*(x-xa)**(-1.d0) !avoid div 0
else
        fx1=0.d0
        if (tmn(mu,1).eq.1) fx1=1.d0         !see the development in the notebook
end if        
if (abs(y-ya).gt.0d-5) then
        fy1=tmn(mu,2)*(y-ya)**(-1.d0)
else
        fy1=0.d0
        if (tmn(mu,2).eq.1) fy1=1.d0
end if
if (abs(z-za).gt.0d-5) then
        fz1=tmn(mu,3)*(z-za)**(-1.d0)
else
        fz1=0.d0
        if (tmn(mu,3).eq.1) fz1=1.d0
end if  
 dpx=fx1-2.d0*Alpha(mu)*(x-xa)
 dpy=fy1-2.d0*Alpha(mu)*(y-ya)
 dpz=fz1-2.d0*Alpha(mu)*(z-za)

!compute derivative of the primitive 
if (((abs(x-xa).eq.0.d0)).and.(tmn(mu,1).eq.1)) then !x=xA and n=1
        dx=dpx*(y-ya)**(tmn(mu,2))*(z-za)**(tmn(mu,3))*expr(x,y,z,xa,ya,za,mu)
else
        dx=dpx*pr
end if
if (((abs(y-ya).eq.0.d0)).and.(tmn(mu,2).eq.1)) then !y=yA and m=1
        dy=dpy*(x-xa)**(tmn(mu,1))*(z-za)**(tmn(mu,3))*expr(x,y,z,xa,ya,za,mu)
else
        dy=dpy*pr  
end if  
if (((abs(z-za).eq.0.d0)).and.(tmn(mu,3).eq.1)) then !z=zA and l=1
        dz=dpz*(x-xa)**(tmn(mu,1))*(y-ya)**(tmn(mu,2))*expr(x,y,z,xa,ya,za,mu)
else
        dz=dpz*pr
end if

end subroutine dprim

subroutine ddprim(mu,x,y,z,xa,ya,za,pval,dx,dy,dz,fxval,fyval,fzval,dxx,dyy,dzz)
use geninfo  
implicit none
double precision, intent(in) :: x,y,z,xa,ya,za
double precision, intent(in) :: pval !value of the primitive
double precision, intent(in) :: dx,dy,dz !value of x y and z derivatives of prim
double precision, intent(in) :: fxval,fyval,fzval !part of the derivative of prim
integer, intent(in) :: mu
double precision, intent(out) :: dxx,dyy,dzz
double precision :: px, py, pz
double precision :: expr
if (abs(x-xa).gt.0.d0) then
        px=-pval*TMN(mu,1)*(x-xa)**(-2.d0) !avoid div 0
        dxx=dx*fxval+pval*(-2.d0*Alpha(mu))+px
else
        if (TMN(mu,1).eq.1) then
           !  dxx=-2.d0*Alpha(mu)*dx
              dxx=0.d0
        else if (TMN(mu,1).eq.2) then
                dxx=2.d0*(y-ya)**(tmn(mu,2))*(z-za)**(tmn(mu,3))*expr(x,y,z,xa,ya,za,mu)
        else if (TMN(mu,1).eq.0) then 
                dxx=pval*(-2.d0*Alpha(mu))  
           
        else
                dxx=0.d0        
        end if
end if
if (abs(y-ya).gt.0.d0) then
        py=-pval*TMN(mu,2)*(y-ya)**(-2.d0)
        dyy=dy*fyval+pval*(-2.d0*Alpha(mu))+py
else
        if (TMN(mu,2).eq.1) then
             !   dyy=-2.d0*Alpha(mu)*dy
                dyy=0.d0
        else if (TMN(mu,2).eq.2) then
                dyy=2.d0*(x-xa)**(tmn(mu,1))*(z-za)**(tmn(mu,3))*expr(x,y,z,xa,ya,za,mu)
               
 
        else if (TMN(mu,2).eq.0) then
                dyy=pval*(-2.d0*Alpha(mu))     
            
        else
               ! dyy=pval*(-2.d0*Alpha(mu))     
                dyy=0.d0
        end if
end if
if (abs(z-za).gt.0.d0) then
        pz=-pval*TMN(mu,3)*(z-za)**(-2.d0)
        dzz=dz*fzval+pval*(-2.d0*Alpha(mu))+pz
else
        if (TMN(mu,3).eq.1) then
               write(*,*) "dzz not zero", TMN(mu,1),TMN(mu,2),TMN(mu,3)
              ! dzz=-2.d0*Alpha(mu)*dz
              ! write(*,*) dzz,pval
               dzz=0.d0
        else if (TMN(mu,3).eq.2) then
               dzz=2.d0*(x-xa)**(tmn(mu,1))*(y-ya)**(tmn(mu,2))*expr(x,y,z,xa,ya,za,mu)      
        else if (TMN(mu,3).eq.0) then
               dzz=pval*(-2.d0*Alpha(mu))
              
        else
              ! dzz=pval*(-2.d0*Alpha(mu))
               dzz=0.d0
        end if
end if        
end subroutine

