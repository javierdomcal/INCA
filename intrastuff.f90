module intrastuff 
use geninfo
use qhermite
implicit none
integer :: i,j,k,l !primitive quartets
!Grid independent variables for the intracule:
double precision :: a_ik, a_jl, a_ijkl 
double precision :: e_ik, e_jl, e_ijkl, sqe 
double precision, dimension(3) :: R_i, R_k, R_j, R_l !X_i, ...
double precision, dimension(3) :: R_ik, R_jl, R_ijkl !X_ik, ...
double precision :: Alf_ijkl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
!!!!!!!!!!!!!Functions for the first integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        function J_ik(i,k,R_i_k_2) !J upper bound for the first integral screening
        implicit none
        double precision :: J_ik       
        integer, intent(in) :: i,k
        double precision, intent(in) :: R_i_k_2
        double precision, dimension(3) :: x_max
        double precision :: r_part
        integer :: n_nodes, ii
    
        do ii=1,3 !loop for J_ik(x_i(*), (y) and (z)
                n_nodes=TMN(i,ii)+TMN(k,ii) + 1          !obtain x_max
                if (n_nodes.eq.1) then
                        x_max(ii)=0.d0
                else if (n_nodes.eq.2) then
                        x_max(ii)=0.707106781186548d0 
                else if (n_nodes.eq.3) then 
                        x_max(ii)=1.224744871391589d0 
                else if (n_nodes.eq.4) then
                        x_max(ii)=1.650680123885785d0
                else if (n_nodes.eq.5) then 
                        x_max(ii)=2.020182870456086d0 
                else if (n_nodes.eq.6) then 
                        x_max(ii)=2.350604973674492d0
                else if (n_nodes.eq.7) then 
                        x_max(ii)=2.651961356835233d0 
                else if (n_nodes.eq.8) then 
                        x_max(ii)=2.930637420257244d0
                else if (n_nodes.eq.9) then 
                        x_max(ii)=3.190993201781528d0  
                else if (n_nodes.eq.10) then 
                        x_max(ii)=3.436159118837738d0
                end if
        end do
                
        J_ik=pi**(1.5d0)*(2.d0*a_ik)**(-(dble(TMN(i,1)+TMN(k,1)+TMN(i,2)+TMN(k,2)+TMN(i,3)+TMN(k,3))+1.5d0))&
        *exp(-2.d0*e_ik*(R_i_k_2))
        
        r_part=1.d0
        do ii=1,3
              r_part=r_part &
              *((x_max(ii)+Alpha(k)*sqrt(2.d0/a_ik)*abs(Cartes(Ra(k),ii)-Cartes(Ra(i),ii)))**(2.d0*dble(TMN(i,ii))))&
              *((x_max(ii)+Alpha(i)*sqrt(2.d0/a_ik)*abs(Cartes(Ra(k),ii)-Cartes(Ra(i),ii)))**(2.d0*dble(TMN(k,ii))))
        end do

        J_ik=J_ik*r_part
        end function
 
        function J_jl(j,l,R_j_l_2) !J upper bound for the first integral screening
        double precision :: J_jl
        integer, intent(in) :: j,l
        double precision, intent(in) :: R_j_l_2
        double precision, dimension(3) :: x_max 
        integer :: n_nodes, ii
        n_nodes=0  
        do ii=1,3 !loop for J_ik(x_i(*), (y) and (z)
                n_nodes=TMN(j,ii)+TMN(l,ii) + 1
                if (n_nodes.eq.1) then 
                        x_max(ii)=0.d0
                else if (n_nodes.eq.2) then
                        x_max(ii)=0.707106781186548d0 
                else if (n_nodes.eq.3) then 
                        x_max(ii)=1.224744871391589d0 
                else if (n_nodes.eq.4) then 
                        x_max(ii)=1.650680123885785d0
                else if (n_nodes.eq.5) then 
                        x_max(ii)=2.020182870456086d0 
                else if (n_nodes.eq.6) then 
                        x_max(ii)=2.350604973674492d0
                else if (n_nodes.eq.7) then 
                        x_max(ii)=2.651961356835233d0 
                else if (n_nodes.eq.8) then 
                        x_max(ii)=2.930637420257244d0
                else if (n_nodes.eq.9) then 
                        x_max(ii)=3.190993201781528d0  
                else if (n_nodes.eq.10) then 
                        x_max(ii)=3.436159118837738d0
                end if
        end do
  
        J_jl=pi**(1.5d0)*(2.d0*a_jl)**(-(dble(TMN(j,1)+TMN(l,1)+TMN(j,2)+TMN(l,2)+TMN(j,3)+TMN(l,3))+1.5d0))&
        *exp(-2.d0*e_jl*(R_j_l_2))&
        *((x_max(1)+Alpha(l)*sqrt(2.d0/a_jl)*abs(Cartes(Ra(l),1)-Cartes(Ra(j),1)))**(2.d0*dble(TMN(j,1))))&
        *((x_max(1)+Alpha(j)*sqrt(2.d0/a_jl)*abs(Cartes(Ra(l),1)-Cartes(Ra(j),1)))**(2.d0*dble(TMN(l,1))))&
        *((x_max(2)+Alpha(l)*sqrt(2.d0/a_jl)*abs(Cartes(Ra(l),2)-Cartes(Ra(j),2)))**(2.d0*dble(TMN(j,2))))&
        *((x_max(2)+Alpha(j)*sqrt(2.d0/a_jl)*abs(Cartes(Ra(l),2)-Cartes(Ra(j),2)))**(2.d0*dble(TMN(l,2))))&
        *((x_max(3)+Alpha(l)*sqrt(2.d0/a_jl)*abs(Cartes(Ra(l),3)-Cartes(Ra(j),3)))**(2.d0*dble(TMN(j,3))))&
        *((x_max(3)+Alpha(j)*sqrt(2.d0/a_jl)*abs(Cartes(Ra(l),3)-Cartes(Ra(j),3)))**(2.d0*dble(TMN(l,3))))
        end function
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

   !     function A_ind(R_i_k_2,R_j_l_2)  !grid independent part of the intracule
   !     implicit none
   !     double precision, intent(in) :: R_i_k_2, R_j_l_2
   !     double precision :: A_ind
   !     A_ind=(a_ijkl)**(-1.5d0)* exp(-e_ik*R_i_k_2-e_jl*R_j_l_2) !eq.18
 ! 
 !       end function A_ind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!eq 15!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!


  function Wr(rhat,r,ax)  !for 1st integral screening
        double precision :: Wr
        double precision, intent(in) :: rhat,r
        integer, intent(in) :: ax
        Wr=((sqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl-0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(i),ax)))**dble(TMN(i,ax))*&
        ((sqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl+0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(j),ax)))**dble(TMN(j,ax))*&
        ((sqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl-0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(k),ax)))**dble(TMN(k,ax))*&
        ((sqrt(a_ijkl)**(-1.d0)*rhat)+(alf_ijkl+0.5d0)*r+(r_ijkl(ax)-Cartes(Ra(l),ax)))**dble(TMN(l,ax))
  end function Wr
      
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!Functions for second integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
 function w_m(m) !eq. A10
 double precision :: w_m
 integer, intent(in) :: m
 if (m.eq.1) then
    w_m=1.d0
 else
    w_m=(dble(m)*0.5d0)**(dble(m)*0.5d0) * exp(-dble(m)*0.5d0)
 end if       
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine polycoef(C_r,np,n,ax)
!Evaluates W function (eqn.15) at xx(np=L+1) points and sets up a system of linear equations (see
!matrix M). 
!Solves this linear equations using Lapack subroutine dgesv.
!As output, we get the polynomial coefficients(C_r(np)) of V_r. This process is repeated for 
!x (ax=1), y (ax=2) and z (ax=3).

!n is the number of nodes of the quadrature
!np is the number of unknowns
!w_n is a matrix with the n weights of the G-Hermite quadrature (in module quadratures)
!rh is a matrix with the n nodes of the G-Hermite quadrature (in module quadratures)
IMPLICIT NONE
!global variables
INTEGER, intent(in) :: np, n, ax !number of points, number of nodes, axis
double precision, allocatable, dimension(:) :: C_r
!local variables
double precision, dimension(np,np) :: M  
double precision, dimension(np) :: xx    !Points where we evaluate W
integer :: Ltot
integer :: i1, j1, l1

!variables for dgbsv!!!!!
integer :: INFO
integer, allocatable, dimension(:) :: ipiv

allocate(C_r(np))

if (np.eq.1) then !only one coefficient
   C_r=w_r(1)     !see equation 15 (in this case we only have 1 node)
else
   Ltot=np-1 !degree of the polynomial    
   !generate 'np' points to evaluate the polynomial (xx(np)) !linearly independent!! 
   do i1=1,np
     xx(i1)=-0.5d0+dble(i1)
   end do   
   !build M matrix (M and C_r form the augmented matrix)
   do i1=1,np
       l1=Ltot !start from the maximum degree
       do j1=1,np
          M(i1,j1)=(sqe*(xx(i1)+r_ik(ax)-r_jl(ax)))**(dble(l1)) !r'^l               
          l1=l1-1
       end do
   end do
   C_r=0.d0
   do i1=1,np
       do j1=1,n        
         C_r(i1)=C_r(i1)+w_r(j1)*Wr(rh(j1),xx(i1),ax) !V_ijkl(X') (~eq.16)
       end do 
   end do 
   !!!Solve linear system with lapack!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate(ipiv(np))
   call dgesv(np,1,M,np,ipiv,C_r,np,INFO)
   if (INFO.ne.0) write(*,*) "ERROR, cannot solve linear equations"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if
END SUBROUTINE Polycoef
      
end module intrastuff
