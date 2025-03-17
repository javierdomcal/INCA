subroutine intracule(dmname,dm2norm)  !Computes vector or radial intracule, or radial integration
                                                            !need .wfx and .dm2p as input, and a grid
                                                            !as output gives a vector with intracule values at each grid point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use geninfo !information about the primitive functions
use intrastuff !subroutines and functions to compute the intracule 
use coeff !C_x,C_y,C_z
use total_integral
use radial_scan
use thresholds
implicit none
!!!!!!!!!!!!!!!!!INPUT file names and normalization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character*40, intent(in) :: dmname !name of the .dm2p file
logical, intent(in) :: dm2norm     !is the dm2p value normalized 
!!!!!!!!!!!!!!!local variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: dfact    !double factorial
integer :: ii,iii,sum, ig, ir, summ, np, sm, smm !integers
double precision :: R_i_k_2, R_j_l_2  !R_ik=(R_i-R_k)^2
double precision :: Aa_ijkl !grid independent part of the intracule
double precision :: n_prim_i, n_prim_j,n_prim_k, n_prim_l !normalization of primitives for DM2
double precision :: n_prim_t !total normalization factor
double precision :: DMval !value of DM2
double precision :: A_ind !eq.18 (grid independent part)
double precision :: screen1, screen2, lim !integral screenings
double precision :: U_x, U_y, U_z !coef. for 2nd integral screening
!Gauss-hermite quadrature
integer :: nn !number of gauss hermite nodes (for x y and z)
!GRID POINTS
!double precision :: r_integral                         !radial integral
!double precision, allocatable, dimension(:) :: r_intra !radial intracule
!check accuracy of calculations
double precision :: trace_DM2prim, trDM2 !normalized and not normalized DM2prim
integer :: npairs      !number of electron pairs
real :: T1,T2
real :: Tread, T1screen, T2screen, Tgrid
 call cpu_time(T1)
 lim=thresh*(dble(nprim)*(dble(nprim)+1.d0)*0.5d0)**(-1.d0) !limit for the 1st integral screening 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 trace_DM2prim=0.d0 
 trDM2=0.d0
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 allocate(I_vec1(np1)); allocate(I_vec2(np2))
 I_vec1=0.d0; I_vec2=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!Start loop over primitive quartets!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!Following Cioslowski and Liu algorithm!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 open(unit=5,file=dmname, form='unformatted') !open binary file
 rewind(5)
 do while (.true.) 
          !call cpu_time(T2)
          read(5,end=100) i,j,k,l,DMval !read a line from binary file .dm2
          if (i.eq.nprim+1) then
                  write(*,*) i,j,k,l,DMval
                  goto 100 !file is finished
          end if        

          if (i.eq.0) goto 100 !file is finished   
          if (abs(DMval).gt.1d-16) then           
           !trDM2=trDM2+DMval
           smm=smm+1 !count primitive quartet
           if (dm2norm) then !dm2 comes from fchk so does not include normalization of prim
              !normalization of DM2--> Normalize the primitives
              N_prim_i=(2.d0*Alpha(i)/pi)**(0.75d0)&
              *sqrt(((4.d0*Alpha(i))**(dble(TMN(i,1)+TMN(i,2)+TMN(i,3))))*&
              dble(dfact(2*TMN(i,1)-1)*dfact(2*TMN(i,2)-1)*dfact(2*TMN(i,3)-1))**(-1.d0))  
   
              N_prim_j=(2.d0*Alpha(j)/pi)**(0.75d0)&
              *sqrt(((4.d0*Alpha(j))**(dble(TMN(j,1)+TMN(j,2)+TMN(j,3))))*&
              dble(dfact(2*TMN(j,1)-1)*dfact(2*TMN(j,2)-1)*dfact(2*TMN(j,3)-1))**(-1.d0)) 
   
              N_prim_k=(2.d0*Alpha(k)/pi)**(0.75d0)&
              *sqrt(((4.d0*Alpha(k))**(dble(TMN(k,1)+TMN(k,2)+TMN(k,3))))*&
              dble(dfact(2*TMN(k,1)-1)*dfact(2*TMN(k,2)-1)*dfact(2*TMN(k,3)-1))**(-1.d0)) 
         
              N_prim_l=(2.d0*Alpha(l)/pi)**(0.75d0)&
              *sqrt(((4.d0*Alpha(l))**(dble(TMN(l,1)+TMN(l,2)+TMN(l,3))))*&
              dble(dfact(2*TMN(l,1)-1)*dfact(2*TMN(l,2)-1)*dfact(2*TMN(l,3)-1))**(-1.d0)) 
  
              !compute DMval with the normalization of primitives
              n_prim_t=N_prim_i*N_prim_j*N_prim_k*N_prim_l
              DMval=n_prim_t*DMval 
              trace_DM2prim=trace_DM2prim+DMval !sum all the DM2 quartets to check accuracy
           end if
           !compute the first variables
           a_ik=Alpha(i)+Alpha(k)  
           a_jl=Alpha(j)+Alpha(l)                  !eqn. 10                         
           e_ik=Alpha(i)*Alpha(k)*a_ik**(-1.d0)
           e_jl=Alpha(j)*Alpha(l)*a_jl**(-1.d0)   

           R_i_k_2=(Cartes(Ra(i),1)-Cartes(Ra(k),1))**2.d0+& !this is needed to compute A_ijkl (eqn. 18)
                   (Cartes(Ra(i),2)-Cartes(Ra(k),2))**2.d0+&      !R_i_k_2=(R_i-R_k)²
                   (Cartes(Ra(i),3)-Cartes(Ra(k),3))**2.d0    
        

           R_j_l_2=(Cartes(Ra(j),1)-Cartes(Ra(l),1))**2.d0+&
                   (Cartes(Ra(j),2)-Cartes(Ra(l),2))**2.d0+&
                   (Cartes(Ra(j),3)-Cartes(Ra(l),3))**2.d0
                    
          !!!!1st integral screening!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                   
          screen1=abs(DMval)*sqrt(J_ik(i,k,R_i_K_2)*J_jl(j,l,R_j_l_2))     
               
          !call cpu_time(TT3)
          !T1screen=T1screen-(TT3-TT2)         
          if (screen1.ge.lim) then    
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                 !compute the other variables (eqn. 12)
                 a_ijkl=a_ik+a_jl
                 e_ijkl=(a_ik*a_jl)*a_ijkl**(-1.d0)
                 sqe=sqrt(e_ijkl)
                 do ii=1,3                   !X_ik,Y_ik,Z_ik,...
                       R_ik(ii)=(Alpha(i)*Cartes(Ra(i),ii)+Alpha(k)*Cartes(Ra(k),ii))*a_ik**(-1.d0)
                       R_jl(ii)=(Alpha(j)*Cartes(Ra(j),ii)+Alpha(l)*Cartes(Ra(l),ii))*a_jl**(-1.d0)
                       R_ijkl(ii)=(a_ik*R_ik(ii)+a_jl*R_jl(ii))*a_ijkl**(-1.d0)
                 end do                      
                 Alf_ijkl=0.5d0*(a_ik-a_jl)*a_ijkl**(-1.d0)
                 !compute Aa_ijkl (the grid independent part)  !eq.18
                 A_ind=(a_ijkl)**(-1.5d0)* exp(-e_ik*R_i_k_2-e_jl*R_j_l_2)
                 Aa_ijkl=DMval*A_ind 
                 !!!!!!!!!!!!Calculate coeficients of V!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 sm=0
                 U_x=0.d0; U_y=0.d0; U_z=0.d0
                 do ii=1,3
                       Lrtot(ii)=TMN(i,ii)+TMN(j,ii)+TMN(k,ii)+TMN(l,ii) !Lrtot=degree of eqn 15 (I don't use Lmax)
                       call gauherm(Lrtot(ii), nn) !obtain Gauss-Hermite nodes(rh) and weights(w_r) (2L+1)   
                       np=Lrtot(ii)+1              !np is the number of coefficients to represent the polynomial V
                       if (ii.eq.1) call polycoef(C_x, np, nn, ii) 
                       if (ii.eq.2) call polycoef(C_y, np, nn, ii)
                       if (ii.eq.3) call polycoef(C_z, np, nn, ii)
                       deallocate(rh)
                       deallocate(w_r)                            
                       !evaluates Vr at Lxtot+1=np points to create the augmented matrix M
                       !diagonalizes M to obtain the coeficients of V --> C 
                        do iii=1,np
                            if (ii.eq.1) U_x=U_x+C_x(iii)*w_m(iii)  !compute U coefficients
                            if (ii.eq.2) U_y=U_y+C_y(iii)*w_m(iii)  !for second integral screening
                            if (ii.eq.3) U_z=U_z+C_z(iii)*w_m(iii)  !eq.43 of the paper  
                        end do
                 end do                                 
                 !!!!!!!!!!2nd Integral Screening!!!!!!!!!!!!!!!!!!!!!                      
                 screen2=abs(Aa_ijkl*U_x*U_y*U_z) !eqn.40
                 !call CPU_time(TT4)
                 !T2screen=T2screen+(TT4-TT3)                
                 if (screen2.ge.lim) then 
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                             
                           !loop over grid points 
                            call intra_eval(np2,r2,Aa_ijkl,I_vec2) !for radial scan
                           
                            call intra_eval(np1,r1,Aa_ijkl,I_vec1) !for total integral
                            
                 else                     
                              summ=summ+1  !count quartets that do not pass 2nd screening 
                 end if 
                 deallocate(C_x)
                 deallocate(C_y)   !deallocate polynomial coefficients
                 deallocate(C_z)                                                                
       else               
               sum=sum+1    !count quartets that do not pass 1st screening                                 
       end if   
      end if   !dmval lt 10-16    
 end do !end loop over quartets 
 100 continue           !it comes here when .dm2 file is finished
 close(5)
 write(*,*) "Ended loop for primitives"
 call CPU_time(T2)
 write(*,*) "Calculation time", T2-T1
 
end subroutine


subroutine intra_eval(ngrid,r,Aa_ijkl,I_vec)
use intrastuff !r_ik,r_jl,sqe
use coeff !C_x,C_y,C_z,Lrtot!Already allocated
implicit none
!Global variables
integer, intent(in) :: ngrid
double precision, intent(in), dimension(3,ngrid) :: r
double precision, intent(in) :: Aa_ijkl
double precision, intent(inout), dimension(ngrid) :: I_vec
!Local variables
double precision, dimension(3) :: rp
double precision :: pot
double precision :: V_x, V_y, V_z !eq 16
integer :: ig, ii, iii
  do ig=1,nGrid
       V_x=0.d0; V_y=0.d0; V_z=0.d0
       do ii=1,3                                       !loop for x y and z  
            rp(ii)=sqe*(r(ii,ig)+r_ik(ii)-r_jl(ii))   !compute R'(eq. 19), 
                                                       !using total grid points
            pot=dble(Lrtot(ii))
            !!!!!!!Compute Vx, Vy and Vz!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do iii=1,Lrtot(ii)+1   !Ltot+1 is the number of coefficients we have
                    if (pot.gt.1d-16) then
                         if (ii.eq.1) V_x=V_x+C_x(iii)*rp(1)**(pot)
                         if (ii.eq.2) V_y=V_y+C_y(iii)*rp(2)**(pot)
                         if (ii.eq.3) V_z=V_z+C_z(iii)*rp(3)**(pot)
                         pot=pot-1.d0
                    else
                         if (ii.eq.1) V_x=V_x+C_x(iii) !when pot is 0
                         if (ii.eq.2) V_y=V_y+C_y(iii)
                         if (ii.eq.3) V_z=V_z+C_z(iii)
                    end if
            end do  !end loop over the polynomial coefficients
       end do       !end loop over x,y,z
       !!!!!!!!!!!!!Calculate intracule at a point!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       I_vec(ig)=I_vec(ig)+Aa_ijkl*exp(-(rp(1)**2.d0+rp(2)**2.d0+rp(3)**2.d0))*V_x*V_y*V_z
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function dfact(a) !computes double factorial
integer :: dfact
integer, intent(in) :: a
dfact=1
if (mod(a,2).eq.0) then !n is even
  do i=1,a/2
      dfact=dfact*(2*i)
  end do
else 
  do i=1,(a+1)/2
      dfact=dfact*(2*i-1)
  end do
end if
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dfact_dp(a) result(res)
    integer, intent(in) :: a
    double precision :: res
    integer :: i

    ! Handle negative input
    if (a < 0) then
        write(*, '(A,I4)') 'WARNING: dfact_dp called with negative input:', a
        res = 1.0d0  ! or whatever makes sense in your physics context
        return
    endif

    res = 1.0d0
    if (mod(a,2) == 0) then !n is even
        do i = 1, a/2
            res = res * (2.0d0*i)
        end do
    else 
        do i = 1, (a+1)/2
            res = res * (2.0d0*i-1.0d0)
        end do
    end if
end function
