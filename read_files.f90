!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filewfx(wfxfilename)                                             !
!  1-Reads wfxfile and stores the important values in wfxinfo module        !
!    (nelec, noccmo, nprim, ...)                                            !
!  2-Determines the wavefunction type we have: RHF, UHF, correlated WF      !
!  +relaxed/unrelaxed density                                               !
!  3-Builds TMN matrix (matrix with primitive exponents)                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use wfxinfo
use geninfo
use located
implicit none
character*40, intent(in) :: wfxfilename
character*80 :: line
character*80 :: zaborra  !coses que no volem llegir
integer :: i,j,k, kk
integer :: mon !molecular orbital number   
double precision :: maxim, minim
double precision, parameter :: trsh=1.d-6 !threeshold for orbital occupancies
integer :: mult
double precision, allocatable, dimension(:) :: ener
character*40 :: spin
integer, allocatable, dimension(:) :: omg
integer :: sma, smb
 !system
 opsh=.false.
 !method
 corr=.false.  
 uhf=.false.
 !relaxed/unrelaxed density
 udens=.false.
   
open(unit=1,file=wfxfilename,status='OLD') 

call locate(1,"Number of Nuclei") 
read(1,'(a80)')line ! llegeix el que hi ha a la linia que hi ha despres de la linia que conte "Number of Nuclei"
read(line(1:4),*) natoms ! del contingut de "line", guarda el que hi ha a les 3 primeres posicions en una variable
                        
rewind 1  

!!!!!ONLY IF DM2PWRITE IS GOING TO BE USED!!!!!!!!!!!!!!!!!!!!
!DETERMINE OPEN/CLOSED-SHELL
call locate(1,"Number of Occupied Molecular Orbitals")
read(1,'(a80)')line
read(line(1:5),*) noccmo

call locate(1,"Electronic Spin Multiplicity")
read(1,*) mult
 
rewind 1 

if (mult.eq.1) then
    write(*,*) "Multiplicity is 1-->", "CLOSED-SHELL"
else
    write(*,*) "OPEN-SHELL"    
    opsh=.true.
    call locate(1,"Number of Alpha Electrons")
    read(1,*) nalfae
    call locate(1,"Number of Beta Electrons")
    read(1,*) nbetae
    allocate(omg(noccmo))    
end if



!!!!!ONLY IF DM2PWRITE IS GOING TO BE USED!!!!!!!!!!!!!!!!!!!!
!DETERMINE CORRELATED-UNCORRELATED METHOD 
allocate(ener(noccmo))
call locate(1,"Molecular Orbital Energies")
do i=1,noccmo
   read(1,*) ener(i)
   write(*,*) "Ener=", ener(i)
end do
if (all(ener.eq.0.d0)) then
   corr=.true.
   write(*,*) "0 energy orbitals-->CORRELATED WAVEFUNCTION"
   deallocate(ener)
   if (opsh) then
   !ASSIGN ALPHA AND BETA ORBITALS
      call locate(1,"Molecular Orbital Spin Types")
      nalfaorb=0; nbetaorb=0
      do i=1,noccmo
          read(1,*) spin
          if (spin.eq."Alpha") then
             omg(i)=1
             nalfaorb=nalfaorb+1
          else
             omg(i)=2 
             nbetaorb=nbetaorb+1
          end if       
      end do
   end if   
end if


rewind 1

allocate(Occ(noccmo))
call locate(1,"Molecular Orbital Occupation Numbers") !check RHF,UHF, or Correlated
do i=1,noccmo                                         !by     
   read(1,*) Occ(i)                                   !reading the MO occupation values
end do

if (.not.corr) then 
!DETERMINE RESTRICTED-UNRESTRICTED HF
    if (minval(Occ).eq.2.d0) then
        write(*,*) "RESTRICTED HARTREE-FOCK"
    else if (minval(Occ).eq.1.d0) then    
        write(*,*) "UNRESTRICTED HARTREE-FOCK"    
        uhf=.true.
    else
       write(*,*) "ERROR", "Uncorrelated wavefunction with fractional occupations"    
    end if
end if
  
!GENERAL STUFF FOR OPEN OR CLOSED-SHELL  
call locate(1,"Net Charge")           
read(1,*) netch

rewind 1

call locate(1,"Number of Electrons")
read(1,'(a80)')line
read(line(1:5),*) nelec


rewind 1

call locate(1,"Number of Primitives")
read(1,'(a80)')line
read(line(1:5),*) nprim

rewind 1

allocate(cartes(natoms,3))
call locate(1,"Nuclear Cartesian Coordinates")
read(1,*) ((cartes(i,j) , j=1,3 ), i=1,natoms)

rewind 1

allocate(Ra(nprim))
call locate(1,"Primitive Centers")
read(1,*) (Ra(i), i=1,nprim)

rewind 1

allocate(Ptyp(nprim))
call locate(1,"Primitive Types")
read(1,*) (Ptyp(i), i=1,nprim)

rewind 1

allocate(Alpha(nprim))
call locate(1,"Primitive Exponents")
read(1,*) (Alpha(i), i=1,nprim)

rewind 1

allocate(an(natoms))
call locate(1,"Atomic Numbers")
do i=1,natoms
 read(1,*) an(i)
end do

allocate(chrg(natoms))
call locate(1,"Nuclear Charges")
do i=1,natoms
  read(1,*) chrg(i)
end do  

allocate(T(noccmo,nprim))                  !read T matrix (MOs in prims)       
call locate(1,"Molecular Orbital Primitive Coefficients")
do kk=1,noccmo
  read(1,'(a80)')line
  read(1,*) mon                          !molecular orbital number "mo"
 read(1,'(a80)') zaborra  
 read(1,*) (T(mon,i), i=1,nprim)
end do

rewind 1

allocate(TMN(nprim,3))      !create a matrix with t m and n exponents for each primitive (see prim equation)
do i=1,nprim                                                                 !t->lx
   if (ptyp(i).eq.1) then   !s prim                                          !m->ly
     do j=1,3                 						!n->lz
       TMN(i,j)=0
     end do
   else if (ptyp(i).eq.2) then  !px
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.3) then  !py
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=0   
   else if (ptyp(i).eq.4) then  !pz
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=1   
   else if (ptyp(i).eq.5) then  !dx²
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=0   
   else if (ptyp(i).eq.6) then  !dy²
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=0  
   else if (ptyp(i).eq.7) then  !dz²
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=2   
   else if (ptyp(i).eq.8) then  !dxy
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=0   
   else if (ptyp(i).eq.9) then  !dxz
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.10) then !dyz 
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=1      
   else if (ptyp(i).eq.11) then !fxxx 
      TMN(i,1)=3
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.12) then !fyyy 
      TMN(i,1)=0
      TMN(i,2)=3
      TMN(i,3)=0
   else if (ptyp(i).eq.13) then !fzzz 
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=3
   else if (ptyp(i).eq.14) then !fxxy 
      TMN(i,1)=2
      TMN(i,2)=1
      TMN(i,3)=0
   else if (ptyp(i).eq.15) then !fxxz 
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.16) then !fyyz 
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=1
   else if (ptyp(i).eq.17) then !fxyy 
      TMN(i,1)=1
      TMN(i,2)=2
      TMN(i,3)=0
   else if (ptyp(i).eq.18) then !fxzz 
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=2
   else if (ptyp(i).eq.19) then !fyzz 
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=2
   else if (ptyp(i).eq.20) then !fxyz 
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=1
   else if (ptyp(i).eq.21) then !gxxxx
      TMN(i,1)=4
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.22) then !gxxxy
      TMN(i,1)=3
      TMN(i,2)=1
      TMN(i,3)=0
   else if (ptyp(i).eq.23) then !gxxxz
      TMN(i,1)=3
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.24) then !gxxyy
      TMN(i,1)=2
      TMN(i,2)=2
      TMN(i,3)=0
   else if (ptyp(i).eq.25) then !gxxyz
      TMN(i,1)=2
      TMN(i,2)=1
      TMN(i,3)=1
   else if (ptyp(i).eq.26) then !gxxzz
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=2
   else if (ptyp(i).eq.27) then !gxyyy
      TMN(i,1)=1
      TMN(i,2)=3
      TMN(i,3)=0
   else if (ptyp(i).eq.28) then !gxyyz
      TMN(i,1)=1
      TMN(i,2)=2
      TMN(i,3)=1
   else if (ptyp(i).eq.29) then !gxyzz
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=2
   else if (ptyp(i).eq.30) then !gxzzz
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=3
   else if (ptyp(i).eq.31) then !gyyyy
      TMN(i,1)=0
      TMN(i,2)=4
      TMN(i,3)=0
   else if (ptyp(i).eq.32) then !gyyyz
      TMN(i,1)=0
      TMN(i,2)=3
      TMN(i,3)=1
   else if (ptyp(i).eq.33) then !gyyzz
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=2
   else if (ptyp(i).eq.34) then !gyzzz
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=3
   else if (ptyp(i).eq.35) then !gzzzz
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=4
   else if (ptyp(i).eq.36) then !hxxxxx
      TMN(i,1)=5
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.37) then !hxxxxy
      TMN(i,1)=4
      TMN(i,2)=1
      TMN(i,3)=0
   else if (ptyp(i).eq.38) then !hxxxxz
      TMN(i,1)=4
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.39) then !hxxxyy
      TMN(i,1)=3
      TMN(i,2)=2
      TMN(i,3)=0
   else if (ptyp(i).eq.40) then !hxxxyz
      TMN(i,1)=3
      TMN(i,2)=1
      TMN(i,3)=1
   else if (ptyp(i).eq.41) then !hxxxzz
      TMN(i,1)=3
      TMN(i,2)=0
      TMN(i,3)=2
   else if (ptyp(i).eq.42) then !hxxyyy
      TMN(i,1)=2
      TMN(i,2)=3
      TMN(i,3)=0
   else if (ptyp(i).eq.43) then !hxxyyz
      TMN(i,1)=2
      TMN(i,2)=2
      TMN(i,3)=1
   else if (ptyp(i).eq.44) then !hxxyzz
      TMN(i,1)=2
      TMN(i,2)=1
      TMN(i,3)=2
   else if (ptyp(i).eq.45) then !hxxzzz
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=3
   else if (ptyp(i).eq.46) then !hxyyyy
      TMN(i,1)=1
      TMN(i,2)=4
      TMN(i,3)=0
   else if (ptyp(i).eq.47) then !hxyyyz
      TMN(i,1)=1
      TMN(i,2)=3
      TMN(i,3)=1
   else if (ptyp(i).eq.48) then !hxyyzz
      TMN(i,1)=1
      TMN(i,2)=2
      TMN(i,3)=2
   else if (ptyp(i).eq.49) then !hxyzzz
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=3
   else if (ptyp(i).eq.50) then !hxzzzz
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=4
   else if (ptyp(i).eq.51) then !hyyyyy
      TMN(i,1)=0
      TMN(i,2)=5
      TMN(i,3)=0
   else if (ptyp(i).eq.52) then !hyyyyz
      TMN(i,1)=0
      TMN(i,2)=4
      TMN(i,3)=1
   else if (ptyp(i).eq.53) then !hyyyzz
      TMN(i,1)=0
      TMN(i,2)=3
      TMN(i,3)=2
   else if (ptyp(i).eq.54) then !hyyzzz
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=3
   else if (ptyp(i).eq.55) then !hyzzzz
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=4
   else if (ptyp(i).eq.56) then !hzzzzz
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=5
   else
     write(*,*) "i orbitals are not implemented"
     stop
   end if
end do

if (opsh) then                !split T matrix into T_a and T_b (alpha and beta)
     write(*,*) "OPEN-SHELL SYSTEM, SPLITTING T MATRIX INTO ALPHA AND BETA"
     allocate(T_a(nalfaorb,nprim)); allocate(occ_a(nalfaorb))
     allocate(T_b(nbetaorb,nprim)); allocate(occ_b(nalfaorb))
     sma=0; smb=0 !sum for alpha and beta orbitals
     do i=1,noccmo
        write(*,*) "SPIN=",omg(i)
        if (omg(i).eq.1) then !alpha
          sma=sma+1
          do j=1,nprim    
                T_a(sma,j)=T(i,j)  !store the prim. coefficients of alpha orb.   
          end do
          
          occ_a(sma)=occ(i) !store occupation numbers
        else if (omg(i).eq.2) then !beta
          smb=smb+1     
          do j=1,nprim
                T_b(smb,j)=T(i,j)     
          end do          
          occ_b(smb)=occ(i)
        else
           write(*,*) "ERROR", "SPIN DIFFERENT FROM ALPHA OR BETA"
        end if   
     end do
     write(*,*) sma, nalfaorb
     write(*,*) smb, nbetaorb
end if 
     
 close(1)

write(*,*) "=== Wavefunction Parameters ===" !JD
write(*,*) "Number of electrons: ", nelec !JD
write(*,*) "Number of occupied MOs: ", noccmo !JD
write(*,*) "Number of primitives: ", nprim !JD
write(*,*) "Open shell flag: ", opsh !JD
write(*,*) "UHF flag: ", uhf !JD
write(*,*) "Correlation flag: ", corr !JD

write(*,*) "First 5 orbital occupations:" !JD
do i=1,min(5,noccmo)
    write(*,*) "Occ(", i, ") = ", Occ(i) !JD
end do

end subroutine filewfx

subroutine wfxhf(wfxhfname)
        use wfxinfo
        use geninfo
        use located
        use wfxinfo_hf !different T matrix and Orbitals than previous calculation
        implicit none
        character*40, intent(in) :: wfxhfname
        integer :: i,j,mon
        character*80 :: line
        character*80 :: zaborra
        character*40 :: spin
        integer, allocatable, dimension(:) :: omg
        integer :: sma, smb
        open(unit=1,file=wfxhfname,status='OLD')
        rewind 1
        call locate(1,"Number of Occupied Molecular Orbitals")
        read(1,'(a80)')line
        read(line(1:5),*) noccmo2
        write(*,*) "Allocating Occ2"
        allocate(Occ2(noccmo2))
  
        rewind 1
        call locate(1,"Molecular Orbital Occupation Numbers") !check RHF,UHF, or Correlated
        do i=1,noccmo2                                         !by
          read(1,*) Occ2(i)                                   !reading the MO occupation values
        end do
        
        if (minval(Occ2).eq.2.d0) then
            write(*,*) "RESTRICTED HARTREE-FOCK" !be careful with the variables!!!!
        else if (minval(Occ).eq.1.d0) then    
            uhf=.true.
        else
           write(*,*) "ERROR", "Uncorrelated wavefunction with fractional occupations"
           write(*,*) "Expected a HF calculation"
           STOP    
        end if      
        
        allocate(T2(noccmo2,nprim))
        call locate(1,"Molecular Orbital Primitive Coefficients")
        do i=1,noccmo2
             read(1,'(a80)')line
             read(1,*) mon                          !molecular orbital number "mo"
             read(1,'(a80)') zaborra
             read(1,*) (T2(mon,j), j=1,nprim)
        end do
        
        if (uhf) then !read the coefficients separatelly 
               call locate(1,"Molecular Orbital Spin Types")
               nalfaorb2=0; nbetaorb2=0
               do i=1,noccmo
                   read(1,*) spin
                   if (spin.eq."Alpha") then
                      omg(i)=1
                       nalfaorb2=nalfaorb+1
                   else
                       omg(i)=2 
                       nbetaorb2=nbetaorb+1
                   end if       
               end do
               do i=1,noccmo
                   if (omg(i).eq.1) then !alpha
                     do j=1,nprim    
                        T_a2(sma,j)=T(i,j)  !store the prim. coefficients of alpha orb.   
                     end do
                     sma=sma+1
                     !noccmo_a(sma)=noccmo(i)
                   else if (omg(i).eq.2) then !beta       
                     do j=1,nprim
                        T_b2(smb,j)=T(i,j)     
                     end do
                     smb=smb+1
                     !noccmo_b(smb)=noccmo(i) !all the orbitals have 1 occ
                   else
                     write(*,*) "ERROR", "SPIN DIFFERENT FROM ALPHA OR BETA"
                     STOP
                   end if   
               end do
        end if
write(*,*) "=== Wavefunction Parameters ===" !JD
write(*,*) "Number of electrons: ", nelec !JD
write(*,*) "Number of occupied MOs: ", noccmo !JD
write(*,*) "Number of primitives: ", nprim !JD
write(*,*) "Open shell flag: ", opsh !JD
write(*,*) "UHF flag: ", uhf !JD
write(*,*) "Correlation flag: ", corr !JD

write(*,*) "First 5 orbital occupations:" !JD
do i=1,min(5,noccmo)
    write(*,*) "Occ(", i, ") = ", Occ(i) !JD
end do
                
end subroutine wfxhf     

!****************************************************************************
subroutine filelog(logfilename)                                             !
!1-reads d (Flg) fixed coeficients of the primitives			       !
!2-Normalizes that coeficient (with respect to the primitives, depending on !
!  primitive type (ptyp)						       !
!3-Builds C matrix (expansion of MOs into AOs)                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use loginfo                    
use geninfo
use wfxinfo
use located
implicit none
integer, parameter :: maxao=100
character*10 :: logfilename
character*80 :: line 
character*80 :: zaborra
character*10 :: tao !atomic orbital type
integer :: i,j,k,l,m
integer :: sm, smm, sm2, sm3, sm4 !sums !smp is the sum of the previous primitives
integer :: smatom
real :: kk !alpha value of the primitives (already got it from wfx) 

open(unit=4 ,file=logfilename, status='OLD') !obre arxiu amb nom name a la unitat 4

allocate(npao(maxao)) !number of primitives of each AO basis. 
allocate(Flg(nprim))  !d coeficient for each primitive
allocate(aotyp(maxao))!atomic orbital type

Flg=0.d0

smm=0 !sum primitives
sm=0  !sum hybrid AOs (S, SP, D, ...)
smatom=1
call locate(4,"AO basis set in the form of general basis input") 
read(4,'(a80)') zaborra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!read fixed coeficients of the AOs in the primitives!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do while (smatom.le.natoms)   
      read(4,'(a80)') line
      read(line(2:3),*) tao 
      if (tao.eq."**") then !new atom finished 
           smatom=smatom+1
           read(4,'(a80)') zaborra   
      else if (tao.eq."S") then
           sm=sm+1  !one more AO
           read(line(5:7),*) npao(sm) !number of primitives of S ao   
           do i=1,npao(sm)
                smm=smm+1 !one more primitive
                read(4,'(a80)') line
                read(line(24:39),*) Flg(smm) !read d coeficient of smm primitive  
           end do
           aotyp(sm)=0  !S type ao
       else if (tao.eq."SP") then 
            sm=sm+1 !one more ao, S type
            aotyp(sm)=0
            read(line(5:7),*) npao(sm)      
            do i=1,npao(sm)
                 smm=smm+1  !one more primitive
                 sm2=sm+1     !px AO
                 aotyp(sm2)=1
                 sm3=sm+2     !py AO
                 aotyp(sm3)=1
                 sm4=sm+3     !pz AO  
                 aotyp(sm4)=1  
                 k=smm+npao(sm) !k is the p type primitives, smm is for the S type primitives     
                 read(4,*) kk, Flg(smm), Flg(k) !read d coeficients for S and for P primitives
                 l=k+npao(sm) !py primitive number                  
                 m=l+npao(sm) !pz primitive number    
                 Flg(l)=Flg(k)     !same value(C_px1=C_py1=C_pz1, 
                                   !(C_px2=C_py2=C_pz2, ...)  
                 Flg(m)=Flg(k)     
            end do   
            smm=smm+3*npao(sm)     !sum all the primitives used in this P type orbital                
            do i=1,3
                 npao(sm+i)=npao(sm) !In an SP hybrid orbital, px, py and pz AOs have the same number of primitives than S
            end do
            sm=sm+3 !sum px, py and pz to the total number of AOs
       else if (tao.eq."D") then
            sm=sm+1  !one more ao
            aotyp=2 !dx, dy, dz  !aotyp is used to the normalization of d coeficient (see norm subroutine)
            read(line(5:7),*) npao(sm)   
            do i=1,npao(sm) !do for the number of primitives of this AO
                 smm=smm+1  
                 read(4,*) kk, Flg(smm)      
                 !check if we have cartesian or spherical coordinates    
                 do j=1,5  !six AOs with same npao and F coeficient (if cartesian)
                        k=smm+j
                        l=sm+j  !current AO 
                        if (l.le.3) then 
                                 aotyp(l)=2  !dx,dy or dz coordinate
                        else
                                 aotyp(l)=3  !square (different normalization)  
                        end if 
                        Flg(k)=Flg(smm)      !same coeficient for each primitive of this AO 
                        npao(sm+j)=npao(sm)    
                 end do
            end do    
            smm=smm+6*npao(sm)  !sum all the primitives and AOs used. 
            sm=sm+5     
       end if
end do
nao=sm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!Normalize for each d coeficient of each primitive!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(N_prim(nprim))  
do i=1,nprim
      if (ptyp(i).eq.1) then 
            N_prim(i)=(pi*(2*Alpha(i))**(-1))**(-3.d0/4.d0)
           write(*,*) N_prim(i)
      else if ((ptyp(i).gt.1).and.(ptyp(i).le.4)) then
            N_prim(i)= (pi*0.5d0)**(-3.d0/4.d0) * (2*(Alpha(i)**(1.25d0)))   
      else if ((ptyp(i).gt.4).and.(ptyp(i).le.7)) then  !dxy,dxz.dyz
            N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))
      else if ((ptyp(i).gt.7).and.(ptyp(i).le.10)) then  !dx²,dy².dz²
            N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))/(sqrt(3.d0))   
      else
            write(*,*) "f type orbital, insert its normalization equation"
      end if    
end do
do i=1, nprim              !Multiply the normalization of the primitives
      Flg(i)=Flg(i)*N_prim(i)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!Create expansion coef. of MOs in AOs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cmatrix() 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!102 format(30(F4.2,1X))
!103 format(30(I4,1X)) 
 close(4)
end subroutine filelog

!*******************************************************************************
subroutine cmatrix()   
!!calculates C with F and T matrices. C are the contracted coefficients of AOs in
!MOs. F are coefficients of AOs in primitives and T are for MOs in Primirives.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use wfxinfo
use loginfo
use geninfo
implicit none
integer :: i, j, k
     
allocate(Ckalk(noccmo,nao))
write(*,*) nao
write(*,*) noccmo
do i=1,noccmo
 k=1
 do j=1,nao
      Ckalk(i,j)=T(i,k)/(Flg(k)) 
      k=k+npao(j)  
  end do
end do   

! 41 format(10(F7.3,1X))
end subroutine cmatrix

