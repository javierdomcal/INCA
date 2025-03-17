subroutine wrthfdm2()
    use wfxinfo_hf
    use geninfo
    implicit none

    ! Local module-level cache variables
    integer :: num_sig_orbs
    double precision, allocatable :: T2_cache(:,:)
    double precision, allocatable :: occ2_cache(:)
    integer, allocatable :: sig_orbitals(:)
    
    ! Local variables
    double precision :: tot_hf, densp_hf, sqrdm1_hf
    double precision :: tot, densp, sqrdm1
    integer :: i, j, k, l, sm, a
    
    ! Arrays for I/O optimization
    integer, parameter :: MAX_BUFFER = 10000
    integer, allocatable :: idx_buffer(:,:)
    double precision, allocatable :: val_buffers(:,:)
    integer :: buffer_count
    
    ! Arrays for primitive screening
    logical, allocatable :: primitive_contributes(:)
    integer :: nonzero_primitives
    
    write(*,*) "Starting DM2HF exchange calculation"
    
    ! Initialize caches
    num_sig_orbs = count(occ2(1:noccmo2) > 1.0d-10)
    allocate(sig_orbitals(num_sig_orbs))
    allocate(occ2_cache(num_sig_orbs))
    
    ! Fill significant orbital arrays
    sm = 0
    do a = 1, noccmo2
        if (occ2(a) > 1.0d-10) then
            sm = sm + 1
            sig_orbitals(sm) = a
            occ2_cache(sm) = occ2(a)
        endif
    end do
    
    ! Cache T2 matrix elements
    allocate(T2_cache(num_sig_orbs, size(T2,2)))
    do a = 1, num_sig_orbs
        T2_cache(a,:) = T2(sig_orbitals(a),:)
    end do
    
    ! Pre-screen primitives
    allocate(primitive_contributes(nprim))
    call prescreen_primitives(primitive_contributes, nonzero_primitives)
    write(*,*) "Number of contributing primitives:", nonzero_primitives
    
    ! Allocate I/O buffers
    allocate(idx_buffer(4,MAX_BUFFER))
    allocate(val_buffers(9,MAX_BUFFER))
    buffer_count = 0
    
    ! Open all output files
    open(unit=10, file="dm2hf.dat", form="UNFORMATTED", buffered='YES')
    open(unit=11, file="dm1hf.dat", form="UNFORMATTED", buffered='YES')
    open(unit=12, file="densphf.dat", form="UNFORMATTED", buffered='YES')
    open(unit=13, file="dm2.dat", form="UNFORMATTED", buffered='YES')
    open(unit=14, file="dm1.dat", form="UNFORMATTED", buffered='YES')
    open(unit=15, file="densp.dat", form="UNFORMATTED", buffered='YES')
    open(unit=16, file="c1tot.dat", form="UNFORMATTED", buffered='YES')
    open(unit=17, file="c1s2.dat", form="UNFORMATTED", buffered='YES')
    open(unit=18, file="hfls2.dat", form="UNFORMATTED", buffered='YES')
    
    sm = 0
    if (opsh) then
        ! Open shell case handling
        do i = 1, nprim
            if (.not. primitive_contributes(i)) cycle
            do k = 1, nprim
                if (.not. primitive_contributes(k)) cycle
                do j = 1, nprim
                    if (.not. primitive_contributes(j)) cycle
                    do l = 1, nprim
                        if (.not. primitive_contributes(l)) cycle
                        
                        call DM2prim_hf_oshell(i,j,k,l,tot_hf,densp_hf,sqrdm1_hf)
                        call DM2prim_oshell(i,j,k,l,tot,densp,sqrdm1)
                        
                        if (any(abs([tot_hf,densp_hf,sqrdm1_hf,tot,densp,sqrdm1]) > 1.d-10)) then
                            call add_to_buffer(i,j,k,l,tot_hf,densp_hf,sqrdm1_hf, &
                                            tot,densp,sqrdm1,buffer_count)
                            sm = sm + 1
                        endif
                    end do
                end do
            end do
        end do
    else
        ! Closed shell case
        do i = 1, nprim
            if (.not. primitive_contributes(i)) cycle
            do k = 1, nprim
                if (.not. primitive_contributes(k)) cycle
                do j = 1, nprim
                    if (.not. primitive_contributes(j)) cycle
                    do l = 1, nprim
                        if (.not. primitive_contributes(l)) cycle
                        
                        call DM2prim_hf(i,j,k,l,tot_hf,densp_hf,sqrdm1_hf)
                        call DM2prim(i,j,k,l,tot,densp,sqrdm1)
                        
                        if (any(abs([tot_hf,densp_hf,sqrdm1_hf,tot,densp,sqrdm1]) > 1.d-10)) then
                            call add_to_buffer(i,j,k,l,tot_hf,densp_hf,sqrdm1_hf, &
                                            tot,densp,sqrdm1,buffer_count)
                            sm = sm + 1
                        endif
                    end do
                end do
            end do
        end do
    end if
    
    ! Flush remaining buffer
    if (buffer_count > 0) then
        call flush_buffers(buffer_count)
    endif
    
    ! Cleanup
    close(10); close(11); close(12); close(13); close(14)
    close(15); close(16); close(17); close(18)
    deallocate(primitive_contributes, idx_buffer, val_buffers)
    deallocate(T2_cache, occ2_cache, sig_orbitals)
    
    write(*,*) "Number of non-zero elements:", sm
    write(*,*) "End writing DM2HF exchange"

contains
    subroutine prescreen_primitives(prim_contrib, count)
        logical, intent(out) :: prim_contrib(:)
        integer, intent(out) :: count
        integer :: i, j
        
        prim_contrib = .false.
        count = 0
        
        !$OMP PARALLEL DO REDUCTION(+:count) PRIVATE(i,j)
        do i = 1, nprim
            do j = 1, num_sig_orbs
                if (abs(T2_cache(j,i)) > 1.d-10) then
                    prim_contrib(i) = .true.
                    count = count + 1
                    exit
                endif
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine

    subroutine add_to_buffer(i,j,k,l,tot_hf,densp_hf,sqrdm1_hf,tot,densp,sqrdm1,count)
        integer, intent(in) :: i,j,k,l
        double precision, intent(in) :: tot_hf,densp_hf,sqrdm1_hf,tot,densp,sqrdm1
        integer, intent(inout) :: count
        
        count = count + 1
        idx_buffer(:,count) = [i,j,k,l]
        val_buffers(:,count) = [tot_hf,densp_hf,sqrdm1_hf,tot,densp,sqrdm1, &
                               tot-tot_hf, tot-tot_hf, tot]
        
        if (count == MAX_BUFFER) then
            call flush_buffers(count)
            count = 0
        endif
    end subroutine

    subroutine flush_buffers(count)
        integer, intent(in) :: count
        integer :: i
        
        do i = 1, count
            write(10) idx_buffer(:,i), val_buffers(1,i)  ! dm2hf
            write(11) idx_buffer(:,i), val_buffers(2,i)  ! dm1hf
            write(12) idx_buffer(:,i), val_buffers(3,i)  ! densphf
            write(13) idx_buffer(:,i), val_buffers(4,i)  ! dm2
            write(14) idx_buffer(:,i), val_buffers(5,i)  ! dm1
            write(15) idx_buffer(:,i), val_buffers(6,i)  ! densp
            write(16) idx_buffer(:,i), val_buffers(7,i)  ! c1tot
            write(17) idx_buffer(:,i), val_buffers(8,i)  ! c1s2
            write(18) idx_buffer(:,i), val_buffers(9,i)  ! hfls2
        end do
    end subroutine

end subroutine wrthfdm2
! Regular DM2 calculation (non-HF)
subroutine DM2prim(i, j, k, l, tot, densp, sqrdm1)
    use wfxinfo
    implicit none
    double precision, intent(out) :: tot, densp, sqrdm1
    integer, intent(in) :: i, j, k, l
    integer :: a, b
    double precision :: coef_prod, occ_prod
    double precision :: Ti_a, Tk_a, Tj_b, Tl_b
    
    ! Initialize outputs
    tot = 0.d0
    densp = 0.d0
    sqrdm1 = 0.d0
    
    ! Quick screening check
    if (all(abs(T(:,i)) < 1.d-10) .or. &
        all(abs(T(:,j)) < 1.d-10) .or. &
        all(abs(T(:,k)) < 1.d-10) .or. &
        all(abs(T(:,l)) < 1.d-10)) return
    
    ! Main computation with optimized loop structure
    do a = 1, noccmo
        if (Occ(a) < 1.d-10) cycle  ! Skip unoccupied orbitals
        
        Ti_a = T(a,i)
        if (abs(Ti_a) < 1.d-10) cycle
        
        Tk_a = T(a,k)
        if (abs(Tk_a) < 1.d-10) cycle
        
        do b = 1, noccmo
            if (Occ(b) < 1.d-10) cycle
            
            Tj_b = T(b,j)
            if (abs(Tj_b) < 1.d-10) cycle
            
            Tl_b = T(b,l)
            if (abs(Tl_b) < 1.d-10) cycle
            
            occ_prod = Occ(a) * Occ(b)
            
            ! Direct term
            coef_prod = Ti_a * Tk_a * Tj_b * Tl_b
            tot = tot + occ_prod * coef_prod
            densp = densp + occ_prod * coef_prod
            
            ! Exchange term
            if (abs(T(a,j)) > 1.d-10 .and. abs(T(b,k)) > 1.d-10) then
                coef_prod = -0.5d0 * Ti_a * T(b,k) * T(a,j) * Tl_b
                tot = tot + occ_prod * coef_prod
                sqrdm1 = sqrdm1 + occ_prod * coef_prod
            endif
        end do
    end do
end subroutine DM2prim

! Open shell DM2 calculation
subroutine DM2prim_oshell(i, j, k, l, tot, densp, sqrdm1)
    use wfxinfo
    implicit none
    double precision, intent(out) :: tot, densp, sqrdm1
    integer, intent(in) :: i, j, k, l
    integer :: a, b
    double precision :: d_aa, d_bb, d_ab
    double precision :: densp_a, densp_b
    double precision :: sqrdm1_a, sqrdm1_b
    double precision :: Ti_a, Tk_a, Tj_b, Tl_b
    
    ! Initialize all components
    tot = 0.d0
    densp = 0.d0
    sqrdm1 = 0.d0
    d_aa = 0.d0
    d_bb = 0.d0
    d_ab = 0.d0
    densp_a = 0.d0
    densp_b = 0.d0
    sqrdm1_a = 0.d0
    sqrdm1_b = 0.d0
    
    ! Alpha-Alpha contributions
    if (nalfaorb > 0) then
        do a = 1, nalfaorb
            Ti_a = T_a(a,i)
            if (abs(Ti_a) < 1.d-10) cycle
            
            Tk_a = T_a(a,k)
            if (abs(Tk_a) < 1.d-10) cycle
            
            do b = 1, nalfaorb
                Tj_b = T_a(b,j)
                if (abs(Tj_b) < 1.d-10) cycle
                
                Tl_b = T_a(b,l)
                if (abs(Tl_b) < 1.d-10) cycle
                
                ! Direct term
                d_aa = d_aa + Ti_a * Tk_a * Tj_b * Tl_b
                densp_a = densp_a + Ti_a * Tk_a * Tj_b * Tl_b
                
                ! Exchange term
                if (abs(T_a(a,j)) > 1.d-10 .and. abs(T_a(b,k)) > 1.d-10) then
                    sqrdm1_a = sqrdm1_a - T_a(a,i) * T_a(b,k) * T_a(a,j) * T_a(b,l)
                    d_aa = d_aa - T_a(a,i) * T_a(b,k) * T_a(a,j) * T_a(b,l)
                endif
            end do
        end do
    endif
    
    ! Beta-Beta contributions
    if (nbetaorb > 0) then
        do a = 1, nbetaorb
            Ti_a = T_b(a,i)
            if (abs(Ti_a) < 1.d-10) cycle
            
            Tk_a = T_b(a,k)
            if (abs(Tk_a) < 1.d-10) cycle
            
            do b = 1, nbetaorb
                Tj_b = T_b(b,j)
                if (abs(Tj_b) < 1.d-10) cycle
                
                Tl_b = T_b(b,l)
                if (abs(Tl_b) < 1.d-10) cycle
                
                ! Direct term
                d_bb = d_bb + Ti_a * Tk_a * Tj_b * Tl_b
                densp_b = densp_b + Ti_a * Tk_a * Tj_b * Tl_b
                
                ! Exchange term
                if (abs(T_b(a,j)) > 1.d-10 .and. abs(T_b(b,k)) > 1.d-10) then
                    sqrdm1_b = sqrdm1_b - T_b(a,i) * T_b(b,k) * T_b(a,j) * T_b(b,l)
                    d_bb = d_bb - T_b(a,i) * T_b(b,k) * T_b(a,j) * T_b(b,l)
                endif
            end do
        end do
    endif
    
    ! Alpha-Beta contributions (no exchange)
    if (nalfaorb > 0 .and. nbetaorb > 0) then
        do a = 1, nalfaorb
            Ti_a = T_a(a,i)
            if (abs(Ti_a) < 1.d-10) cycle
            
            Tk_a = T_a(a,k)
            if (abs(Tk_a) < 1.d-10) cycle
            
            do b = 1, nbetaorb
                Tj_b = T_b(b,j)
                if (abs(Tj_b) < 1.d-10) cycle
                
                Tl_b = T_b(b,l)
                if (abs(Tl_b) < 1.d-10) cycle
                
                d_ab = d_ab + Ti_a * Tk_a * Tj_b * Tl_b
            end do
        end do
    endif
    
    ! Combine all components
    tot = d_aa + d_bb + d_ab
    densp = densp_a + densp_b + d_ab
    sqrdm1 = sqrdm1_a + sqrdm1_b
end subroutine DM2prim_oshell

! HF DM2 calculation for closed shell
subroutine DM2prim_hf(i, j, k, l, tot_hf, densp_hf, sqrdm1_hf)
    use wfxinfo_hf
    use geninfo
    implicit none
    double precision, intent(out) :: tot_hf, densp_hf, sqrdm1_hf
    integer, intent(in) :: i, j, k, l
    integer :: a, b
    double precision :: coef_prod, occ_prod
    double precision :: Ti_a, Tk_a, Tj_b, Tl_b
    
    ! Initialize outputs
    tot_hf = 0.d0
    densp_hf = 0.d0
    sqrdm1_hf = 0.d0
    
    ! Quick screening check
    if (all(abs(T2(:,i)) < 1.d-10) .or. &
        all(abs(T2(:,j)) < 1.d-10) .or. &
        all(abs(T2(:,k)) < 1.d-10) .or. &
        all(abs(T2(:,l)) < 1.d-10)) return
    
    ! Main computation with optimized loop structure
    do a = 1, noccmo2
        if (Occ2(a) < 1.d-10) cycle  ! Skip unoccupied orbitals
        
        Ti_a = T2(a,i)
        if (abs(Ti_a) < 1.d-10) cycle
        
        Tk_a = T2(a,k)
        if (abs(Tk_a) < 1.d-10) cycle
        
        do b = 1, noccmo2
            if (Occ2(b) < 1.d-10) cycle
            
            Tj_b = T2(b,j)
            if (abs(Tj_b) < 1.d-10) cycle
            
            Tl_b = T2(b,l)
            if (abs(Tl_b) < 1.d-10) cycle
            
            occ_prod = Occ2(a) * Occ2(b)
            
            ! Direct term
            coef_prod = Ti_a * Tk_a * Tj_b * Tl_b
            tot_hf = tot_hf + occ_prod * coef_prod
            densp_hf = densp_hf + occ_prod * coef_prod
            
            ! Exchange term
            if (abs(T2(b,k)) > 1.d-10 .and. abs(T2(a,j)) > 1.d-10) then
                coef_prod = -0.5d0 * Ti_a * T2(b,k) * T2(a,j) * Tl_b
                tot_hf = tot_hf + occ_prod * coef_prod
                sqrdm1_hf = sqrdm1_hf + occ_prod * coef_prod
            endif
        end do
    end do
    
    ! Normalize values
    tot_hf = tot_hf * 0.5d0
    densp_hf = densp_hf * 0.5d0
    sqrdm1_hf = sqrdm1_hf * 0.5d0
end subroutine DM2prim_hf

! HF DM2 calculation for open shell
subroutine DM2prim_hf_oshell(i, j, k, l, tot_hf, densp_hf, sqrdm1_hf)
    use wfxinfo_hf
    use geninfo
    implicit none
    double precision, intent(out) :: tot_hf, densp_hf, sqrdm1_hf
    integer, intent(in) :: i, j, k, l
    integer :: a, b
    double precision :: hf_aa, hf_bb, hf_ab
    double precision :: densp_hf_a, densp_hf_b
    double precision :: sqrdm1_a, sqrdm1_b
    double precision :: Ti_a, Tk_a, Tj_b, Tl_b
    
    ! Initialize all components
    tot_hf = 0.d0
    densp_hf = 0.d0
    sqrdm1_hf = 0.d0
    hf_aa = 0.d0
    hf_bb = 0.d0
    hf_ab = 0.d0
    densp_hf_a = 0.d0
    densp_hf_b = 0.d0
    sqrdm1_a = 0.d0
    sqrdm1_b = 0.d0
    
    ! Alpha-Alpha contributions
    if (nalfaorb2 > 0) then
        do a = 1, nalfaorb2
            Ti_a = T_a2(a,i)
            if (abs(Ti_a) < 1.d-10) cycle
            
            Tk_a = T_a2(a,k)
            if (abs(Tk_a) < 1.d-10) cycle
            
            do b = 1, nalfaorb2
                Tj_b = T_a2(b,j)
                if (abs(Tj_b) < 1.d-10) cycle
                
                Tl_b = T_a2(b,l)
                if (abs(Tl_b) < 1.d-10) cycle
                
                ! Direct term
                hf_aa = hf_aa + Ti_a * Tk_a * Tj_b * Tl_b
                densp_hf_a = densp_hf_a + Ti_a * Tk_a * Tj_b * Tl_b
                
                ! Exchange term
                if (abs(T_a2(a,j)) > 1.d-10 .and. abs(T_a2(b,k)) > 1.d-10) then
                    sqrdm1_a = sqrdm1_a - Ti_a * T_a2(b,k) * T_a2(a,j) * Tl_b
                    hf_aa = hf_aa - Ti_a * T_a2(b,k) * T_a2(a,j) * Tl_b
                endif
            end do
        end do
    endif
    
    ! Beta-Beta contributions
    if (nbetaorb2 > 0) then
        do a = 1, nbetaorb2
            Ti_a = T_b2(a,i)
            if (abs(Ti_a) < 1.d-10) cycle
            
            Tk_a = T_b2(a,k)
            if (abs(Tk_a) < 1.d-10) cycle
            
            do b = 1, nbetaorb2
                Tj_b = T_b2(b,j)
                if (abs(Tj_b) < 1.d-10) cycle
                
                Tl_b = T_b2(b,l)
                if (abs(Tl_b) < 1.d-10) cycle
                
                ! Direct term
                hf_bb = hf_bb + Ti_a * Tk_a * Tj_b * Tl_b
                densp_hf_b = densp_hf_b + Ti_a * Tk_a * Tj_b * Tl_b
                
                ! Exchange term
                if (abs(T_b2(a,j)) > 1.d-10 .and. abs(T_b2(b,k)) > 1.d-10) then
                    sqrdm1_b = sqrdm1_b - Ti_a * T_b2(b,k) * T_b2(a,j) * Tl_b
                    hf_bb = hf_bb - Ti_a * T_b2(b,k) * T_b2(a,j) * Tl_b
                endif
            end do
        end do
    endif
    
    ! Alpha-Beta contributions (no exchange)
    if (nalfaorb2 > 0 .and. nbetaorb2 > 0) then
        do a = 1, nalfaorb2
            Ti_a = T_a2(a,i)
            if (abs(Ti_a) < 1.d-10) cycle
            
            Tk_a = T_a2(a,k)
            if (abs(Tk_a) < 1.d-10) cycle
            
            do b = 1, nbetaorb2
                Tj_b = T_b2(b,j)
                if (abs(Tj_b) < 1.d-10) cycle
                
                Tl_b = T_b2(b,l)
                if (abs(Tl_b) < 1.d-10) cycle
                
                hf_ab = hf_ab + Ti_a * Tk_a * Tj_b * Tl_b
            end do
        end do
    endif
    
    ! Combine all components
    tot_hf = hf_aa + hf_bb + hf_ab
    densp_hf = densp_hf_a + densp_hf_b + hf_ab
    sqrdm1_hf = sqrdm1_a + sqrdm1_b
end subroutine DM2prim_hf_oshell
