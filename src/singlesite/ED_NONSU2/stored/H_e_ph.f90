!Electron-Phono hamiltonian H_eph =  \sum_sigma \sum_m,n g_m,n cdagger_m,sigma c_n,sigma ( bdg + b ) 
do i=MpiIstart,MpiIend
    i_el = mod(i-1,DimEl) + 1
    iph = (i-1)/DimEl + 1
    m  = Hsector%H(1)%map(i_el)
    ib = bdecomp(m,2*Ns)
    do iorb=1,Norb
       nup(iorb)=dble(ib(iorb))
       ndw(iorb)=dble(ib(iorb+Ns))
    enddo

     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb
     !
     if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|     
        htmp = zero
        do iorb=1,Norb
           htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
        enddo
        htmp = htmp*sqrt(dble(iph))
        j = i_el + (iph)*DimEl
        select case(MpiStatus)
            case (.true.)
                call sp_insert_element(MpiComm,spH0,htmp,i,j)
            case (.false.)
                call sp_insert_element(spH0,htmp,i,j)
        end select
    endif
    if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|  
        htmp = zero
        do iorb=1,Norb
           htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
        enddo
        htmp = htmp*sqrt(dble(iph-1))
        j = i_el + (iph-2)*DimEl
        select case(MpiStatus)
            case (.true.)
                call sp_insert_element(MpiComm,spH0,htmp,i,j)
            case (.false.)
                call sp_insert_element(spH0,htmp,i,j)
        end select
    endif
    !
    ! UP
    do iorb=1,Norb
        do jorb=1,Norb
            ! g_ij cdag_i,up c_j,up (bdag + b)
            if( g_ph(iorb,jorb)/= (0.d0,0.d0) .and. nup(jorb)==1 .and. nup(iorb)==0 )then
                call c(  jorb, m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j_el = binary_search(Hsector%H(1)%map,k2)
                ! N.B.here iph = n+1
                if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
                    htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph))
                    j = j_el + (iph)*DimEl
                    select case(MpiStatus)
                        case (.true.)
                            call sp_insert_element(MpiComm,spH0,htmp,i,j)
                        case (.false.)
                            call sp_insert_element(spH0,htmp,i,j)
                    end select
                endif
                if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
                    htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph-1))
                    j = j_el + (iph-2)*DimEl
                    select case(MpiStatus)
                        case (.true.)
                            call sp_insert_element(MpiComm,spH0,htmp,i,j)
                        case (.false.)
                            call sp_insert_element(spH0,htmp,i,j)
                    end select
                endif
            endif
        end do
    end do
    !
    ! DW
    do iorb=1,Norb
        do jorb=1,Norb
            ! g_ij cdag_i,dw c_j,dw (bdag + b)
            if( g_ph(iorb,jorb)/= (0.d0,0.d0) .and. ndw(jorb)==1 .and. ndw(iorb)==0 )then
                call c(  jorb+Ns, m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j_el = binary_search(Hsector%H(1)%map,k2)
                ! N.B.here iph = n+1
                if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
                    htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph))
                    j = j_el + (iph)*DimEl
                    select case(MpiStatus)
                        case (.true.)
                            call sp_insert_element(MpiComm,spH0,htmp,i,j)
                        case (.false.)
                            call sp_insert_element(spH0,htmp,i,j)
                    end select
                endif
                if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
                    htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph-1))
                    j = j_el + (iph-2)*DimEl
                    select case(MpiStatus)
                        case (.true.)
                            call sp_insert_element(MpiComm,spH0,htmp,i,j)
                        case (.false.)
                            call sp_insert_element(spH0,htmp,i,j)
                    end select
                endif
            endif
        end do
    end do
    !
enddo
