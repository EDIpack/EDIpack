!Electron-Phonon hamiltonian H_eph =  \sum_sigma \sum_m,n g_m,n cdagger_m,sigma c_n,sigma ( bdg + b ) 
!
! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb*(bdg + b)
htmp=zero
do iorb=1,Norb
    htmp = htmp + g_ph(iorb,iorb)*( nup(iorb)+ndw(iorb) ) !electron part
enddo
!
if(iph<DimPh)then !bdg
    i = i_el + (iph  )*DimEl
    Hv(j-MpiIshift) = Hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph))
endif
if(iph>1)then !b
    i = i_el + (iph-2)*DimEl
    Hv(j-MpiIshift) = Hv(j-MpiIshift) + htmp*vin(i)*sqrt(dble(iph-1))
endif
! UP
do iorb=1,Norb
    do jorb=1,Norb
        ! g_ij cdag_i,up c_j,up (bdag + b)
        if( g_ph(iorb,jorb)/= (0.d0,0.d0) .and. nup(jorb)==1 .and. nup(iorb)==0 )then
            call c(  jorb, m,k1,sg1)
            call cdg(iorb,k1,k2,sg2)
            i_el = binary_search(Hsector%H(1)%map,k2)
            ! N.B.here iph = n+1
            if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
                htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph))
                i = i_el + (iph)*DimEl
                hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
            endif
            if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
                htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph-1))
                i = i_el + (iph-2)*DimEl
                hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
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
            i_el = binary_search(Hsector%H(1)%map,k2)
            ! N.B.here iph = n+1
            if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
                htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph))
                i = i_el + (iph)*DimEl
                hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
            endif
            if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
                htmp = g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph-1))
                i = i_el + (iph-2)*DimEl
                hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
            endif
        endif
    end do
end do
