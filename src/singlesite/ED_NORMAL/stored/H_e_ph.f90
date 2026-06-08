  !We build the electronic part of the electron-phonon interaction:
  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0e_eph,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0e_eph,htmp,i,i)
     end select
     !
     ! Off-Diagonal terms: Sum_iorb,jorb g_iorb,jorb cdg_iorb*c_jorb
     ! UP spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (g_ph(iorb,jorb)/=zero) .AND. &
                (Nup(jorb)==1) .AND. (Nup(iorb)==0)
           if(Jcondition)then
              call c(jorb,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup  = binary_search(Hsector%H(1)%map,k2)
              jdw  = idw
              j    = jup + (jdw-1)*DimUp
              htmp = g_ph(iorb,jorb)*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0e_eph,htmp,j,i)
              case (.false.)
                 call sp_insert_element(spH0e_eph,htmp,j,i)
              end select
           endif
        enddo
     enddo
     ! DW spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (g_ph(iorb,jorb)/=zero) .AND. &
                (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
           if (Jcondition) then
              call c(jorb,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup  = iup
              jdw  = binary_search(Hsector%H(2)%map,k2)
              j    = jup + (jdw-1)*DimUp
              htmp = g_ph(iorb,jorb)*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0e_eph,htmp,j,i)
              case (.false.)
                 call sp_insert_element(spH0e_eph,htmp,j,i)
              end select
              !
           endif
              
        enddo
     enddo
     !
  enddo

  ! Here we build the phononic part of the electron-phonon interaction: (b^+ + b)
  htmp = zero
  do iph=1,DimPh
     ! N.B. here iph = n+1
     if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
        htmp = sqrt(dble(iph))
        call sp_insert_element(spH0ph_eph,htmp,iph+1,iph)
     end if
     if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
        htmp = sqrt(dble(iph-1))
        call sp_insert_element(spH0ph_eph,htmp,iph-1,iph)
     end if
  end do
