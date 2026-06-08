  !We build the electronic part of the electron-phonon interaction: Sum_iorb g_iorb n_iorb
  do i=MpiIstart,MpiIend
     call state2indices(i,[DimUps,DimDws],Indices)
     do iud=1,Ns_Ud
        mup = Hsector%H(iud)%map(Indices(iud))
        mdw = Hsector%H(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb)
     enddo
     Nup = Breorder(Nups)
     Ndw = Breorder(Ndws)
     !
     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0e_eph,one*htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0e_eph,one*htmp,i,i)
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
              j   = jup + (jdw-1)*DimUp
              htmp = g_ph(iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0e_eph,one*htmp,j,i)
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
              j   = jup + (jdw-1)*DimUp
              htmp = g_ph(iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0e_eph,one*htmp,j,i)
              !
           endif
              
        enddo
     enddo
     !
  enddo

  !Here we build the phononc part of the electron-phonon interaction: (b^+ + b)
  htmp = zero
  do iph=1,DimPh
     i = iph + 1
     if(i <= DimPh) then
        htmp = sqrt(dble(iph))
        call sp_insert_element(spH0ph_eph,htmp,iph,i)
     end if
     i = iph - 1
     if(i>0) then
        htmp = sqrt(dble(iph - 1))
        call sp_insert_element(spH0ph_eph,one*one*htmp,iph,i)
     end if
  end do













