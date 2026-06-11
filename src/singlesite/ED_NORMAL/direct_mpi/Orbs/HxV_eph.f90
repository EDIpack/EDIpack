  !We build the electronic part of the electron-phonon interaction:
  do i=1,Nloc
     i_el = mod(i-1,DimUp*MpiQdw) + 1
     iph = (i-1)/(DimUp*MpiQdw) + 1
     !
     call state2indices(i_el+mpiIshift,[DimUps,DimDws],Indices)
     do iud=1,Ns_Ud
        mup = Hsector%H(iud)%map(Indices(iud))
        mdw = Hsector%H(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![1+Nbath]*Norb
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb) ![1+Nbath]*Norb
     enddo
     Nup = Breorder(Nups)
     Ndw = Breorder(Ndws)
     !
     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb*(bdg + b)
     htmp=zero
     do iorb=1,Norb
        htmp = htmp + one*g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
     enddo
     !
     do jj = 1,DimPh
        if(jj .eq. iph+1) then
           j = i_el + (jj-1)*DimUp*MpiQdw 
           Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
        endif
        !
        if(jj .eq. iph-1) then
           j = i_el + (jj-1)*DimUp*MpiQdw
           Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
        endif
     enddo
     !
     ! Off-Diagonal terms: Sum_iorb,jorb g_iorb,jorb cdg_iorb*c_jorb*(bdg+b)
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
              htmp = one*g_ph(iorb,jorb)*sg1*sg2
              !
              if(iph<Nph)then !bdg
                 j     = jup + (jdw-1)*dimUp + (iph)*DimUp*MpiQdw
                 hv(j) = hv(j) + htmp*vin(i)
              elseif(iph>1)then !b
                 j     = jup + (jdw-1)*dimUp + (iph-2)*DimUp*MpiQdw
                 hv(j) = hv(j) + htmp*vin(i)
              endif
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
           if(Jcondition)then
              call c(jorb,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw = binary_search(Hsector%H(2)%map,k2)
              jup  = iup
              htmp = one*g_ph(iorb,jorb)*sg1*sg2
              !
              if(iph<Nph)then !bdg
                 j     = jup + (jdw-1)*dimUp + (iph)*DimUp*MpiQdw
                 hv(j) = hv(j) + htmp*vin(i)
              elseif(iph>1)then !b
                 j     = jup + (jdw-1)*dimUp + (iph-2)*DimUp*MpiQdw
                 hv(j) = hv(j) + htmp*vin(i)
              endif
           endif
           !              
        enddo
     enddo
     
  enddo
