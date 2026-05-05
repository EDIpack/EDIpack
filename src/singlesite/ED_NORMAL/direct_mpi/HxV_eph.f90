  !We build the electronic part of the electron-phonon interaction:
  do i=1,Nloc
     i_el = mod(i-1,DimUp*MpiQdw) + 1
     iph = (i-1)/(DimUp*MpiQdw) + 1
     !
     iup = iup_index(i_el+mpiIshift,DimUp)
     idw = idw_index(i_el+mpiIshift,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb*(bdg + b)
     htmp=zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
     enddo
     !
     if( iph<DimPh )then !bdg
        j     = i_el + (iph  )*DimUp*MpiQdw
        Hv(j) = Hv(j) + htmp*sqrt(dble(iph))*vin(i)
     endif
     if(iph>1) then  !b
        j     = i_el + (iph-2)*DimUp*MpiQdw
        Hv(j) = Hv(j) + htmp*sqrt(dble(iph-1))*vin(i)
     endif
     !
     ! Off-Diagonal terms: Sum_iorb,jorb g_iorb,jorb cdg_iorb*c_jorb*(bdg+b)
     ! UP spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = (g_ph(iorb,jorb)/=zero) .AND. &
                (nup(jorb)==1) .AND. (nup(iorb)==0)
           if(Jcondition)then
              call c(jorb,mup,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup  = binary_search(Hsector%H(1)%map,k2)
              jdw  = idw
              htmp = g_ph(iorb,jorb)*sg1*sg2
              !
              if(iph<Nph)then !bdg
                 j     = jup + (jdw-1)*dimUp + (iph)*DimUp*MpiQdw
                 hv(j) = hv(j) + sqrt(dble(iph))*htmp*vin(i)
              end if
              if(iph>1)then !b
                 j     = jup + (jdw-1)*dimUp + (iph-2)*DimUp*MpiQdw
                 hv(j) = hv(j) + sqrt(dble(iph-1))*htmp*vin(i)
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
              call c(jorb,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw = binary_search(Hsector%H(2)%map,k2)
              jup  = iup
              htmp = g_ph(iorb,jorb)*sg1*sg2
              !
              if(iph<Nph)then !bdg
                 j     = jup + (jdw-1)*dimUp + (iph)*DimUp*MpiQdw
                 hv(j) = hv(j) + sqrt(dble(iph))*htmp*vin(i)
              end if
              if(iph>1  )then !b
                 j     = jup + (jdw-1)*dimUp + (iph-2)*DimUp*MpiQdw
                 hv(j) = hv(j) + sqrt(dble(iph-1))*htmp*vin(i)
              endif
           endif
           !              
        enddo
     enddo
  enddo
