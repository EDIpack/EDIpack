  do j=1,Nloc                   !DimUp*MpiQdw*DimPh.
     j_el = mod(j-1,DimUp*MpiQdw) + 1 !strip the phonon first
     iph = (j-1)/(DimUp*MpiQdw) + 1   !get phonon index
     !
     jup = iup_index(j_el+mpiIshift,DimUp)
     jdw = idw_index(j_el+mpiIshift,DimUp)
     !
     mup = Hsector%H(1)%map(jup)
     mdw = Hsector%H(2)%map(jdw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! SPIN-EXCHANGE (S-E) TERMS
     !    S-E: J c^+_a_up c^+_b_dw c_a_dw c_b_up
     !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
     if(Norb>1.AND.any((Jx_internal/=0d0)))then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (nup(jorb)==1).AND.&
                   (ndw(iorb)==1).AND.&
                   (ndw(jorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(iorb,mdw,k1,sg1)  !DW
                 call cdg(jorb,k1,k2,sg2) !DW
                 idw=binary_search(Hsector%H(2)%map,k2)
                 call c(jorb,mup,k3,sg3)  !UP
                 call cdg(iorb,k3,k4,sg4) !UP
                 iup=binary_search(Hsector%H(1)%map,k4)
                 htmp = Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4
                 i = iup + (idw-1)*DimUp + (iph-1)*DimUp*MpiQdw
                 !
                 Hv(j) = Hv(j) + htmp*vt(i)
                 !
              endif
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
     !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
     if(any((Jp_internal/=0d0)))then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (nup(jorb)==1).AND.&
                   (ndw(jorb)==1).AND.&
                   (ndw(iorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(jorb,mdw,k1,sg1)       !c_jorb_dw
                 call cdg(iorb,k1,k2,sg2)      !c^+_iorb_dw
                 idw = binary_search(Hsector%H(2)%map,k2)
                 call c(jorb,mup,k3,sg3)       !c_jorb_up
                 call cdg(iorb,k3,k4,sg4)      !c^+_iorb_up
                 iup = binary_search(Hsector%H(1)%map,k4)
                 htmp = Jp_internal(iorb,jorb)*sg1*sg2*sg3*sg4
                 i = iup + (idw-1)*dimup + (iph-1)*DimUp*DimDw
                 !
                 Hv(j) = Hv(j) + htmp*vt(i)
                 !
              endif
           enddo
        enddo
     endif
     !     
  enddo
