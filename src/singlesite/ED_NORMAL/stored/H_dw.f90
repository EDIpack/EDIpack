  do jdw=1,DimDw
     mdw  = Hsector%H(2)%map(jdw)
     Ndw  = bdecomp(mdw,Ns)
     !
     !
     !> H_imp: Off-diagonal elements, i.e. non-local part. 
     !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (impHloc(Nspin,Nspin,iorb,jorb)+mfHloc(2,2,iorb,jorb)/=zero) .AND. &
                (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
           if (Jcondition) then
              call c(jorb,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              idw = binary_search(Hsector%H(2)%map,k2)
              htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
              htmp = htmp + mfHloc(2,2,iorb,jorb)*sg1*sg2
              !
              call sp_insert_element(spH0dws(1),htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !
     !
     !> H_Bath: inter-orbital bath hopping contribution.
     if(bath_type=="replica".or.bath_type=="general") then
        do kp=1,Nbath
           do iorb=1,Norb
              do jorb=1,Norb
                 !
                 ialfa = getBathStride(iorb,kp)
                 ibeta = getBathStride(jorb,kp)
                 Jcondition = &
                      (hbath_tmp(Nspin,Nspin,iorb,jorb,kp)/=zero) &
                      .AND. (Ndw(ibeta)==1) .AND. (Ndw(ialfa)==0)
                 !
                 if (Jcondition)then
                    call c(ibeta,mdw,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    idw = binary_search(Hsector%H(2)%map,k2)
                    htmp = hbath_tmp(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                    !
                    call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                    !
                 endif
              enddo
           enddo
        enddo
        !
     endif
     !
     !
     !>H_hyb: hopping terms for a given spin (imp <--> bath)
     do iorb=1,Norb
        do kp=1,Nbath
           ialfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                .AND. (ndw(iorb)==1) .AND. (ndw(ialfa)==0) )then
              call c(iorb,mdw,k1,sg1)
              call cdg(ialfa,k1,k2,sg2)
              idw=binary_search(Hsector%H(2)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(1),htmp,idw,jdw)
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                .AND. (ndw(iorb)==0) .AND. (ndw(ialfa)==1) )then
              call c(ialfa,mdw,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              idw=binary_search(Hsector%H(2)%map,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(1),htmp,idw,jdw)
              !
           endif
        enddo
     enddo
     !


     !F_0. T^0_ab :=   F_0 . (+ C^+_{a,dw}C_{b,dw})
     !F_z. T^z_ab :=   F_z . (- C^+_{a,dw}C_{b,dw})
     if(any(exc_field/=0d0))then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 idw = binary_search(Hsector%H(2)%map,k2)
                 !
                 htmp = exc_field(1)*sg1*sg2
                 call sp_insert_element(spH0dws(1),htmp,idw,jdw)
                 !
                 htmp = -exc_field(4)*sg1*sg2
                 call sp_insert_element(spH0dws(1),htmp,idw,jdw)
              endif
           enddo
        enddo
     endif



  enddo

