  do iph=1,DimPh
     do jdw=1,MpiQdw
        do jup=1,DimUp
           mup  = Hsector%H(1)%map(jup)
           nup  = bdecomp(mup,Ns)
           !
           j    = jup + (jdw-1)*dimUp + (iph-1)*DimUp*MpiQdw
           !
           !> H_imp: Off-diagonal elements, i.e. non-local part. 
           !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
           do iorb=1,Norb
              do jorb=1,Norb
                 Jcondition = &
                      ((impHloc(1,1,iorb,jorb)+mfHloc(1,1,iorb,jorb))/=zero) .AND. &
                      (nup(jorb)==1) .AND. (nup(iorb)==0)
                 if (Jcondition) then
                    call c(jorb,mup,k1,sg1)
                    call cdg(iorb,k1,k2,sg2)
                    iup  = binary_search(Hsector%H(1)%map,k2)
                    idw  = jdw
                    i    = iup + (idw-1)*dimUp + (iph-1)*DimUp*MpiQdw
                    htmp = one*(impHloc(1,1,iorb,jorb)+mfHloc(1,1,iorb,jorb))*sg1*sg2
                    !
                    Hv(j) = Hv(j) + htmp*vin(i)
                    ! Hv(i) = Hv(i) + htmp*vin(j)
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
                            (hbath_tmp(1,1,iorb,jorb,kp)/=zero) .AND.&
                            (nup(ibeta)==1) .AND. (nup(ialfa)==0)
                       !
                       if (Jcondition)then
                          call c(ibeta,mup,k1,sg1)
                          call cdg(ialfa,k1,k2,sg2)
                          iup  = binary_search(Hsector%H(1)%map,k2)
                          idw  = jdw
                          i    = iup + (idw-1)*dimUp + (iph-1)*DimUp*MpiQdw
                          htmp = one*hbath_tmp(1,1,iorb,jorb,kp)*sg1*sg2
                          !
                          Hv(j) = Hv(j) + htmp*vin(i)
                          ! Hv(i) = Hv(i) + htmp*vin(j)
                          !
                       endif
                    enddo
                 enddo
              enddo
           end if
           !
           !
           !>H_hyb: hopping terms for a given spin (imp <--> bath)
           do iorb=1,Norb
              do kp=1,Nbath
                 ialfa=getBathStride(iorb,kp)
                 !
                 if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                      (nup(iorb)==1) .AND. (nup(ialfa)==0) )then              
                    call c(iorb,mup,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    iup  = binary_search(Hsector%H(1)%map,k2)
                    idw  = jdw
                    i    = iup + (idw-1)*dimUp + (iph-1)*DimUp*MpiQdw
                    htmp = one*diag_hybr(1,iorb,kp)*sg1*sg2
                    !
                    Hv(j) = Hv(j) + htmp*vin(i)
                    ! Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
                 !
                 if( (diag_hybr(1,iorb,kp)/=0d0) .AND. &
                      (nup(iorb)==0) .AND. (nup(ialfa)==1) )then
                    call c(ialfa,mup,k1,sg1)
                    call cdg(iorb,k1,k2,sg2)
                    iup  = binary_search(Hsector%H(1)%map,k2)
                    idw  = jdw
                    i    = iup + (idw-1)*dimUp + (iph-1)*DimUp*MpiQdw
                    htmp = one*diag_hybr(1,iorb,kp)*sg1*sg2
                    !
                    Hv(j) = Hv(j) + htmp*vin(i)
                    ! Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
           !
           !
           !F_0. T^0_ab :=   F_0 . C^+_{a,up}C_{b,up}
           !F_z. T^z_ab :=   F_z . C^+_{a,up}C_{b,up}
           if(any(exc_field/=0d0))then
              do iorb=1,Norb
                 do jorb=iorb+1,Norb
                    Jcondition = (Nup(jorb)==1) .AND. (Nup(iorb)==0)
                    if (Jcondition) then
                       call c(jorb,mup,k1,sg1)
                       call cdg(iorb,k1,k2,sg2)
                       iup = binary_search(Hsector%H(1)%map,k2)
                       idw = jdw 
                       i    = iup + (idw-1)*dimUp + (iph-1)*DimUp*MpiQdw
                       !
                       htmp = one*exc_field(1)*sg1*sg2
                       Hv(j) = Hv(j) + htmp*vin(i)
                       ! Hv(i) = Hv(i) + htmp*vin(j)
                       !
                       htmp = one*exc_field(4)*sg1*sg2
                       Hv(j) = Hv(j) + htmp*vin(i)
                       ! Hv(i) = Hv(i) + htmp*vin(j)
                    endif
                 enddo
              enddo
           endif

        enddo
     enddo
  enddo


