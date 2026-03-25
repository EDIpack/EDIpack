  do i=MpiIstart,MpiIend
     i_el = mod(i-1,DimEl) + 1
     iph = (i-1)/DimEl + 1
     m  = Hsector%H(1)%map(i_el)
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo

     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     select case(bath_type)
     case default 
        htmp=zero
        do iorb=1,size(dmft_bath%e,2)
           do kp=1,Nbath
              ialfa=getBathStride(iorb,kp)
              htmp =htmp + dmft_bath%e(1,iorb,kp)*ib(ialfa)        !UP
              htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ib(ialfa+Ns) !DW
           enddo
        enddo
        !
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0,htmp,i,i)
        end select
        !
     case("replica","general")
        htmp=zero
        do kp=1,Nbath
           do iorb=1,Norb
              ialfa = getBathStride(iorb,kp)
              htmp = htmp + bath_diag(1    ,iorb,kp)*ib(ialfa)    !UP
              htmp = htmp + bath_diag(Nspin,iorb,kp)*ib(ialfa+Ns) !DW
           enddo
        enddo
        !
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0,htmp,i,i)
        end select
        !
        !off-diagonal elements
        !1. same spin:
        do kp=1,Nbath
           do iorb=1,Norb
              do jorb=1,Norb
                 !UP
                 ialfa = getBathStride(iorb,kp)
                 ibeta = getBathStride(jorb,kp)
                 Jcondition = &
                      (hbath_tmp(1,1,iorb,jorb,kp)/=zero)   .AND. &
                      (ib(ibeta)==1) .AND. (ib(ialfa)==0)
                 if (Jcondition)then
                    call c(ibeta,m,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    j_el=binary_search(Hsector%H(1)%map,k2)
                    j = j_el * (iph-1)*DimEl
                    htmp = conjg(hbath_tmp(1,1,iorb,jorb,kp))*sg1*sg2
                    !
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0,htmp,i,j)
                    end select
                    !
                 endif
                 !DW
                 ialfa = getBathStride(iorb,kp) + Ns
                 ibeta = getBathStride(jorb,kp) + Ns
                 Jcondition = &
                      (hbath_tmp(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                      (ib(ibeta)==1).AND. (ib(ialfa)==0)
                 if (Jcondition)then
                    call c(ibeta,m,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    j_el=binary_search(Hsector%H(1)%map,k2)
                    j = j_el * (iph-1)*DimEl
                    htmp =conjg(hbath_tmp(Nspin,Nspin,iorb,jorb,kp))*sg1*sg2
                    !
                    select case(MpiStatus)
                    case (.true.)
                       call sp_insert_element(MpiComm,spH0,htmp,i,j)
                    case (.false.)
                       call sp_insert_element(spH0,htmp,i,j)
                    end select
                    !
                 endif
              enddo
           enddo
        enddo
        !
        !2. spin-flip part (only for the nonSU2 channel!)
        do kp=1,Nbath
           do ispin=1,Nspin
              jspin = 3-ispin !ispin=1,jspin=2, ispin=2,jspin=1
              do iorb=1,Norb
                 do jorb=1,Norb
                    ialfa = getBathStride(iorb,kp) + (ispin-1)*Ns
                    ibeta = getBathStride(jorb,kp) + (jspin-1)*Ns
                    Jcondition=&
                         (hbath_tmp(ispin,jspin,iorb,jorb,kp)/=zero) .AND. &
                         (ib(ibeta)==1) .AND. (ib(ialfa)==0)
                    if(Jcondition)then
                       call c(ibeta,m,k1,sg1)
                       call cdg(ialfa,k1,k2,sg2)
                       j_el=binary_search(Hsector%H(1)%map,k2)
                       j = j_el * (iph-1)*DimEl
                       htmp = conjg(hbath_tmp(ispin,jspin,iorb,jorb,kp))*sg1*sg2
                       !
                       select case(MpiStatus)
                       case (.true.)
                          call sp_insert_element(MpiComm,spH0,htmp,i,j)
                       case (.false.)
                          call sp_insert_element(spH0,htmp,i,j)
                       end select
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo

     end select

  enddo
