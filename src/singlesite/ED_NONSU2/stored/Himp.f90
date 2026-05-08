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

     !> H_Imp: Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + (impHloc(1,1,iorb,iorb) + mfHloc(1,1,iorb,iorb))*nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
        htmp = htmp + mfHloc(2,2,iorb,iorb)*ndw(iorb)
        htmp = htmp - xmu*(nup(iorb)+ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0,htmp,i,i)
     end select


     !Off-diagonal elements, i.e. non-local part
     !1. same spin:
     do iorb=1,Norb
        do jorb=1,Norb
           !this loop considers only the orbital off-diagonal terms
           !because iorb=jorb can not have simultaneously
           !occupation 0 and 1, as required by this if Jcondition:
           !UP
           Jcondition = &
                (impHloc(1,1,iorb,jorb)/=zero .OR. &
                 mfHloc(1,1,iorb,jorb)/=zero).AND. &
                (ib(jorb)==1) .AND. (ib(iorb)==0)
           if (Jcondition) then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j_el = binary_search(Hsector%H(1)%map,k2)
              j = j_el + (iph-1)*DimEl
              htmp = conjg(impHloc(1,1,iorb,jorb)+mfHloc(1,1,iorb,jorb))*sg1*sg2
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
           Jcondition = &
                (impHloc(Nspin,Nspin,iorb,jorb)/=zero .OR. &
                 mfHloc(2,2,iorb,jorb)/=zero).AND. &
                (ib(jorb+Ns)==1) .AND. (ib(iorb+Ns)==0)
           if (Jcondition) then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j_el = binary_search(Hsector%H(1)%map,k2)
              j = j_el + (iph-1)*DimEl
              htmp = conjg(impHloc(Nspin,Nspin,iorb,jorb))*sg1*sg2
              htmp = htmp + conjg(mfHloc(2,2,iorb,jorb))*sg1*sg2
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


     !2. spin-flip part (only for the nonSU2 channel!)
     do ispin=1,Nspin
        jspin = 3-ispin !ispin=1,jspin=2, ispin=2,jspin=1
        do iorb=1,Norb
           do jorb=1,Norb
              ialfa = iorb + (ispin-1)*Ns
              ibeta = jorb + (jspin-1)*Ns
              Jcondition=&
                   ((impHloc(ispin,jspin,iorb,jorb)/=zero .OR. &
                     mfHloc(ispin,jspin,iorb,jorb)/=zero) .AND. &
                   (ib(ibeta)==1) .AND. (ib(ialfa)==0))
              if(Jcondition)then
                 call c(ibeta,m,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 j_el = binary_search(Hsector%H(1)%map,k2)
                 j = j_el + (iph-1)*DimEl
                 htmp = conjg(impHloc(ispin,jspin,iorb,jorb)+mfHloc(ispin,jspin,iorb,jorb))*sg1*sg2
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
              endif
              !
           enddo
        enddo
     enddo



     !Evaluate: F.dot.T = sum_a={0,x,y,z} F_a.T^a
     if(any(exc_field/=0d0))then
        do iorb=1,Norb
           do jorb=1,Norb
              !
              !F_0. T^0_ab :=   F_0 . (C^+_{a,up}C_{b,up} + C^+_{a,dw}C_{b,dw})
              !F_z. T^z_ab :=   F_z . (C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw})
              !
              !C^+_{a,dw}C_{b,dw}
              Jcondition = (ib(jorb+Ns)==1) .AND. (ib(iorb+Ns)==0)
              if (Jcondition) then
                 call c(jorb+Ns,m,k1,sg1)
                 call cdg(iorb+Ns,k1,k2,sg2)
                 j_el = binary_search(Hsector%H(1)%map,k2)
                 j = j_el + (iph-1)*DimEl
                 !
                 htmp = exc_field(1)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
                 htmp = -exc_field(4)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              !
              !C^+_{a,up}C_{b,up}
              Jcondition = (ib(jorb)==1) .AND. (ib(iorb)==0)
              if (Jcondition) then
                 call c(jorb,m,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 j_el = binary_search(Hsector%H(1)%map,k2)
                 j = j_el + (iph-1)*DimEl
                 !
                 htmp = exc_field(1)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
                 htmp = exc_field(4)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              !
              !
              !F_x . T^x_ab:= F_x .   C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}
              !F_y . T^y_ab:= F_y .-i(C^+_{a,up}C_{b,dw} - C^+_{a,dw}C_{b,up})
              !
              !C^+_{a,dw}C_{b,up}
              Jcondition = (ib(jorb)==1) .AND. (ib(iorb+Ns)==0)
              if (Jcondition) then
                 call c(jorb,m,k1,sg1)
                 call cdg(iorb+Ns,k1,k2,sg2)
                 j_el = binary_search(Hsector%H(1)%map,k2)
                 j = j_el + (iph-1)*DimEl
                 !
                 htmp = exc_field(2)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
                 htmp = -xi*exc_field(3)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              !
              !C^+_{a,up}C_{b,dw}
              Jcondition = (ib(jorb+Ns)==1) .AND. (ib(iorb)==0)
              if (Jcondition) then
                 call c(jorb+Ns,m,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 j_el = binary_search(Hsector%H(1)%map,k2)
                 j = j_el + (iph-1)*DimEl
                 !
                 htmp = exc_field(2)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
                 !
                 htmp = xi*exc_field(3)*sg1*sg2
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              endif
              !
           enddo
        enddo
     endif





     !Evaluate: F.dot.S = sum_i={x,y,z} F_i.S^i     
     if(any(spin_field/=0d0))then
        !
        !F_z.S^z:= F_z.(C^+_{a,up}C_{a,up}-C^+_{a,dw}C_{a,dw})=F_z.n_up-n_dw
        do iorb=1,Norb
           htmp = htmp + spin_field(iorb,3)*(nup(iorb)-ndw(iorb))
        enddo
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0,htmp,i,i)
        end select
        !
        !
        do iorb=1,Norb
           !F_x . S^x:= F_x .   C^+_{a,up}C_{a,dw}  + C^+_{a,dw}C_{a,up}
           !F_y . S^y:= F_y . -i(C^+_{a,up}C_{a,dw} - C^+_{a,dw}C_{a,up})
           Jcondition = (ib(iorb)==1) .AND. (ib(iorb+Ns)==0)
           if (Jcondition) then
              call c(iorb,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j = binary_search(Hsector%H(1)%map,k2)
              !
              htmp = spin_field(iorb,1)*sg1*sg2
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
              !
              htmp = -xi*spin_field(iorb,2)*sg1*sg2
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
           endif
           !
           Jcondition = (ib(iorb+Ns)==1) .AND. (ib(iorb)==0)
           if (Jcondition) then
              call c(iorb+Ns,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j = binary_search(Hsector%H(1)%map,k2)
              !
              htmp = spin_field(iorb,1)*sg1*sg2
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
              !
              htmp = xi*spin_field(iorb,2)*sg1*sg2
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
           endif
           !
        enddo
     endif

     !
  enddo
