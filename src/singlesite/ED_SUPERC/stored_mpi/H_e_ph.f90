 !We build the electronic part of the electron-phonon interaction:
  do i=MpiIstart,MpiIend
     i_el = mod(i-1,DimEl) + 1
     iph = (i-1)/DimEl + 1
     !
     m  = Hsector%H(1)%map(i_el)
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     ! Diagonal terms: Sum_iorb g_iorb,iorb n_iorb
     htmp_el = zero
     do iorb=1,Norb
        htmp_el = htmp_el + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
     enddo
     !
     if(iph < DimPh) then ! N* bdg|iph>
        htmp = sqrt(dble(iph))*htmp_el
        j = i_el + (iph)*DimEl
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0,htmp,i,j)
        case (.false.)
           call sp_insert_element(spH0,htmp,i,j)
        end select
     end if
     if(iph>1) then ! N* b|iph>
        htmp = sqrt(dble(iph - 1))*htmp_el
        j = i_el + (iph-2)*DimEl
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0,htmp,i,j)
        case (.false.)
           call sp_insert_element(spH0,htmp,i,j)
        end select
     end if
     !
     ! Off-Diagonal terms: Sum_iorb,jorb g_iorb,jorb cdg_iorb*c_jorb
     ! UP spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (g_ph(iorb,jorb)/=zero) .AND. &
                (ib(jorb)==1) .AND. (ib(iorb)==0)
           if(Jcondition)then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j_el = binary_search(Hsector%H(1)%map,k2)
              htmp_el = conjg(g_ph(iorb,jorb))*sg1*sg2
              if(iph < DimPh) then ! bdg|iph>
                 htmp = sqrt(dble(iph))*htmp_el
                 j = j_el + (iph)*DimEl
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              end if
              if(iph>1) then ! b|iph>
                 htmp = sqrt(dble(iph - 1))*htmp_el
                 j = j_el + (iph-2)*DimEl
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              end if
              !
           endif
        enddo
     enddo
     ! DW spin
     ! remark: iorb=jorb can't have simultaneously n=0 and n=1 (Jcondition)
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition = &
                (g_ph(iorb,jorb)/=zero) .AND. &
                (ib(jorb+Ns)==1) .AND. (ib(iorb+Ns)==0)
           if (Jcondition) then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j_el  = binary_search(Hsector%H(1)%map,k2)
              j     = j_el + (iph-1)*DimEl
              htmp_el = conjg(g_ph(iorb,jorb))*sg1*sg2
              !
              if(iph < DimPh) then ! bdg|iph>
                 htmp = sqrt(dble(iph))*htmp_el
                 j = j_el + (iph)*DimEl
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              end if
              if(iph>1) then ! b|iph>
                 htmp = sqrt(dble(iph - 1))*htmp_el
                 j = j_el + (iph-2)*DimEl
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0,htmp,i,j)
                 end select
              end if
              !
           endif
              
        enddo
     enddo
     !
  enddo
