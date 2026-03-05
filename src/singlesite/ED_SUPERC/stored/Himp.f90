  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
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
        htmp = htmp + (impHloc(Nspin,Nspin,iorb,iorb))*ndw(iorb)
        htmp = htmp + (mfHloc(2,2,iorb,iorb))*ndw(iorb) !Needed because these terms come from the anticommutator
        htmp = htmp - xmu*(nup(iorb)+ndw(iorb))
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0,htmp,i,i)
     end select
     !

     !Off-diagonal elements, i.e. non-local part
     do iorb=1,Norb
        do jorb=1,Norb
           !this loop considers only the orbital off-diagonal terms
           !because iorb=jorb can not have simultaneously
           !occupation 0 and 1, as required by this if Jcondition:
           !UP
           Jcondition = &
                (impHloc(1,1,iorb,jorb)/=zero .OR. &
                 mfHloc(1,1,iorb,jorb)/=zero) .AND. &
                (ib(jorb)==1)                 .AND. &
                (ib(iorb)==0)
           if (Jcondition) then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j    = binary_search(Hsector%H(1)%map,k2)
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
           !DW: 
           Jcondition = &
                (impHloc(Nspin,Nspin,iorb,jorb)/=zero  .OR. &
                 mfHloc(2,2,iorb,jorb)/=zero) .AND. &
                (ib(jorb+Ns)==1)                      .AND. &
                (ib(iorb+Ns)==0)
           if (Jcondition) then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j    = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(impHloc(Nspin,Nspin,iorb,jorb))*sg1*sg2
              htmp = htmp + conjg(mfHloc(2,2,iorb,jorb))*sg1*sg2  !Needed because these terms come from the anticommutator
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
     !



     !Evaluate: Fd . D = Fd . (C^+_{a,up}C^+_{b,dw} + C_{b,dw}C_{a,up})
     if(any(pair_field/=0d0).or.any(abs(impHloc_anomalous)/=0d0))then
        do iorb=1,Norb
          do jorb = 1,Norb
            !
            Jcondition = (ib(iorb)==1) .AND. (ib(jorb+Ns)==1)
            if(Jcondition)then
               call c(iorb,m,k1,sg1)
               call c(jorb+Ns,k1,k2,sg2)
               j=binary_search(Hsector%H(1)%map,k2)
               htmp=one*impHloc_anomalous(1,1,iorb,jorb)*sg1*sg2
               if(iorb == jorb) htmp= htmp + one*pair_field(iorb)*sg1*sg2
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
            Jcondition = (ib(iorb)==0) .AND. (ib(jorb+Ns)==0)
            if(Jcondition)then
               call cdg(jorb+Ns,m,k1,sg1)
               call cdg(iorb,k1,k2,sg2)
               j=binary_search(Hsector%H(1)%map,k2)
               htmp=one*conjg(impHloc_anomalous(1,1,iorb,jorb))*sg1*sg2
               if(iorb == jorb) htmp = htmp + one*pair_field(iorb)*sg1*sg2
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
     endif

  enddo

