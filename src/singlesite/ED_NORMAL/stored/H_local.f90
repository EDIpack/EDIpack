  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     !
     !> H_Imp: Diagonal Elements, i.e. local part
     htmp = zero
     if(Nspin > 1 .and. (any(mfHloc(1,2,:,:)/=0d0) .or. any(mfHloc(2,1,:,:)/=0d0)))STOP "mfHloc cannot have spin-mixing in NORMAL mode"
     do iorb=1,Norb
        htmp = htmp + (impHloc(1,1,iorb,iorb) + mfHloc(1,1,iorb,iorb))*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb) * Ndw(iorb)
        htmp = htmp + mfHloc(2,2,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*( Nup(iorb)+Ndw(iorb) )
     enddo
     if(any(spin_field(:,3)/=0d0))then
        !F_z.S^z:= F_z.(n_up-n_dw)
        do iorb=1,Norb
           htmp = htmp + spin_field(iorb,3)*(nup(iorb)-ndw(iorb))
        enddo
     endif
     !
     !
     !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
     !
     !density-density interaction: same orbital, opposite spins:
     ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     do iorb=1,Norb
        htmp = htmp + Uloc_internal(iorb)*nup(iorb)*ndw(iorb)
     enddo
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc_internal-2*Jh_internal(iorb,jorb))*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust_internal(iorb,jorb)*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc_internal-3*Jh_internal(iorb,jorb))*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + (Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do iorb=1,Norb
           htmp = htmp - 0.5d0*Uloc_internal(iorb)*(nup(iorb)+ndw(iorb)) + 0.25d0*Uloc_internal(iorb)
        enddo
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 htmp=htmp-0.5d0*Ust_internal(iorb,jorb)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.5d0*Ust_internal(iorb,jorb)
                 htmp=htmp-0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))
              enddo
           enddo
        endif
     endif
     !
     !
     !> H_Bath: local bath energy contribution.
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     do iorb=1,size(bath_diag,2)
        do kp=1,Nbath
           ialfa = getBathStride(iorb,kp)
           htmp =htmp + bath_diag(1    ,iorb,kp)*Nup(ialfa) !UP
           htmp =htmp + bath_diag(Nspin,iorb,kp)*Ndw(ialfa) !DW
        enddo
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0d,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0d,htmp,i,i)
     end select
     !
  enddo

















