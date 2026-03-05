  !density-density interaction: same orbital, opposite spins:
  ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
  htmp = zero
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
  hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
  !


  ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
  !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
  !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
  if(Norb>1.AND.any((Jx_internal/=0d0)))then
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition=(&
                (iorb/=jorb).AND.&
                (ib(jorb)==1).AND.&
                (ib(iorb+Ns)==1).AND.&
                (ib(jorb+Ns)==0).AND.&
                (ib(iorb)==0))
           if(Jcondition)then
              call c(jorb,m,k1,sg1)
              call c(iorb+Ns,k1,k2,sg2)
              call cdg(jorb+Ns,k2,k3,sg3)
              call cdg(iorb,k3,k4,sg4)
              j_el = binary_search(Hsector%H(1)%map,k4)
              j    = j_el + (iph-1)*DimEl
              htmp = one*Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4
              !
              hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
              !
           endif
        enddo
     enddo
  endif
  !
  ! PAIR-HOPPING (P-H) TERMS
  !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
  !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
  if(Norb>1.AND.any((Jp_internal/=0d0)))then
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition=(&
                (iorb/=jorb).AND.&
                (ib(jorb)==1).AND.&
                (ib(jorb+Ns)==1).AND.&
                (ib(iorb+Ns)==0).AND.&
                (ib(iorb)==0))
           if(Jcondition)then
              call c(jorb,m,k1,sg1)
              call c(jorb+Ns,k1,k2,sg2)
              call cdg(iorb+Ns,k2,k3,sg3)
              call cdg(iorb,k3,k4,sg4)
              j_el = binary_search(Hsector%H(1)%map,k4)
              j    = j_el + (iph-1)*DimEl
              htmp = one*Jp_internal(iorb,jorb)*sg1*sg2*sg3*sg4
              !
              hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
              !
           endif
        enddo
     enddo
  endif
  
 !Sundry Coulomb terms
 if(allocated(coulomb_sundry))then
    do iline=1,size(coulomb_sundry)
       orbvec_dag  = [coulomb_sundry(iline)%cd_i(1), coulomb_sundry(iline)%cd_j(1)]
       orbvec      = [coulomb_sundry(iline)%c_k(1),  coulomb_sundry(iline)%c_l(1) ]
       spinvec_dag = [coulomb_sundry(iline)%cd_i(2), coulomb_sundry(iline)%cd_j(2)]
       spinvec     = [coulomb_sundry(iline)%c_k(2),  coulomb_sundry(iline)%c_l(2) ]
       
      
       !- we need to operate on either m_up or m_dw depending on the spins
       !- the operators need to be applied from right to left as c -> cd -> c -> cd to be consistent with the
       !  density terms originating from the anticommutation relations that are included in mfHloc
       
       Jcondition=.true. !Start applying operators
       !
       !last annihilation operator
       if(Jcondition)then 
         call c(orbvec(2) + Ns * (spinvec(2)-1), m ,k1 ,sg1 ,Jcondition)  !last annihilation operator
         if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
       endif
       !
       !last creation operator
       if(Jcondition)then 
           call cdg(orbvec_dag(2) + Ns * (spinvec_dag(2)-1), k1, k2, sg2, Jcondition)   !last annihilation operator
         if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
       endif
       !
       !first annihilation operator
       if(Jcondition)then 
         call c(orbvec(1) + Ns * (spinvec(1)-1), k2 ,k3 ,sg3 ,Jcondition)  !last annihilation operator
         if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
       endif
       !
       !last creation operator
       if(Jcondition)then 
           call cdg(orbvec_dag(1) + Ns * (spinvec_dag(1)-1), k3, k4, sg4, Jcondition)   !last annihilation operator
         if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
       endif
       !
       j_el = binary_search(Hsector%H(1)%map,k4)
       if (j_el == 0)then
         STOP "H_sundry: impossible operator"
       endif
       j    = j_el + (iph-1)*DimEl
       htmp = one*coulomb_sundry(iline)%U *sg1*sg2*sg3*sg4
       !
       hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
       !
      !
     enddo
  endif
