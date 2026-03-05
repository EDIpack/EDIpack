  do j=1,Nloc
     j_el = mod(j-1,DimUp*DimDw) + 1 !electronic index
     iph = (j-1)/(DimUp*DimDw) + 1   !phononic index
     !
     jup = iup_index(j_el,DimUp)
     jdw = idw_index(j_el,DimUp)
     !
     mup = Hsector%H(1)%map(jup)
     mdw = Hsector%H(2)%map(jdw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     !
     if(allocated(coulomb_sundry))then
        do iline=1,size(coulomb_sundry)
           orbvec_dag  = [coulomb_sundry(iline)%cd_i(1), coulomb_sundry(iline)%cd_j(1)]
           orbvec      = [coulomb_sundry(iline)%c_k(1),  coulomb_sundry(iline)%c_l(1) ]
           spinvec_dag = [coulomb_sundry(iline)%cd_i(2), coulomb_sundry(iline)%cd_j(2)]
           spinvec     = [coulomb_sundry(iline)%c_k(2),  coulomb_sundry(iline)%c_l(2) ]
           
           !This needs to happen:
           !- this is normal, so no terms that unbalance the spin count can exist (e.g. no c_up c_up cd_dw cd_dw)
           spinchange = 0
           !from right to left:
           !apply second annihilation operator
           spinchange = spinchange + (-1)**(spinvec(2)-1)
           !apply second creation operator
           spinchange = spinchange - (-1)**(spinvec_dag(2)-1)
           !apply first annihilation operator
           spinchange = spinchange + (-1)**(spinvec(1)-1)
           !apply second creation operator
           spinchange = spinchange - (-1)**(spinvec_dag(1)-1)
           !
           if(spinchange /= 0) STOP "In NORMAL mode, operators that change the total spin are forbidden. Check your umatrix file"
           
           !- we need to operate on either m_up or m_dw depending on the spins
           !- the operators need to be applied from right to left as c -> cd -> c -> cd to be consistent with the
           !  density terms originating from the anticommutation relations that are included in mfHloc
           
           Jcondition=.true. !Start applying operators
           p_up_new = mup
           p_dw_new = mdw
           !
           !last annihilation operator
           if(Jcondition)then 
             p_up_old = p_up_new
             p_dw_old = p_dw_new
             if (spinvec(2)==1)then
               call c(orbvec(2), p_up_old ,p_up_new ,sg1 ,Jcondition)  !last annihilation operator
             elseif (spinvec(2)==2)then
               call c(orbvec(2), p_dw_old ,p_dw_new ,sg1 ,Jcondition)  !last annihilation operator
             endif
             if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
           endif
           !
           !last creation operator
           if(Jcondition)then 
             p_up_old = p_up_new
             p_dw_old = p_dw_new
             if (spinvec_dag(2)==1)then
               call cdg(orbvec_dag(2), p_up_old, p_up_new, sg2, Jcondition)   !last annihilation operator
             elseif (spinvec_dag(2)==2)then
               call cdg(orbvec_dag(2), p_dw_old, p_dw_new, sg2, Jcondition)   !last annihilation operator
             endif
             if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
           endif
           !
           !first annihilation operator
           if(Jcondition)then 
             p_up_old = p_up_new
             p_dw_old = p_dw_new
             if (spinvec(1)==1)then
               call c(orbvec(1), p_up_old ,p_up_new ,sg3 ,Jcondition)  !last annihilation operator
             elseif (spinvec(1)==2)then
               call c(orbvec(1), p_dw_old ,p_dw_new ,sg3 ,Jcondition)  !last annihilation operator
             endif
             if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
           endif
           !
           !first creation operator
           if(Jcondition)then 
             p_up_old = p_up_new
             p_dw_old = p_dw_new
             if (spinvec_dag(1)==1)then
               call cdg(orbvec_dag(1), p_up_old ,p_up_new ,sg4 ,Jcondition)  !last annihilation operator
             elseif (spinvec_dag(1)==2)then
               call cdg(orbvec_dag(1), p_dw_old ,p_dw_new ,sg4 ,Jcondition)  !last annihilation operator
             endif
             if (.not. Jcondition) cycle                 !this gives zero, no hamiltonian element added
           endif
           !
           idw=binary_search(Hsector%H(2)%map,p_dw_new)
           iup=binary_search(Hsector%H(1)%map,p_up_new)
           
           if (idw == 0 .or. iup == 0)then
              STOP "H_sundry: impossible operator"
           endif
           
           htmp = coulomb_sundry(iline)%U * sg1 * sg2 * sg3 * sg4
           !
           i = iup + (idw-1)*dimup + (iph-1)*DimUp*DimDw
           !
           Hv(j) = Hv(j) + htmp*vin(i)
          !
         enddo
      endif
     !
  enddo
