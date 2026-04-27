MODULE ED_OBSERVABLES_NORMAL
  !:synopsis: Routines for observable calculation, :code:`NORMAL` case
  !This module calculates a series of observables, and stores them in aptly named plain-text files. :f:var:`ed_mode` = :code:`normal`
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_BATH
  USE ED_HAMILTONIAN_NORMAL
  !
  implicit none
  private
  !
  public :: observables_normal
  public :: local_energy_normal

  real(8),dimension(:),allocatable   :: dens ! orbital-resolved charge density
  real(8),dimension(:),allocatable   :: dens_up ! orbital-resolved spin-:math:`\uparrow` electron density
  real(8),dimension(:),allocatable   :: dens_dw ! orbital-resolved spin-:math:`\downarrow` electron density
  real(8),dimension(:),allocatable   :: docc ! orbital-resolved double occupation
  real(8),dimension(:),allocatable   :: magz ! orbital-resolved magnetization ( :code:`z` component )
  real(8),dimension(:,:),allocatable :: n2 ! :math:`\langle n_{i} n_{j} \rangle` for i,j orbitals
  real(8),dimension(:,:),allocatable :: sz2! :math:`\langle S^{z}_{i} S^{z}_{j} \rangle` for i,j orbitals
  real(8),dimension(:,:),allocatable :: exct_s0 !excitonic order parameter :math:`\langle c^{\dagger}_{is}\sigma^{0}c_{js^{'}} \rangle`
  real(8),dimension(:,:),allocatable :: exct_tz !excitonic order parameter :math:`\langle c^{\dagger}_{is}\sigma^{z}c_{js^{'}} \rangle`
  real(8)                            :: dens_ph !phonon density
  real(8)                            :: X_ph ! :math:`\langle X \rangle` with :math:`X = \frac{b + b^{dagger}}{2}`
  real(8)                            :: X2_ph ! :math:`\langle X^{2} \rangle` with :math:`X = \frac{b + b^{dagger}}{2}`
  real(8)                            :: s2tot ! :math:`\langle S_{z}^{2} \rangle`
  real(8)                            :: Egs ! Ground-state energy
  real(8)                            :: Ei
  real(8),dimension(:),allocatable   :: Prob ! Probability distribution
  real(8),dimension(:),allocatable   :: prob_ph ! Phonon probability
  real(8),dimension(:),allocatable   :: pdf_ph ! Phonon probability distribution :f:func:`prob_distr_ph`
  real(8),dimension(:,:),allocatable :: pdf_part ! Lattice probability distribution as obtained by :f:func:`prob_distr_ph`
  !
  integer                            :: iorb,jorb,iorb1,jorb1
  integer                            :: ispin,jspin
  integer                            :: isite,jsite
  integer                            :: ibath
  integer                            :: r,m,k,k1,k2,k3,k4
  integer                            :: iup,idw
  integer                            :: jup,jdw
  integer                            :: mup,mdw
  integer                            :: iph,i_el,j_el,isectorDim
  real(8)                            :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                            :: gs_weight
  !
  real(8)                            :: peso
  real(8)                            :: norm
  !
  integer                            :: i,j,ii
  integer                            :: isector,jsector
  !
  real(8),dimension(:),allocatable   :: vvinit,v1,v2
  real(8),dimension(:),allocatable   :: v_state
  logical                            :: Jcondition
  !
  type(sector)                       :: sectorI,sectorJ


contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine observables_normal()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Calculate the values of the local observables
    integer                         :: iprob,istate,Nud(2,Ns),iud(2),jud(2),val
    integer,dimension(2*Ns_Ud)      :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb) :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)           :: IbUp,IbDw  ![Ns]
    real(8),dimension(Norb)         :: nup,ndw,Sz,nt
    real(8),dimension(Norb,Norb)    :: theta_upup,theta_dwdw
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG observables_normal"
#endif
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(exct_S0(Norb,Norb),exct_Tz(Norb,Norb))
    allocate(Prob(3**Norb))
    allocate(prob_ph(DimPh))
    allocate(pdf_ph(Lpos))
    allocate(pdf_part(Lpos,3**Norb)) ! each orbital can be empty '0', single occupied '1' and doubly occupied '2'
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    exct_s0 = 0d0
    exct_tz = 0d0
    theta_upup = 0d0
    theta_dwdw = 0d0
    Prob    = 0.d0
    prob_ph = 0.d0
    dens_ph = 0.d0
    X_ph = 0.d0
    X2_ph= 0.d0
    pdf_ph  = 0.d0
    pdf_part= 0.d0

    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_normal: get local observables"
#endif
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state =  es_return_dvec(state_list,istate)
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_normal: get contribution from state:"//str(istate)
#endif
       !
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             gs_weight=peso*abs(v_state(i))**2
             call build_op_Ns(i,IbUp,IbDw,sectorI)
             nup = IbUp(1:Norb)
             ndw = IbDw(1:Norb)
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !Configuration probability
             iprob=1
             do iorb=1,Norb
                iprob=iprob+nint(nt(iorb))*3**(iorb-1)
             end do
             Prob(iprob) = Prob(iprob) + gs_weight
             !
             !Evaluate averages of observables:
             do iorb=1,Norb
                dens(iorb)     = dens(iorb)      +  nt(iorb)*gs_weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*gs_weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*gs_weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*gs_weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*gs_weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight
                do jorb=iorb+1,Norb
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
             !
             iph = (i-1)/(sectorI%DimEl) + 1
             i_el = mod(i-1,sectorI%DimEl) + 1
             prob_ph(iph) = prob_ph(iph) + gs_weight
             dens_ph = dens_ph + (iph-1)*gs_weight
             !
             !<X> and <X^2> with X=(b+bdg)/sqrt(2)
             if(iph<DimPh)then
               ! only bdag
                j= i_el + (iph)*sectorI%DimEl
                X_ph = X_ph + sqrt(2.d0*dble(iph))*(v_state(i)*v_state(j))*peso
             end if
             ! bdag*b+b*bdag = 1+2*bdag*b  contribution to X^2
             X2_ph = X2_ph + 0.5d0*(1+2*(iph-1))*gs_weight
             if(iph<DimPh-1)then
                ! only bdag*bdag contribution to X^2
                j= i_el + (iph+1)*sectorI%DimEl
                X2_ph = X2_ph + sqrt(dble((iph)*(iph+1)))*(v_state(i)*v_state(j))*peso
             end if
             !
             !compute the lattice probability distribution function
             if(Dimph>1 .AND. iph==1) then
                val = 1
                ! val encodes the occupation of each orbital in basis 3 (since each orbital can have occupation 0,1,2)
                ! the ternary representation of 'val-1' gives the occupation of each orbital
                ! val-1 mod 3 is the first orbital, (val-1)/3 mod 3 is the second and so on...
                ! 
                ! OLD : val= 1 + Nr. of polarized orbitals (full or empty) makes sense only for 2 orbs
                do iorb=1,Norb
                   !val = val + abs(nint(sign((nt(iorb) - 1.d0),real(g_ph(iorb,iorb)))))
                   val = val + nt(iorb)*( 3**(iorb-1) )
                enddo
                call prob_distr_ph(v_state,val)
             end if
          enddo
          !
          call delete_sector(sectorI)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
    !
    !EVALUATE EXCITON OP <S_ab> AND <T^z_ab>
    !<S_ab>  :=   <C^+_{a,up}C_{b,up} + C^+_{a,dw}C_{b,dw}>
    !<T^z_ab>:=   <C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw}>
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_normal: eval exciton OP Singlet, Triplet_Z"
#endif
    if (ed_total_ud) then
      do istate=1,state_list%size
         isector = es_return_sector(state_list,istate)
         Ei      = es_return_energy(state_list,istate)
         v_state  =  es_return_dvec(state_list,istate)
#ifdef _DEBUG
         if(ed_verbose>3)write(Logfile,"(A)")&
              "DEBUG observables_normal: get contribution from state:"//str(istate)
#endif

         !
         peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
         peso = peso/zeta_function
         !
         do iorb=1,Norb
            do jorb=iorb+1,Norb
               !
               !\Theta_upup = <v|v>, |v> = (C_aup + C_bup)|>
               jsector = getCsector(1,1,isector)
               if(jsector/=0)then
                  vvinit =  apply_Cops(v_state,[1d0,1d0],[-1,-1],[iorb,jorb],[1,1],isector,jsector)
                  if(Mpimaster)then
                     theta_upup(iorb,jorb) = theta_upup(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                  endif
                  if(allocated(vvinit))deallocate(vvinit)
               endif
               !
               !\Theta_dwdw = <v|v>, |v> = (C_adw + C_bdw)|>
               jsector = getCsector(1,2,isector)
               if(jsector/=0)then
                  vvinit =  apply_Cops(v_state,[1d0,1d0],[-1,-1],[iorb,jorb],[2,2],isector,jsector)
                  if(Mpimaster)then
                     theta_dwdw(iorb,jorb) = theta_dwdw(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                  endif
                  if(allocated(vvinit))deallocate(vvinit)
               endif
               !
            enddo
         enddo
         !
         if(allocated(v_state))deallocate(v_state)
         !
      enddo
      !
      do iorb=1,Norb
         do jorb=iorb+1,Norb
            exct_s0(iorb,jorb) = 0.5d0*(theta_upup(iorb,jorb) + theta_dwdw(iorb,jorb) - dens(iorb) - dens(jorb))
            exct_tz(iorb,jorb) = 0.5d0*(theta_upup(iorb,jorb) - theta_dwdw(iorb,jorb) - magZ(iorb) - magZ(jorb))
         enddo
      enddo
#ifdef _DEBUG
      if(ed_verbose>2)write(Logfile,"(A)")""
#endif
  else
    write(Logfile,"(A)")"Cannot calculate exciton operators with ed_total_ud=F"
  endif

    !
    !SINGLE PARTICLE IMPURITY DENSITY MATRIX
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_normal: eval single particle density matrix <C^+_a C_b>"
#endif
    if(allocated(single_particle_density_matrix)) deallocate(single_particle_density_matrix)
    allocate(single_particle_density_matrix(Nspin,Nspin,Norb,Norb));single_particle_density_matrix=zero
    if(allocated(full_denmat)) deallocate(full_denmat)
    allocate(full_denmat(Nspin,Nspin,Ns,Ns));full_denmat=zero
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state =  es_return_dvec(state_list,istate)              
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_normal: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             iph = (i-1)/(sectorI%DimEl) + 1
             i_el = mod(i-1,sectorI%DimEl) + 1
             call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
             !
             call build_op_Ns(i,IbUp,IbDw,sectorI)
             Nud(1,:)=IbUp
             Nud(2,:)=IbDw
             !
             !Diagonal densities
             do ispin=1,Nspin
                do iorb=1,Norb
                   single_particle_density_matrix(ispin,ispin,iorb,iorb) = &
                        single_particle_density_matrix(ispin,ispin,iorb,iorb) + &
                        peso*nud(ispin,iorb)*(v_state(i))*v_state(i)
                enddo
             enddo
             !
             !Off-diagonal
             if(ed_total_ud)then
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         !
                         if((Nud(ispin,jorb)==1).and.(Nud(ispin,iorb)==0))then
                            iud(1) = sectorI%H(1)%map(Indices(1))
                            iud(2) = sectorI%H(2)%map(Indices(2))
                            call c(jorb,iud(ispin),r,sgn1)
                            call cdg(iorb,r,k,sgn2)
                            Jndices = Indices
                            Jndices(1+(ispin-1)*Ns_Ud) = &
                                 binary_search(sectorI%H(1+(ispin-1)*Ns_Ud)%map,k)
                            call indices2state(Jndices,[sectorI%DimUps,sectorI%DimDws],j)
                            !
                            j = j + (iph-1)*sectorI%DimEl
                            !
                            single_particle_density_matrix(ispin,ispin,iorb,jorb) = &
                                 single_particle_density_matrix(ispin,ispin,iorb,jorb) + &
                                 peso*sgn1*v_state(i)*sgn2*(v_state(j))
                         endif
                      enddo
                   enddo
                enddo
             endif
             do ispin=1,Nspin
               do iorb=1,Ns
                  do jorb=1,Ns
                  ! diagonal
                  if(iorb==jorb )then
                     full_denmat(ispin,ispin,iorb,jorb) = &
                     full_denmat(ispin,ispin,iorb,jorb) + Nud(ispin,iorb)*peso*(v_state(i))*v_state(i)
                  elseif((Nud(ispin,jorb)==1).and.(Nud(ispin,iorb)==0))then
                     iud(1) = sectorI%H(1)%map(Indices(1))
                     iud(2) = sectorI%H(2)%map(Indices(2))
                     call c(jorb,iud(ispin),r,sgn1)
                     call cdg(iorb,r,k,sgn2)
                     Jndices = Indices
                     Jndices(1+(ispin-1)*Ns_Ud) = &
                           binary_search(sectorI%H(1+(ispin-1)*Ns_Ud)%map,k)
                     call indices2state(Jndices,[sectorI%DimUps,sectorI%DimDws],j)
                     !
                     j = j + (iph-1)*sectorI%DimEl
                     !
                     full_denmat(ispin,ispin,iorb,jorb) = &
                        full_denmat(ispin,ispin,iorb,jorb) + peso*sgn1*v_state(i)*sgn2*(v_state(j))
                  endif
                  enddo
               enddo
             enddo
             !
             !
          enddo
          call delete_sector(sectorI)         
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
    !
    !
    write(LOGfile,"(A,10f18.12,f18.12)")&
         " dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    write(LOGfile,"(A,10f18.12)")&
         " docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    if(Nspin==2)then
       write(LOGfile,"(A,10f18.12,A)")&
            " magZ"//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
    endif
    if(any(abs(exct_S0)>1d-9))write(LOGfile,"(A,10f18.12)")&
         "excS0"//reg(ed_file_suffix)//"=",((exct_S0(iorb,jorb),iorb=1,Norb),jorb=1,Norb)
    if(any(abs(exct_tz)>1d-9))write(LOGfile,"(A,10f18.12)")&
         "excTZ"//reg(ed_file_suffix)//"=",((exct_Tz(iorb,jorb),iorb=1,Norb),jorb=1,Norb)
    !
    if(DimPh>1)then
       call write_pdf()
    endif
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_mag(3,iorb)  =magZ(iorb)
       ed_docc(iorb)   =docc(iorb)
       do jorb=1,Norb
          ed_exct(1,iorb,jorb) = exct_S0(iorb,jorb)
          ed_exct(4,iorb,jorb) = exct_Tz(iorb,jorb)
       enddo
    enddo
    !
    ed_imp_info=[s2tot,egs]
    ed_dens_ph = dens_ph
    ed_X_ph = X_ph
    ed_X2_ph = X2_ph
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_mag)
       call Bcast_MPI(MpiComm,ed_exct)
       call Bcast_MPI(MpiComm,ed_imp_info)
       if(allocated(single_particle_density_matrix))call Bcast_MPI(MpiComm,single_particle_density_matrix)
       if(allocated(full_denmat))call Bcast_MPI(MpiComm,full_denmat)
    endif
#endif
    !
    if(MPIMASTER)then
       call write_observables()
    endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2,Prob)
    deallocate(exct_S0,exct_Tz)
    deallocate(prob_ph,pdf_ph,pdf_part)
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
  end subroutine observables_normal





  !+------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_normal()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Calculate the value of the local energy components
    integer                             :: istate,iud(2),jud(2),iline, p_up_new, p_up_old, p_dw_new, p_dw_old
    integer,dimension(2)                :: orbvec, orbvec_dag, spinvec, spinvec_dag
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    real(8),dimension(Ns)               :: Nup,Ndw
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG local_energy_normal"
#endif
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    ed_Eeph    = 0.d0
    ed_Eph    = 0.d0
    !
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state =  es_return_dvec(state_list,istate)
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG local_energy_normal: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       !Master:
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             iph = (i-1)/(sectorI%DimEl) + 1
             i_el = mod(i-1,sectorI%DimEl) + 1
             !
             call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
             do ii=1,Ns_Ud
                mup = sectorI%H(ii)%map(Indices(ii))
                mdw = sectorI%H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
                Nups(ii,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
                Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
             enddo
             Nup = Breorder(Nups)
             Ndw = Breorder(Ndws)
             !
             gs_weight=peso*abs(v_state(i))**2
             !
             !> H_Imp: Diagonal Elements, i.e. local part
             do iorb=1,Norb
                ed_Eknot = ed_Eknot + (impHloc(1,1,iorb,iorb))*Nup(iorb)*gs_weight
                ed_Eknot = ed_Eknot + (impHloc(Nspin,Nspin,iorb,iorb))*Ndw(iorb)*gs_weight
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part. 
             if(ed_total_ud)then
                iup = Indices(1)
                idw = Indices(2)
                mup = sectorI%H(1)%map(iup)
                mdw = sectorI%H(2)%map(idw)
                do iorb=1,Norb
                   do jorb=1,Norb
                      !UP
                      Jcondition = &
                           (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                           (Nup(jorb)==1) .AND. (Nup(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mup,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jup = binary_search(sectorI%H(1)%map,k2)
                         j_el   = jup + (idw-1)*sectorI%DimUp
                         j = j_el + (iph-1)*sectorI%DimEl
                         ed_Eknot = ed_Eknot + &
                              (impHloc(1,1,iorb,jorb))*sg1*sg2*v_state(i)*(v_state(j))*peso
                      endif
                      !
                      !DW
                      Jcondition = &
                           (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                           (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mdw,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jdw = binary_search(sectorI%H(2)%map,k2)
                         j_el   = iup + (jdw-1)*sectorI%DimUp
                         j = j_el + (iph-1)*sectorI%DimEl
                         ed_Eknot = ed_Eknot + &
                              (impHloc(Nspin,Nspin,iorb,jorb))*sg1*sg2*v_state(i)*(v_state(j))*peso
                      endif
                   enddo
                enddo
                ! 
                !SPIN-EXCHANGE Jx
                if(Norb>1 .AND. any((Jx_internal /= 0.0d0)))then
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
                            jdw=binary_search(sectorI%H(2)%map,k2)
                            call c(jorb,mup,k3,sg3)  !UP
                            call cdg(iorb,k3,k4,sg4) !UP
                            jup=binary_search(sectorI%H(1)%map,k4)
                            j_el   = jup + (idw-1)*sectorI%DimUp
                            j = j_el + (iph-1)*sectorI%DimEl
                            !
                            ed_Epot = ed_Epot + Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4*v_state(i)*v_state(j)*peso
                            ed_Dse = ed_Dse + sg1*sg2*sg3*sg4*v_state(i)*v_state(j)*peso
                            !
                         endif
                      enddo
                   enddo
                endif
                !
                ! PAIR-HOPPING Jp
                if(Norb>1 .AND. any((Jp_internal /= 0.0d0)))then
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
                            jdw = binary_search(sectorI%H(2)%map,k2)
                            call c(jorb,mup,k3,sg3)       !c_jorb_up
                            call cdg(iorb,k3,k4,sg4)      !c^+_iorb_up
                            jup = binary_search(sectorI%H(1)%map,k4)
                            j_el   = jup + (jdw-1)*sectorI%DimUp
                            j = j_el + (iph-1)*sectorI%DimEl
                            !
                            ed_Epot = ed_Epot + Jp_internal(iorb,jorb)*sg1*sg2*sg3*sg4*v_state(i)*v_state(j)*peso
                            ed_Dph = ed_Dph + sg1*sg2*sg3*sg4*v_state(i)*v_state(j)*peso
                            !
                         endif
                      enddo
                   enddo
                endif
                !
                ! Mean-field from sundry terms
                do iorb=1,Norb
                  ed_Epot = ed_Epot + (mfHloc(1,1,iorb,iorb))*Nup(iorb)*gs_weight
                  ed_Epot = ed_Epot + (mfHloc(2,2,iorb,iorb))*Ndw(iorb)*gs_weight
                  do jorb=1,Norb
                      !UP
                      Jcondition = &
                           (mfHloc(1,1,iorb,jorb)/=zero) .AND. &
                           (Nup(jorb)==1) .AND. (Nup(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mup,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jup = binary_search(sectorI%H(1)%map,k2)
                         j_el   = jup + (idw-1)*sectorI%DimUp
                         j = j_el + (iph-1)*sectorI%DimEl
                         ed_Epot = ed_Epot + mfHloc(1,1,iorb,jorb)*sg1*sg2*v_state(i)*(v_state(j))*peso
                      endif
                      !
                      !DW
                      Jcondition = &
                           (mfHloc(2,2,iorb,jorb)/=zero) .AND. &
                           (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                      if (Jcondition) then
                         call c(jorb,mdw,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jdw = binary_search(sectorI%H(2)%map,k2)
                         j_el   = iup + (jdw-1)*sectorI%DimUp
                         j = j_el + (iph-1)*sectorI%DimEl
                         ed_Epot = ed_Epot + mfHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*v_state(i)*(v_state(j))*peso
                      endif
                   enddo
                enddo                
                
                ! SUNDRY TERMS
                if(allocated(coulomb_sundry))then
                  do iline=1,size(coulomb_sundry)
                       orbvec_dag  = [coulomb_sundry(iline)%cd_i(1), coulomb_sundry(iline)%cd_j(1)]
                       orbvec      = [coulomb_sundry(iline)%c_k(1),  coulomb_sundry(iline)%c_l(1) ]
                       spinvec_dag = [coulomb_sundry(iline)%cd_i(2), coulomb_sundry(iline)%cd_j(2)]
                       spinvec     = [coulomb_sundry(iline)%c_k(2),  coulomb_sundry(iline)%c_l(2) ]
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
                       !
                       jdw=binary_search(sectorI%H(2)%map, p_dw_new)
                       jup=binary_search(sectorI%H(1)%map, p_up_new)
                       
                       if (jup == 0 .or. jdw == 0)then
                          STOP "H_sundry: impossible operator"
                       endif
                       
                       j_el = jup + (jdw-1)*sectorI%DimUp
                       j = j_el + (iph-1)*sectorI%DimEl
                       !
                       ed_Epot = ed_Epot + coulomb_sundry(iline)%U*sg1*sg2*sg3*sg4*v_state(i)*v_state(j)*peso
                       !
                  enddo
                endif
             endif
             !
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
             do iorb=1,Norb
                ed_Epot = ed_Epot + Uloc_internal(iorb)*nup(iorb)*ndw(iorb)*gs_weight
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + Ust_internal(iorb,jorb)*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                      ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + (Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                      ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
                do iorb=1,Norb
                   ed_Ehartree=ed_Ehartree - 0.5d0*Uloc_internal(iorb)*(nup(iorb)+ndw(iorb))*gs_weight + 0.25d0*Uloc_internal(iorb)*gs_weight
                enddo
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         ed_Ehartree=ed_Ehartree - 0.5d0*Ust_internal(iorb,jorb)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.5d0*Ust_internal(iorb,jorb)*gs_weight
                         ed_Ehartree=ed_Ehartree - 0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*gs_weight
                      enddo
                   enddo
                endif
             endif
             !
             ! Electron-Phonon coupling term
             if(DimPh>1)then 
                ed_Eph = ed_Eph + (iph-1)*w0_ph*gs_weight
                !Diagonal part
                do iorb=1,Norb
                   ed_Eeph = ed_Eeph + g_ph(iorb,iorb)*(nup(iorb)+ndw(iorb))
                enddo
                ! UP
                do iorb=1,Norb
                   do jorb=1,Norb
                      ! g_ij cdag_i,up c_j,up (bdag + b)
                      if( g_ph(iorb,jorb)/= (0.d0,0.d0) .and. nup(jorb)==1 .and. nup(iorb)==0 )then
                         call c( jorb,mup,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         jup = binary_search(sectorI%H(1)%map,k2)
                         j_el = jup + (idw-1)*sectorI%DimUp
                         ! N.B.here iph = n+1
                         if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
                            j = j_el + (iph)*sectorI%DimEl
                            ed_Eeph = ed_Eeph + g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph))*v_state(i)*(v_state(j))*peso
                         endif
                         if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
                            j = j_el + (iph-2)*sectorI%DimEl
                            ed_Eeph = ed_Eeph + g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph-1))*v_state(i)*(v_state(j))*peso
                         endif
                      endif
                   end do
                end do
               !
               ! DW
               do iorb=1,Norb
                  do jorb=1,Norb
                     ! g_ij cdag_i,dw c_j,dw (bdag + b)
                     if( g_ph(iorb,jorb)/= (0.d0,0.d0) .and. ndw(jorb)==1 .and. ndw(iorb)==0 )then
                        call c( jorb+Ns,mdw,k1,sg1)
                        call cdg(iorb+Ns,k1,k2,sg2)
                        jdw = binary_search(sectorI%H(2)%map,k2)
                        j_el = iup + (jdw-1)*sectorI%DimUp
                        ! N.B.here iph = n+1
                        if(iph < DimPh) then !bdg = sum_n |n+1> sqrt(n+1) <n|
                           j = j_el + (iph)*sectorI%DimEl
                           ed_Eeph = ed_Eeph + g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph))*v_state(i)*(v_state(j))*peso
                        endif
                        if(iph > 1) then !b = sum_n |n-1> sqrt(n) <n|
                           j = j_el + (iph-2)*sectorI%DimEl
                           ed_Eeph = ed_Eeph + g_ph(iorb,jorb)*sg1*sg2*sqrt(dble(iph-1))*v_state(i)*(v_state(j))*peso
                        endif
                     endif
                  end do
               end do
             endif
             !
          enddo
          call delete_sector(sectorI)         
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
#ifdef _DEBUG
    write(Logfile,"(A)")""
#endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
       call Bcast_MPI(MpiComm,ed_Eeph)
       call Bcast_MPI(MpiComm,ed_Eph)
    endif
#endif
    !
    ed_Eint = ed_Epot !without Hartree. This variable was previously unassigned
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose>=3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A,10f18.12)")"Eeph    =",ed_Eeph
       write(LOGfile,"(A,10f18.12)")"Eph     =",ed_Eph
    endif
    !
    if(MPIMASTER)then
       call write_energy()
    endif
    !
    !
  end subroutine local_energy_normal








  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################




  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    !Write a plain-text file called :code:`observables_info.ed` detailing the names and contents of the observable output files.
    !Write the observable output files. Filenames with suffix :code:`_all` contain values for all DMFT interations, those with suffix :code:`_last` 
    !only values for the last iteration
    call write_obs_info()
    call write_obs_last()
    if(ed_obs_all)call write_obs_all()
  end subroutine write_observables




  subroutine write_obs_info()
    integer :: unit,iorb,jorb,ispin
    !
    !Parameters used:
    unit = free_unit()
      if(.not.ed_read_umatrix)then
      open(unit,file="parameters_info.ed")
      write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
           (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
           reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh",reg(txtfy(2+Norb+3))//"Jx",reg(txtfy(2+Norb+4))//"Jp"
      close(unit)
    endif
    !
    !Generic observables 
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,*(A10,6X))")"#",&
         (str(iorb)//"dens_"//str(iorb),iorb=1,Norb),&
         (str(Norb+iorb)//"docc_"//str(iorb),iorb=1,Norb),&
         (str(2*Norb+iorb)//"nup_"//str(iorb),iorb=1,Norb),&
         (str(3*Norb+iorb)//"ndw_"//str(iorb),iorb=1,Norb),&
         (str(4*Norb+iorb)//"mag_"//str(iorb),iorb=1,Norb),&
         str(5*Norb+1)//"s2tot",str(5*Norb+2)//"egs"
    close(unit)
    !
    !Spin-Spin correlation
    unit = free_unit()
    open(unit,file="Sz2_info.ed")
    write(unit,"(A1,2A6,A15)")"#","a","b","Sz.Sz(a,b)"
    close(unit)
    !
    !Density-Density correlation
    unit = free_unit()
    open(unit,file="N2_info.ed")
    write(unit,"(A1,2A6,A15)")"#","a","b","N.N(a,b)"
    close(unit)
    !
    !Exciton Singlet-Triplet Z
    if(any(abs(exct_s0)>1d-9).OR.any(abs(exct_tz)>1d-9))then
       open(unit,file="Exct_info.ed")
       write(unit,"(A1,*(A10,6X))") "# ",(("Exct_"//str(iorb)//str(jorb),jorb=iorb+1,Norb),iorb=1,Norb)
       close(unit)
       !
    endif
    !
    !Phonons info
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="nph_info.ed")
       write(unit,"(A1,*(A10,6X))") "#","1nph","2X_ph", "3X2_ph"
       close(unit)
       !
       !N_ph probability:
       unit = free_unit()
       open(unit,file="Nph_probability_info.ed")
       write(unit,"(A1,90(A10,6X))")"#",&
            (reg(txtfy(i+1))//"Nph="//reg(txtfy(i)),i=0,DimPh-1)
       close(unit)
    endif
  end subroutine write_obs_info


  subroutine write_obs_last()
    integer :: unit
    integer :: iorb,jorb,ispin
    !
    !Parameters used:
    if(.not.ed_read_umatrix)then
      unit = free_unit()
      open(unit,file="parameters.ed")
      write(unit,"(90F15.9)")xmu,beta,(Uloc_internal(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
      close(unit)
    endif
    !
    !Generic observables 
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs
    close(unit)
    !
    !Spin-Spin correlation
    unit = free_unit()
    open(unit,file="Sz2_last"//reg(ed_file_suffix)//".ed")
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,sz2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Density-Density correlation
    unit = free_unit()
    open(unit,file="N2_last"//reg(ed_file_suffix)//".ed")
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,n2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !
    !Exciton Singlet-Triplet Z
    if(any(abs(exct_s0)>1d-9).OR.any(abs(exct_tz)>1d-9))then
       unit = free_unit()
       open(unit,file="ExctS0_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(*(F15.9,1X))") ((exct_s0(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="ExctTz_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(*(F15.9,1X))") ((exct_tz(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
       close(unit)
    endif
    !
    !Phonons info
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="nph_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90(F15.9,1X))") dens_ph, X_ph, X2_ph
       close(unit)
       !
       unit = free_unit()
       open(unit,file="Occupation_prob"//reg(ed_file_suffix)//".ed")
       write(unit,"(125F15.9)")Prob,sum(Prob)
       close(unit)
       !
       !N_ph probability:
       unit = free_unit()
       open(unit,file="Nph_probability"//reg(ed_file_suffix)//".ed")
       write(unit,"(90(F15.9,1X))") (prob_ph(i),i=1,DimPh)
       close(unit)
    endif
  end subroutine write_obs_last





  subroutine write_obs_all()
    integer :: unit,iorb,jorb,ispin
    !
    !Generic observables 
    unit = free_unit()
    open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(*(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs
    close(unit)
    !
    !Spin-Spin correlation
    unit = free_unit()
    open(unit,file="Sz2_all"//reg(ed_file_suffix)//".ed",position='append')
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,sz2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Density-Density correlation
    unit = free_unit()
    open(unit,file="N2_all"//reg(ed_file_suffix)//".ed",position='append')
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,n2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Exciton Singlet-Triplet Z
    if(any(abs(exct_s0)>1d-9).OR.any(abs(exct_tz)>1d-9))then
       unit = free_unit()
       open(unit,file="ExctS0_all"//reg(ed_file_suffix)//".ed",position='append')
       write(unit,"(*(F15.9,1X))") ((exct_s0(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="ExctTz_all"//reg(ed_file_suffix)//".ed",position='append')
       write(unit,"(*(F15.9,1X))") ((exct_tz(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
       close(unit)
    endif
    !
    !Phonons info
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="nph_all"//reg(ed_file_suffix)//".ed",position='append')
       write(unit,"(*(F15.9,1X))") dens_ph,X_ph, X2_ph
       close(unit)
    endif
  end subroutine write_obs_all







  subroutine write_energy()
    !Write the latest iteration values of energy observables
    integer :: unit
    !
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>",&
         reg(txtfy(7))//"<Dse>",&
         reg(txtfy(8))//"<Dph>",&
         reg(txtfy(9))//"<Eeph>",&
         reg(txtfy(10))//"<Eph>"
    close(unit)
    !
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph,ed_Eeph,ed_Eph
    close(unit)
  end subroutine write_energy








  subroutine write_pdf()
    !Write the lattice probability distribution function
    integer :: unit,i
    real(8) :: x,dx
    unit = free_unit()
    open(unit,file="lattice_prob"//reg(ed_file_suffix)//".ed")
    dx = (xmax-xmin)/dble(Lpos)
    x = xmin
    do i=1,Lpos
       write(unit,"(50F15.9)") x,pdf_ph(i),pdf_part(i,:)
       x = x + dx
    enddo
    close(unit)
  end subroutine write_pdf


  subroutine prob_distr_ph(vec,val)
    ! Compute the contribution of the currect state to the local lattice probability distribution function (PDF)
    ! ,i.e. the local probability of displacement operator X_ph=(bdag+b) as a function of the displacement itself
    !
    ! P(X) = <x\| rho_ph \|x> where rho_ph = Tr_fermions rho_aim
    ! and \|x> = \sum_n psi*_n(x)\|n> with psi_n(x) eigenstates of the quantum harmonic oscillator
    real(8),dimension(:) :: vec
    real(8)              :: psi(0:DimPh-1)
    real(8)              :: x,dx
    integer              :: i,j,i_ph,j_ph,val
    integer              :: istart,jstart,iend,jend
    !
    dx = (xmax-xmin)/dble(Lpos)
    !
    x = xmin
    do i=1,Lpos
       call Hermite(x,psi)
       !
       do i_ph=1,DimPh
          istart = i_el + (i_ph-1)*sectorI%DimEl
          !
          do j_ph=1,DimPh
             jstart = i_el + (j_ph-1)*sectorI%DimEl
             ! Global PDF
             pdf_ph(i) = pdf_ph(i) + peso*psi(i_ph-1)*psi(j_ph-1)*vec(istart)*vec(jstart)
             ! PDF restricted to a sector of orbital occupation
             pdf_part(i,val) = pdf_part(i,val) + peso*psi(i_ph-1)*psi(j_ph-1)*vec(istart)*vec(jstart)
             !
          enddo
       enddo
       !
       x = x + dx
    enddo
  end subroutine prob_distr_ph


  subroutine Hermite(x,psi)
    !Compute the Hermite functions (i.e. harmonic oscillator eigenfunctions)
    !the output is a vector with the functions up to order Dimph-1 evaluated at position x
    real(8),intent(in)  ::  x
    real(8),intent(out) ::  psi(0:DimPh-1)
    integer             ::  i
    real(8)             ::  den
    !
    den=1.331335373062335d0!pigr**(0.25d0)
    !
    psi(0)=exp(-0.5d0*x*x)/den
    psi(1)=exp(-0.5d0*x*x)*sqrt(2d0)*x/den
    !
    do i=2,DimPh-1
       psi(i)=2*x*psi(i-1)/sqrt(dble(2*i))-psi(i-2)*sqrt(dble(i-1)/dble(i))
    enddo
  end subroutine Hermite


end MODULE ED_OBSERVABLES_NORMAL


















