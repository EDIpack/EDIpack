MODULE ED_OBSERVABLES_SUPERC
  !:synopsis: Routines for direct observables calculation, :code:`SUPERC` case
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
  USE ED_HAMILTONIAN_SUPERC
  !
  implicit none
  private
  !
  public :: observables_superc
  public :: local_energy_superc


  real(8),dimension(:),allocatable      :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable      :: docc
  real(8),dimension(:),allocatable      :: magZ,magX,magY
  complex(8),dimension(:,:),allocatable :: phiscAB
  real(8),dimension(:,:),allocatable    :: rePhi,imPhi
  real(8),dimension(:),allocatable      :: phisc
  real(8),dimension(:,:),allocatable    :: sz2,n2
  real(8)                               :: s2tot
  real(8)                               :: Egs
  real(8)                               :: Ei
  real(8),dimension(:),allocatable      :: Prob
  real(8),dimension(:),allocatable      :: prob_ph
  real(8),dimension(:),allocatable      :: pdf_ph
  real(8),dimension(:,:),allocatable    :: pdf_part
  real(8)                               :: dens_ph
  real(8)                               :: X_ph, X2_ph

  !
  integer                               :: iorb,jorb,istate
  integer                               :: ispin,jspin
  integer                               :: isite,jsite
  integer                               :: ibath
  integer                               :: r,m,k,k1,k2,k3,k4
  integer                               :: iup,idw
  integer                               :: jup,jdw
  integer                               :: mup,mdw
  integer                               :: iph,i_el,j_el,isz
  real(8)                               :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                               :: gs_weight
  real(8)                               :: peso
  real(8)                               :: norm

  !
  integer                               :: i,j,ii,iprob
  integer                               :: isector,jsector
  !
  complex(8),dimension(:),allocatable   :: vvinit
  complex(8),dimension(:),allocatable   :: veta,vkappa
  complex(8),dimension(:),allocatable   :: v_state
  logical                               :: Jcondition
  !
  type(sector)                          :: sectorI,sectorJ


contains 



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_superc()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Calculate the values of the local observables
    integer                 :: val
    integer,dimension(2*Ns) :: ib
    integer,dimension(2,Ns) :: Nud
    integer,dimension(Ns)   :: IbUp,IbDw
    real(8),dimension(Norb) :: nup,ndw,Sz,nt
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG observables_superc"
#endif
    !
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(phisc(Norb),phiscAB(Norb,Norb),RePhi(Norb,Norb),ImPhi(Norb,Norb))
    allocate(Prob(3**Norb))
    allocate(prob_ph(DimPh))
    allocate(pdf_ph(Lpos))
    allocate(pdf_part(Lpos,3))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    phisc   = 0.d0
    phiscAB = zero
    rePhi   = 0.d0
    imPhi   = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    Prob    = 0.d0
    prob_ph = 0.d0
    X_ph = 0.d0
    X2_ph = 0.d0
    dens_ph = 0.d0
    pdf_ph  = 0.d0
    pdf_part= 0.d0
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_superc: get local observables"
#endif
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)       
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_superc: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(Mpimaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             gs_weight=peso*abs(v_state(i))**2
             i_el = mod(i-1,sectorI%DimEl)+1
             m    = sectorI%H(1)%map(i_el)
             ib   = bdecomp(m,2*Ns)
             do iorb=1,Norb
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
             enddo
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
                j= i_el + (iph)*sectorI%DimEl
                X_ph = X_ph + sqrt(2.d0*dble(iph))*real(v_state(i)*conjg(v_state(j)))*peso
             end if
             X2_ph = X2_ph + 0.5d0*(1+2*(iph-1))*gs_weight
             if(iph<DimPh-1)then
                j= i_el + (iph+1)*sectorI%DimEl
                X2_ph = X2_ph + sqrt(dble((iph)*(iph+1)))*real(v_state(i)*conjg(v_state(j)))*peso
             end if
             !compute the lattice probability distribution function
             if(Dimph>1 .AND. iph==1) then
                val = 1
                !val = 1 + Nr. of polarized orbitals (full or empty) makes sense only for 2 orbs
                do iorb=1,Norb
                   val = val + abs(nint(sign((nt(iorb) - 1.d0),real(g_ph(iorb,iorb)))))
                enddo
                call prob_distr_ph(v_state,val)
             end if
          enddo
          !
          !
          !SUPERCONDUCTING ORDER PARAMETER
#ifdef _DEBUG
          if(ed_verbose>2)write(Logfile,"(A)")"DEBUG observables_superc: get OP"
#endif
          do ispin=1,Nspin 
             !GET: <(b_up + adg_dw)(bdg_up + a_dw)> = 
             !<b_up*bdg_up> + <adg_dw*a_dw> + <b_up*a_dw> + <adg_dw*bdg_up> = 
             !<n_a,dw> + < 1 - n_b,up> + 2*<PHI>_ab
             !EVALUATE [a_dw + bdg_up]|gs> = [1,1].[C_{-1},C_{+1}].[iorb,jorb].[dw,up]

             !Get:
             !\eta   = <(B_up + A^+_dw) (B^+_up + A_dw)>  = <n_Adw> + <1-n_Bup> + 2RePhi_AB
             !\kappa = <(B_up +i.A^+_dw)(B^+_up -i.A_dw)> = <n_Adw> + <1-n_Bup> + 2ImPhi_AB

             do iorb=1,Norb !A
                do jorb=1,Norb !B 
                   isz = getsz(isector)
                   if(isz<Ns)then
                      jsector = getsector(isz+1,1)
                      veta    = apply_Cops(v_state,[one,one],[-1,1],[iorb,jorb],[2,1],isector,jsector)
                      vkappa  = apply_Cops(v_state,[one,xi], [-1,1],[iorb,jorb],[2,1],isector,jsector)
                      RePhi(iorb,jorb) = RePhi(iorb,jorb) + dot_product(veta,veta)*peso
                      ImPhi(iorb,jorb) = ImPhi(iorb,jorb) + dot_product(vkappa,vkappa)*peso
                      deallocate(veta,vkappa)
                   endif
                enddo
             enddo
          enddo
          !
          call delete_sector(sectorI)
          !
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
    !
    do iorb=1,Norb
       do jorb=1,Norb
          RePhi(iorb,jorb) = 0.5d0*(RePhi(iorb,jorb) - dens_dw(iorb) - (1.d0-dens_up(jorb)))
          ImPhi(iorb,jorb) = 0.5d0*(ImPhi(iorb,jorb) - dens_dw(iorb) - (1.d0-dens_up(jorb)))
          ! phiscAB(iorb,jorb) = 0.5d0*(phiscAB(iorb,jorb) - dens_dw(iorb) - (1.d0-dens_up(jorb)))
          phiscAB(iorb,jorb) = dcmplx(RePhi(iorb,jorb),ImPhi(iorb,jorb))
       enddo
       phisc(iorb)=abs(phiscAB(iorb,iorb)) !sqrt(rePhi**2 + imPhi**2)
    enddo
    !
    !STATUS: IMPORTED FROM NORMAL, TO BE UPDATE TO NAMBU BASIS <\psi^+_a Psi_a>
    !     !
    !     !SINGLE PARTICLE IMPURITY DENSITY MATRIX
    ! #ifdef _DEBUG
    !     if(ed_verbose>2)write(Logfile,"(A)")&
    !          "DEBUG observables_superc: eval single particle density matrix <C^+_a C_b>"
    ! #endif
    !     if(allocated(single_particle_density_matrix)) deallocate(single_particle_density_matrix)
    !     allocate(single_particle_density_matrix(Nspin,Nspin,Norb,Norb));single_particle_density_matrix=zero
    !     do istate=1,state_list%size
    !        isector = es_return_sector(state_list,istate)
    !        Ei      = es_return_energy(state_list,istate)
    ! #ifdef _DEBUG
    !        if(ed_verbose>3)write(Logfile,"(A)")&
    !             "DEBUG observables_normal: get contribution from state:"//str(istate)
    ! #endif
    ! #ifdef _MPI
    !        if(MpiStatus)then
    !           call es_return_dvector(MpiComm,state_list,istate,state_dvec) 
    !        else
    !           call es_return_dvector(state_list,istate,state_dvec) 
    !        endif
    ! #else
    !        call es_return_dvector(state_list,istate,state_dvec) 
    ! #endif
    !        !
    !        peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
    !        peso = peso/zeta_function
    !        !
    !        if(MpiMaster)then
    !           call build_sector(isector,sectorI)
    !           do i=1,sectorI%Dim
    !              iph = (i-1)/(sectorI%DimEl) + 1
    !              i_el = mod(i-1,sectorI%DimEl) + 1
    !              call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
    !              !
    !              call build_op_Ns(i,IbUp,IbDw,sectorI)
    !              Nud(1,:)=IbUp
    !              Nud(2,:)=IbDw
    !              !
    !              !Diagonal densities
    !              do ispin=1,Nspin
    !                 do iorb=1,Norb
    !                    single_particle_density_matrix(ispin,ispin,iorb,iorb) = &
    !                         single_particle_density_matrix(ispin,ispin,iorb,iorb) + &
    !                         peso*nud(ispin,iorb)*(state_dvec(i))*state_dvec(i)
    !                 enddo
    !              enddo
    !              !
    !              !Off-diagonal
    !              if(ed_total_ud)then
    !                 do ispin=1,Nspin
    !                    do iorb=1,Norb
    !                       do jorb=1,Norb
    !                          !
    !                          if((Nud(ispin,jorb)==1).and.(Nud(ispin,iorb)==0))then
    !                             iud(1) = sectorI%H(1)%map(Indices(1))
    !                             iud(2) = sectorI%H(2)%map(Indices(2))
    !                             call c(jorb,iud(ispin),r,sgn1)
    !                             call cdg(iorb,r,k,sgn2)
    !                             Jndices = Indices
    !                             Jndices(1+(ispin-1)*Ns_Ud) = &
    !                                  binary_search(sectorI%H(1+(ispin-1)*Ns_Ud)%map,k)
    !                             call indices2state(Jndices,[sectorI%DimUps,sectorI%DimDws],j)
    !                             !
    !                             j = j + (iph-1)*sectorI%DimEl
    !                             !
    !                             single_particle_density_matrix(ispin,ispin,iorb,jorb) = &
    !                                  single_particle_density_matrix(ispin,ispin,iorb,jorb) + &
    !                                  peso*sgn1*state_dvec(i)*sgn2*(state_dvec(j))
    !                          endif
    !                       enddo
    !                    enddo
    !                 enddo
    !              endif
    !              !
    !              !
    !           enddo
    !           call delete_sector(sectorI)         
    !        endif
    !        !
    !        if(allocated(state_dvec))deallocate(state_dvec)
    !        !
    !     enddo
    ! #ifdef _DEBUG
    !     if(ed_verbose>2)write(Logfile,"(A)")""
    ! #endif
    !     !
    write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    write(LOGfile,"(A,10f18.12,A)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    write(LOGfile,"(A,20f18.12,A)")    "phiAB "//reg(ed_file_suffix)//"=",((phiscAB(iorb,jorb),iorb=1,Norb),jorb=1,Norb)
    write(LOGfile,"(A,20f18.12,A)")     " | phiAA*Uloc ",(abs(Uloc_internal(iorb))*phisc(iorb),iorb=1,Norb)
    if(Nspin==2)then
       write(LOGfile,"(A,10f18.12,A)")    "magZ"//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
    endif
    !
    if(DimPh>1)then
       call write_pdf()
    endif
    !
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
       ed_mag(3,iorb)  =magZ(iorb)
       do jorb=1,Norb
          ed_phisc(iorb,jorb)  =phiscAB(iorb,jorb)
       enddo
    enddo
    !
    ed_imp_info=[s2tot,egs]
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_phisc)
       call Bcast_MPI(MpiComm,ed_mag)
       call Bcast_MPI(MpiComm,ed_imp_info)
    endif
#endif
    !
    if(MPIMASTER)then
       call write_observables()
    endif
    !
    deallocate(dens,docc,phiscAB,phisc,dens_up,dens_dw,magz,sz2,n2,Prob)
    deallocate(prob_ph,pdf_ph,pdf_part)
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
  end subroutine observables_superc







  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_superc()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Calculate the values of the local observables
    integer,dimension(2*Ns) :: ib
    integer,dimension(2,Ns) :: Nud
    integer,dimension(Ns)   :: IbUp,IbDw
    real(8),dimension(Norb)         :: nup,ndw
    real(8),dimension(Nspin,Norb)   :: eloc
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG local_energy_superc"
#endif
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG local_energy_superc: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(Mpimaster)then
          !
          call build_sector(isector,sectorI)
          do i=1,sectorI%Dim
             iph  = (i-1)/(sectorI%DimEl)+1
             i_el = mod(i-1,sectorI%DimEl)+1
             m    = sectorI%H(1)%map(i_el)
             ib   = bdecomp(m,2*Ns)
             do iorb=1,Norb
                nup(iorb)=dble(ib(iorb))
                ndw(iorb)=dble(ib(iorb+Ns))
             enddo
             !
             gs_weight=peso*abs(v_state(i))**2
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !> H_Imp: Diagonal Elements, i.e. local part
             do iorb=1,Norb
                ed_Eknot = ed_Eknot + impHloc(1,1,iorb,iorb)*Nup(iorb)*gs_weight
                ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)*gs_weight
             enddo
             ! !> H_imp: Off-diagonal elements, i.e. non-local part. 
             do iorb=1,Norb
                do jorb=1,Norb
                   !SPIN UP
                   Jcondition = &
                        (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                        (ib(jorb)==1)                  .AND. &
                        (ib(iorb)==0)
                   if (Jcondition) then
                      call c(jorb,m,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      j_el=binary_search(sectorI%H(1)%map,k2)
                      j   = j_el + (iph-1)*sectorI%DimEl
                      ed_Eknot = ed_Eknot + impHloc(1,1,iorb,jorb)*sg1*sg2*v_state(i)*conjg(v_state(j))*peso
                   endif
                   !SPIN DW
                   Jcondition = &
                        (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                        (ib(jorb+Ns)==1)                       .AND. &
                        (ib(iorb+Ns)==0)
                   if (Jcondition) then
                      call c(jorb+Ns,m,k1,sg1)
                      call cdg(iorb+Ns,k1,k2,sg2)
                      j_el=binary_search(sectorI%H(1)%map,k2)
                      j   = j_el + (iph-1)*sectorI%DimEl
                      ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*v_state(i)*conjg(v_state(j))*peso
                   endif
                enddo
             enddo
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
             !SPIN-EXCHANGE (S-E) TERMS
             !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
             if(Norb>1.AND.(any((Jx_internal/=0d0)).OR.any((Jp_internal/=0d0))))then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(iorb+Ns)==1)   .AND.&
                           (ib(jorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(iorb+Ns,k1,k2,sg2)
                         call cdg(jorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j_el=binary_search(sectorI%H(1)%map,k4)
                         j   = j_el + (iph-1)*sectorI%DimEl
                         ed_Epot = ed_Epot + Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                         ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                      endif
                   enddo
                enddo
             endif
             !
             !
             !PAIR-HOPPING (P-H) TERMS
             !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
             !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
             if(Norb>1.AND.(any((Jx_internal/=0d0)).OR.any((Jp_internal/=0d0))))then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(jorb+Ns)==1)   .AND.&
                           (ib(iorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(jorb+Ns,k1,k2,sg2)
                         call cdg(iorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j_el=binary_search(sectorI%H(1)%map,k4)
                         j   = j_el + (iph-1)*sectorI%DimEl
                         ed_Epot = ed_Epot + Jp_internal(iorb,jorb)*sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                         ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                      endif
                   enddo
                enddo
             endif
             !
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then               
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
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
       call Bcast_MPI(MpiComm,ed_Dse)
       call Bcast_MPI(MpiComm,ed_Dph)
    endif
#endif
    !
    ed_Eint = ed_Epot
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose>=3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
       write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
    endif
    if(MPIMASTER)then
       call write_energy()
    endif
    !
    !
  end subroutine local_energy_superc



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
    !Parameters used:
    if(.not.ed_read_umatrix)then
       unit = free_unit()
       open(unit,file="parameters_info.ed")
       write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
            (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
            reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
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
    !SC order parameters
    unit = free_unit()
    open(unit,file="phi_info.ed")
    write(unit,"(A1,*(A10,6X))")"#",(("phi_"//str(iorb)//str(jorb),jorb=1,Norb),iorb=1,Norb)
    close(unit)  
    !
    !Phonons info
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="nph_info.ed")
       write(unit,"(A1,*(A10,6X))") "#","1nph", "2X_ph", "3X2_ph"
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
    integer :: unit,iorb,jorb,ispin
    !
    !Parameters used:
    if(.not.ed_read_umatrix)then
       unit = free_unit()
       open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
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
    !SC order parameters
    unit = free_unit()
    open(unit,file="phi_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(1X,*(F15.9,1x))")((phiscAB(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
    close(unit)
    !
    !Phonons info
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="nph_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(90(F15.9,1X))") dens_ph,X_ph, X2_ph
       close(unit)
       !
       !
       if(.not.ed_read_umatrix)then
          unit = free_unit()
          open(unit,file="Occupation_prob"//reg(ed_file_suffix)//".ed")
          write(unit,"(125F15.9)")Uloc(1),Prob,sum(Prob)
          close(unit)
       endif
       !
       !N_ph probability:
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
    !SC order parameters
    unit = free_unit()
    open(unit,file="phi_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(1X,*(F15.9,1x))")((phiscAB(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
    close(unit)
    !
    !Phonons info
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="nph_all"//reg(ed_file_suffix)//".ed",position='append')
       write(unit,"(*(F15.9,1X))") dens_ph, X_ph, X2_ph
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
         reg(txtfy(8))//"<Dph>"
    close(unit)
    !
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy






  subroutine write_pdf()
    integer :: unit,i
    real(8) :: x,dx
    unit = free_unit()
    open(unit,file="lattice_prob"//reg(ed_file_suffix)//".ed")
    dx = (xmax-xmin)/dble(Lpos)
    x = xmin
    do i=1,Lpos
       write(unit,"(5F15.9)") x,pdf_ph(i),pdf_part(i,:)
       x = x + dx
    enddo
    close(unit)
  end subroutine write_pdf


  subroutine prob_distr_ph(vec,val)
    !Compute the local lattice probability distribution function (PDF), i.e. the local probability of displacement
    !as a function of the displacement itself
    complex(8),dimension(:) :: vec
    real(8)              :: psi(0:DimPh-1)
    real(8)              :: x,dx
    integer              :: i,j,j_ph,val
    integer              :: istart,jstart,iend,jend
    !
    dx = (xmax-xmin)/dble(Lpos)
    !
    x = xmin
    do i=1,Lpos !cycle over x
       call Hermite(x,psi)
       !
       istart = i_el + (iph-1)*sectorI%DimEl !subroutine already inside a sectorI cycle
       !
       !phonon diagonal part
       pdf_ph(i) = pdf_ph(i) + peso*psi(iph-1)*psi(iph-1)*abs(vec(istart))**2
       pdf_part(i,val) = pdf_part(i,val) + peso*psi(iph-1)*psi(iph-1)*abs(vec(istart))**2
       !
       !phonon off-diagonal part
       do j_ph=iph+1,DimPh
          jstart = i_el + (j_ph-1)*sectorI%DimEl
          pdf_ph(i)       = pdf_ph(i)       + peso*psi(iph-1)*psi(j_ph-1)*2.d0*real( vec(istart)*conjg(vec(jstart)) )
          pdf_part(i,val) = pdf_part(i,val) + peso*psi(iph-1)*psi(j_ph-1)*2.d0*real( vec(istart)*conjg(vec(jstart)) )
       enddo
       !
       x = x + dx
    enddo
  end subroutine prob_distr_ph

  !Compute the Hermite functions (i.e. harmonic oscillator eigenfunctions)
  !the output is a vector with the functions up to order Dimph-1 evaluated at position x
  subroutine Hermite(x,psi)
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



end MODULE ED_OBSERVABLES_SUPERC
