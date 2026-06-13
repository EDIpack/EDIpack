!KSC=Krylov-state complexity
!KOC=Krylov-operator complexity
MODULE ED_KRYLOV
  !
  !:synopsis: A Module to evaluate the Krylov-state and Krylov-operator complexities
  ! for the impurity problem.
  !
  USE SF_CONSTANTS, only: xi
  USE SF_IOTOOLS,  only: str,reg,to_lower,free_unit
#ifdef _MPI
  USE SF_MPI
#endif
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SECTOR
  USE ED_AUX_FUNX, only: allocate_grids
  USE ED_BATH,     only: allocate_dmft_bath,set_dmft_bath,deallocate_dmft_bath,check_bath_dimension
  USE ED_KRYLOV_NORMAL
  USE ED_KRYLOV_SUPERC
  USE ED_KRYLOV_NONSU2
  implicit none
  private

  interface KSC_Krylov_Basis
     module procedure :: KSC_Krylov_Basis_d
     module procedure :: KSC_Krylov_Basis_c
  end interface KSC_Krylov_Basis

  interface KSC_AddFromState
     module procedure :: KSC_AddFromState_d
     module procedure :: KSC_AddFromState_c
  end interface KSC_AddFromState 
  
  

  public :: krylov_state_complexity

  integer :: KSCmax=0


contains




      
  subroutine krylov_state_complexity(bath)
    !
    ! Build the Krylov-state complexity for the impurity problem for the states obtained 
    ! adding one electrons with a defined spin and orbital index AND for the states with 
    ! given orbital resolved charge density. 
    !
    real(8),dimension(:),intent(in),optional          :: bath
    real(8),dimension(:),allocatable                   :: Kt,Sk,Nt
    real(8),dimension(:,:),allocatable                 :: Pnk
    real(8),dimension(1)                               :: bath_dummy
    integer                                            :: iorb,ispin
    logical                                            :: bath_allocated_here,check
    !
    if(.not.state_list%status.OR.state_list%size==0)&
      stop "krylov_state_complexity error: state_list is empty."
    bath_allocated_here = .false.
    if(present(bath))then
      if(Nbath>0)then
        check = check_bath_dimension(bath)
        if(.not.check)stop "krylov_state_complexity error: wrong bath dimensions"
      endif
      call allocate_dmft_bath()
      call set_dmft_bath(bath)
      bath_allocated_here = .true.
    elseif(.not.dmft_bath%status)then
      if(Nbath>0)then
        stop "krylov_state_complexity error: dmft_bath is not allocated. Pass the bath array."
      else
        bath_dummy = zero
        call allocate_dmft_bath()
        call set_dmft_bath(bath_dummy)
        bath_allocated_here = .true.
      endif
    endif
    !
    call allocate_grids()
    !
    KSCmax = max_krylov_size()
    !
    if(allocated(ed_Kcdg))   deallocate(ed_Kcdg)
    if(allocated(ed_Scdg))   deallocate(ed_Scdg)
    if(allocated(ed_Normcdg))deallocate(ed_Normcdg)
    if(allocated(ed_Pcdg))   deallocate(ed_Pcdg)
    if(allocated(ed_Kn))     deallocate(ed_Kn)
    if(allocated(ed_Sn))     deallocate(ed_Sn)
    if(allocated(ed_Normn))  deallocate(ed_Normn)
    if(allocated(ed_Pn))     deallocate(ed_Pn)
    allocate(ed_Kcdg(Nspin,Norb,Ltimes))
    allocate(ed_Scdg(Nspin,Norb,Ltimes))
    allocate(ed_Normcdg(Nspin,Norb,Ltimes))
    allocate(ed_Pcdg(Nspin,Norb,KSCmax,Ltimes))
    allocate(ed_Kn(Nspin,Norb,Ltimes))
    allocate(ed_Sn(Nspin,Norb,Ltimes))
    allocate(ed_Normn(Nspin,Norb,Ltimes))
    allocate(ed_Pn(Nspin,Norb,KSCmax,Ltimes))
    !
    ed_Kcdg   =0d0
    ed_Scdg   =0d0
    ed_Normcdg=0d0
    ed_Pcdg   =0d0
    ed_Kn     =0d0
    ed_Sn     =0d0
    ed_Normn  =0d0
    ed_Pn     =0d0
    !
    !States with one electron with a given spin and orbital index C^+_{iorb,ispin}
    do ispin=1,Nspin
      do iorb=1,Norb
          call KSC_build_OpStateList("cdg",iorb,ispin,Kt,Sk,Pnk,Nt)
          ed_Kcdg(ispin,iorb,:)    = Kt
          ed_Scdg(ispin,iorb,:)    = Sk
          ed_Normcdg(ispin,iorb,:) = Nt
          ed_Pcdg(ispin,iorb,:,:)  = Pnk
          if(allocated(Kt)) deallocate(Kt)
          if(allocated(Sk)) deallocate(Sk)
          if(allocated(Pnk))deallocate(Pnk)
          if(allocated(Nt)) deallocate(Nt)
      enddo
    enddo
    !States with orbital resolved density N_{iorb}
    do iorb=1,Norb
      call KSC_build_OpStateList("N",iorb,1,Kt,Sk,Pnk,Nt)
      ed_Kn(1,iorb,:)    = Kt
      ed_Sn(1,iorb,:)    = Sk
      ed_Normn(1,iorb,:) = Nt
      ed_Pn(1,iorb,:,:)  = Pnk
      if(allocated(Kt)) deallocate(Kt)
      if(allocated(Sk)) deallocate(Sk)
      if(allocated(Pnk))deallocate(Pnk)
      if(allocated(Nt)) deallocate(Nt)
    enddo
    !
    if(MPIMASTER)call KSC_write()
    if(bath_allocated_here)call deallocate_dmft_bath()
    !
  end subroutine krylov_state_complexity      







  !Build the Krylov-complexity for the states in the states_list with an applied Op.  
  subroutine KSC_build_OpStateList(op,iorb,ispin,Kt,Sk,Pnk,Nt)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),allocatable,intent(out)   :: Kt
    real(8),dimension(:),allocatable,intent(out)   :: Sk
    real(8),dimension(:,:),allocatable,intent(out) :: Pnk
    real(8),dimension(:),allocatable,intent(out)   :: Nt
    real(8)                                        :: spectral_weight
    !
    if(KSCmax==0)stop "build_ksc_OpStateList error: KSCmax=0"
    !
    allocate(Kt(Ltimes),Sk(Ltimes))
    allocate(Pnk(KSCmax,Ltimes))
    allocate(Nt(Ltimes))
    Kt              = 0d0
    Sk              = 0d0
    Pnk             = 0d0
    Nt              = 0d0
    spectral_weight = 0d0
    !
    select case(ed_mode)
      case default  ;call KSC_ApplyOp_StateList_normal(op,iorb,ispin,Kt,Sk,Pnk,Nt,spectral_weight)   
      case("superc");call KSC_ApplyOp_StateList_superc(op,iorb,ispin,Kt,Sk,Pnk,Nt,spectral_weight)
      case("nonsu2");call KSC_ApplyOp_StateList_nonsu2(op,iorb,ispin,Kt,Sk,Pnk,Nt,spectral_weight)
    end select
    !
    if(spectral_weight==0d0)then
      if(MPIMASTER)write(LOGfile,"(A)")&
      "WARNING KSC_build_OpStateList WARNING: zero spectral weight for op="//&
      trim(op)//", iorb="//str(iorb)//", ispin="//str(ispin)
    endif
    !
    if(spectral_weight>0d0)then
      Kt  = Kt/spectral_weight
      Sk  = Sk/spectral_weight
      Pnk = Pnk/spectral_weight
      Nt  = Nt/spectral_weight
    endif
    !
  end subroutine KSC_build_OpStateList



  subroutine KSC_ApplyOp_StateList_normal(op,iorb,ispin,Kt,Sk,Pnk,Nt,spectral_weight)
    character(len=*),intent(in)          :: op
    integer,intent(in)                   :: iorb,ispin
    real(8),dimension(:),intent(inout)   :: Kt,Sk
    real(8),dimension(:,:),intent(inout) :: Pnk
    real(8),dimension(:),intent(inout)   :: Nt
    real(8),intent(inout)                :: spectral_weight
    real(8),dimension(:),allocatable     :: v_state,vvinit
    integer                              :: istate,isector,jsector,ialfa
    !
    ialfa = 1
    if(.not.ed_total_ud)ialfa = iorb
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       v_state = es_return_dvec(state_list,istate)
       select case(to_lower(trim(op)))
       case("cdg")
          jsector = getCDGsector(ialfa,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("n")
          jsector = isector
          vvinit = apply_op_N(v_state,iorb,isector)
          call KSC_AddFromState(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
          if(allocated(vvinit))deallocate(vvinit)
       case default
          stop "ed_build_krylov_operator: op must be 'cdg' or 'N'"
       end select
       if(allocated(v_state))deallocate(v_state)
    enddo
  end subroutine KSC_ApplyOp_StateList_normal


  subroutine KSC_ApplyOp_StateList_superc(op,iorb,ispin,Kt,Sk,Pnk,Nt,spectral_weight)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),intent(inout)             :: Kt,Sk
    real(8),dimension(:,:),intent(inout)           :: Pnk
    real(8),dimension(:),intent(inout),optional    :: Nt
    real(8),intent(inout)                          :: spectral_weight
    complex(8),dimension(:),allocatable            :: v_state,vvinit
    integer                                        :: istate,isector,jsector
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       v_state = es_return_cvec(state_list,istate)
       select case(to_lower(trim(op)))
       case("cdg")
          jsector = getCDGsector(1,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("n")
          jsector = isector
          vvinit = apply_op_N(v_state,iorb,isector)
          call KSC_AddFromState(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
          if(allocated(vvinit))deallocate(vvinit)
       case default
          stop "ed_build_krylov_operator: op must be 'cdg' or 'N'"
       end select
       if(allocated(v_state))deallocate(v_state)
    enddo
  end subroutine KSC_ApplyOp_StateList_superc


  subroutine KSC_ApplyOp_StateList_nonsu2(op,iorb,ispin,Kt,Sk,Pnk,Nt,spectral_weight)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),intent(inout)             :: Kt,Sk
    real(8),dimension(:,:),intent(inout)           :: Pnk
    real(8),dimension(:),intent(inout),optional    :: Nt
    real(8),intent(inout)                          :: spectral_weight
    complex(8),dimension(:),allocatable            :: v_state,vvinit
    integer                                        :: istate,isector,jsector
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       v_state = es_return_cvec(state_list,istate)
       select case(to_lower(trim(op)))
       case("cdg")
          jsector = getCDGsector(1,ispin,isector)
          if(Jz_basis)jsector = getCDGsector_Jz(iorb,ispin,isector)
          if(getN(isector)/=Nlevels.and.jsector>=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("n")
          jsector = isector
          vvinit = apply_op_N(v_state,iorb,isector)
          call KSC_AddFromState(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
          if(allocated(vvinit))deallocate(vvinit)
       case default
          stop "ed_build_krylov_operator: op must be 'cdg' or 'N'"
       end select
       if(allocated(v_state))deallocate(v_state)
    enddo
  end subroutine KSC_ApplyOp_StateList_nonsu2









  !INTERFACED ABOVE TO SINGLE FUNCTION: KSC_AddState
  subroutine KSC_AddFromState_d(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
    integer,intent(in)                   :: istate,jsector
    real(8),dimension(:),intent(inout)   :: vvinit
    real(8),dimension(:),intent(inout)   :: Kt,Sk
    real(8),dimension(:,:),intent(inout) :: Pnk
    real(8),dimension(:),intent(inout)   :: Nt
    real(8),intent(inout)                :: spectral_weight
    real(8),dimension(:),allocatable     :: alanc,blanc,Ktmp,Stmp,Ntmp
    real(8),dimension(:,:),allocatable   :: Ptmp
    real(8)                              :: norm2,weight
    !
    norm2 = 0d0
    if(MpiMaster)norm2 = dot_product(vvinit,vvinit)
#ifdef _MPI
    if(MpiStatus)call Bcast_MPI(MpiComm_Global,norm2)
#endif
    if(norm2 <= 0d0)return
    !
    call KSC_Krylov_Basis(jsector,vvinit,alanc,blanc,norm2)
    if(norm2 <= 0d0)return
    !
    call KSC_Build_Complexity(alanc,blanc,Ktmp,Stmp,Ptmp,Ntmp)
    weight = thermal_weight(istate)*norm2
    call KSC_add_weighted_krylov(weight,Ktmp,Stmp,Ptmp,Ntmp,Kt,Sk,Pnk,Nt,spectral_weight)
    !
    if(allocated(alanc))deallocate(alanc)
    if(allocated(blanc))deallocate(blanc)
    if(allocated(Ktmp)) deallocate(Ktmp)
    if(allocated(Stmp)) deallocate(Stmp)
    if(allocated(Ptmp)) deallocate(Ptmp)
    if(allocated(Ntmp)) deallocate(Ntmp)
  end subroutine KSC_AddFromState_d


  subroutine KSC_AddFromState_c(istate,jsector,vvinit,Kt,Sk,Pnk,Nt,spectral_weight)
    integer,intent(in)                    :: istate,jsector
    complex(8),dimension(:),intent(inout) :: vvinit
    real(8),dimension(:),intent(inout)    :: Kt,Sk
    real(8),dimension(:,:),intent(inout)  :: Pnk
    real(8),dimension(:),intent(inout)    :: Nt
    real(8),intent(inout)                 :: spectral_weight
    real(8),dimension(:),allocatable      :: alanc,blanc,Ktmp,Stmp,Ntmp
    real(8),dimension(:,:),allocatable    :: Ptmp
    real(8)                               :: norm2,weight
    !
    norm2 = 0d0
    if(MpiMaster)norm2 = dot_product(vvinit,vvinit)
#ifdef _MPI
    if(MpiStatus)call Bcast_MPI(MpiComm_Global,norm2)
#endif
    if(norm2 <= 0d0)return
    !
    !Get the Krylov basis of alancs/blancs: tridiag H(jsector)
    call KSC_Krylov_Basis(jsector,vvinit,alanc,blanc,norm2)
    if(norm2 <= 0d0)return
    !
    call KSC_Build_Complexity(alanc,blanc,Ktmp,Stmp,Ptmp,Ntmp)
    weight = thermal_weight(istate)*norm2
    call KSC_add_weighted_krylov(weight,Ktmp,Stmp,Ptmp,Ntmp,Kt,Sk,Pnk,Nt,spectral_weight)
    !
    if(allocated(alanc))deallocate(alanc)
    if(allocated(blanc))deallocate(blanc)
    if(allocated(Ktmp)) deallocate(Ktmp)
    if(allocated(Stmp)) deallocate(Stmp)
    if(allocated(Ptmp)) deallocate(Ptmp)
    if(allocated(Ntmp)) deallocate(Ntmp)
  end subroutine KSC_AddFromState_c






  !INTERFACE TO GENERAL PROCEDURE KSC_Krylov_Basis
  subroutine KSC_Krylov_Basis_d(isector,vvinit,alanc,blanc,norm2)
    integer,intent(in)                 :: isector
    real(8),dimension(:),intent(inout) :: vvinit
    real(8),dimension(:),allocatable   :: alanc
    real(8),dimension(:),allocatable   :: blanc
    real(8)                            :: norm2
    select case(ed_mode)
    case default; call KSC_Krylov_Basis_normal(isector,vvinit,alanc,blanc,norm2)
    case("superc","nonsu2");stop "ed_build_krylov_state: Dble states only for ed_mode=normal"
    end select
  end subroutine KSC_Krylov_Basis_d
  !
  subroutine KSC_Krylov_Basis_c(isector,vvinit,alanc,blanc,norm2)
    integer,intent(in)                    :: isector
    complex(8),dimension(:),intent(inout) :: vvinit
    real(8),dimension(:),allocatable      :: alanc
    real(8),dimension(:),allocatable      :: blanc
    real(8)                               :: norm2
    select case(ed_mode)
    case default;stop "ed_build_krylov_state: Cmplx states only for ed_mode = superc/nonsu2"      
    case("superc");call KSC_Krylov_Basis_superc(isector,vvinit,alanc,blanc,norm2)
    case("nonsu2");call KSC_Krylov_Basis_nonsu2(isector,vvinit,alanc,blanc,norm2)
    end select
  end subroutine KSC_Krylov_Basis_c












  !Build the Krylov coefficients dyanamics \phi_n(t) so that \Psi(t) = \sum_n \phi_n(t) |k_n>
  !and return the Krylov complexity Kt, Sk, Pnk
  subroutine KSC_Build_Complexity(alanc,blanc,Kt,Sk,Pnk,Nt)
    real(8),dimension(:),intent(in)                :: alanc
    real(8),dimension(size(alanc)),intent(in)      :: blanc
    real(8),dimension(:),allocatable,intent(out)   :: Kt
    real(8),dimension(:),allocatable,intent(out)   :: Sk
    real(8),dimension(:,:),allocatable,intent(out) :: Pnk
    real(8),dimension(:),allocatable,intent(out)   :: Nt
    complex(8),dimension(:,:),allocatable          :: phi
    real(8)                                        :: P
    integer                                        :: it,in,Nlanc
    !
    !Build the Krylov coefficients dyanamics \phi_n(t)
    call KSC_Evolve_Phi(alanc,blanc,phi)  !<- allocate Phi in here:
    !
    Nlanc=size(alanc)
    if(allocated(Kt)) deallocate(Kt)
    if(allocated(Sk)) deallocate(Sk)
    if(allocated(Pnk))deallocate(Pnk)
    allocate(Kt(Ltimes))       ;Kt =0d0
    allocate(Sk(Ltimes))       ;Sk =0d0
    allocate(Pnk(Nlanc,Ltimes));Pnk=0d0
    allocate(Nt(Ltimes))       ;Nt =0d0
    !
    do it=1,Ltimes
      do in=1,Nlanc
          P          = abs(phi(in,it))**2
          Pnk(in,it) = P
          Kt(it)     = Kt(it) + dble(in-1)*P
          if(P>0d0)Sk(it)   = Sk(it) - P*log(P)
          Nt(it) = Nt(it) + P
      enddo
    enddo
    !
    deallocate(phi)
    !
  end subroutine KSC_Build_Complexity       






  !Build the Krylov coefficients dyanamics \phi_n(t)
  !solving the coupled equations:
  !\phi_n(t) = \sum_m a_{n,m} \phi_m(t) + b_{n,m} \phi_m(t+\tau)
  subroutine KSC_Evolve_Phi(alanc,blanc,phi)
    real(8),dimension(:),intent(in)                   :: alanc
    real(8),dimension(size(alanc)),intent(in)         :: blanc
    complex(8),dimension(:,:),allocatable,intent(out) :: phi
    complex(8),dimension(size(alanc))                 :: state
    real(8)                                           :: tnow,target,dt_step
    integer                                           :: it,nsteps,istep,Nlanc

    Nlanc = size(alanc)
    if(Nlanc==0)stop "ed_get_krylov_phi error: empty Lanczos basis Nlanc=0"
    !
    call allocate_grids
    !
    allocate(phi(Nlanc,Ltimes))
    phi      = zero
    state    = zero
    state(1) = one
    !
    tnow = 0d0
    do it=1,Ltimes
       target = times(it)
       if(target > tnow)then
          nsteps = max(1,ceiling((target-tnow)/dtimes))
          dt_step = (target-tnow)/dble(nsteps)
          do istep=1,nsteps
             call rk4_krylov_step(state,alanc,blanc,dt_step)
          enddo
          tnow = target
       endif
       phi(:,it) = state
    enddo
      !
  end subroutine KSC_Evolve_Phi


  






  






  !#############################################################################
  !
  !                       RUNGE-KUTTA 4 STEP
  !
  !#############################################################################
  subroutine rk4_krylov_step(phi,alanc,blanc,dt_loc)
    complex(8),dimension(:),intent(inout)       :: phi
    real(8),dimension(:),intent(in)             :: alanc
    real(8),dimension(size(alanc)),intent(in)   :: blanc
    real(8),intent(in)                          :: dt_loc
    complex(8),dimension(size(phi))             :: k1,k2,k3,k4,tmp
    !
    call krylov_rhs(phi,alanc,blanc,k1)
    tmp = phi + 0.5d0*dt_loc*k1
    call krylov_rhs(tmp,alanc,blanc,k2)
    tmp = phi + 0.5d0*dt_loc*k2
    call krylov_rhs(tmp,alanc,blanc,k3)
    tmp = phi + dt_loc*k3
    call krylov_rhs(tmp,alanc,blanc,k4)
    phi = phi + dt_loc*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
  end subroutine rk4_krylov_step
  !
  subroutine krylov_rhs(phi,alanc,blanc,dphi)
    complex(8),dimension(:),intent(in)             :: phi
    real(8),dimension(:),intent(in)                :: alanc
    real(8),dimension(size(alanc)),intent(in)      :: blanc
    complex(8),dimension(size(phi)),intent(out)    :: dphi
    complex(8)                                     :: hphi
    integer                                        :: n,Nlanc
    !
    Nlanc = size(alanc)
    dphi  = zero
    do n=1,Nlanc
       hphi = alanc(n)*phi(n)
       if(n>1)hphi = hphi + blanc(n)*phi(n-1)
       if(n<Nlanc)hphi = hphi + blanc(n+1)*phi(n+1)
       dphi(n) = -xi*hphi
    enddo
  end subroutine krylov_rhs











  !#########################################################################
  !
  !                            AUXILIARY ROUTINES
  !
  !#########################################################################

  subroutine KSC_write()
    integer :: iorb,ispin
    !
    if(.not.allocated(times))call allocate_grids()
    if(.not.allocated(ed_Kcdg))return
    call KSC_write_info()
    do ispin=1,Nspin
      do iorb=1,Norb
        call KSC_write_trace("KSC_cdg",iorb,ispin,ed_Kcdg(ispin,iorb,:),ed_Scdg(ispin,iorb,:),ed_Normcdg(ispin,iorb,:))
        call KSC_write_prob("KSC_Pcdg",iorb,ispin,ed_Pcdg(ispin,iorb,:,:))
      enddo
    enddo
    do iorb=1,Norb
      call KSC_write_trace("KSC_N",iorb,1,ed_Kn(1,iorb,:),ed_Sn(1,iorb,:),ed_Normn(1,iorb,:))
      call KSC_write_prob("KSC_PN",iorb,1,ed_Pn(1,iorb,:,:))
    enddo
  end subroutine KSC_write


  subroutine KSC_write_info()
    integer :: unit

    open(free_unit(unit),file="KSC_info"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"# Krylov-state complexity output"
    write(unit,"(A)")"# Trace files: columns are time, K(t), S_K(t), norm(t)."
    write(unit,"(A)")"# Probability files: columns are time, P_1(t), ..., P_Nk(t)."
    write(unit,"(A)")"# cdg files are spin/orbital resolved d^+_{a,sigma}|state> channels."
    write(unit,"(A)")"# N files are total orbital density N_a|state> channels; spin index is fixed to 1 in phase 1."
    write(unit,"(A,I8)")"# Ltimes = ",Ltimes
    write(unit,"(A,I8)")"# KSCmax = ",KSCmax
    close(unit)
  end subroutine KSC_write_info


  subroutine KSC_write_trace(prefix,iorb,ispin,Kt,Sk,Nt)
    character(len=*),intent(in)      :: prefix
    integer,intent(in)               :: iorb,ispin
    real(8),dimension(:),intent(in)  :: Kt,Sk,Nt
    character(len=256)               :: file
    integer                          :: unit,it
    file = trim(prefix)//"_l"//reg(str(iorb))//"_s"//reg(str(ispin))//reg(ed_file_suffix)//".ed"
    open(free_unit(unit),file=file)
    write(unit,"(A)")"# time K(t) S_K(t) norm(t)"
    do it=1,Ltimes
      write(unit,"(4ES24.16)")times(it),Kt(it),Sk(it),Nt(it)
    enddo
    close(unit)
  end subroutine KSC_write_trace


  subroutine KSC_write_prob(prefix,iorb,ispin,Pnk)
    character(len=*),intent(in)        :: prefix
    integer,intent(in)                 :: iorb,ispin
    real(8),dimension(:,:),intent(in)  :: Pnk
    character(len=256)                 :: file
    integer                            :: unit,it,in,Nk
    Nk = size(Pnk,1)
    file = trim(prefix)//"_l"//reg(str(iorb))//"_s"//reg(str(ispin))//reg(ed_file_suffix)//".ed"
    open(free_unit(unit),file=file)
    write(unit,"(A)")"# time P_1(t) ... P_Nk(t)"
    do it=1,Ltimes
      write(unit,"(*(ES24.16,1X))")times(it),(Pnk(in,it),in=1,Nk)
    enddo
    close(unit)
  end subroutine KSC_write_prob



  subroutine KSC_add_weighted_krylov(weight,Ktmp,Stmp,Ptmp,Ntmp,Kt,Sk,Pnk,Nt,spectral_weight)
    real(8),intent(in)                   :: weight
    real(8),dimension(:),intent(in)      :: Ktmp,Stmp,Ntmp
    real(8),dimension(:,:),intent(in)    :: Ptmp
    real(8),dimension(:),intent(inout)   :: Kt,Sk
    real(8),dimension(:,:),intent(inout) :: Pnk
    real(8),dimension(:),intent(inout)   :: Nt
    real(8),intent(inout)                :: spectral_weight
    if(weight<=0d0)return
    !
    Kt = Kt  + weight*Ktmp
    Sk = Sk  + weight*Stmp
    Pnk= Pnk + weight*Ptmp
    Nt = Nt  + weight*Ntmp
    spectral_weight = spectral_weight + weight
    !
  end subroutine KSC_add_weighted_krylov





  function max_krylov_size() result(nmax)
    integer :: nmax
    nmax = max(1,lanc_nGFiter)
    if(allocated(getDim))nmax = max(1,min(lanc_nGFiter,maxval(getDim)))
  end function max_krylov_size



  function thermal_weight(istate) result(peso)
    integer,intent(in) :: istate
    real(8)            :: peso,Ei,Egs,Z
    !
    Ei  = es_return_energy(state_list,istate)
    Egs = state_list%emin
    Z   = zeta_function
    peso= 1d0/zeta_function
    if(finiteT)then
       if(beta*(Ei-Egs) < 200d0)peso = exp(-beta*(Ei-Egs))/zeta_function
    endif
  end function thermal_weight





end MODULE ED_KRYLOV
