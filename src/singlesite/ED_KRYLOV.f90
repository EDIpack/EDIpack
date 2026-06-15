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
  public :: krylov_operator_complexity

  integer :: KSCmax=0
  integer :: KOCmax=0

  character(len=32),dimension(:),allocatable :: KOC_Ops,KSC_Ops    
contains




  !#############################################################################
  !
  !                       Krylov-operator complexity 
  !
  !#############################################################################
  subroutine krylov_operator_complexity(bath,Ops)
    !
    ! Build the T=0 Krylov-operator complexity for d^+_{a,sigma},
    ! d_{a,sigma}, and gamma1=d^+_{a,sigma}+d_{a,sigma}.
    ! Phase 1 supports only the NORMAL channel.
    !
    real(8),dimension(:),intent(in),optional          :: bath
    character(len=*),dimension(:),intent(in),optional :: Ops
    real(8),dimension(:),allocatable                  :: Kt,Sk
    real(8),dimension(:,:),allocatable                :: Pnk
    real(8),dimension(1)                              :: bath_dummy
    integer                                           :: ic,iorb,ispin,Nops
    logical                                           :: bath_allocated_here,check
    !
    !
    if(MpiMaster)write(LOGfile,*)'Krylov operator complexity'
    !
    !Checks:
    if(finiteT)stop "krylov_operator_complexity error: T=0 only."
    if(.not.state_list%status.OR.state_list%size==0)&
      stop "krylov_operator_complexity error: state_list is empty."
    !
    !Check and allocate dmft bath  
    bath_allocated_here = .false.
    if(present(bath))then
      if(Nbath>0)then
        check = check_bath_dimension(bath)
        if(.not.check)stop "krylov_operator_complexity error: wrong bath dimensions"
      endif
      call allocate_dmft_bath()
      call set_dmft_bath(bath)
      bath_allocated_here = .true.
    elseif(.not.dmft_bath%status)then
      if(Nbath>0)stop "krylov_operator_complexity error: no bath given && dmft_bath.status=F && Nbath>0."
      bath_dummy = zero
      call allocate_dmft_bath()
      call set_dmft_bath(bath_dummy)
      bath_allocated_here = .true.
    endif
    !
    call allocate_grids()
    KOCmax = max_krylov_size()
    !
    if(allocated(KOC_Ops))deallocate(KOC_Ops)
    if(present(Ops))then
      Nops    = size(Ops)
      allocate(KOC_Ops(Nops), source=Ops)
    else
      Nops = 3
      allocate(KOC_Ops(3))
      KOC_Ops = [character(len=32) :: "cdg","c","g1"]
    endif 
    !
    !Allocate global KOC variables:
    if(allocated(ed_KOC))   deallocate(ed_KOC)
    if(allocated(ed_SOC))   deallocate(ed_SOC)
    if(allocated(ed_POC))   deallocate(ed_POC)
    allocate(ed_KOC(Nops,Nspin,Norb,Ltimes))
    allocate(ed_SOC(Nops,Nspin,Norb,Ltimes))
    allocate(ed_POC(Nops,Nspin,Norb,KOCmax,Ltimes))
    ed_KOC=0d0
    ed_SOC=0d0
    ed_POC=0d0
    !
    do ispin=1,Nspin
      do iorb=1,Norb
        do ic=1,size(KOC_Ops)
          call KOC_build_Complexity(str(KOC_Ops(ic)),iorb,ispin,Kt,Sk,Pnk)
          ed_KOC(ic,ispin,iorb,:)    = Kt
          ed_SOC(ic,ispin,iorb,:)    = Sk
          ed_POC(ic,ispin,iorb,:,:)  = Pnk
        enddo
      enddo
    enddo
    !
    if(allocated(Kt)) deallocate(Kt)
    if(allocated(Sk)) deallocate(Sk)
    if(allocated(Pnk))deallocate(Pnk)
    !
    if(MPIMASTER)call KOC_write()
    if(bath_allocated_here)call deallocate_dmft_bath()
    !
  end subroutine krylov_operator_complexity

  !The actual KOC workhorse:
  subroutine KOC_build_Complexity(op,iorb,ispin,Kt,Sk,Pnk)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),allocatable,intent(out)   :: Kt
    real(8),dimension(:),allocatable,intent(out)   :: Sk
    real(8),dimension(:,:),allocatable,intent(out) :: Pnk
    real(8),dimension(:),allocatable               :: alanc,blanc
    real(8),dimension(:,:),allocatable             :: Ptmp
    real(8)                                        :: norm2
    integer                                        :: Nk
    !
    if(allocated(Kt)) deallocate(Kt)
    if(allocated(Sk)) deallocate(Sk)
    if(allocated(Pnk))deallocate(Pnk)
    !
    select case(str(ed_mode))
      case default  ;call KOC_Krylov_Basis_normal(op,iorb,ispin,alanc,blanc,norm2)
      case("superc");stop "KOC_build_Complexity error: ed_mode=normal only."
      case("nonsu2");stop "KOC_build_Complexity error: ed_mode=normal only."
    end select
    call Krylov_write_coeff("KOC_coeff_"//str(to_lower(op)),iorb,ispin,alanc,blanc)
    call Krylov_Build_Complexity(alanc,blanc,Kt,Sk,Ptmp)
    allocate(Pnk(KOCmax,Ltimes))
    Pnk= 0d0
    Nk = min(size(Ptmp,1),KOCmax)
    Pnk(1:Nk,:) = Ptmp(1:Nk,:)
    if(allocated(Ptmp)) deallocate(Ptmp)
    if(allocated(alanc))deallocate(alanc)
    if(allocated(blanc))deallocate(blanc)
  end subroutine KOC_build_Complexity





      




  
  !#############################################################################
  !
  !                       Krylov-state complexity 
  !
  !#############################################################################
  subroutine krylov_state_complexity(bath,Ops)
    !
    ! Build the Krylov-state complexity for the impurity problem for the states obtained 
    ! applying d^+_{a,sigma}, d_{a,sigma}, and gamma1=d^+_{a,sigma}+d_{a,sigma}.
    !
    real(8),dimension(:),intent(in),optional          :: bath
    character(len=*),dimension(:),intent(in),optional :: Ops
    real(8),dimension(:),allocatable                  :: Kt,Sk
    real(8),dimension(:,:),allocatable                :: Pnk
    real(8),dimension(1)                              :: bath_dummy
    integer                                           :: ic,iorb,ispin,Nops
    logical                                           :: bath_allocated_here,check
    !
    if(MpiMaster)write(LOGfile,*)'Krylov state complexity'
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
    if(allocated(KSC_Ops))deallocate(KSC_Ops)
    if(present(Ops))then
      Nops    = size(Ops)
      allocate(KSC_Ops(Nops), source=Ops)
    else
      Nops = 3
      allocate(KSC_Ops(3))
      KSC_Ops = [character(len=32) :: "cdg","c","g1"]
    endif 
    !
    if(allocated(ed_KSC))   deallocate(ed_KSC)
    if(allocated(ed_SSC))   deallocate(ed_SSC)
    if(allocated(ed_PSC))   deallocate(ed_PSC)
    allocate(ed_KSC(Nops,Nspin,Norb,Ltimes))
    allocate(ed_SSC(Nops,Nspin,Norb,Ltimes))
    allocate(ed_PSC(Nops,Nspin,Norb,KSCmax,Ltimes))
    !
    ed_KSC   =0d0
    ed_SSC   =0d0
    ed_PSC   =0d0
    !
    do ispin=1,Nspin
      do iorb=1,Norb
        do ic=1,Nops
          call KSC_Build_Complexity(str(KSC_Ops(ic)),iorb,ispin,Kt,Sk,Pnk)
          ed_KSC(ic,ispin,iorb,:)    = Kt
          ed_SSC(ic,ispin,iorb,:)    = Sk
          ed_PSC(ic,ispin,iorb,:,:)  = Pnk
        enddo
      enddo
    enddo
    !
    if(allocated(Kt)) deallocate(Kt)
    if(allocated(Sk)) deallocate(Sk)
    if(allocated(Pnk))deallocate(Pnk)          
    !
    if(MPIMASTER)call KSC_write()
    if(bath_allocated_here)call deallocate_dmft_bath()
    !
    if(MpiMaster)write(LOGfile,*)'...Done'
  end subroutine krylov_state_complexity      


  !THe actual workhorse of Krylov-state complexity
  subroutine KSC_Build_Complexity(op,iorb,ispin,Kt,Sk,Pnk)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),allocatable,intent(out)   :: Kt
    real(8),dimension(:),allocatable,intent(out)   :: Sk
    real(8),dimension(:,:),allocatable,intent(out) :: Pnk
    real(8)                                        :: spectral_weight
    !
    if(KSCmax==0)stop "build_ksc_OpStateList error: KSCmax=0"
    !
    if(allocated(Kt)) deallocate(Kt)
    if(allocated(Sk)) deallocate(Sk)
    if(allocated(Pnk))deallocate(Pnk)  
    !
    allocate(Kt(Ltimes))
    allocate(Sk(Ltimes))
    allocate(Pnk(KSCmax,Ltimes))
    Kt              = 0d0
    Sk              = 0d0
    Pnk             = 0d0
    spectral_weight = 0d0
    if(MPIMASTER)call KSC_init_state_coeff_file(op,iorb,ispin)
    !
    select case(ed_mode)
      case default  ;call KSC_ApplyOp_StateList_normal(op,iorb,ispin,Kt,Sk,Pnk,spectral_weight)
      case("superc");call KSC_ApplyOp_StateList_superc(op,iorb,ispin,Kt,Sk,Pnk,spectral_weight)
      case("nonsu2");call KSC_ApplyOp_StateList_nonsu2(op,iorb,ispin,Kt,Sk,Pnk,spectral_weight)
    end select
    !
    if(spectral_weight==0d0)then
      if(MPIMASTER)write(LOGfile,"(A)")&
      "WARNING KSC_build_OpStateList WARNING: zero spectral weight for op="//&
      trim(op)//", iorb="//str(iorb)//", ispin="//str(ispin)
      return
    endif
    !
    Kt  = Kt/spectral_weight
    Sk  = Sk/spectral_weight
    Pnk = Pnk/spectral_weight
    !
  end subroutine KSC_Build_Complexity








  subroutine KSC_ApplyOp_StateList_normal(op,iorb,ispin,Kt,Sk,Pnk,spectral_weight)
    character(len=*),intent(in)          :: op
    integer,intent(in)                   :: iorb,ispin
    real(8),dimension(:),intent(inout)   :: Kt,Sk
    real(8),dimension(:,:),intent(inout) :: Pnk
    real(8),intent(inout)                :: spectral_weight
    real(8),dimension(:),allocatable     :: v_state,vvinit
    integer                              :: istate,isector,jsector,ialfa
    !
    ialfa = 1
    if(.not.ed_total_ud)ialfa = iorb
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       v_state = es_return_dvec(state_list,istate)
       select case(to_lower(str(op)))
       case("cdg")
          jsector = getCDGsector(ialfa,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("c")
          jsector = getCsector(ialfa,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_C(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("g1")
          jsector = getCDGsector(ialfa,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
          jsector = getCsector(ialfa,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_C(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case default
          stop "ed_build_krylov_operator: op must be 'cdg', 'c', or 'g1'"
       end select
       if(allocated(v_state))deallocate(v_state)
    enddo
  end subroutine KSC_ApplyOp_StateList_normal


  subroutine KSC_ApplyOp_StateList_superc(op,iorb,ispin,Kt,Sk,Pnk,spectral_weight)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),intent(inout)             :: Kt,Sk
    real(8),dimension(:,:),intent(inout)           :: Pnk
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
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("c")
          jsector = getCsector(1,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_C(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("g1")
          jsector = getCDGsector(1,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
          jsector = getCsector(1,ispin,isector)
          if(jsector/=0)then
             vvinit = apply_op_C(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case default
          stop "ed_build_krylov_operator: op must be 'cdg', 'c', or 'g1'"
       end select
       if(allocated(v_state))deallocate(v_state)
    enddo
  end subroutine KSC_ApplyOp_StateList_superc


  subroutine KSC_ApplyOp_StateList_nonsu2(op,iorb,ispin,Kt,Sk,Pnk,spectral_weight)
    character(len=*),intent(in)                    :: op
    integer,intent(in)                             :: iorb,ispin
    real(8),dimension(:),intent(inout)             :: Kt,Sk
    real(8),dimension(:,:),intent(inout)           :: Pnk
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
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("c")
          jsector = getCsector(1,ispin,isector)
          if(Jz_basis)jsector = getCsector_Jz(iorb,ispin,isector)
          if(getN(isector)/=0.and.jsector>=0)then
             vvinit = apply_op_C(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case("g1")
          jsector = getCDGsector(1,ispin,isector)
          if(Jz_basis)jsector = getCDGsector_Jz(iorb,ispin,isector)
          if(getN(isector)/=Nlevels.and.jsector>=0)then
             vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
          jsector = getCsector(1,ispin,isector)
          if(Jz_basis)jsector = getCsector_Jz(iorb,ispin,isector)
          if(getN(isector)/=0.and.jsector>=0)then
             vvinit = apply_op_C(v_state,iorb,ispin,isector,jsector)
             call KSC_AddFromState(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       case default
          stop "ed_build_krylov_operator: op must be 'cdg', 'c', or 'g1'"
       end select
       if(allocated(v_state))deallocate(v_state)
    enddo
  end subroutine KSC_ApplyOp_StateList_nonsu2


  !INTERFACED ABOVE TO SINGLE FUNCTION: KSC_AddState
  subroutine KSC_AddFromState_d(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
    character(len=*),intent(in)          :: op
    integer,intent(in)                   :: iorb,ispin,istate,jsector
    real(8),dimension(:),intent(inout)   :: vvinit
    real(8),dimension(:),intent(inout)   :: Kt,Sk
    real(8),dimension(:,:),intent(inout) :: Pnk
    real(8),intent(inout)                :: spectral_weight
    real(8),dimension(:),allocatable     :: alanc,blanc,Ktmp,Stmp
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
    call Krylov_Build_Complexity(alanc,blanc,Ktmp,Stmp,Ptmp)
    weight = thermal_weight(istate)*norm2
    if(MPIMASTER)call KSC_append_state_coeff(op,iorb,ispin,istate,jsector,weight,alanc,blanc)
    if(weight>0d0)then
      Kt = Kt  + weight*Ktmp
      Sk = Sk  + weight*Stmp
      Pnk= Pnk + weight*Ptmp
      spectral_weight = spectral_weight + weight
    endif 
    ! call KSC_add_weighted_krylov(weight,Ktmp,Stmp,Ptmp,Kt,Sk,Pnk,spectral_weight)
    !
    if(allocated(alanc))deallocate(alanc)
    if(allocated(blanc))deallocate(blanc)
    if(allocated(Ktmp)) deallocate(Ktmp)
    if(allocated(Stmp)) deallocate(Stmp)
    if(allocated(Ptmp)) deallocate(Ptmp)
  end subroutine KSC_AddFromState_d


  subroutine KSC_AddFromState_c(op,iorb,ispin,istate,jsector,vvinit,Kt,Sk,Pnk,spectral_weight)
    character(len=*),intent(in)           :: op
    integer,intent(in)                    :: iorb,ispin,istate,jsector
    complex(8),dimension(:),intent(inout) :: vvinit
    real(8),dimension(:),intent(inout)    :: Kt,Sk
    real(8),dimension(:,:),intent(inout)  :: Pnk
    real(8),intent(inout)                 :: spectral_weight
    real(8),dimension(:),allocatable      :: alanc,blanc,Ktmp,Stmp
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
    call Krylov_Build_Complexity(alanc,blanc,Ktmp,Stmp,Ptmp)
    weight = thermal_weight(istate)*norm2
    if(MPIMASTER)call KSC_append_state_coeff(op,iorb,ispin,istate,jsector,weight,alanc,blanc)
    if(weight>0d0)then
      Kt = Kt  + weight*Ktmp
      Sk = Sk  + weight*Stmp
      Pnk= Pnk + weight*Ptmp
      spectral_weight = spectral_weight + weight
    endif
    ! call KSC_add_weighted_krylov(weight,Ktmp,Stmp,Ptmp,Kt,Sk,Pnk,spectral_weight)
    !
    if(allocated(alanc))deallocate(alanc)
    if(allocated(blanc))deallocate(blanc)
    if(allocated(Ktmp)) deallocate(Ktmp)
    if(allocated(Stmp)) deallocate(Stmp)
    if(allocated(Ptmp)) deallocate(Ptmp)
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








  

  !#############################################################################
  !
  !                    KRYLOV COMPLEXITY (AGNOSTIC)
  !
  !#############################################################################

  !Build the Krylov coefficients dyanamics \phi_n(t) so that \Psi(t) = \sum_n \phi_n(t) |k_n>
  !and return the Krylov complexity Kt, Sk, Pnk
  subroutine Krylov_Build_Complexity(alanc,blanc,Kt,Sk,Pnk)
    real(8),dimension(:),intent(in)                :: alanc
    real(8),dimension(size(alanc)),intent(in)      :: blanc
    real(8),dimension(:),allocatable,intent(out)   :: Kt
    real(8),dimension(:),allocatable,intent(out)   :: Sk
    real(8),dimension(:,:),allocatable,intent(out) :: Pnk
    complex(8),dimension(:,:),allocatable          :: phi
    real(8)                                        :: P
    integer                                        :: it,in,Nlanc
    !
    !Build the Krylov coefficients dyanamics \phi_n(t)
    call Krylov_Evolve_Phi(alanc,blanc,phi)  !<- allocate Phi in here:
    !
    Nlanc=size(alanc)
    if(allocated(Kt)) deallocate(Kt)
    if(allocated(Sk)) deallocate(Sk)
    if(allocated(Pnk))deallocate(Pnk)
    allocate(Kt(Ltimes))       ;Kt =0d0
    allocate(Sk(Ltimes))       ;Sk =0d0
    allocate(Pnk(Nlanc,Ltimes));Pnk=0d0
    !
    do it=1,Ltimes
      do in=1,Nlanc
          P          = abs(phi(in,it))**2
          Pnk(in,it) = P
          Kt(it)     = Kt(it) + dble(in-1)*P
          if(P>0d0)Sk(it)   = Sk(it) - P*log(P)
      enddo
    enddo
    !
    deallocate(phi)
    !
  end subroutine Krylov_Build_Complexity       






  !Build the Krylov coefficients dyanamics \phi_n(t)
  !solving the coupled equations:
  !\phi_n(t) = \sum_m a_{n,m} \phi_m(t) + b_{n,m} \phi_m(t+\tau)
  subroutine Krylov_Evolve_Phi(alanc,blanc,phi)
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
  end subroutine Krylov_Evolve_Phi


  






  






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
    character(len=32) :: opname
    integer           :: ic,iorb,ispin
    !
    if(.not.allocated(times))call allocate_grids()
    if(.not.allocated(ed_KSC))return
    if(.not.allocated(KSC_Ops))return
    call KSC_write_info()
    do ic=1,size(KSC_Ops)
      opname = trim(adjustl(to_lower(KSC_Ops(ic))))
      do ispin=1,Nspin
        do iorb=1,Norb
          call Krylov_write_trace("KSC_"//trim(opname),iorb,ispin,ed_KSC(ic,ispin,iorb,:),ed_SSC(ic,ispin,iorb,:))
          call Krylov_write_prob("KSC_P"//trim(opname),iorb,ispin,ed_PSC(ic,ispin,iorb,:,:))
        enddo
      enddo
    enddo
  end subroutine KSC_write


  subroutine KOC_write()
    character(len=32) :: opname
    integer           :: ic,iorb,ispin
    !
    if(.not.allocated(times))call allocate_grids()
    if(.not.allocated(ed_KOC))return
    if(.not.allocated(KOC_Ops))return
    call KOC_write_info()
    do ic=1,size(KOC_Ops)
      opname = trim(adjustl(to_lower(KOC_Ops(ic))))
      do ispin=1,Nspin
        do iorb=1,Norb
          call Krylov_write_trace("KOC_"//trim(opname),iorb,ispin,ed_KOC(ic,ispin,iorb,:),ed_SOC(ic,ispin,iorb,:))
          call Krylov_write_prob("KOC_P"//trim(opname),iorb,ispin,ed_POC(ic,ispin,iorb,:,:))
        enddo
      enddo
    enddo
  end subroutine KOC_write


  subroutine KOC_write_info()
    integer :: ic,unit

    open(free_unit(unit),file="KOC_info"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"# Krylov-operator complexity output"
    write(unit,"(A)")"# Phase 1: T=0, ed_mode=normal."
    write(unit,"(A)")"# Trace files: columns are time, K_OC(t), S_OC(t)."
    write(unit,"(A)")"# Probability files: columns are time, P_1(t), ..., P_Nk(t)."
    write(unit,"(A)")"# Coefficient files: columns are n, alpha_n, beta_n."
    if(allocated(KOC_Ops))then
      do ic=1,size(KOC_Ops)
        write(unit,"(A,I8,A,A)")"# op(",ic,") = ",trim(adjustl(KOC_Ops(ic)))
      enddo
    endif
    write(unit,"(A,I8)")"# Ltimes = ",Ltimes
    write(unit,"(A,I8)")"# KOCmax = ",KOCmax
    close(unit)
  end subroutine KOC_write_info


  subroutine KSC_write_info()
    integer :: ic,unit

    open(free_unit(unit),file="KSC_info"//reg(ed_file_suffix)//".ed")
    write(unit,"(A)")"# Krylov-state complexity output"
    write(unit,"(A)")"# Trace files: columns are time, K(t), S_K(t)."
    write(unit,"(A)")"# Probability files: columns are time, P_1(t), ..., P_Nk(t)."
    write(unit,"(A)")"# Coefficient files: columns are state, sector, weight, n, alpha_n, beta_n."
    if(allocated(KSC_Ops))then
      do ic=1,size(KSC_Ops)
        write(unit,"(A,I8,A,A)")"# op(",ic,") = ",trim(adjustl(KSC_Ops(ic)))
      enddo
    endif
    write(unit,"(A,I8)")"# Ltimes = ",Ltimes
    write(unit,"(A,I8)")"# KSCmax = ",KSCmax
    close(unit)
  end subroutine KSC_write_info


  subroutine Krylov_write_trace(prefix,iorb,ispin,Kt,Sk)
    character(len=*),intent(in)      :: prefix
    integer,intent(in)               :: iorb,ispin
    real(8),dimension(:),intent(in)  :: Kt,Sk
    character(len=256)               :: file
    integer                          :: unit,it
    file = trim(prefix)//"_l"//reg(str(iorb))//"_s"//reg(str(ispin))//reg(ed_file_suffix)//".ed"
    open(free_unit(unit),file=file)
    write(unit,"(A)")"# time K(t) S_K(t)"
    do it=1,Ltimes
      write(unit,"(3ES24.16)")times(it),Kt(it),Sk(it)
    enddo
    close(unit)
  end subroutine Krylov_write_trace


  subroutine Krylov_write_coeff(prefix,iorb,ispin,alanc,blanc)
    character(len=*),intent(in)     :: prefix
    integer,intent(in)              :: iorb,ispin
    real(8),dimension(:),intent(in) :: alanc,blanc
    character(len=256)              :: file
    integer                         :: unit,in
    file = trim(prefix)//"_l"//reg(str(iorb))//"_s"//reg(str(ispin))//reg(ed_file_suffix)//".ed"
    open(free_unit(unit),file=file)
    write(unit,"(A)")"# n alpha_n beta_n"
    do in=1,size(alanc)
      write(unit,"(I8,2ES24.16)")in,alanc(in),blanc(in)
    enddo
    close(unit)
  end subroutine Krylov_write_coeff


  subroutine KSC_init_state_coeff_file(op,iorb,ispin)
    character(len=*),intent(in) :: op
    integer,intent(in)          :: iorb,ispin
    character(len=256)          :: file
    integer                     :: unit
    file = "KSC_coeff_"//trim(to_lower(op))//"_l"//reg(str(iorb))//"_s"//reg(str(ispin))//reg(ed_file_suffix)//".ed"
    open(free_unit(unit),file=file)
    write(unit,"(A)")"# state sector weight n alpha_n beta_n"
    close(unit)
  end subroutine KSC_init_state_coeff_file


  subroutine KSC_append_state_coeff(op,iorb,ispin,istate,isector,weight,alanc,blanc)
    character(len=*),intent(in)     :: op
    integer,intent(in)              :: iorb,ispin,istate,isector
    real(8),intent(in)              :: weight
    real(8),dimension(:),intent(in) :: alanc,blanc
    character(len=256)              :: file
    integer                         :: unit,in
    file = "KSC_coeff_"//trim(to_lower(op))//"_l"//reg(str(iorb))//"_s"//reg(str(ispin))//reg(ed_file_suffix)//".ed"
    open(free_unit(unit),file=file,position="append")
    do in=1,size(alanc)
      write(unit,"(2I8,ES24.16,I8,2ES24.16)")istate,isector,weight,in,alanc(in),blanc(in)
    enddo
    close(unit)
  end subroutine KSC_append_state_coeff


  subroutine Krylov_write_prob(prefix,iorb,ispin,Pnk)
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
  end subroutine Krylov_write_prob



  ! subroutine KSC_add_weighted_krylov(weight,Ktmp,Stmp,Ptmp,Kt,Sk,Pnk,spectral_weight)
  !   real(8),intent(in)                   :: weight
  !   real(8),dimension(:),intent(in)      :: Ktmp,Stmp
  !   real(8),dimension(:,:),intent(in)    :: Ptmp
  !   real(8),dimension(:),intent(inout)   :: Kt,Sk
  !   real(8),dimension(:,:),intent(inout) :: Pnk
  !   real(8),intent(inout)                :: spectral_weight
  !   if(weight<=0d0)return
  !   !
  !   Kt = Kt  + weight*Ktmp
  !   Sk = Sk  + weight*Stmp
  !   Pnk= Pnk + weight*Ptmp
  !   spectral_weight = spectral_weight + weight
  !   !
  ! end subroutine KSC_add_weighted_krylov





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
