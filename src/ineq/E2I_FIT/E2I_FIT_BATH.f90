MODULE E2I_BATH_FIT
  !:synopsis: Bath fitting routines: real-space DMFT extension
  !Contains routines that fit the Impurity model bath
  USE EDIPACK
  !
  USE E2I_VARS_GLOBAL
  USE E2I_AUX_FUNX
  !
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgplus,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,operator(.x.)
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape
  !
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif


  implicit none
  private

  interface ed_chi2_fitgf
     !This subroutine realizes the :math:`\chi^2` fit of the Weiss field or hybridization function via
     !an impurity model non-interacting Green's function. The bath levels (levels/internal structure
     !and hybridization strength) are supplied by the user in the :f:var:`bath` array
     !and are the parameters of the fit.
     !The function(s) to fit can have different shapes:
     !
     !  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
     !    :math:`\cdot` :f:var:`norb`, :f:var:`lfit`  ]  
     !  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lfit` ] 
     !  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lfit` ]
     !
     !where :f:var:`nlat` is the number of impurity sites in real-space DMFT. Accordingly, the bath array or arrays have rank 2 or 3.
     !Some global variables directly influence the way the fit is performed and can be modified in the input file. See :f:mod:`ed_input_vars`
     !for the description of :f:var:`lfit`, :f:var:`cg_method` , :f:var:`cg_grad`, :f:var:`cg_ftol`, :f:var:`cg_stop` , :f:var:`cg_niter` ,
     !:f:var:`cg_weight` , :f:var:`cg_scheme` , :f:var:`cg_pow` , :f:var:`cg_minimize_ver` , :f:var:`cg_minimize_hh` .      
     !
     module procedure chi2_fitgf_lattice_normal_n3
     module procedure chi2_fitgf_lattice_normal_n4
     module procedure chi2_fitgf_lattice_normal_n6
     module procedure chi2_fitgf_lattice_superc_n3
     module procedure chi2_fitgf_lattice_superc_n4
     module procedure chi2_fitgf_lattice_superc_n6
  end interface ed_chi2_fitgf

  public :: ed_chi2_fitgf


contains




  subroutine chi2_fitgf_lattice_normal_n3(g,bath,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8),intent(inout)       :: bath(:,:)
    complex(8),dimension(:,:,:) :: g
    integer,optional            :: ispin
    !
    complex(8)                  :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                     :: bath_tmp(size(bath,1),size(bath,2))
    integer                     :: i,ispin_
    integer                     :: ilat
    integer                     :: iorb,is,io
    integer                     :: jorb,js,jo
    integer                     :: Nsites
    logical                     :: check_dim
    character(len=5)            :: tmp_suffix
    integer                     :: MPI_ID=0
    integer                     :: MPI_SIZE=1
    logical                     :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    ! do ilat=1+MPI_ID,Nsites,MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    fg = zero
    !
    if(ed_verbose>1)write(Logfile,"(A)")"Chi**2 get G with rank:"//str(rank(g))
    call assert_shape(g,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_normal','g')
    fg = lso2nnn_reshape(g(1:Nsites*Nspin*Norb,1:Nsites*Nspin*Norb,1:Lmats),Nsites,Nspin,Norb,Lmats)
    !
    !
    bath_tmp=0d0
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
       if(ed_verbose>1)write(Logfile,"(A)")"Start Chi**2 fit for site:"//str(ilat)
       !
       call ed_set_suffix(ilat)
       !
       bath_tmp(ilat,:)=bath(ilat,:)
       !
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal_n3

  subroutine chi2_fitgf_lattice_normal_n4(g,bath,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8),intent(inout)         :: bath(:,:)
    complex(8),dimension(:,:,:,:) :: g
    integer,optional              :: ispin
    !
    complex(8)                    :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                       :: bath_tmp(size(bath,1),size(bath,2))
    integer                       :: i,ispin_
    integer                       :: ilat
    integer                       :: iorb,is,io
    integer                       :: jorb,js,jo
    integer                       :: Nsites
    logical                       :: check_dim
    character(len=5)              :: tmp_suffix
    integer                       :: MPI_ID=0
    integer                       :: MPI_SIZE=1
    logical                       :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    ! do ilat=1+MPI_ID,Nsites,MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    fg = zero
    !
    if(ed_verbose>1)write(Logfile,"(A)")"Chi**2 get G with rank:"//str(rank(g))
    call assert_shape(g,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_normal','g')
    do ilat=1,Nsites
       fg(ilat,:,:,:,:,:)  = so2nn_reshape(g(ilat,1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
    enddo
    !
    !
    bath_tmp=0d0
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
       if(ed_verbose>1)write(Logfile,"(A)")"Start Chi**2 fit for site:"//str(ilat)
       call ed_set_suffix(ilat)
       !
       bath_tmp(ilat,:)=bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal_n4

  subroutine chi2_fitgf_lattice_normal_n6(g,bath,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8),intent(inout)             :: bath(:,:)
    complex(8),dimension(:,:,:,:,:,:) :: g
    integer,optional                  :: ispin
    !
    complex(8)                        :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                           :: bath_tmp(size(bath,1),size(bath,2))
    integer                           :: i,ispin_
    integer                           :: ilat
    integer                           :: iorb,is,io
    integer                           :: jorb,js,jo
    integer                           :: Nsites
    logical                           :: check_dim
    character(len=5)                  :: tmp_suffix
    integer                           :: MPI_ID=0
    integer                           :: MPI_SIZE=1
    logical                           :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    ! do ilat=1+MPI_ID,Nsites,MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    fg = zero
    !
    if(ed_verbose>1)write(Logfile,"(A)")"Chi**2 get G with rank:"//str(rank(g))
    call assert_shape(g,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_normal','g')
    fg = g(1:Nsites,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    !
    bath_tmp=0d0
    do ilat = 1+MPI_ID,Nsites,MPI_SIZE
       if(ed_verbose>1)write(Logfile,"(A)")"Start Chi**2 fit for site:"//str(ilat)
       call ed_set_suffix(ilat)
       !
       bath_tmp(ilat,:)=bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_,fmpi=.false.)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),fmpi=.false.)
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_normal_n6









  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################






  subroutine chi2_fitgf_lattice_superc_n3(g,f,bath,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8),intent(inout)       :: bath(:,:)
    complex(8),dimension(:,:,:) :: g,f
    integer,optional            :: ispin
    !
    complex(8)                  :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                     :: bath_tmp(size(bath,1),size(bath,2))
    integer                     :: ilat,i,iorb,ispin_
    integer                     :: Nsites
    logical                     :: check_dim
    character(len=5)            :: tmp_suffix
    integer                     :: MPI_ID=0
    integer                     :: MPI_SIZE=1
    logical                     :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    ! do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    call assert_shape(g,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','g')
    fg(1,:,:,:,:,:,:) = lso2nnn_reshape(g(1:Nsites*Nspin*Norb,1:Nsites*Nspin*Norb,1:Lmats),Nsites,Nspin,Norb,Lmats)
    !
    call assert_shape(f,[Nsites*Nspin*Norb,Nsites*Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','f')
    fg(2,:,:,:,:,:,:) = lso2nnn_reshape(f(1:Nsites*Nspin*Norb,1:Nsites*Nspin*Norb,1:Lmats),Nsites,Nspin,Norb,Lmats)
    !
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
       write(Logfile,"(A)")"ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
       call ed_set_suffix(ilat)
       !
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc_n3






  subroutine chi2_fitgf_lattice_superc_n4(g,f,bath,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8),intent(inout)         :: bath(:,:)
    complex(8),dimension(:,:,:,:) :: g,f
    integer,optional              :: ispin
    !
    complex(8)                    :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                       :: bath_tmp(size(bath,1),size(bath,2))
    integer                       :: ilat,i,iorb,ispin_
    integer                       :: Nsites
    logical                       :: check_dim
    character(len=5)              :: tmp_suffix
    integer                       :: MPI_ID=0
    integer                       :: MPI_SIZE=1
    logical                       :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    ! do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    call assert_shape(g,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','g')
    do ilat=1,Nsites
       fg(1,ilat,:,:,:,:,:)  = so2nn_reshape(g(ilat,1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
    enddo
    call assert_shape(f,[Nsites,Nspin*Norb,Nspin*Norb,Lmats],'chi2_fitgf_generic_superc','f')
    do ilat=1,Nsites
       fg(2,ilat,:,:,:,:,:)  = so2nn_reshape(f(ilat,1:Nspin*Norb,1:Nspin*Norb,1:Lmats),Nspin,Norb,Lmats)
    enddo
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
       write(Logfile,"(A)")"ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
       !
       call ed_set_suffix(ilat)
       !
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc_n4







  subroutine chi2_fitgf_lattice_superc_n6(g,f,bath,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8),intent(inout)             :: bath(:,:)
    complex(8),dimension(:,:,:,:,:,:) :: g,f
    integer,optional                  :: ispin
    !
    complex(8)                        :: fg(2,size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                           :: bath_tmp(size(bath,1),size(bath,2))
    integer                           :: ilat,i,iorb,ispin_
    integer                           :: Nsites
    logical                           :: check_dim
    character(len=5)                  :: tmp_suffix
    integer                           :: MPI_ID=0
    integer                           :: MPI_SIZE=1
    logical                           :: MPI_MASTER=.true.
    !
    write(Logfile,"(A)")""
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    ! do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    call assert_shape(g,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_superc','g')
    fg(1,:,:,:,:,:,:) = g(1:Nsites,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)
    !
    call assert_shape(f,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],'chi2_fitgf_generic_superc','f')
    fg(2,:,:,:,:,:,:) = f(1:Nsites,1:Nspin,1:Nspin,1:Norb,1:Norb,1:Lmats)    
    !
    bath_tmp=0.d0
    !
    do ilat= 1 + MPI_ID, Nsites, MPI_SIZE
       write(Logfile,"(A)")"ed_fit_bath_sites_superc: Start Chi**2 fit for site:"//str(ilat)
       call ed_set_suffix(ilat)
       !
       bath_tmp(ilat,:) = bath(ilat,:)
       call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
       !
       if(present(ispin))then
          call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
       else
          do ispin_=1,Nspin
             call ed_chi2_fitgf(fg(1,ilat,:,:,:,:,:),fg(2,ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin_,fmpi=.false.)
          enddo
       endif
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call AllReduce_MPI(MPI_COMM_WORLD,bath_tmp,bath)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix
  end subroutine chi2_fitgf_lattice_superc_n6




end MODULE E2I_BATH_FIT
















