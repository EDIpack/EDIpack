MODULE ED_BATH_REPLICA
  !:synopsis: Routines for replica bath manipulation
  !In these module we implement the functions to set the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` and the initial variational parameters :math:`\vec{\lambda}` used to decompose each local bath hamiltonian for the  :f:var:`replica` and :f:var:`general` bath types. Note that the two bath types share the very same matrix basis, yet for consistency we provide distinct interfaces to the two for each procedure.    
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  !
  USE ED_BATH_AUX
  implicit none

  private


  !This part deals with the matrix structure of the bath components in the replica/general channel.
  !The parameters \lambda are handled by the effective_bath structure in ED_BATH_DMFT. 
  type Omatrix
     ! The matrix of rank [:f:var:`nspin` , :f:var:`nspin` , :f:var:`norb` , :f:var:`norb` ] corresponding to a single element of the matrix basis decomposing the replica/general bath Hamiltonian
     ! :math:`H_p=\sum_{i=1}^{N_{basis}} \lambda_i(p) O_i`, where :math:`N_{basis}` is the dimension of the user defined basis.  
     complex(8),dimension(:,:,:,:),allocatable :: O !Replica/General hamiltonian
  end type Omatrix

  type Hreplica
     ! The object storing all the elements :math:`O_i` and the initial :math:`\lambda` parameters determining each bath Hamiltonian
     ! decomposition: :math:`H_p=\sum_{i=1}^{N_{basis}} \lambda_i(p) O_i`.     
     type(Omatrix),dimension(:),allocatable    :: basis   ![Nsym]
     real(8),dimension(:,:),allocatable        :: linit   ![Nbath,Nsym]
     integer                                   :: nsym=0
     logical                                   :: status=.false.
  end type Hreplica

  type(Hreplica),public                        :: Hb !a default, global and shared instance of :f:var:`Hreplica` to store matrix decomposition






  interface set_Hreplica
     !
     ! This function sets the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` used to decompose the single bath hamiltonian :math:`h^p`. It also sets the initial values of the variational parameters :math:`\vec{\lambda}` in the :f:var:`replica` ath type.
     !
     ! Input: :f:var:`Hvec` 
     !  * rank-5, dimensions: [ |Nns| , |Nns|, |Norb| , |Norb| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nnso| , |Nnso|, |Nsym| ]
     !
     ! Input: :f:var:`lvec`
     !  * rank-2, dimensions: [ |Nbath| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nlat| , |Nbath| , |Nsym| ]
     !
     module procedure :: init_Hreplica_symmetries_d5
     module procedure :: init_Hreplica_symmetries_d3
     module procedure :: init_Hreplica_symmetries_legacy
  end interface set_Hreplica



  interface allocate_Hgeneral
     ! A clone of :f:var:`allocate_hreplica` for the :code:`GENERAL` bath
     module procedure :: allocate_Hreplica
  end interface allocate_Hgeneral

  interface deallocate_Hgeneral
     ! A clone of :f:var:`deallocate_hreplica` for the :code:`GENERAL` bath
     module procedure :: deallocate_Hreplica
  end interface deallocate_Hgeneral


  interface save_Hgeneral
     ! A clone of :f:var:`save_hreplica` for the :code:`GENERAL` bath
     module procedure :: save_Hreplica
  end interface save_Hgeneral

  interface read_Hgeneral
     ! A clone of :f:var:`read_hreplica` for the :code:`GENERAL` bath
     module procedure :: read_Hreplica
  end interface read_Hgeneral

  interface set_Hgeneral
     ! A clone of :f:var:`set_hreplica` for the :code:`GENERAL` bath
     module procedure :: init_Hreplica_symmetries_d5
     module procedure :: init_Hreplica_symmetries_d3
     module procedure :: init_Hreplica_symmetries_legacy
  end interface set_Hgeneral


  interface set_linit_Hgeneral
     ! A clone of :f:var:`set_linit_hreplica` for the :code:`GENERAL` bath
     module procedure :: set_linit_Hreplica
  end interface set_linit_Hgeneral


  interface set_hsym_Hgeneral
     ! A clone of :f:var:`set_hsym_hreplica` for the :code:`GENERAL` bath
     module procedure :: set_hsym_Hreplica
  end interface set_hsym_Hgeneral


  interface build_Hgeneral
     ! A clone of :f:var:`build_hreplica` for the :code:`GENERAL` bath
     module procedure :: build_Hreplica
  end interface build_Hgeneral


  interface print_Hgeneral
     ! A clone of :f:var:`print_hreplica` for the :code:`GENERAL` bath
     module procedure :: print_Hreplica
  end interface print_Hgeneral


  interface Hgeneral_mask
     ! A clone of :f:func:`hreplica_mask` for the :code:`GENERAL` bath
     module procedure :: Hreplica_mask
  end interface Hgeneral_mask




  public :: allocate_Hreplica
  public :: deallocate_Hreplica
  public :: set_Hreplica
  public :: save_Hreplica
  public :: read_Hreplica
  public :: set_linit_Hreplica
  public :: set_hsym_Hreplica
  public :: build_Hreplica
  public :: print_Hreplica
  public :: Hreplica_mask


  public :: allocate_Hgeneral
  public :: deallocate_Hgeneral
  public :: set_Hgeneral
  public :: save_Hgeneral
  public :: read_Hgeneral
  public :: set_linit_Hgeneral
  public :: set_hsym_Hgeneral
  public :: build_Hgeneral
  public :: print_Hgeneral
  public :: Hgeneral_mask


  integer :: ibath

contains




  !##################################################################
  !
  !     H_REPLICA ROUTINES:
  !
  !##################################################################

  subroutine allocate_Hreplica(Nsym)
    ! Allocate global basis for :f:var:`Hb` and the initial vector coefficient.
    integer          :: Nsym    !Number of elements in the matrix decomposition
    integer          :: isym
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG allocate_Hreplica"
#endif
    if(allocated(Hb%basis))deallocate(Hb%basis)
    if(allocated(Hb%linit))deallocate(Hb%linit)
    !
    Hb%Nsym = Nsym
    allocate(Hb%basis(Nsym))
    allocate(Hb%linit(Nbath,Nsym))
    do isym=1,Hb%Nsym
       allocate(Hb%basis(isym)%O(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb))
       Hb%basis(isym)%O=zero
       Hb%linit(:,isym)=1d0
    enddo
    Hb%status=.true.
  end subroutine allocate_Hreplica


  subroutine deallocate_Hreplica()
    ! Deallocate the global :f:var:`Hb` as well as the initial vector coefficient.
    integer              :: isym
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG deallocate_Hreplica"
#endif
    if(.not.Hb%status)return
    if(allocated(Hb%basis))then
       do isym=1,size(Hb%basis)
          if(allocated(Hb%basis(isym)%O))deallocate(Hb%basis(isym)%O)
       enddo
       deallocate(Hb%basis)
    endif
    if(allocated(Hb%linit))deallocate(Hb%linit)
    Hb%Nsym=0
    Hb%status=.false.
  end subroutine deallocate_Hreplica



  subroutine save_Hreplica(file)
    ! Save :f:var:`Hb` to file. 
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*)                                          :: file !Name of the file to store :f:var:`Hb`
    integer                                                   :: unit
    integer                                                   :: iorb,jorb
    integer                                                   :: ispin,jspin
    integer                                                   :: isym
    !
    open(free_unit(unit),file=reg(file))
    write(LOGfile,"(A)")"save Hb.basis to file :"//reg(file)
    !
    write(unit,*)Hb%Nsym
    !
    do isym=1,Hb%Nsym
       !
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             write(unit,*)((&
                  Hb%basis(isym)%O(ispin,jspin,iorb,jorb),&
                  jorb =1,Norb),jspin=1,Nnambu*Nspin)
          enddo
       enddo
       write(unit,*)""
       !
    enddo
    !
    close(unit)
  end subroutine save_Hreplica


  subroutine read_Hreplica(file)
    ! Read :f:var:`Hb` from specified file. 
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*)                                          :: file !Name of the file to store :f:var:`Hb`
    integer                                                   :: unit
    integer                                                   :: iorb,jorb
    integer                                                   :: ispin,jspin
    integer                                                   :: isym,Nsym
    !
    open(free_unit(unit),file=reg(file))
    write(LOGfile,"(A)")"read Hb.basis from file :"//reg(file)
    !
    read(unit,*)Nsym
    if(.not.Hb%status)call allocate_Hreplica(Nsym)
    if(Hb%Nsym/=Nsym)stop "read_Hreplica error: Nsym != Hb.Nsym"
    !
    do isym=1,Hb%Nsym
       !
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             read(unit,*)((&
                  Hb%basis(isym)%O(ispin,jspin,iorb,jorb),&
                  jorb =1,Norb),jspin=1,Nnambu*Nspin)
          enddo
       enddo
       read(unit,*)
       !
    enddo
    !
    close(unit)
  end subroutine read_Hreplica




  !
  subroutine init_Hreplica_symmetries_d5(Hvec,lvec)
    complex(8),dimension(:,:,:,:,:) :: Hvec ![Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym]
    real(8),dimension(:,:)          :: lvec ![Nbath,Nsym]
    integer                         :: isym
    integer                         :: Nsym
    logical                         :: bool
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_site: from {[Hs,Lam]}_b"
#endif
    !
    if(size(lvec,1)/=Nbath)then
       write(LOGfile,*) "                                                                               "
       write(LOGfile,*) "ERROR: if you are trying to init Hreplica for inequivalent sites please note   "
       write(LOGfile,*) "       that the lvec array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(LOGfile,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(LOGfile,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(LOGfile,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    endif
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym = size(lvec,2)
    !
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries_d5","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default  ; bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc"); bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_d5 ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       Hb%basis(isym)%O = Hvec(:,:,:,:,isym)
       Hb%linit(:,isym) = Lvec(:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(LOGfile,*) "Hreplica #"//str(ibath)//":"
          call print_Hreplica( build_Hreplica( Hb%linit(ibath,:) ) )
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_d5

  subroutine init_Hreplica_symmetries_d3(Hvec,lvec)
    complex(8),dimension(:,:,:) :: Hvec      ![Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym]
    real(8),dimension(:,:)      :: lvec ![Nbath,Nsym]
    integer                     :: isym
    integer                     :: Nsym
    logical                     :: bool
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_d3: from {[Hs,Lam]}_b"
#endif
    !
    if(size(lvec,1)/=Nbath)then
       write(LOGfile,*) "                                                                               "
       write(LOGfile,*) "ERROR: if you are trying to init Hreplica for inequivalent sites please note   "
       write(LOGfile,*) "       that the lvec array /MUST/ have [Nineq]x[Nbath]x[Nsym] shape.      "
       write(LOGfile,*) "       The legacy [Nineq]x[Nsym] is not supported anymore, for it would shadow "
       write(LOGfile,*) "       the new recommended [Nbath]x[Nsym] shape for the single impurity case.  "
       write(LOGfile,*) "                                                                               "
       stop ! This unfortunately still leaves room for nasty problems if Nbath==Nineq, but that's it...
    endif
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym = size(lvec,2)
    !
    call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hreplica_symmetries_d3","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default  ; bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
       case("superc"); bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    do isym=1,Nsym
       Hb%basis(isym)%O = so2nn_reshape(Hvec(:,:,isym), Nnambu*Nspin,Norb)
       Hb%linit(:,isym) = lvec(:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(LOGfile,"(A)") "Hreplica #"//str(ibath)//":"
          call print_Hreplica( build_Hreplica( Hb%linit(ibath,:) ) )
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_d3


  subroutine init_Hreplica_symmetries_legacy(Hvec,lvec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:)            :: lvec ![Nsym]
    integer                         :: isym,Nsym
    logical                         :: bool
    !
    if(ed_mode=="superc")Nnambu=2
    Nsym=size(lvec)
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_legacy: from {[Hs,Lam]}_b"
#endif
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default  ;bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc");bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call allocate_hreplica(Nsym)
    !
    !> BACK-COMPATIBILITY PATCH (cfr init_dmft_bath)
    do isym=1,Nsym
       Hb%basis(isym)%O = Hvec(:,:,:,:,isym)
       do ibath=1,Nbath
          Hb%linit(ibath,isym) = lvec(isym)
       enddo
    enddo
    !
    ! PRINT DEPRECATION MESSAGE TO LOG
    write(LOGfile,*) "                                                                               "
    write(LOGfile,*) "WARNING: Passing a single lambdasym vector to ed_set_Hreplica is /deprecated/. "
    write(LOGfile,*) "         You should instead define a different lambda for each bath component, "
    write(LOGfile,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
    write(LOGfile,*) "         Your single lambda vector has been internally copied into the required"
    write(LOGfile,*) "         higher-rank array, so giving each replica the same set of lambdas.    "
    write(LOGfile,*) "         >>> This back-compatibility patch might be removed in a future update."
    write(LOGfile,*) "                                                                               "
    !
    if(ed_verbose>2)then
       do ibath=1,Nbath
          write(LOGfile,*) "Hreplica #"//str(ibath)//":"
          call print_Hreplica( build_Hreplica( Hb%linit(ibath,:)) )
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_legacy








  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  subroutine set_linit_Hreplica(lvec)
    !
    !This function is used to set the initial value of :math:`\lambda` parameters in the bath Hamiltonian matrix decomposition
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    real(8),dimension(:,:) :: lvec   !the input vector of bath parameters [Nsym,Nbath]
    !
    if(.not.Hb%status)stop "set_linit_Hreplica error: Hb.status=F"
    call assert_shape(lvec,[Nbath,Hb%Nsym],"set_linit_Hreplica","lvec")
    Hb%linit = lvec       
  end subroutine set_linit_Hreplica


  subroutine set_hsym_Hreplica(isym,Hsym)
    !
    !This function is used to set a matrix element in the basis of the bath Hamiltonian matrix decomposition
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    integer                       :: isym
    complex(8),dimension(:,:,:,:) :: Hsym      ![Nambu*Nspin,Nambu*Nspin,Norb,Norb]
    !
    if(.not.Hb%status)stop "set_basis_Hreplica error: Hb.status=F"
    !
    if(ed_mode=="superc")Nnambu=2
    call assert_shape(Hsym,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb],"set_hsym_Hreplica","Hsym")
    Hb%basis(isym)%O  = Hsym
  end subroutine set_hsym_Hreplica



  function build_Hreplica(lvec) result(H)
    !
    !This function is used to reconstruct the local bath Hamiltonian from basis expansion given the vector of :math:`\vec{\lambda}` parameters :math:`h^p=\sum_i \lambda^p_i O_i`. The resulting Hamiltonian has dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| ]
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    real(8),dimension(:),optional                             :: lvec   !The input vector of bath parameters
    real(8),dimension(:),allocatable                          :: lambda
    integer                                                   :: isym,nsym
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H !The replica matrix
    !
    if(.not.Hb%status)stop "ERROR build_Hreplica: Hb.basis is not setup"
    !
    Nsym = Hb%Nsym ! == size(Hb%basis)
    allocate(lambda(Nsym))
    lambda=1d0
    if(present(lvec))then
       if(size(lvec)/=Nsym) stop "ERROR build_Hreplica: size(lvec) != Hb.Nsym"
       lambda = lvec
    endif
    !
    H=zero
    do isym=1,size(lambda)
       H = H + lambda(isym)*Hb%basis(isym)%O
    enddo
  end function build_Hreplica




  subroutine print_Hreplica(H,file)
  !Prints the replica Hamiltonian to a file or to standard output
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H !The matrix to print
    character(len=*),optional                                 :: file !Filename
    integer                                                   :: iorb,jorb,ispin,jspin,Nso,unit
    !
    unit=LOGfile
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hbuild to file :"//reg(file)
    endif
    !
    do ispin=1,Nnambu*Nspin
       do iorb=1,Norb
          write(unit,"(*(A1,F8.4,A1,F8.4,A1,2x))")&
               (&
               (&
               '(',dreal(H(ispin,jspin,iorb,jorb)),',',dimag(H(ispin,jspin,iorb,jorb)),')',&
               jorb =1,Norb),&
               jspin=1,Nnambu*Nspin)
       enddo
    enddo
    write(unit,*)""
    !
    if(present(file))close(unit)
  end subroutine print_Hreplica






  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################



  function Hreplica_mask(wdiag,uplo) result(Hmask)
! Creates a logical matrix :f:var:`Hmask` which is :code:`True` where the bath replica matrix
! is nonzero. If :f:var:`wdiag` is :code:`True` the major diagonal is set to :code:`True` -
!If :f:var:`uplo` is :code:`True` only the upper triangular part is kept, the rest
! is set to :code:`False`
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    logical,optional                                          :: wdiag !Flag to set the major diagonal to :code:`True`
    logical,optional                                          :: uplo !Flag to only set the upper triangular matrix
    logical                                                   :: wdiag_,uplo_
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: H
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)    :: Hmask !The resulting logical matrix
    integer                                                   :: iorb,jorb,ispin,jspin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    H = build_Hreplica( Hb%linit(Nbath,:) )
    Hmask=.false.
    where(abs(H)>1d-6)Hmask=.true.
    !
    if(ED_MODE=="normal" .and. Nspin==2)then
      Hmask(1,2,:,:)=.false.
      Hmask(2,1,:,:)=.false.
    endif
    !
    !
    if(wdiag_)then
       do ispin=1,Nnambu*Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nnambu*Nspin
          do jspin=1,Nnambu*Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = index_stride_so(ispin,iorb)
                   jo = index_stride_so(jspin,jorb)
                   if(io>jo)Hmask(ispin,jspin,iorb,jorb)=.false.
                enddo
             enddo
          enddo
       enddo
    endif
    !
  end function Hreplica_mask




END MODULE ED_BATH_REPLICA















! subroutine Hreplica_site(site)
! #if __INTEL_COMPILER
!   use ED_INPUT_VARS, only: Nspin,Norb,Nbath
! #endif
!   integer :: site
!   if(site<1.OR.site>size(Hreplica_lambda_ineq,1))stop "ERROR Hreplica_site: site not in [1,Nlat]"
!   if(.not.allocated(Hreplica_lambda_ineq))stop "ERROR Hreplica_site: Hreplica_lambda_ineq not allocated"
!   Hreplica_lambda(:,:)  = Hreplica_lambda_ineq(site,:,:)
! end subroutine Hreplica_site




! subroutine init_Hreplica_symmetries_lattice_d5(Hvec,lvec)
!   complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
!   real(8),dimension(:,:,:)        :: lvec ![Nlat,Nbath,Nsym]
!   integer                         :: ilat,Nlat
!   integer                         :: isym,Nsym
!   logical                         :: bool
!   !
! #ifdef _DEBUG
!   if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
! #endif
!   if(ed_mode=="superc")Nnambu=2
!   Nlat=size(lvec,1)
!   Nsym=size(lvec,3)
!   !
!   call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
!   !
!   !CHECK NAMBU and HERMITICTY of each Hvec
!   do isym=1,Nsym
!      select case(ed_mode)
!      case default
!         bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
!      case("superc")
!         bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
!      end select
!      if(.not.bool)then
!         write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
!         stop
!      endif
!   enddo
!   !
!   if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
!   allocate(Hreplica_lambda_ineq(Nlat,Nbath,Nsym))
!   call allocate_hreplica(Nsym)
!   !
!   do isym=1,Nsym
!      Hreplica_lambda_ineq(:,:,isym)  = lvec(:,:,isym)
!      Hreplica_basis(isym)%O          = Hvec(:,:,:,:,isym)
!   enddo
!   !
!   if(ed_verbose>2)then
!      do ilat=1,Nlat
!         write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
!         do ibath=1,Nbath
!            write(LOGfile,"(A)")"> Hreplica #"//str(ibath)//":"
!            call print_Hbuild(build_Hreplica(Hreplica_lambda_ineq(ilat,ibath,:)))
!         enddo
!      enddo
!   endif
!   !
! end subroutine init_Hreplica_symmetries_lattice_d5

! subroutine init_Hreplica_symmetries_lattice_d3(Hvec,lvec)
!   complex(8),dimension(:,:,:) :: Hvec      ![size(H),Nsym]
!   real(8),dimension(:,:,:)    :: lvec ![Nlat,Nbath,Nsym]
!   integer                     :: ilat,Nlat
!   integer                     :: isym,Nsym
!   logical                     :: bool
!   !
! #ifdef _DEBUG
!   if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
! #endif
!   if(ed_mode=="superc")Nnambu=2
!   Nlat=size(lvec,1)
!   Nsym=size(lvec,3)
!   !
!   call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hreplica_symmetries","Hvec")
!   !
!   !CHECK NAMBU and HERMITICTY of each Hvec
!   do isym=1,Nsym
!      select case(ed_mode)
!      case default
!         bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
!      case("superc")
!         bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
!      end select
!      if(.not.bool)then
!         write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
!         stop
!      endif
!   enddo
!   !
!   if(allocated(Hreplica_lambda_ineq))deallocate(Hreplica_lambda_ineq)
!   allocate(Hreplica_lambda_ineq(Nlat,Nbath,Nsym))
!   call allocate_hreplica(Nsym)
!   !
!   do isym=1,Nsym
!      Hreplica_lambda_ineq(:,:,isym)  = lvec(:,:,isym)
!      Hreplica_basis(isym)%O          = so2nn_reshape(Hvec(:,:,isym),Nnambu*Nspin,Norb)
!   enddo
!   !
!   if(ed_verbose>2)then
!      do ilat=1,Nlat
!         write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
!         do ibath=1,Nbath
!            write(LOGfile,"(A)")"> Hreplica #"//str(ibath)//":"
!            call print_Hbuild(build_Hreplica(Hreplica_lambda_ineq(ilat,ibath,:)))
!         enddo
!      enddo
!   endif
!   !
! end subroutine init_Hreplica_symmetries_lattice_d3


! subroutine init_Hreplica_symmetries_legacy(Hvec,lvec)
!   complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
!   real(8),dimension(:)            :: lvec ![Nsym]
!   integer                         :: isym,Nsym
!   logical                         :: bool
!   !
!   if(ed_mode=="superc")Nnambu=2
!   Nsym=size(lvec)
! #ifdef _DEBUG
!   if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_legacy: from {[Hs,Lam]}_b"
! #endif
!   call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
!   !
!   !CHECK NAMBU and HERMITICTY of each Hvec
!   do isym=1,Nsym
!      select case(ed_mode)
!      case default
!         bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
!      case("superc")
!         bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
!      end select
!      if(.not.bool)then
!         write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
!         stop
!      endif
!   enddo
!   !
!   call allocate_hreplica(Nsym)
!   !
!   do isym=1,Nsym
!      do ibath=1,Nbath
!         !> BACK-COMPATIBILITY PATCH (cfr init_dmft_bath)
!         Hreplica_lambda(ibath,isym) = lvec(isym)
!      enddo
!      Hreplica_basis(isym)%O = Hvec(:,:,:,:,isym)
!   enddo
!   !
!   ! PRINT DEPRECATION MESSAGE TO LOG
!   write(LOGfile,*) "                                                                               "
!   write(LOGfile,*) "WARNING: Passing a single lambdasym vector to ed_set_Hreplica is /deprecated/. "
!   write(LOGfile,*) "         You should instead define a different lambda for each bath component, "
!   write(LOGfile,*) "         namely passing a [Nbath]x[Nsym] array instead of a [Nsym] vector.     "
!   write(LOGfile,*) "         Your single lambda vector has been internally copied into the required"
!   write(LOGfile,*) "         higher-rank array, so giving each replica the same set of lambdas.    "
!   write(LOGfile,*) "         >>> This back-compatibility patch might be removed in a future update."
!   write(LOGfile,*) "                                                                               "
!   !
!   if(ed_verbose>2)then
!      do ibath=1,Nbath
!         write(LOGfile,*) "Hreplica #"//str(ibath)//":"
!         call print_Hbuild(build_Hreplica(Hreplica_lambda(ibath,:)))
!      enddo
!   endif
!   !
! end subroutine init_Hreplica_symmetries_legacy
