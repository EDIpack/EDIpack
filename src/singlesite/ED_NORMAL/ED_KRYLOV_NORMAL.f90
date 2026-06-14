!KSC=Krylov-state complexity
!KOC=Krylov-operator complexity
MODULE ED_KRYLOV_NORMAL
  !:synopsis: Low-level Krylov routines, :code:`NORMAL` case
  USE SF_IOTOOLS, only: str,to_lower
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SECTOR
  USE ED_SPARSE_MATRIX
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private

  public :: KSC_Krylov_Basis_normal
  public :: KOC_Krylov_Basis_normal

  !KOC symmetry blocks of a given operator O to be used in Liouvillian evolution:
  ! [L,O] = O*L - L*O ==> O^{qLeft}.H^{qRight} - H^{qLeft}.O^{qRight}
  ! where qLeft and qRight are the sector indices associate to the action 
  ! of the operator O. 
  type koc_block_normal
     integer                 :: left_sector=0
     integer                 :: right_sector=0
     type(sparse_matrix_csr) :: mat
  end type koc_block_normal
  !
  ! KOC as a vector
  type koc_vector_normal
     type(koc_block_normal),dimension(:),allocatable :: block
  end type koc_vector_normal


  real(8),parameter :: OpTol=1d-14

contains

  !Krylov State Complexity: Build the Krylov basis |k_n> from given |v>=O|\psi>
  subroutine KSC_Krylov_Basis_normal(isector,vvinit,alanc,blanc,norm2)
    integer,intent(in)                  :: isector
    real(8),dimension(:),intent(inout)    :: vvinit
    real(8),dimension(:),allocatable     :: alanc
    real(8),dimension(:),allocatable     :: blanc
    real(8)                             :: norm2
    call tridiag_Hv_sector_normal(isector,vvinit,alanc,blanc,norm2)
  end subroutine KSC_Krylov_Basis_normal

  !Kryloc Operator Complexity: Build the Krylov basis |k_n> for a given op vector |O)
  subroutine KOC_Krylov_Basis_normal(op,iorb,ispin,alanc,blanc,norm2)
    character(len=*),intent(in)        :: op
    integer,intent(in)                 :: iorb,ispin
    real(8),dimension(:),allocatable   :: alanc
    real(8),dimension(:),allocatable   :: blanc
    real(8)                            :: norm2
    type(koc_vector_normal)            :: qprev,qcurr,w
    real(8),dimension(:),allocatable   :: atmp,btmp
    real(8)                            :: beta,beta_prev,wnorm
    integer                            :: n,Nlanc,Nmax
    !
    Nmax = max(1,lanc_nGFiter)
    if(allocated(atmp))deallocate(atmp)
    if(allocated(btmp))deallocate(btmp)
    allocate(atmp(Nmax))
    allocate(btmp(Nmax))
    atmp=0d0
    btmp=0d0
    !
    if(MpiMaster)write(LOGfile,*)'Evaluating Krylov basis for operator: ',op
    !
    !Build the initial vector |Op)
    call KOC_build_seed_normal(op,iorb,ispin,qcurr)
    !
    norm2 = KOC_inner_product_normal(qcurr,qcurr)
    if(norm2<=0d0)then
       allocate(alanc(1));alanc=0d0
       allocate(blanc(1));blanc=0d0
       call KOC_delete_vector(qcurr)
       return
    endif
    call KOC_scale_vector(1d0/sqrt(norm2),qcurr)
    !
    !Start the Krylov basis construction
    beta_prev = 0d0
    Nlanc     = 1
    do n=1,Nmax
      print*,"DEBUG",n
      !L|Op_n)=|Op_{n+1}>
      call KOC_apply_liouvillian_normal(qcurr,w)
      !a_n=(Op_n,Op_{n+1})
      atmp(n) = KOC_inner_product_normal(qcurr,w)
      !|Op_{n+1})= -a_n|Op_n)+|Op_{n+1})
      call KOC_axpy_vector(-atmp(n),qcurr,w)
      !b_n:-b_{n-1}|Op_{n-1})+|Op_{n+1})
      if(n>1)call koc_axpy_vector(-beta_prev,qprev,w)
      !(Op_{n+1},Op_{n+1})
      wnorm = koc_inner_product_normal(w,w)
      if(wnorm<0d0.AND.abs(wnorm)<1d-12)wnorm=0d0
      if(wnorm<=1d-24.OR.n==Nmax)then
        Nlanc = n
        exit
      endif
      !
      beta = sqrt(wnorm)
      btmp(n+1) = beta
      call koc_copy_vector(qcurr,qprev)
      call koc_copy_vector(w,qcurr)
      call koc_scale_vector(1d0/beta,qcurr)
      call koc_delete_vector(w)
      beta_prev = beta
      Nlanc = n+1
    enddo
    !
    allocate(alanc(Nlanc),blanc(Nlanc))
    alanc = atmp(1:Nlanc)
    blanc = btmp(1:Nlanc)
    !
    call koc_delete_vector(qprev)
    call koc_delete_vector(qcurr)
    call koc_delete_vector(w)
    deallocate(atmp,btmp)
  end subroutine KOC_Krylov_Basis_normal








  !############################################################
  !
  !        KRYLOV BUILD COMPUTATIONAL ROUTINES
  !
  !############################################################
  !
  !Build the seed vector |Op) to start the Krylov basis construction:
  subroutine KOC_build_seed_normal(op,iorb,ispin,O)
    character(len=*),intent(in) :: op
    integer,intent(in)          :: iorb,ispin
    type(koc_vector_normal)     :: O
    integer                     :: istate,isector,lsector,rsector,ialfa
    !
    call KOC_delete_vector(O)
    !
    ialfa=1
    if(.not.ed_total_ud)ialfa=iorb
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       select case(to_lower(trim(op)))
       case default;stop "KOC_Krylov_Basis_normal error: op must be 'cdg' or 'N'"
       case("cdg")
          lsector = getCDGsector(ialfa,ispin,isector)
          rsector = isector
          if(lsector/=0)call koc_append_seed_block("cdg",iorb,ispin,lsector,rsector,O)
          lsector = isector
          rsector = getCsector(ialfa,ispin,isector)
          if(rsector/=0)call koc_append_seed_block("cdg",iorb,ispin,lsector,rsector,O)
       case("n")
          call koc_append_seed_block("n",iorb,ispin,isector,isector,O)
       end select
    enddo
  end subroutine KOC_build_seed_normal


  !Apply the liouvillian operator to a op-vector |A):
  !|LA) =  H^{qLeft} @ A - A@ H^{qRight}
  subroutine KOC_apply_liouvillian_normal(A,LA)
    type(koc_vector_normal),intent(in) :: A
    type(koc_vector_normal)            :: LA
    type(sparse_matrix_csr)            :: Hl,Hr,HlA,AHr
    integer                            :: ib
    !
    call koc_delete_vector(LA)
    !
    if(.not.allocated(A%block))return
    !
    allocate(LA%block(size(A%block)))
    do ib=1,size(A%block)
      LA%block(ib)%left_sector  = A%block(ib)%left_sector
      LA%block(ib)%right_sector = A%block(ib)%right_sector
      !
      !Build the H^{qLeft/qRight} operator matrices
      call KOC_build_hamiltonian_block(A%block(ib)%left_sector,Hl)
      call KOC_build_hamiltonian_block(A%block(ib)%right_sector,Hr)
      !H^{qLeft} @ A
      call sp_matmul_matrix(Hl,A%block(ib)%mat,HlA,OpTol)
      !A @ H^{qRight}
      call sp_matmul_matrix(A%block(ib)%mat,Hr,AHr,OpTol)
      !H^{qLeft} @ A - A@ H^{qRight}
      call sp_axpy_matrix(1d0,HlA,-1d0,AHr,LA%block(ib)%mat,OpTol)
      !
      !Free memory
      call sp_delete_matrix(Hl)
      call sp_delete_matrix(Hr)
      call sp_delete_matrix(HlA)
      call sp_delete_matrix(AHr)
    enddo
  end subroutine koc_apply_liouvillian_normal


  subroutine koc_build_hamiltonian_block(isector,H)
    integer,intent(in)                    :: isector
    type(sparse_matrix_csr),intent(inout) :: H
    real(8),dimension(:),allocatable      :: vin,Hv
    integer                               :: i,j,Dim
    !
    Dim=getDim(isector)
    !
    if(H%status)call sp_delete_matrix(H)
    !
    call sp_init_matrix(H,Dim,Dim)
    call build_Hv_sector_normal(isector)
    allocate(vin(Dim),Hv(Dim))
    do j=1,Dim
       vin=0d0;vin(j)=1d0;Hv=0d0
       call spHtimesV_p(Dim,vin,Hv)
       do i=1,Dim
          if(abs(Hv(i))>OpTol)call sp_insert_element(H,Hv(i),i,j)
       enddo
    enddo
    deallocate(vin,Hv)
    call delete_Hv_sector_normal()
  end subroutine koc_build_hamiltonian_block






  subroutine koc_append_seed_block(op,iorb,ispin,lsector,rsector,O)
    character(len=*),intent(in)        :: op
    integer,intent(in)                 :: iorb,ispin,lsector,rsector
    type(koc_vector_normal)            :: O
    type(koc_block_normal),allocatable :: tmp(:)
    integer                            :: ib,Nold
    !
    if(lsector<=0.OR.rsector<=0)return
    !
    if(allocated(O%block))then
      !Check if this block already exists:
      do ib=1,size(O%block)
        if(O%block(ib)%left_sector==lsector.AND.&
        O%block(ib)%right_sector==rsector)return
      enddo
      !Extend the block array: (this is not the most efficient way)
      Nold=size(O%block)
      allocate(tmp(Nold+1))
      do ib=1,Nold
        tmp(ib)%left_sector  = O%block(ib)%left_sector
        tmp(ib)%right_sector = O%block(ib)%right_sector
        call sp_copy_matrix(O%block(ib)%mat,tmp(ib)%mat)
        call sp_delete_matrix(O%block(ib)%mat)
      enddo
      deallocate(O%block)
      call move_alloc(tmp,O%block)
    else
      allocate(O%block(1))
      Nold=0
    endif
    !Actual append:
    ib=Nold+1
    O%block(ib)%left_sector  = lsector
    O%block(ib)%right_sector = rsector
    call koc_build_operator_block(op,iorb,ispin,lsector,rsector,O%block(ib)%mat)
    !
  end subroutine koc_append_seed_block

  !Build the operator matrix for a given block
  subroutine KOC_build_operator_block(op,iorb,ispin,lsector,rsector,Omat)
    character(len=*),intent(in)          :: op
    integer,intent(in)                   :: iorb,ispin,lsector,rsector
    type(sparse_matrix_csr),intent(inout):: Omat
    real(8),dimension(:),allocatable     :: vin,vout
    integer                              :: i,j
    !
    if(Omat%status)call sp_delete_matrix(Omat)
    !
    call sp_init_matrix(Omat,getDim(lsector),getDim(rsector))
    !
    allocate(vin(getDim(rsector)))
    do j=1,getDim(rsector)
      vin   =0d0
      vin(j)=1d0
      select case(to_lower(trim(op)))
      case("cdg");vout = apply_op_CDG(vin,iorb,ispin,rsector,lsector)
      case("n")  ;vout = apply_op_N(vin,iorb,rsector)
      end select
      do i=1,size(vout)
        if(abs(vout(i))>OpTol)call sp_insert_element(Omat,vout(i),i,j)
      enddo
      if(allocated(vout))deallocate(vout)
    enddo
    deallocate(vin)
    !
  end subroutine koc_build_operator_block



  !################################################################
  !################################################################
  !
  !      AUXILIARY ROUTINES TO HANDLE KOC VECTORS 
  !
  !################################################################
  !################################################################
  




  !This is the inner product between two koc vectors
  ! it is the implementation of scalar products
  ! between two KOC vectors. at the heart of KOC. 
  function koc_inner_product_normal(A,B) result(dot)
    type(koc_vector_normal),intent(in) :: A,B
    real(8)                            :: dot
    real(8),dimension(:),allocatable   :: gs,Av,Bv
    integer                            :: istate,isector,ib,jb,ngs
    !
    dot=0d0
    if(.not.allocated(A%block).OR..not.allocated(B%block))return
    !
    ngs=0
    do istate=1,state_list%size
       ngs=ngs+1
       isector = es_return_sector(state_list,istate)
       gs      = es_return_dvec(state_list,istate)
       !There should be a better way to match the blocks but hey...
       do ib=1,size(A%block)
          do jb=1,size(B%block)
            !Check if the blocks match:
            if(A%block(ib)%left_sector/=B%block(jb)%left_sector)cycle
            if(A%block(ib)%right_sector/=B%block(jb)%right_sector)cycle
            !
            if(A%block(ib)%right_sector==isector)then
              allocate(Av(A%block(ib)%mat%Nrow),Bv(B%block(jb)%mat%Nrow))
              call sp_matvec_matrix(A%block(ib)%mat,gs,Av)
              call sp_matvec_matrix(B%block(jb)%mat,gs,Bv)
              dot = dot + 0.5d0*dot_product(Av,Bv)
              deallocate(Av,Bv)
            endif
            if(A%block(ib)%left_sector==isector)then
              allocate(Av(A%block(ib)%mat%Ncol),Bv(B%block(jb)%mat%Ncol))
              call sp_Tmatvec_matrix(A%block(ib)%mat,gs,Av)
              call sp_Tmatvec_matrix(B%block(jb)%mat,gs,Bv)
              dot = dot + 0.5d0*dot_product(Av,Bv)
              deallocate(Av,Bv)
            endif
          enddo
       enddo
       if(allocated(gs))deallocate(gs)
    enddo
    if(ngs>0)dot=dot/dble(ngs)
  end function koc_inner_product_normal




  subroutine koc_copy_vector(A,B)
    type(koc_vector_normal),intent(in) :: A
    type(koc_vector_normal)            :: B
    integer                            :: ib
    !
    call koc_delete_vector(B)
    if(.not.allocated(A%block))return
    allocate(B%block(size(A%block)))
    do ib=1,size(A%block)
       B%block(ib)%left_sector  = A%block(ib)%left_sector
       B%block(ib)%right_sector = A%block(ib)%right_sector
       call sp_copy_matrix(A%block(ib)%mat,B%block(ib)%mat)
    enddo
  end subroutine koc_copy_vector


  subroutine koc_scale_vector(alpha,A)
    real(8),intent(in)      :: alpha
    type(koc_vector_normal) :: A
    integer                 :: ib,i
    !
    if(.not.allocated(A%block))return
    !
    do ib=1,size(A%block)
       do i=1,A%block(ib)%mat%Nrow
          if(A%block(ib)%mat%row(i)%size>0)&
             A%block(ib)%mat%row(i)%dvals = alpha*A%block(ib)%mat%row(i)%dvals
       enddo
    enddo
  end subroutine koc_scale_vector

  !Returns: alpha*X + Y
  subroutine koc_axpy_vector(alpha,X,Y)
    real(8),intent(in)              :: alpha
    type(koc_vector_normal),intent(in) :: X
    type(koc_vector_normal)         :: Y
    type(sparse_matrix_csr)         :: tmp
    integer                         :: ib
    !
    if(.not.allocated(X%block))return
    if(.not.allocated(Y%block))then
       call koc_copy_vector(X,Y)
       call koc_scale_vector(alpha,Y)
       return
    endif
    if(size(X%block)/=size(Y%block))stop "koc_axpy_vector error: block mismatch"
    do ib=1,size(Y%block)
       if(X%block(ib)%left_sector/=Y%block(ib)%left_sector.OR.&
          X%block(ib)%right_sector/=Y%block(ib)%right_sector)&
          stop "koc_axpy_vector error: sector mismatch"
       call sp_axpy_matrix(alpha,X%block(ib)%mat,1d0,Y%block(ib)%mat,tmp,OpTol)
       call sp_delete_matrix(Y%block(ib)%mat)
       call sp_copy_matrix(tmp,Y%block(ib)%mat)
       call sp_delete_matrix(tmp)
    enddo
  end subroutine koc_axpy_vector


  subroutine koc_delete_vector(A)
    type(koc_vector_normal) :: A
    integer                 :: ib
    if(.not.allocated(A%block))return
    do ib=1,size(A%block)
       call sp_delete_matrix(A%block(ib)%mat)
    enddo
    deallocate(A%block)
  end subroutine koc_delete_vector



end MODULE ED_KRYLOV_NORMAL
