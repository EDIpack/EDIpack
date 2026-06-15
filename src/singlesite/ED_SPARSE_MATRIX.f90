MODULE ED_SPARSE_MATRIX
  !:synopsis: Types and routines for sparse matrices
  !This class defines a data structure to efficiently store sparse matrices into dedicated CSR matrices, featuring support to MPI parallel storage, so that each matrix is spread across the threads.
  !
  !
  USE ED_INPUT_VARS
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private

  complex(8),parameter :: zero=(0d0,0d0)

  type sparse_row_csr
     !
     ! The sparse row data structure containing the non-zero elements on any given row of sparse matrix.
     !
     integer                                   :: size !the current size of the row, corresponding to number of non-zero elements
     real(8),dimension(:),allocatable          :: dvals !rank-1 array for double precision values
     complex(8),dimension(:),allocatable       :: cvals !rank-1 array for double complex values
     integer,dimension(:),allocatable          :: cols  !rank-1 array for column indices
  end type sparse_row_csr

  type sparse_matrix_csr
     !
     ! The sparse matrix data structure realized as an allocatable array of of :f:var:`sparse_row_csr` types.  
     !
     type(sparse_row_csr),dimension(:),pointer :: row !the array of :f:var:`sparse_row_csr` 
     integer                                   :: Nrow !the total number of rows
     integer                                   :: Ncol !the total number of columns
     logical                                   :: status=.false. !Allocation status
#ifdef _MPI
     type(sparse_row_csr),dimension(:),pointer :: loc !the array of :f:var:`sparse_row_csr` for the diagonal blocks
     integer                                   :: istart=0 !starting index for the MPI decomposition of the sparse matrix
     integer                                   :: iend=0   !ending index for the MPI decomposition of the sparse matrix
     integer                                   :: ishift=0 !integer shift index for the MPI decomposition of the sparse matrix
     logical                                   :: mpi=.false.
#endif
  end type sparse_matrix_csr



  !INIT SPARSE MATRICES 
  interface sp_init_matrix
     !
     ! Initialization of the :f:var:`sparse_matrix_csr` via memory allocation. An empty matrix is returned 
     !
     module procedure :: sp_init_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_init_matrix_csr
#endif
  end interface sp_init_matrix



  !DELETE SPARSE MATRIX 
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_delete_matrix_csr
#endif
  end interface sp_delete_matrix


  !INSERT ELEMENTS
  interface sp_insert_element
     !
     ! Matrix element insertion at a given row and columns. If an active MPI communicator is passed as input the element is stored in the matrix chunk of the corresponding thread. Double precision and double complex values are supported.
     !
     !
     module procedure :: sp_insert_element_csr_d
     module procedure :: sp_insert_element_csr_c
     module procedure :: mpi_sp_insert_element_csr_d
     module procedure :: mpi_sp_insert_element_csr_c
  end interface sp_insert_element


  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     !
     ! Dump the :f:var:`sparse_matrix_csr` into a dense matrix. 
     !
     module procedure :: sp_dump_matrix_csr_d
     module procedure :: sp_dump_matrix_csr_c
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_csr_d
     module procedure :: mpi_sp_dump_matrix_csr_c
#endif
  end interface sp_dump_matrix

  interface sp_copy_matrix
     module procedure :: sp_copy_matrix_csr_d
  end interface sp_copy_matrix

  interface sp_axpy_matrix
     module procedure :: sp_axpy_matrix_csr_d
  end interface sp_axpy_matrix

  interface sp_axpy_matrix_inplace
     module procedure :: sp_axpy_matrix_inplace_csr_d
  end interface sp_axpy_matrix_inplace

  interface sp_scale_matrix
     module procedure :: sp_scale_matrix_csr_d
  end interface sp_scale_matrix

  interface sp_matmul_matrix
     module procedure :: sp_matmul_matrix_csr_d
  end interface sp_matmul_matrix

  interface sp_matvec_matrix
     module procedure :: sp_matvec_matrix_csr_d
  end interface sp_matvec_matrix


  interface sp_Tmatvec_matrix
      module procedure :: sp_transpose_matvec_matrix_csr_d
  end interface sp_Tmatvec_matrix

#ifdef _MPI  
  interface sp_set_mpi_matrix
     !
     !  Set up the MPI parameters in the :f:var:`sparse_matrix_csr` for automatic spread of the values across the threads
     !
     module procedure :: sp_set_mpi_matrix_csr
  end interface sp_set_mpi_matrix
#endif

  !Linked-List Sparse Matrix
  public :: sparse_matrix_csr

  public :: sp_init_matrix 
  public :: sp_delete_matrix
  public :: sp_insert_element
  public :: sp_dump_matrix
  public :: sp_copy_matrix
  public :: sp_axpy_matrix
  public :: sp_axpy_matrix_inplace
  public :: sp_scale_matrix
  public :: sp_matmul_matrix
  public :: sp_matvec_matrix
  public :: sp_Tmatvec_matrix
#ifdef _MPI
  public :: sp_set_mpi_matrix
#endif





  interface add_to
     module procedure :: add_to_I
     module procedure :: add_to_D
     module procedure :: add_to_Z
  end interface add_to



  integer :: MpiIerr





contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_csr(sparse,N,N1)
    type(sparse_matrix_csr),intent(inout) :: sparse !sparse matrix to be initialized
    integer                               :: N      !Number of rows 
    integer,optional                      :: N1     !Number of columns [Optional]. If not present :code:`N1=N`
    integer                               :: i
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG sp_init_matrix_csr: allocate sparse"
#endif
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix: already allocated can not init"
    !
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       sparse%row(i)%size=0
       allocate(sparse%row(i)%dvals(0))
       allocate(sparse%row(i)%cvals(0))
       allocate(sparse%row(i)%cols(0))
    end do
    !
    sparse%status=.true.
    !
  end subroutine sp_init_matrix_csr



  subroutine mpi_sp_init_matrix_csr(MpiComm,sparse,N,N1)
    integer                               :: MpiComm !MPI global communicator
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                               :: N
    integer,optional                      :: N1
    integer                               :: i,Ncol,Nloc
    !
#ifdef _MPI
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG MPI_sp_init_matrix_csr: allocate sparse"
#endif
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse,"mpi_sp_init_matrix_csr")
    !
    Ncol = N
    if(present(N1))Ncol=N1
    !
    Nloc = sparse%iend-sparse%istart+1
    !
    call sp_init_matrix_csr(sparse,Nloc,Ncol)
    !
    allocate(sparse%loc(Nloc))
    do i=1,Nloc
       sparse%loc(i)%size=0
       allocate(sparse%loc(i)%dvals(0)) !empty array
       allocate(sparse%loc(i)%cvals(0)) !empty array
       allocate(sparse%loc(i)%cols(0)) !empty array
    end do
    !
#else
  STOP "mpi_sp_init_matrix_csr called in a non-mpi build"
#endif
  end subroutine mpi_sp_init_matrix_csr








  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_csr(sparse)    
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row_csr),pointer          :: row
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG sp_delete_matrix_csr: delete sparse"
#endif
    if(.not.sparse%status)return !stop "Error SPARSE/sp_delete_matrix: sparse is not allocated."
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%dvals)
       deallocate(sparse%row(i)%cvals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
  end subroutine sp_delete_matrix_csr



  subroutine mpi_sp_delete_matrix_csr(MpiComm,sparse)
    integer                              :: MpiComm
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                              :: i
    type(sparse_row_csr),pointer          :: row
    !
#ifdef _MPI
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG MPI_sp_delete_matrix_csr: delete sparse"
#endif
    if(.not.sparse%status)return !stop "Error SPARSE/mpi_sp_delete_matrix: sparse is not allocated."
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%dvals)
       deallocate(sparse%row(i)%cvals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
       !
       deallocate(sparse%loc(i)%dvals)
       deallocate(sparse%loc(i)%cvals)
       deallocate(sparse%loc(i)%cols)
       sparse%loc(i)%Size  = 0
    enddo
    deallocate(sparse%row)
    deallocate(sparse%loc)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
    !
    sparse%istart=0
    sparse%iend=0
    sparse%ishift=0
    sparse%mpi=.false.
    !
#else
    STOP "mpi_sp_delete_matrix_csr called in a non-mpi build"
#endif
  end subroutine mpi_sp_delete_matrix_csr













  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_csr_d(sparse,value,i,j)
    type(sparse_matrix_csr),intent(inout) :: sparse !
    real(8),intent(in)                    :: value  !matrix value to be inserted 
    integer,intent(in)                    :: i      !row index of the matrix value to be inserted
    integer,intent(in)                    :: j      !column index of the matrix value to be inserted
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,2I8)")"DEBUG sp_insert_element_csr_d: insert element in sparse @",i,j
#endif
    !
    column = j
    !
    row => sparse%row(i)
    !
    iadd = .false.
    if(any(row%cols == column))then
       pos = binary_search_spmat(row%cols,column)
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then     
       row%dvals(pos)=row%dvals(pos) + value  
    else                                    
       call add_to(row%dvals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element_csr_d

  subroutine sp_insert_element_csr_c(sparse,value,i,j)
    type(sparse_matrix_csr),intent(inout) :: sparse
    complex(8),intent(in)                 :: value
    integer,intent(in)                    :: i,j
    !
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,2I8)")"DEBUG sp_insert_element_csr_c: insert element in sparse @",i,j
#endif
    !
    column = j
    !
    row => sparse%row(i)
    !
    iadd = .false.                   
    if(any(row%cols == column))then   
       pos = binary_search_spmat(row%cols,column)
       iadd=.true.                          
    endif
    !
    if(iadd)then                           
       row%cvals(pos)=row%cvals(pos) + value  
    else                                    
       call add_to(row%cvals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element_csr_c


  subroutine mpi_sp_insert_element_csr_d(MpiComm,sparse,value,i,j)
    integer                               :: MpiComm !MPI global communicator
    type(sparse_matrix_csr),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
#ifdef _MPI
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,2I8)")"DEBUG MPI_sp_insert_element_csr_d: insert element in sparse @",i,j
#endif
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_insert_element_csr")
    !
    column = j
    !
    row => sparse%row(i-sparse%Ishift)
    !
    iadd = .false.                          
    if(any(row%cols == column))then         
       pos = binary_search_spmat(row%cols,column) 
       iadd=.true.                          
    endif
    !
    if(iadd)then                            
       row%dvals(pos)=row%dvals(pos) + value
    else                                    
       call add_to(row%dvals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "mpi_sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
#else
  STOP "mpi_sp_insert_element_csr_d called in a non-mpi build"
#endif
  end subroutine mpi_sp_insert_element_csr_d

  subroutine mpi_sp_insert_element_csr_c(MpiComm,sparse,value,i,j)
    integer                               :: MpiComm
    type(sparse_matrix_csr),intent(inout) :: sparse
    complex(8),intent(in)                 :: value
    integer,intent(in)                    :: i,j
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
#ifdef _MPI
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,2I8)")"DEBUG MPI_sp_insert_element_csr_c: insert element in sparse @",i,j
#endif
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_insert_element_csr")
    !
    column = j
    !
    if (ED_MODE == "normal") then
      row => sparse%row(i-sparse%Ishift)
    else
      if(column>=sparse%Istart.AND.column<=sparse%Iend)then
         row => sparse%loc(i-sparse%Ishift)
      else
         row => sparse%row(i-sparse%Ishift)
      endif
    endif
    !
    iadd = .false.                          
    if(any(row%cols == column))then         
       pos = binary_search_spmat(row%cols,column) 
       iadd=.true.                          
    endif
    !
    if(iadd)then                            
       row%cvals(pos)=row%cvals(pos) + value
    else                                    
       call add_to(row%cvals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "mpi_sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
#else
  STOP "mpi_sp_insert_element_csr_c called in a non-mpi build"
#endif
  end subroutine mpi_sp_insert_element_csr_c







  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_csr_d(sparse,matrix)
    type(sparse_matrix_csr),intent(in)   :: sparse !sparse matrix to dump
    real(8),dimension(:,:),intent(inout) :: matrix !dense matrix having the same type of the sparse one
    integer                              :: i,j,Ndim1,Ndim2
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG sp_dump_matrix_csr_d: dump sparse"
#endif
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%dvals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix_csr_d

  subroutine sp_dump_matrix_csr_c(sparse,matrix)
    type(sparse_matrix_csr),intent(in)      :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    integer                                 :: i,j,Ndim1,Ndim2
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG sp_dump_matrix_csr_c: dump sparse"
#endif
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    matrix=0.d0
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%cvals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix_csr_c




  !+------------------------------------------------------------------+
  !PURPOSE: copy a real sparse matrix without aliasing pointer rows
  !+------------------------------------------------------------------+
  subroutine sp_copy_matrix_csr_d(A,B)
    type(sparse_matrix_csr),intent(in)    :: A
    type(sparse_matrix_csr),intent(inout) :: B
    integer                               :: i,nn
    !
    if(.not.A%status)stop "sp_copy_matrix_csr_d ERROR: A is not allocated"
    if(B%status)call sp_delete_matrix(B)
    call sp_init_matrix(B,A%Nrow,A%Ncol)
    do i=1,A%Nrow
       nn = A%row(i)%size
       if(nn==0)cycle
       B%row(i)%size = nn
       deallocate(B%row(i)%dvals,B%row(i)%cvals,B%row(i)%cols)
       allocate(B%row(i)%dvals(nn))
       allocate(B%row(i)%cvals(0))
       allocate(B%row(i)%cols(nn))
       B%row(i)%dvals = A%row(i)%dvals
       B%row(i)%cols  = A%row(i)%cols
    enddo
  end subroutine sp_copy_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: real sparse matrix linear combination C = alpha*A + beta*B
  !+------------------------------------------------------------------+
  subroutine sp_axpy_matrix_csr_d(alpha,A,beta,B,C,tol)
    real(8),intent(in)                    :: alpha,beta
    type(sparse_matrix_csr),intent(in)    :: A,B
    type(sparse_matrix_csr),intent(inout) :: C
    real(8),intent(in),optional           :: tol
    real(8)                               :: threshold,value
    real(8),dimension(:),allocatable      :: work
    integer,dimension(:),allocatable      :: used_cols
    logical,dimension(:),allocatable      :: used
    integer                               :: i,j,col,nused
    !
    if(.not.A%status)stop "sp_axpy_matrix_csr_d ERROR: A is not allocated"
    if(.not.B%status)stop "sp_axpy_matrix_csr_d ERROR: B is not allocated"
    if(A%Nrow/=B%Nrow)stop "sp_axpy_matrix_csr_d ERROR: Nrow mismatch"
    if(A%Ncol/=B%Ncol)stop "sp_axpy_matrix_csr_d ERROR: Ncol mismatch"
    threshold = 0d0;if(present(tol))threshold=tol
    if(C%status)call sp_delete_matrix(C)
    call sp_init_matrix(C,A%Nrow,A%Ncol)
    allocate(work(A%Ncol),used_cols(A%Ncol),used(A%Ncol))
    work=0d0
    used=.false.
    do i=1,A%Nrow
       nused=0
       do j=1,A%row(i)%size
          col = A%row(i)%cols(j)
          if(.not.used(col))then
             nused=nused+1
             used_cols(nused)=col
             used(col)=.true.
          endif
          work(col) = work(col) + alpha*A%row(i)%dvals(j)
       enddo
       do j=1,B%row(i)%size
          col = B%row(i)%cols(j)
          if(.not.used(col))then
             nused=nused+1
             used_cols(nused)=col
             used(col)=.true.
          endif
          work(col) = work(col) + beta*B%row(i)%dvals(j)
       enddo
       call sp_set_row_csr_d(C,i,work,used_cols,nused,threshold)
       do j=1,nused
          col=used_cols(j)
          work(col)=0d0
          used(col)=.false.
       enddo
    enddo
    deallocate(work,used_cols,used)
  end subroutine sp_axpy_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: real sparse matrix in-place update Y = alpha*A + Y
  !+------------------------------------------------------------------+
  subroutine sp_axpy_matrix_inplace_csr_d(alpha,A,Y,tol)
    real(8),intent(in)                    :: alpha
    type(sparse_matrix_csr),intent(in)    :: A
    type(sparse_matrix_csr),intent(inout) :: Y
    real(8),intent(in),optional           :: tol
    real(8)                               :: threshold
    real(8),dimension(:),allocatable      :: work
    integer,dimension(:),allocatable      :: used_cols
    logical,dimension(:),allocatable      :: used
    integer                               :: i,j,col,nused
    !
    if(.not.A%status)stop "sp_axpy_matrix_inplace_csr_d ERROR: A is not allocated"
    if(.not.Y%status)stop "sp_axpy_matrix_inplace_csr_d ERROR: Y is not allocated"
    if(A%Nrow/=Y%Nrow)stop "sp_axpy_matrix_inplace_csr_d ERROR: Nrow mismatch"
    if(A%Ncol/=Y%Ncol)stop "sp_axpy_matrix_inplace_csr_d ERROR: Ncol mismatch"
    threshold = 0d0;if(present(tol))threshold=tol
    allocate(work(Y%Ncol),used_cols(Y%Ncol),used(Y%Ncol))
    work=0d0
    used=.false.
    do i=1,Y%Nrow
       nused=0
       do j=1,Y%row(i)%size
          col = Y%row(i)%cols(j)
          if(.not.used(col))then
             nused=nused+1
             used_cols(nused)=col
             used(col)=.true.
          endif
          work(col) = work(col) + Y%row(i)%dvals(j)
       enddo
       do j=1,A%row(i)%size
          col = A%row(i)%cols(j)
          if(.not.used(col))then
             nused=nused+1
             used_cols(nused)=col
             used(col)=.true.
          endif
          work(col) = work(col) + alpha*A%row(i)%dvals(j)
       enddo
       call sp_set_row_csr_d(Y,i,work,used_cols,nused,threshold)
       do j=1,nused
          col=used_cols(j)
          work(col)=0d0
          used(col)=.false.
       enddo
    enddo
    deallocate(work,used_cols,used)
  end subroutine sp_axpy_matrix_inplace_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: real sparse matrix in-place scale A = alpha*A
  !+------------------------------------------------------------------+
  subroutine sp_scale_matrix_csr_d(alpha,A)
    real(8),intent(in)                    :: alpha
    type(sparse_matrix_csr),intent(inout) :: A
    integer                               :: i
    !
    if(.not.A%status)stop "sp_scale_matrix_csr_d ERROR: A is not allocated"
    do i=1,A%Nrow
       if(A%row(i)%size>0)A%row(i)%dvals = alpha*A%row(i)%dvals
    enddo
  end subroutine sp_scale_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: real sparse matrix product C = A*B
  !+------------------------------------------------------------------+
  subroutine sp_matmul_matrix_csr_d(A,B,C,tol)
    type(sparse_matrix_csr),intent(in)    :: A,B
    type(sparse_matrix_csr),intent(inout) :: C
    real(8),intent(in),optional           :: tol
    real(8)                               :: threshold,value
    real(8),dimension(:),allocatable      :: work
    integer,dimension(:),allocatable      :: used_cols
    logical,dimension(:),allocatable      :: used
    integer                               :: i,j,k,icol,col,nused
    !
    if(.not.A%status)stop "sp_matmul_matrix_csr_d ERROR: A is not allocated"
    if(.not.B%status)stop "sp_matmul_matrix_csr_d ERROR: B is not allocated"
    if(A%Ncol/=B%Nrow)stop "sp_matmul_matrix_csr_d ERROR: matrix dimensions mismatch"
    threshold = 0d0;if(present(tol))threshold=tol
    if(C%status)call sp_delete_matrix(C)
    call sp_init_matrix(C,A%Nrow,B%Ncol)
    allocate(work(B%Ncol),used_cols(B%Ncol),used(B%Ncol))
    work=0d0
    used=.false.
    do i=1,A%Nrow
       nused=0
       do j=1,A%row(i)%size
          k = A%row(i)%cols(j)
          do icol=1,B%row(k)%size
             col   = B%row(k)%cols(icol)
             value = A%row(i)%dvals(j)*B%row(k)%dvals(icol)
             if(.not.used(col))then
                nused=nused+1
                used_cols(nused)=col
                used(col)=.true.
             endif
             work(col) = work(col) + value
          enddo
       enddo
       call sp_set_row_csr_d(C,i,work,used_cols,nused,threshold)
       do j=1,nused
          col=used_cols(j)
          work(col)=0d0
          used(col)=.false.
       enddo
    enddo
    deallocate(work,used_cols,used)
  end subroutine sp_matmul_matrix_csr_d


  subroutine sp_set_row_csr_d(A,i,work,used_cols,nused,threshold)
    type(sparse_matrix_csr),intent(inout) :: A
    integer,intent(in)                    :: i,nused
    real(8),dimension(:),intent(in)       :: work
    integer,dimension(:),intent(in)       :: used_cols
    real(8),intent(in)                    :: threshold
    integer                               :: j,col,nn
    !
    nn=0
    do j=1,nused
       col=used_cols(j)
       if(abs(work(col))>threshold)nn=nn+1
    enddo
    if(allocated(A%row(i)%dvals))deallocate(A%row(i)%dvals)
    if(allocated(A%row(i)%cvals))deallocate(A%row(i)%cvals)
    if(allocated(A%row(i)%cols))deallocate(A%row(i)%cols)
    A%row(i)%size=nn
    allocate(A%row(i)%dvals(nn))
    allocate(A%row(i)%cvals(0))
    allocate(A%row(i)%cols(nn))
    nn=0
    do j=1,nused
       col=used_cols(j)
       if(abs(work(col))<=threshold)cycle
       nn=nn+1
       A%row(i)%dvals(nn)=work(col)
       A%row(i)%cols(nn)=col
    enddo
  end subroutine sp_set_row_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: real sparse matrix-vector product Hv = H*v
  !+------------------------------------------------------------------+
  subroutine sp_matvec_matrix_csr_d(H,v,Hv)
    type(sparse_matrix_csr),intent(in) :: H
    real(8),dimension(:),intent(in)    :: v
    real(8),dimension(:),intent(out)   :: Hv
    integer                            :: i,j
    !
    if(.not.H%status)stop "sp_matvec_matrix_csr_d ERROR: H is not allocated"
    if(size(v)/=H%Ncol)stop "sp_matvec_matrix_csr_d ERROR: input vector dimension mismatch"
    if(size(Hv)/=H%Nrow)stop "sp_matvec_matrix_csr_d ERROR: output vector dimension mismatch"
    Hv = 0d0
    do i=1,H%Nrow
       do j=1,H%row(i)%size
          Hv(i) = Hv(i) + H%row(i)%dvals(j)*v(H%row(i)%cols(j))
       enddo
    enddo
  end subroutine sp_matvec_matrix_csr_d


  subroutine sp_transpose_matvec_matrix_csr_d(A,v,Atv)
   type(sparse_matrix_csr),intent(in) :: A
   real(8),dimension(:),intent(in)    :: v
   real(8),dimension(:),intent(out)   :: Atv
   integer                            :: i,j
   !
   if(size(v)/=A%Nrow)stop "sp_transpose_matvec_matrix_csr_d error: input vector dimension mismatch"
   if(size(Atv)/=A%Ncol)stop "sp_transpose_matvec_matrix_csr_d error: output vector dimension mismatch"
   Atv=0d0
   do i=1,A%Nrow
      do j=1,A%row(i)%size
         Atv(A%row(i)%cols(j)) = Atv(A%row(i)%cols(j)) + A%row(i)%dvals(j)*v(i)
      enddo
   enddo
 end subroutine   sp_transpose_matvec_matrix_csr_d



#ifdef _MPI
  subroutine mpi_sp_dump_matrix_csr_d(MpiComm,sparse,matrix)
    integer                              :: MpiComm !MPI global communicator
    type(sparse_matrix_csr),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    real(8),dimension(:,:),allocatable   :: matrix_tmp
    integer                              :: i,impi,j,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG MPI_sp_dump_matrix_csr_d: dump sparse"
#endif
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_dump_matrix_csr")
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    N1_  = sparse%Nrow
    N2_  = sparse%Ncol
    Nrow = 0
    Ncol = 0
    call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
    !
    allocate(matrix_tmp(Ndim1,Ndim2)) ; matrix_tmp=0d0
    do i=sparse%Istart,sparse%Iend
       impi = i - sparse%Ishift
       do j=1,sparse%row(impi)%Size
          matrix_tmp(i,sparse%row(impi)%cols(j))=matrix_tmp(i,sparse%row(impi)%cols(j))+sparse%row(impi)%dvals(j)
       enddo
    enddo
    !
    call AllReduce_MPI(MpiCOmm,Matrix_tmp,Matrix)
    !
  end subroutine mpi_sp_dump_matrix_csr_d

  subroutine mpi_sp_dump_matrix_csr_c(MpiComm,sparse,matrix)
    integer                                 :: MpiComm
    type(sparse_matrix_csr),intent(in)      :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    complex(8),dimension(:,:),allocatable   :: matrix_tmp
    integer                                 :: i,impi,j,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG MPI_sp_dump_matrix_csr_c: dump sparse"
#endif
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_dump_matrix_csr")
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    N1_  = sparse%Nrow
    N2_  = sparse%Ncol
    Nrow = 0
    Ncol = 0
    call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
    !
    allocate(matrix_tmp(Ndim1,Ndim2)) ; matrix_tmp=zero
    do i=sparse%Istart,sparse%Iend
       impi = i - sparse%Ishift
       !Local part:
       do j=1,sparse%loc(impi)%Size
          matrix_tmp(i,sparse%loc(impi)%cols(j))=matrix_tmp(i,sparse%loc(impi)%cols(j))+sparse%loc(impi)%cvals(j)
       enddo
       !
       !Non-local part:
       do j=1,sparse%row(impi)%Size
          matrix_tmp(i,sparse%row(impi)%cols(j))=matrix_tmp(i,sparse%row(impi)%cols(j))+sparse%row(impi)%cvals(j)
       enddo
    enddo
    !
    call AllReduce_MPI(MpiCOmm,Matrix_tmp,Matrix)
    !
  end subroutine mpi_sp_dump_matrix_csr_c
#endif




#ifdef _MPI
  subroutine sp_set_mpi_matrix_csr(MpiComm,sparse,istart,iend,ishift)
    integer                               :: MpiComm !MPI global communicator
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                               :: istart !starting index for the MPI decomposition of the sparse matrix
    integer                               :: iend   !ending index for the MPI decomposition of the sparse matrix
    integer                               :: ishift !shift index for the MPI decomposition of the sparse matrix
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    sparse%istart = istart
    sparse%iend   = iend
    sparse%ishift = ishift
    sparse%mpi    = .true.
  end subroutine sp_set_mpi_matrix_csr

  subroutine sp_test_matrix_mpi(MpiComm,sparse,text)
    integer                              :: MpiComm
    type(sparse_matrix_csr),intent(in)    :: sparse
    character(len=*)                     :: text
    integer                              :: MpiRank
    !
    if(MpiComm==Mpi_Comm_Null)stop "sp_test_matrix_mpi ERROR: called in with MpiComm = Mpi_Comm_Null"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    if(.not.sparse%mpi)then
       print*,"Rank, Error in "//trim(text)//": mpi no set"
       stop
    endif
  end subroutine sp_test_matrix_mpi
#endif





  !##################################################################
  !##################################################################
  !              AUXILIARY COMPUTATIONAL ROUTINES
  !##################################################################
  !##################################################################
  recursive function binary_search_spmat(Ain,value) result(bsresult)
    integer,intent(in)           :: Ain(:), value
    integer                      :: bsresult, mid
    integer,dimension(size(Ain)) :: A,Order
    !
    a = ain
    call sort_array(a,Order)
    !
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search_spmat(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search_spmat(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
    !
    bsresult = Order(bsresult)
    !
  end function binary_search_spmat




  subroutine add_to_I(vec,val)
    integer,dimension(:),allocatable,intent(inout) :: vec
    integer,intent(in)                             :: val  
    integer,dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_I

  subroutine add_to_D(vec,val)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),intent(in)                             :: val  
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_D

  subroutine add_to_Z(vec,val)
    complex(8),dimension(:),allocatable,intent(inout) :: vec
    complex(8),intent(in)                             :: val  
    complex(8),dimension(:),allocatable               :: tmp
    integer                                           :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_Z





  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order)
    implicit none
    integer,dimension(:)                    :: array
    integer,dimension(size(array))          :: order
    integer,dimension(size(array))          :: backup
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort(array, order,1, size(array))
    do i=1,size(array)
       backup(i)=array(order(i))
    enddo
    array=backup
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:) :: array
      integer, dimension(:) :: order
      integer               :: left
      integer               :: right
      integer               :: i
      integer               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:) :: order
      integer               :: first, second
      integer               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer function qsort_rand( lower, upper )
      integer               :: lower, upper
      real(8)               :: r
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      integer               :: f,g
      integer               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array



end module ED_SPARSE_MATRIX



