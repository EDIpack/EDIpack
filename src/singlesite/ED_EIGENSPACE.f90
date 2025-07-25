module ED_EIGENSPACE
  !:synopsis: Data types for the eigenspace
  !A class implementing a data structure to efficiently store the low part of the Fock space spectrum, automatically spreading/retrieving the eigenstates among/from MPI threads. 
  !
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SECTOR
  implicit none
  private

  type full_espace
     real(8),dimension(:),allocatable   :: e
     real(8),dimension(:,:),allocatable :: M
  end type full_espace


  type sparse_estate
     !
     !A single element of the linked list :f:var:`sparse_espace`. The :f:var:`sparse_estate` gather all the information relative to a single eigen-state of the Hamiltonian spectrum. It includes a logical :f:var:`itwin` and pointer :f:var:`twin` to identiy  the equivalent state in the twin sector, i.e. a degenerate state with opposite quantum numbers, without actually storing the eigenvector. 
     !
     integer                             :: sector        !Symmetry sector index
     real(8)                             :: e             !energy of the eigen-state, used to order the list
     real(8),dimension(:),allocatable    :: dvec          !double precision eigen-vector
     complex(8),dimension(:),allocatable :: cvec          !double complex eigen-vector
     logical                             :: itwin=.false. !twin sector logical label
     type(sparse_estate),pointer         :: twin=>null()  !link to twin :f:var:`sparse_estate` in the list 
     type(sparse_estate),pointer         :: next=>null()  !link to next :f:var:`sparse_estate` used to construct the ordered linked list
  end type sparse_estate

  type sparse_espace
     !
     !Ordered single linked list storing the lower part of the Hamiltonian spectrum state by state. 
     !
     integer                     :: size               !The current size of the list
     real(8)                     :: emax               !The maximum energy of the list, fixed by the condition :math:`e^{-\beta {\rm emax}}<` :f:var:`cutoff`
     real(8)                     :: emin               !The minimum energy of the list, i.e. the groundstate energy
     logical                     :: status=.false.     !Allocation status of the list
     type(sparse_estate),pointer :: root=>null()       !Root of the linked list
  end type sparse_espace


  interface es_insert_state
     !
     !Insert a :f:var:`sparse_estate` into the :f:var:`sparse_espace` using :f:var:`e`, :f:var:`vec`, :f:var:`sector` and optionally a :f:var:`itwin` label
     !
     module procedure :: es_insert_state_d
     module procedure :: es_insert_state_c
  end interface es_insert_state

  interface es_add_state
     !
     !Insert a :f:var:`sparse_estate` into the :f:var:`sparse_espace` using :f:var:`e`, :f:var:`vec`, :f:var:`sector` and optionally a :f:var:`itwin` label optionally  filling the list  respecting the energy threshold limit. 
     !
     module procedure :: es_add_state_d
     module procedure :: es_add_state_c
  end interface es_add_state


  interface es_return_dvector
     !
     !Returns the double precision vector of a given :f:var:`sparse_estate` indicated by the position :f:var:`n` in the list. If MPI execution is active the vector is returned already split in chunks assigned to each thread according to the current parallel algorithm. Should the :f:var:`sparse_estate` correspond to a twin state, the vector is reconstructed on the fly using :f:var:`twin` pointer and suitable reordering (as the basis of the two symmetry sector with opposite quantum numbers do have different ordering in general).     
     !
     module procedure :: es_return_dvector_default
#ifdef _MPI
     module procedure :: es_return_dvector_mpi
#endif
  end interface es_return_dvector



  interface es_return_cvector
     !
     !Returns the double complex vector of a given :f:var:`sparse_estate` indicated by the position :f:var:`n` in the list. If MPI execution is active the vector is returned already split in chunks assigned to each thread according to the current parallel algorithm. Should the :f:var:`sparse_estate` correspond to a twin state, the vector is reconstructed on the fly using :f:var:`twin` pointer and suitable reordering (as the basis of the two symmetry sector with opposite quantum numbers do have different ordering in general).     
     !
     module procedure :: es_return_cvector_default
#ifdef _MPI
     module procedure :: es_return_cvector_mpi
#endif
  end interface es_return_cvector


  public :: sparse_estate
  public :: sparse_espace
  !
  public :: es_init_espace      !init the espace                 !checked
  public :: es_delete_espace    !del the espace                  !checked
  public :: es_free_espace      !free the espace                 !checked
  !
  public :: es_insert_state     !insert a state                  !checked
  public :: es_add_state        !add a state w/ costraint        !checked
  public :: es_pop_state        !pop a state                     !checked
  !
  public :: es_return_sector       !get the sector of a state       !checked
  public :: es_return_energy       !get the energy of a state       !checked
  public :: es_return_evals        !get the eigenvalues list
  !subroutine
  public :: es_return_dvector      !get the vector of a state       !checked
  public :: es_return_cvector      !get the vector of a state       !checked
  !functions
  public :: es_return_dvec
  public :: es_return_cvec
  public :: es_return_gs_degeneracy!get the number of degenerate GS !checked
  !
  type(sparse_espace)                        :: state_list
  type(full_espace),dimension(:),allocatable :: espace
  public :: state_list          !The shared instance of the :f:var:`sparse_espace` used in the :code:`EDIpack` library
  public :: espace




contains        !some routine to perform simple operation on the lists



  !+------------------------------------------------------------------+
  !PURPOSE  : Setting up the Full ED eigen-Space
  !+------------------------------------------------------------------+
  subroutine setup_eigenspace
    integer :: isector,dim,jsector
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG setup_eigenspace"
#endif
    if(allocated(espace)) deallocate(espace)
    allocate(espace(1:Nsectors))
    do isector=1,Nsectors
       dim=GetDim(isector);if(dim==0)stop "setup_eigenspace: dim==0!"
       allocate(espace(isector)%e(dim))
       allocate(espace(isector)%M(dim,dim))
    enddo
  end subroutine setup_eigenspace





  !+------------------------------------------------------------------+
  !PURPOSE  : Deleting the Full ED eigen-Space (free the memory)
  !+------------------------------------------------------------------+
  subroutine delete_eigenspace
    integer :: isector
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG delete_eigenspace"
#endif
    if(allocated(espace))then
       do isector=1,size(espace)
          deallocate(espace(isector)%e)
          deallocate(espace(isector)%M)
       end do
       deallocate(espace)
    endif
  end subroutine delete_eigenspace





  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################



  !+------------------------------------------------------------------+
  !PURPOSE  : initialize the list of states
  !+------------------------------------------------------------------+
  function es_init_espace() result(space)
    type(sparse_espace) :: space
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG es_init_espace"
#endif
    allocate(space%root)
    space%status=.true.
    space%root%next => null()
    space%size=0
    space%emax=-huge(1d0)
    space%emin= huge(1d0)
  end function es_init_espace



  !+------------------------------------------------------------------+
  !PURPOSE  : destroy the list of states
  !+------------------------------------------------------------------+
  subroutine es_delete_espace(space)
    !Destroys the list of states.
    type(sparse_espace),intent(inout) :: space !the eigenspace to clear
    type(sparse_estate),pointer       :: p,c
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG es_delete_espace"
#endif
    if(.not.space%status)return
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       if(allocated(c%dvec))deallocate(c%dvec)
       if(allocated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
    end do
    deallocate(space%root)
    space%status=.false.
    p=>null()
    c=>null()
  end subroutine es_delete_espace




  !+------------------------------------------------------------------+
  !PURPOSE  : empty the list of states
  !+------------------------------------------------------------------+
  subroutine es_free_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG es_free_espace"
#endif
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next            !
       c%next=>null()
       if(allocated(c%dvec))deallocate(c%dvec)
       if(allocated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
    end do
    space%size=0
    space%emax=-huge(1.d0)
    space%emin=huge(1.d0)
    p=>null()
    c=>null()
  end subroutine es_free_espace







  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_add_state_d(espace,e,vec,sector,twin,size,verbose)
    type(sparse_espace),intent(inout) :: espace
    real(8),intent(in)                :: e        !The eigenenergy of the state to be added
    real(8),dimension(:),intent(in)   :: vec      !The eigenvector of the state to be added
    integer,intent(in)                :: sector   !The symetry sector index of the state to be added 
    integer,intent(in),optional       :: size     !The size threshold of the list [optional]
    logical,intent(in),optional       :: verbose
    logical,intent(in),optional       :: twin     !The twin state lable [optional]
    logical                           :: twin_
    logical                           :: verb_    
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG es_add_state_d"
#endif
    twin_=.false.;if(present(twin))twin_=twin
    verb_=.false.;if(present(verbose))verb_=verbose
    if(present(size))then !if present size add respecting the size constraint.
       if(espace%size<size)then
          call es_insert_state_d(espace,e,vec,sector,twin_)
       else
          if(e < es_return_energy(espace))then
             if(verb_)print*,"found a new state:"
             call es_pop_state(espace)
             call es_insert_state_d(espace,e,vec,sector,twin_)
          endif
       endif
    else                      !else add normally
       call es_insert_state_d(espace,e,vec,sector,twin_)
    endif
  end subroutine es_add_state_d

  subroutine es_add_state_c(espace,e,vec,sector,twin,size,verbose)
    type(sparse_espace),intent(inout) :: espace
    real(8),intent(in)                :: e
    complex(8),dimension(:),intent(in):: vec
    integer,intent(in)                :: sector
    integer,intent(in),optional       :: size
    logical,intent(in),optional       :: verbose
    logical,intent(in),optional       :: twin
    logical                           :: twin_
    logical                           :: verb_
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG es_add_state_c"
#endif
    twin_=.false.;if(present(twin))twin_=twin
    verb_=.false.;if(present(verbose))verb_=verbose
    if(present(size))then !if present size add respecting the size constraint.
       if(espace%size<size)then
          call es_insert_state_c(espace,e,vec,sector,twin_)
       else
          if(e < es_return_energy(espace))then
             if(verb_)print*,"found a new state:"
             call es_pop_state(espace)
             call es_insert_state_c(espace,e,vec,sector,twin_)
          endif
       endif
    else                      !else add normally
       call es_insert_state_c(espace,e,vec,sector,twin_)
    endif
  end subroutine es_add_state_c







  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_insert_state_d(space,e,vec,sector,twin)
    type(sparse_espace),intent(inout) :: space 
    real(8),intent(in)                :: e      !The eigenenergy of the state to be added
    real(8),dimension(:),intent(in)   :: vec    !The eigenvector of the state to be added
    integer,intent(in)                :: sector!The symmetry sector of the state to be added
    logical                           :: twin  !The twin label of the state to be added
    type(sparse_estate),pointer       :: p,c
    p => space%root
    c => p%next
    do
       if(.not.associated(c))exit
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)
    p%next%e = e
    if(e > space%emax)space%emax=e
    if(e < space%emin)space%emin=e
    allocate(p%next%dvec(size(vec)))
    p%next%dvec = vec
    p%next%itwin=.false.
    p%next%sector=sector
    space%size = space%size+1
    if(twin)then
       allocate(p%next%next)
       p%next%next%e = e    
       p%next%next%itwin=.true.
       p%next%next%sector=get_twin_sector(sector)
       p%next%next%twin => p%next
       p%next%twin      => p%next%next !wiggled arrow pointing to the twin :f:var:`sparse_estate`
       space%size = space%size+1
    endif
    if(.not.associated(c))then
       if(.not.twin)then
          p%next%next  => null()
       else
          p%next%next%next  => null()
       endif
    else
       if(.not.twin)then
          p%next%next  => c
       else
          p%next%next%next  => c
       endif
    end if
    if(associated(p))nullify(p)
    if(associated(c))nullify(c)
  end subroutine es_insert_state_d

  subroutine es_insert_state_c(space,e,vec,sector,twin)
    type(sparse_espace),intent(inout)  :: space
    real(8),intent(in)                 :: e
    complex(8),dimension(:),intent(in) :: vec
    integer,intent(in)                 :: sector
    logical                            :: twin
    type(sparse_estate),pointer        :: p,c
    p => space%root
    c => p%next
    do                         
       if(.not.associated(c))exit
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)              
    p%next%e = e
    if(e > space%emax)space%emax=e
    if(e < space%emin)space%emin=e
    allocate(p%next%cvec(size(vec)))
    p%next%cvec = vec
    p%next%itwin=.false.
    p%next%sector=sector
    space%size = space%size+1
    if(twin)then 
       allocate(p%next%next)
       p%next%next%e = e    
       p%next%next%itwin=.true.
       p%next%next%sector=get_twin_sector(sector)
       p%next%next%twin => p%next
       p%next%twin      => p%next%next !wiggled arrow pointing to the twin :f:var:`sparse_estate`
       space%size = space%size+1
    endif
    if(.not.associated(c))then
       if(.not.twin)then
          p%next%next  => null()
       else
          p%next%next%next  => null()
       endif
    else
       if(.not.twin)then
          p%next%next  => c 
       else
          p%next%next%next  => c
       endif
    end if
    p=>null()
    c=>null()
  end subroutine es_insert_state_c






  !+------------------------------------------------------------------+
  !PURPOSE  : remove last element from the list, if +n is given remove 
  ! the n-th element, if +e is given remove the state with state%e=e
  ! hint: CIRCLE = twin state (a state flagged with itwin=T bearing no vector)
  !       SQUARE = normal state (itwin=F bearing vector)
  !+------------------------------------------------------------------+
  subroutine es_pop_state(space,n)
    type(sparse_espace),intent(inout) :: space
    integer,optional,intent(in)       :: n
    integer                           :: i,pos
    type(sparse_estate),pointer       :: pp,p,c
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")"DEBUG es_pop_state"
#endif
    pos= space%size ; if(present(n))pos=n
    !
    if(pos>space%size)stop "es_pop_state: pos > espace.size"
    if(space%size==0)stop "es_pop_state: empty list"
    pp => null()
    p  => null()
    c  => space%root
    do i=1,pos
       pp => p
       p  => c
       c  => c%next
       if(.not.associated(c))return !empty or end of the list
    end do
    !c is a circle, so the prev/next are necessarily squares: remove c and p (twins)
    if(c%itwin)then             
       pp%next => c%next
       !delete C
       if(allocated(c%dvec))deallocate(c%dvec)
       if(allocated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
       !delete P
       if(allocated(p%dvec))deallocate(p%dvec)
       if(allocated(p%cvec))deallocate(p%cvec)
       if(associated(p%twin))p%twin=>null()
       deallocate(p)
       p => pp
       space%size=space%size-2
    else
       !c is a square:
       !if c%next is associated:
       ! if it is a circle: delete c and c%next (twins)
       ! if it is a square: delete c
       !else c%next is not associated: delete c
       if(associated(c%next))then
          if(c%next%itwin)then
             p%next => c%next%next 
             !delete C
             if(allocated(c%dvec))deallocate(c%dvec)
             if(allocated(c%cvec))deallocate(c%cvec)
             if(associated(c%twin))c%twin=>null()
             deallocate(c)
             !delete C%NEXT
             if(allocated(c%next%dvec))deallocate(c%next%dvec)
             if(allocated(c%next%cvec))deallocate(c%next%cvec)
             if(associated(c%next%twin))c%next%twin=>null()
             deallocate(c%next)
             space%size=space%size-2
          else
             p%next => c%next
             if(allocated(c%dvec))deallocate(c%dvec)
             if(allocated(c%cvec))deallocate(c%cvec)
             if(associated(c%twin))c%twin=>null()
             deallocate(c)
             space%size=space%size-1
          endif
       else
          p%next => c%next
          if(allocated(c%dvec))deallocate(c%dvec)
          if(allocated(c%cvec))deallocate(c%cvec)
          if(associated(c%twin))c%twin=>null()
          deallocate(c)
          space%size=space%size-1
       endif
    endif
    if(space%size>0)then
       !    space%root%next => null()
       !    space%size=0
       !    space%emax=-huge(1.d0)
       !    space%emin=huge(1.d0)
       ! else
       ! if(pos==space%size)then     !pop last term carrying e=emax, update emax
       !    space%emax = p%e
       ! elseif(pos==1)then          !pop first term carrying e=emin, update emin
       !    space%emin = p%e
       ! endif
       !
       space%emax = es_return_energy(space,space%size)
       space%emin = es_return_energy(space,1)
    endif
    if(associated(pp))nullify(pp)
    if(associated(p))nullify(p)
    if(associated(c))nullify(c)
  end subroutine es_pop_state






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_gs_degeneracy(space,gsthreshold) result(numzero)
    type(sparse_espace),intent(in) :: space
    real(8),optional               :: gsthreshold
    real(8)                        :: gsthreshold_
    integer                        :: numzero,pos
    type(sparse_estate),pointer    :: c
    real(8)                        :: oldzero,enemin
    gsthreshold_=1.d-9;if(present(gsthreshold))gsthreshold_=gsthreshold
    if(.not.space%status) stop "es_return_gs_degeneracy: espace not allocated"
    oldzero=1000.d0
    numzero=0
    c => space%root
    pos=0
    do 
       c => c%next
       pos=pos+1
       if(.not.associated(c))exit !end of the list
       enemin = c%e
       if (enemin < oldzero-10.d0*gsthreshold_) then
          numzero=1
          oldzero=enemin
       elseif(abs(enemin-oldzero) <= gsthreshold_)then
          numzero=numzero+1
          oldzero=min(oldzero,enemin)
       endif
    end do
    if(associated(c))nullify(c)
  end function es_return_gs_degeneracy




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_sector(space,n) result(sector)
    !
    ! Returns the symmetry sector index of the indicated :f:var:`sparse_estate` in the list. If no position is indicated it returns the energy of the last state.
    !
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n !positions of the indicated :f:var:`sparse_estate` in the list     
    integer                        :: sector !the symmetry sector index of the indicated :f:var:`sparse_estate` in the list
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_sector: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_sector: n > espace.size"
    sector=0
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(space%size==0)return
    sector = c%sector
    if(associated(c))nullify(c)
  end function es_return_sector






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_energy(space,n) result(egs)
    !
    ! Returns the eigen-energy of the indicated :f:var:`sparse_estate` in the list. If no position is indicated it returns the energy of the last state. 
    !
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n !positions of the indicated :f:var:`sparse_estate` in the list 
    real(8)                        :: egs !eigen-energy of the indicated  :f:var:`sparse_estate` in the list
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_energy: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)    stop "es_return_energy: n > espace.size"
    c => space%root
    egs=space%emax
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(space%size==0)return
    if(.not.c%itwin)then
       egs = c%e
    else
       egs = c%twin%e
    endif
    if(associated(c))nullify(c)
  end function es_return_energy





  subroutine es_return_evals(space,evals)
    !
    ! Returns the list of eigen-energies of the indicated :f:var:`sparse_estate` in the list. 
    !
    type(sparse_espace),intent(in)   :: space
    real(8),dimension(:),allocatable :: evals
    integer                          :: i
    if(allocated(evals))deallocate(evals)
    if(.not.space%status) stop "es_return_evals: space%stauts=F, not allocated"
    if(space%size==0)stop "es_return_evals: space%size=0"
    allocate(evals(space%size))
    do i=1,space%size
       evals(i) = es_return_energy(space,i)
    enddo
  end subroutine es_return_evals
  


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine es_return_dvector_default(space,n,vector)
    type(sparse_espace),intent(in)   :: space
    integer,optional,intent(in)      :: n      !The position in the list corresponding to the vector to be returned
    real(8),dimension(:),allocatable :: vector !The selected eigen-vector 
    type(sparse_estate),pointer      :: c
    integer                          :: i,pos
    integer                          :: dim
    integer,dimension(:),allocatable :: order
    !
    if(.not.space%status) stop "es_return_dvector ERRROR: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_dvector ERRROR: n > espace.size"
    if(space%size==0)stop "es_return_cvector ERRROR: espace emtpy"
    !
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    !
    Dim = getdim(c%sector)
    allocate(vector(Dim)) ; vector = 0d0
    if(.not.c%itwin)then
       vector = c%dvec
    else
       allocate(Order(dim))
       call twin_sector_order(c%twin%sector,Order)
       do i=1,dim
          vector(i) = c%twin%dvec(Order(i))
       enddo
       deallocate(order)
    endif
    if(associated(c))nullify(c)
  end subroutine es_return_dvector_default

  subroutine es_return_cvector_default(space,n,vector)
    type(sparse_espace),intent(in)      :: space
    integer,optional,intent(in)         :: n
    complex(8),dimension(:),allocatable :: vector
    type(sparse_estate),pointer         :: c
    integer                             :: i,pos
    integer                             :: dim
    integer,dimension(:),allocatable    :: order
    !
    if(.not.space%status) stop "es_return_cvector ERRROR: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_cvector ERRROR: n > espace.size"
    if(space%size==0)stop "es_return_cvector ERRROR: espace emtpy"
    !
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    !
    Dim = getdim(c%sector)
    allocate(vector(Dim)) ; vector = zero
    if(.not.c%itwin)then
       vector = c%cvec
    else
       allocate(Order(dim))
       call twin_sector_order(c%twin%sector,Order)
       do i=1,dim
          vector(i) = c%twin%cvec(Order(i))
       enddo
       deallocate(order)
    endif
    if(associated(c))nullify(c)
  end subroutine  es_return_cvector_default

#ifdef _MPI
  subroutine es_return_dvector_mpi(MpiComm,space,n,vector)
    integer                          :: MpiComm
    type(sparse_espace),intent(in)   :: space
    integer,optional,intent(in)      :: n
    real(8),dimension(:),allocatable :: vtmp
    real(8),dimension(:),allocatable :: vector
    type(sparse_estate),pointer      :: c
    integer                          :: i,pos,Nloc,Ndim
    integer                          :: dim,ierr
    logical                          :: MpiMaster
    integer,dimension(:),allocatable :: order
    !
    if(MpiComm==MPI_COMM_NULL)return
    if(MpiComm==MPI_UNDEFINED)stop "es_return_dvector_MPI ERRROR: MpiComm = MPI_UNDEFINED"
    !
    if(.not.space%status) stop "es_return_dvector_MPI ERRROR: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_dvector_MPI ERRROR: n > espace.size"
    if(space%size==0)stop "es_return_dvector_MPI ERRROR: espace emtpy"
    !
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    !
    if(.not.c%itwin)then
       Nloc = size(c%dvec)
    else
       Nloc = size(c%twin%dvec)
    endif
    !
    !Ensure that the sum of the dimension of all vector chunks equals the sector dimension.
    Dim  = getdim(c%sector)
    Ndim = 0
    call Allreduce_MPI(MpiComm,Nloc,Ndim)
    if(Dim/=Ndim)stop "es_return_dvector_MPI ERROR: Dim != Ndim from v chunks"
    !
    MpiMaster = get_master_MPI(MpiComm)
    !
    if(.not.c%itwin)then
       if(MpiMaster)then
          allocate(Vector(Ndim))
       else
          allocate(Vector(1))
       endif
       Vector = 0d0
       call gather_vector_MPI(MpiComm,c%dvec,Vector)
    else
       !
       if(MpiMaster)then
          allocate(Vtmp(Ndim))
          allocate(Order(Dim))
          call twin_sector_order(c%twin%sector,Order)
       else
          allocate(Vtmp(1))
       endif
       Vtmp = 0d0
       call gather_vector_MPI(MpiComm,c%twin%dvec,Vtmp)
       if(MpiMaster)then
          allocate(Vector(Ndim))
          forall(i=1:Dim)Vector(i) = Vtmp(Order(i))
          deallocate(Order)
       else
          allocate(Vector(1))
          Vector = 0d0
       endif
       deallocate(Vtmp)
    endif
    if(associated(c))nullify(c)
  end subroutine  es_return_dvector_mpi

  subroutine es_return_cvector_mpi(MpiComm,space,n,vector)
    integer                             :: MpiComm
    type(sparse_espace),intent(in)      :: space
    integer,optional,intent(in)         :: n
    complex(8),dimension(:),allocatable :: vtmp
    complex(8),dimension(:),allocatable :: vector
    type(sparse_estate),pointer         :: c
    integer                             :: i,pos,Nloc,Ndim
    integer                             :: dim,ierr
    logical                             :: MpiMaster
    integer,dimension(:),allocatable    :: order
    !
    if(MpiComm==MPI_UNDEFINED)stop "es_return_cvector ERRROR: MpiComm = MPI_UNDEFINED"
    !
    if(.not.space%status) stop "es_return_cvector ERRROR: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_cvector ERRROR: n > espace.size"
    if(space%size==0)stop "es_return_cvector ERRROR: espace emtpy"
    !
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    !
    if(.not.c%itwin)then
       Nloc = size(c%cvec)
    else
       Nloc = size(c%twin%cvec)
    endif
    !
    !Ensure that the sum of the dimension of all vector chunks equals the sector dimension.
    Dim  = getdim(c%sector)
    Ndim = 0
    call Allreduce_MPI(MpiComm,Nloc,Ndim)
    if(Dim/=Ndim)stop "es_return_cvector_MPI ERROR: Dim != Ndim from v chunks"
    !
    MpiMaster = get_master_MPI(MpiComm)
    !
    if(.not.c%itwin)then
       if(MpiMaster)then
          allocate(Vector(Ndim))
       else
          allocate(Vector(1))
       endif
       Vector = zero
       call gather_vector_MPI(MpiComm,c%cvec,Vector)
    else
       !
       if(MpiMaster)then
          allocate(Vtmp(Ndim))
          allocate(Order(Dim))
          call twin_sector_order(c%twin%sector,Order)
       else
          allocate(Vtmp(1))
       endif
       Vtmp = zero
       call gather_vector_MPI(MpiComm,c%twin%cvec,Vtmp)
       if(MpiMaster)then
          allocate(Vector(Ndim))
          forall(i=1:Dim)Vector(i) = Vtmp(Order(i))
          deallocate(Order)
       else
          allocate(Vector(1))
          Vector = zero
       endif
       deallocate(Vtmp)
    endif
    if(associated(c))nullify(c)
  end subroutine   es_return_cvector_mpi
#endif




  function es_return_dvec(space,n) result(vector)
    type(sparse_espace),intent(in)   :: space
    integer,optional,intent(in)      :: n
    real(8),dimension(:),allocatable :: vector
    if(allocated(vector))deallocate(vector)
#ifdef _MPI    
    if(MpiStatus)then
       call es_return_dvector_mpi(MpiComm,space,n,vector)
    else
       call es_return_dvector_default(space,n,vector)
    endif
#else
    call es_return_dvector_default(space,n,vector)
#endif
  end function es_return_dvec


  function es_return_cvec(space,n) result(vector)
    type(sparse_espace),intent(in)      :: space
    integer,optional,intent(in)         :: n
    complex(8),dimension(:),allocatable :: vector
    if(allocated(vector))deallocate(vector)
#ifdef _MPI    
    if(MpiStatus)then
       call es_return_cvector_mpi(MpiComm,space,n,vector)
    else
       call es_return_cvector_default(space,n,vector)
    endif
#else
    call es_return_cvector_default(space,n,vector)
#endif
  end function es_return_cvec







end module ED_EIGENSPACE
