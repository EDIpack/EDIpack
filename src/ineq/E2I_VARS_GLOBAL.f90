MODULE E2I_VARS_GLOBAL
  !:synopsis: Global variables, inequivalent sites version
  !Contains all variables, arrays and derived types instances shared throughout the code.
  !Specifically, it contains definitions of the :f:var:`effective_bath`, the :f:var:`gfmatrix` and the :f:var:`sector` data structures. 
  !
  USE EDIPACK
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:free_unit,reg,str
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  integer :: Nnambu


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                   :: wm,tau,wr,vm,vr





  !--------------- LATTICE WRAP VARIABLES -----------------!
  complex(8),dimension(:,:,:,:,:),allocatable        :: Hloc_ineq
  complex(8),dimension(:,:,:,:,:),allocatable,save   :: single_particle_density_matrix_ineq
  complex(8),dimension(:,:,:),allocatable,save       :: impurity_density_matrix_ineq
  real(8),dimension(:,:),allocatable,save            :: dens_ineq 
  real(8),dimension(:,:),allocatable,save            :: docc_ineq
  real(8),dimension(:,:,:),allocatable,save          :: mag_ineq
  real(8),dimension(:,:,:),allocatable,save          :: phisc_ineq
  real(8),dimension(:,:,:),allocatable,save          :: argsc_ineq
  real(8),dimension(:,:,:,:),allocatable             :: exct_ineq ![1:4,Norb,Norb]
  real(8),dimension(:,:),allocatable,save            :: dd_ineq,e_ineq
  integer,allocatable,dimension(:,:)                 :: neigen_sector_ineq
  integer,allocatable,dimension(:)                   :: neigen_total_ineq
  real(8),dimension(:,:,:),allocatable               :: Hreplica_lambda_ineq ![Nineq,Nbath,Nsym]
  real(8),dimension(:,:,:),allocatable               :: Hgeneral_lambda_ineq ![Nineq,Nbath,Nsym]



  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ed_file_suffix=""       !suffix string attached to the output files.
  character(len=10)                                  :: ineq_site_suffix="_ineq"
  integer                                            :: site_indx_padding=4


  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                            :: MpiComm_Global=MPI_COMM_NULL
  integer                                            :: MpiComm=MPI_COMM_NULL
  integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                            :: MpiGroup=MPI_GROUP_NULL
#else
  integer                                            :: MpiComm_Global=0
  integer                                            :: MpiComm=0
  integer                                            :: MpiGroup_Global=0
  integer                                            :: MpiGroup=0
#endif
  logical                                            :: MpiStatus=.false.
  logical                                            :: MpiMaster=.true.
  integer                                            :: MpiRank=0
  integer                                            :: MpiSize=1
  logical                                            :: mpiAllThreads=.true.
  !

contains






  !=========================================================
  subroutine ed_set_MpiComm()
#ifdef _MPI
    integer :: ierr
    ! call MPI_Comm_dup(Comm,MpiComm_Global,ierr)
    ! call MPI_Comm_dup(Comm,MpiComm,ierr)
    MpiComm_Global = MPI_COMM_WORLD
    MpiComm        = MPI_COMM_WORLD
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_MpiComm: setting MPI comm"
#endif
#endif
  end subroutine ed_set_MpiComm

  subroutine ed_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_del_MpiComm: deleting MPI comm"
#endif
#endif
  end subroutine ed_del_MpiComm
  !=========================================================




END MODULE E2I_VARS_GLOBAL
