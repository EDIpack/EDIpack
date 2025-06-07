module ED_MAIN
  !:synopsis: Routines to initialize, run and finalize the solver
  !Contains routine that initialize, run and finalize the Impurity model solver
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC,only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_PARSE_UMATRIX
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_CHI_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_RDM
  USE ED_DIAG
  USE ED_IO
  implicit none
  private

  !>INIT ED SOLVER
  interface ed_init_solver
     !
     !Initialize the Exact Diagonalization solver of `EDIpack`. This procedure reserves and allocates all the  
     !memory required by the solver, performs all the consistency check and initializes the bath instance guessing or reading from a file.      
     !It requires as an input a double precision array of rank-1 [ :f:var:`nb` ]. 
     !where :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .
     !If the simulation is to be run without a bath, no input is required.
     !
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_single_nobath
  end interface ed_init_solver


  !> ED SOLVER
  interface ed_solve
     !
     !Launch the Exact Diagonalizaton solver for the single-site and multiple-site (R-DMFT) quantum impurity problem.
     !It requires as an input a double precision array of rank-1 [ :f:var:`nb` ]. 
     !where :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .
     !If the simulation is to be run without a bath, no input is required.
     !
     ! The solution is achieved in this sequence:
     !
     !  #. setup the MPI environment, if any 
     !  #. Set the internal bath instance :f:var:`dmft_bath` copying from the user provided input :f:var:`bath` , if existing
     !  #. Get the low energy spectrum: call :f:func:`diagonalize_impurity`
     !  #. Get the impurity Green's functions: call :f:func:`buildgf_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity susceptibilities, if any: call :f:func:`buildchi_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity observables: call :f:func:`observables_impurity`
     !  #. Get the impurity local energy: call :f:func:`local_energy_impurity`
     !  #. Get the impurity reduced density matric: call :f:func:`rdm_impurity`
     !  #. Delete MPI environment and deallocate used structures :f:var:`state_list` and :f:var:`dmft_bath`
     !
     module procedure :: ed_solve_single
     module procedure :: ed_solve_single_nobath
  end interface ed_solve



  !> FINALIZE SOLVER AND CLEANUP ENVIRONMENT
  interface ed_finalize_solver
     ! 
     ! Finalize the Exact Diagonalization solver, clean up all the allocated memory. 
     !
     module procedure :: ed_finalize_solver_single
  end interface ed_finalize_solver



  public :: ed_init_solver
  public :: ed_solve
  public :: ed_finalize_solver



  logical,save :: isetup=.true. !Allow :f:func:`init_ed_structure` and :f:func:`setup_global` to be called. Set to :f:code:`.false.` by :f:func:`ed_init_solver`, reset by :f:func:`ed_finalize_solver` . Default :code:`.true.`



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  subroutine ed_init_solver_single(bath)
    real(8),dimension(:),intent(inout)        :: bath !user bath input array
    logical                                   :: check
    integer                                   :: i
    !
    !SET THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    if (ed_file_suffix == "") then
      write(LOGfile,"(A)")"INIT SOLVER"
    else
      write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    end if
    
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure() 
    !
    if(Nbath>0)then
      check = check_bath_dimension(bath)
      if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
      !
      bath = 0d0
    endif
    !
    call allocate_dmft_bath()
    call init_dmft_bath()
    call get_dmft_bath(bath)
    !
    if(isetup)then
       call setup_global
    endif
    call deallocate_dmft_bath()
    isetup=.false.
    !
    !DELETE THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_del_MpiComm()
#endif
    !
  end subroutine ed_init_solver_single
  
  subroutine ed_init_solver_single_nobath()
    real(8),dimension(1)     :: bath_dummy
    !
    bath_dummy=zero
    call ed_init_solver_single(bath_dummy)
  end subroutine ed_init_solver_single_nobath








  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,flag_gf,flag_mpi)
    real(8),dimension(:),intent(in)              :: bath  !user bath input array
    logical,optional                             :: flag_gf !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` . 
    logical,optional                             :: flag_mpi  !flag to solve the impurity problem parallely ( :code:`.true.` ) or not ( :code:`.false.` ). Default :code:`.true.` . 
    logical                                      :: flag_mpi_, flag_gf_
    logical                                      :: check,iflag
    !
    flag_mpi_=.true.;if(present(flag_mpi))flag_mpi_=flag_mpi
    flag_gf_=.true.;if(present(flag_gf))flag_gf_=flag_gf
    !
#ifdef _MPI    
    if(check_MPI().AND.flag_mpi_)call ed_set_MpiComm()
#endif
    !
    if(.not.allocated(impHloc))stop "ED_SOLVE ERROR: impHloc not allocated. Please call ed_set_Hloc first."
    !
    if(Nbath>0)then
      check   = check_bath_dimension(bath)
      if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    else
      write(LOGfile,"(A)")"NBATH=0: Solving isolated impurity (no bath)"
    endif
    !  
    if(MpiMaster.and.flag_mpi_)call save_input_file(str(ed_input_file))
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath)
    call write_dmft_bath()
    call save_dmft_bath(used=.true.)
    !
    call set_umatrix()
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()
    if(flag_gf_)then
       call buildgf_impurity()
       call buildchi_impurity()
    endif
    call observables_impurity()
    call local_energy_impurity()
    !
    if(Nbath>0)then
      call rdm_impurity()
    else
      print*,"**WARNING** RDM calculation not available without a bath, for now."
    endif
    !
    call deallocate_dmft_bath()
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
#ifdef _MPI    
    if(check_MPI().AND.flag_mpi_)call ed_del_MpiComm()
#endif    
    !
    nullify(spHtimesV_cc)
    nullify(spHtimesV_p)
    write(Logfile,"(A)")""
  end subroutine ed_solve_single

subroutine ed_solve_single_nobath(flag_gf,flag_mpi)
    real(8),dimension(1) :: bath_dummy  !user bath input array
    logical,optional     :: flag_gf     !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` . 
    logical,optional     :: flag_mpi    !flag to solve the impurity problem parallely ( :code:`.true.` ) or not ( :code:`.false.` ). 
    logical              :: flag_gf_,flag_mpi_
    flag_mpi_=.true.;if(present(flag_mpi))flag_mpi_=flag_mpi
    flag_gf_ =.true.;if(present(flag_gf))flag_gf_=flag_gf
    bath_dummy=zero
    call ed_solve_single(bath_dummy,flag_gf_,flag_mpi_)
  end subroutine ed_solve_single_nobath


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: deallocate and finalize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_finalize_solver_single()
    logical                            :: check
    integer                            :: i
    !
    !SET THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_set_MpiComm()
#endif
    !
    write(LOGfile,"(A)")"FINALIZE SOLVER "
    !
    !just in case deallocate some objects
    call deallocate_dmft_bath()
    call es_delete_espace(state_list)
    call deallocate_grids
    nullify(spHtimesV_cc)
    nullify(spHtimesV_p)
    if(Hb%status)call deallocate_Hreplica()
    !
    !Delete ED Structure & memory
    call delete_ed_structure()
    !
    !Ready to be setup again
    isetup=.true.
    !
    !DELETE THE MPI FRAMEWORK:
#ifdef _MPI
    if(check_MPI())call ed_del_MpiComm()
#endif
    !
  end subroutine ed_finalize_solver_single










end module ED_MAIN




