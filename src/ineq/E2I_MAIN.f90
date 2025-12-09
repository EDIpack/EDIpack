module E2I_MAIN
  !:synopsis: Main solver routines, inequivalent sites version
  !Contains routine that initialize, run and finalize the Impurity model solver
  USE EDIPACK
  USE ED_INPUT_VARS
  !
  USE E2I_VARS_GLOBAL
  USE E2I_AUX_FUNX
  USE E2I_BATH
  !
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC,only: assert_shape
  implicit none
  private

  !>INIT ED SOLVER
  interface ed_init_solver
     !
     !Initialize the Exact Diagonalization solver of `EDIpack`. This procedure reserves and allocates all the  
     !memory required by the solver, performs all the consistency check and initializes the bath instance guessing or reading from a file.      
     !It requires as an input a double precision array of rank-2 [ :f:var:`nb` , :f:var:`nlat` ] for the Real space DMFT case. :f:var:`nlat` is the number of inequivalent impurity sites,
     !while :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` .  If the simulation is to be run without a bath, an integer corresponding
     !to :f:var:`nlat` has to be passed instead.
     !
     module procedure :: ed_init_solver_lattice
     module procedure :: ed_init_solver_lattice_nobath
  end interface ed_init_solver


  !> ED SOLVER
  interface ed_solve
     !
     !Launch the Exact Diagonalizaton solver for the single-site and multiple-site (R-DMFT) quantum impurity problem.
     !It requires as an input a double precision array of  rank-2 [ :f:var:`nb` , :f:var:`nlat` ] for the Real space DMFT case. :f:var:`nlat` is the number of inequivalent impurity sites,
     !while :f:var:`nb` depends on the bath size and geometry and can be obtained from :f:func:`get_bath_dimension` . If the simulation is to be run without a bath, an integer corresponding
     !to :f:var:`nlat` has to be passed instead.
     !
     ! The solution is achieved in this sequence:
     !
     !  #. setup the MPI environment, if any 
     !  #. Set the internal bath instance :f:var:`dmft_bath` copying from the user provided input :f:var:`bath`, if existing
     !  #. Get the low energy spectrum: call :f:func:`diagonalize_impurity`
     !  #. Get the impurity Green's functions: call :f:func:`buildgf_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity susceptibilities, if any: call :f:func:`buildchi_impurity` (if :f:var:`sflag` = :code:`.true.` )
     !  #. Get the impurity observables: call :f:func:`observables_impurity`
     !  #. Get the impurity local energy: call :f:func:`local_energy_impurity`
     !  #. Get the impurity reduced density matric: call :f:func:`rdm_impurity`
     !  #. Delete MPI environment and deallocate used structures :f:var:`state_list` and :f:var:`dmft_bath`
     !
     module procedure :: ed_solve_lattice
     module procedure :: ed_solve_lattice_nobath
  end interface ed_solve



  !> FINALIZE SOLVER AND CLEANUP ENVIRONMENT
  interface ed_finalize_solver
     ! 
     ! Finalize the Exact Diagonalization solver, clean up all the allocated memory. 
     !
     module procedure :: ed_finalize_solver_lattice
  end interface ed_finalize_solver



  public :: ed_init_solver
  public :: ed_solve
  public :: ed_finalize_solver



  logical,save :: isetup=.true. !Allow :f:func:`init_ed_structure` and :f:func:`setup_global` to be called. Set to :f:code:`.false.` by :f:func:`ed_init_solver`, reset by :f:func:`ed_finalize_solver` . Default :code:`.true.`



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  subroutine ed_init_solver_lattice(bath)
    real(8),dimension(:,:),intent(inout) :: bath !user bath input array
    integer                              :: ilat,Nineq,Nsectors
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(phisc_ineq))deallocate(phisc_ineq)
    if(allocated(argsc_ineq))deallocate(argsc_ineq)
    if(allocated(exct_ineq))deallocate(exct_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(single_particle_density_matrix_ineq))deallocate(single_particle_density_matrix_ineq)
    if(allocated(impurity_density_matrix_ineq))deallocate(impurity_density_matrix_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    Nineq = size(bath,1)
    if(Nbath>0)then
       select case(bath_type)
       case default
       case('replica','general')
          if(.not.allocated(Hlambda_ineq))&
               stop "ERROR ed_init_solver: initial lambda not defined for all sites"
       end select
    endif
    !
    allocate(dens_ineq(Nineq,Norb))
    allocate(docc_ineq(Nineq,Norb))
    allocate(mag_ineq(Nineq,3,Norb))
    allocate(phisc_ineq(Nineq,Norb,Norb))
    allocate(argsc_ineq(Nineq,Norb,Norb))
    allocate(exct_ineq(Nineq,4,Norb,Norb))
    allocate(e_ineq(Nineq,4))
    allocate(dd_ineq(Nineq,4))
    allocate(single_particle_density_matrix_ineq(Nineq,Nspin,Nspin,Norb,Norb))
    allocate(impurity_density_matrix_ineq(Nineq,4**Norb,4**Norb))
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       if(Nbath>0)then
          select case(bath_type)
          case default;
          case ('replica','general');call Hreplica_site(ilat)
          end select
       endif
       !set the ilat-th lambda vector basis for the replica bath
       call ed_init_solver(bath(ilat,:))
    enddo
#ifdef _MPI

    if(check_MPI())call Barrier_MPI()
#endif
    call ed_reset_suffix
    !
    !This follows because Nsectors is defined after ED is initialized
    Nsectors = ed_get_nsectors()
    allocate(neigen_sector_ineq(Nineq,Nsectors))
    allocate(neigen_total_ineq(Nineq))
    do ilat=1,Nineq       
       call ed_get_neigen_sector( neigen_sector_ineq(ilat,:) )
       neigen_total_ineq(ilat)  = lanc_nstates_total
    end do
    !
  end subroutine ed_init_solver_lattice


  subroutine ed_init_solver_lattice_nobath(Nlat)
    integer                :: Nlat
    real(8),dimension(Nlat,1) :: bath_dummy !user bath input array
    !
    bath_dummy=zero
    call ed_init_solver_lattice(bath_dummy)
  end subroutine ed_init_solver_lattice_nobath




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_lattice(bath,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii,flag_gf)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8)                                       :: bath(:,:) !user bath input array
    logical,optional                              :: mpi_lanc  !parallelization strategy flag: if :code:`.false.` each core serially solves an inequivalent site, if :code:`.true.` all cores parallely solve each site in sequence. Default :code:`.false.` .
    real(8),optional,dimension(size(bath,1),Norb) :: Uloc_ii !site-dependent values for :f:var:`uloc` , overriding the ones in the input file. It has dimension [ :f:var:`nlat` , :f:var:`norb` ].
    real(8),optional,dimension(size(bath,1))      :: Ust_ii  !site-dependent values for :f:var:`ust` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(size(bath,1))      :: Jh_ii   !site-dependent values for :f:var:`jh` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(size(bath,1))      :: Jp_ii   !site-dependent values for :f:var:`jp` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(size(bath,1))      :: Jx_ii   !site-dependent values for :f:var:`jx` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    logical,optional                              :: flag_gf   !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` . 
    !
    !MPI  auxiliary vars
    real(8)                                       :: dens_tmp(size(bath,1),Norb)
    real(8)                                       :: docc_tmp(size(bath,1),Norb)
    real(8)                                       :: mag_tmp(size(bath,1),3,Norb)
    real(8)                                       :: phisc_tmp(size(bath,1),Norb,Norb)
    real(8)                                       :: argsc_tmp(size(bath,1),Norb,Norb)
    real(8)                                       :: exct_tmp(size(bath,1),4,Norb,Norb)
    real(8)                                       :: e_tmp(size(bath,1),4)
    real(8)                                       :: dd_tmp(size(bath,1),4)
    !    
    complex(8)                                    :: single_particle_density_matrix_tmp(size(bath,1),Nspin,Nspin,Norb,Norb)
    complex(8)                                    :: impurity_density_matrix_tmp(size(bath,1),4**Norb,4**Norb)
    !
    integer,allocatable,dimension(:,:)            :: neigen_sector_tmp !(size(bath,1),Nsectors)
    integer                                       :: neigen_total_tmp(size(bath,1))
    ! 
    integer                                       :: i,j,ilat,iorb,jorb,ispin,jspin,Nsectors
    integer                                       :: Nineq
    logical                                       :: check_dim,mpi_lanc_,flag_gf_
    character(len=5)                              :: tmp_suffix
    !
    integer                                       :: MPI_ID=0
    integer                                       :: MPI_SIZE=1
    logical                                       :: MPI_MASTER=.true.
    !
    integer                                       :: mpi_err 
    !
#ifdef _MPI    
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI()
       MPI_SIZE   = get_Size_MPI()
       MPI_MASTER = get_Master_MPI()
    endif
#endif
    !
    mpi_lanc_=.false.;if(present(mpi_lanc))mpi_lanc_=mpi_lanc
    !
    flag_gf_=.true.;if(present(flag_gf))flag_gf=flag_gf_
    !
    ! Check dimensions
    Nineq=size(bath,1)
    Nsectors = ed_get_nsectors()
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nineq"
    !
    if(.not.allocated(Hloc_ineq))stop "ed_solve_lattice error: Hloc_ineq not allocated. Call ed_set_Hloc first."
    !Check the dimensions of the bath are ok.
    !This can always be done in parallel no issues with mpi_lanc
    ! do ilat=1+MPI_ID,Nineq,MPI_SIZE
    !    check_dim = check_bath_dimension(bath(ilat,:))
    !    if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    ! end do
    !
    allocate(neigen_sector_tmp(Nineq,Nsectors))
    !
    dens_ineq     = 0d0  ; docc_ineq     = 0d0
    mag_ineq      = 0d0  ; phisc_ineq    = 0d0  ; argsc_ineq = 0d0 
    e_ineq        = 0d0  ; dd_ineq       = 0d0  ; exct_ineq  = 0d0
    single_particle_density_matrix_ineq = zero
    impurity_density_matrix_ineq = zero
    !
    dens_tmp   = 0d0  ; docc_tmp   = 0d0
    mag_tmp    = 0d0  ; phisc_tmp  = 0d0 ; argsc_tmp = 0d0
    e_tmp      = 0d0  ; dd_tmp     = 0d0 ; exct_tmp  = 0d0
    neigen_sector_tmp = 0
    neigen_total_tmp  = 0
    single_particle_density_matrix_tmp = zero
    impurity_density_matrix_tmp = zero
    !
    !
    select case(mpi_lanc_)
    case default              !mpi_lanc=False => solve sites with MPI
       if(MPI_MASTER)call start_timer(unit=LOGfile)
       do ilat = 1 + MPI_ID, Nineq, MPI_SIZE
          write(LOGfile,*)"CPU: "//str(MPI_ID)//" SOLVES INEQ SITE: "//str(ilat,Npad=4)
          !
          call ed_set_suffix(ilat)
          !
          !If required set the local value of U per each site
          if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
          if(present(Ust_ii)) Ust = Ust_ii(ilat)
          if(present(Jh_ii))  Jh  = Jh_ii(ilat)
          if(present(Jp_ii))  Jp  = Jp_ii(ilat)
          if(present(Jx_ii))  Jx  = Jx_ii(ilat)
          !
          !Solve the impurity problem for the ilat-th site, this are set at init time
          lanc_nstates_total = neigen_total_ineq(ilat)
          call ed_set_neigen_sector( neigen_sector_ineq(ilat,:) )
          !
          if(allocated(Hloc_anomalous_ineq))then
            call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:),Hloc_anomalous_ineq(ilat,:,:,:,:))
          else
            call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
          endif
          !
          if(MPI_MASTER)call save_input_file(str(ed_input_file))
          !
          call ed_solve(bath(ilat,:),flag_gf=flag_gf_,flag_mpi=mpi_lanc_)
          !
          neigen_total_tmp(ilat)      = lanc_nstates_total
          call ed_get_neigen_sector ( neigen_sector_tmp(ilat,:) )
          call ed_get_dens( dens_tmp(ilat,1:Norb) )
          call ed_get_docc( docc_tmp(ilat,1:Norb) ) 
          call ed_get_mag( mag_tmp(ilat,:,1:Norb) )
          call ed_get_phi( phisc_tmp(ilat,1:Norb,1:Norb) )
          call ed_get_argphi( argsc_tmp(ilat,1:Norb,1:Norb) )
          call ed_get_exct( exct_tmp(ilat,1:4,1:Norb,1:Norb) )
          call ed_get_eimp( e_tmp(ilat,:) )
          call ed_get_doubles( dd_tmp(ilat,:) )
          call ed_get_sp_dm( single_particle_density_matrix_tmp(ilat,:,:,:,:) )
          if(RDM_FLAG)call ed_get_impurity_rdm( impurity_density_matrix_tmp(ilat,:,:) )
       enddo
#ifdef _MPI
       call Barrier_MPI()
#endif
       if(MPI_MASTER)call stop_timer
       call ed_reset_suffix
       !
#ifdef _MPI
       if(check_MPI())then
          call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,docc_tmp,docc_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,mag_tmp,mag_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,phisc_tmp,phisc_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,argsc_tmp,argsc_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,exct_tmp,exct_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,e_tmp,e_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,dd_tmp,dd_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,impurity_density_matrix_tmp,impurity_density_matrix_ineq)
          neigen_sector_ineq=0
          neigen_total_ineq=0
          call AllReduce_MPI(MPI_COMM_WORLD,neigen_sector_tmp,neigen_sector_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,neigen_total_tmp,neigen_total_ineq)
          call Barrier_MPI()
       else
          dens_ineq               = dens_tmp
          docc_ineq               = docc_tmp
          mag_ineq                = mag_tmp
          phisc_ineq              = phisc_tmp
          argsc_ineq              = argsc_tmp
          exct_ineq               = exct_tmp
          e_ineq                  = e_tmp
          dd_ineq                 = dd_tmp
          neigen_sector_ineq      = neigen_sector_tmp
          neigen_total_ineq       = neigen_total_tmp
          impurity_density_matrix_ineq = impurity_density_matrix_tmp
       endif
#else
       dens_ineq               = dens_tmp
       docc_ineq               = docc_tmp
       mag_ineq                = mag_tmp
       phisc_ineq              = phisc_tmp
       argsc_ineq              = argsc_tmp
       exct_ineq               = exct_tmp
       e_ineq                  = e_tmp
       dd_ineq                 = dd_tmp
       neigen_sector_ineq      = neigen_sector_tmp
       neigen_total_ineq       = neigen_total_tmp
       impurity_density_matrix_ineq = impurity_density_matrix_tmp
#endif
       !
       !
       !
    case(.true.)                !solve sites serial, Lanczos with MPI
       if(MPI_MASTER)call start_timer(unit=LOGfile)
       do ilat = 1, Nineq
          write(LOGfile,*)" SOLVES INEQ SITE w/ MPI LANCZOS: "//str(ilat,Npad=4)
          !
          call ed_set_suffix(ilat)
          !          
          !If required set the local value of U per each site
          if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
          if(present(Ust_ii)) Ust = Ust_ii(ilat)
          if(present(Jh_ii))  Jh  = Jh_ii(ilat)
          if(present(Jp_ii))  Jp  = Jp_ii(ilat)
          if(present(Jx_ii))  Jx  = Jx_ii(ilat)
          !
          !Solve the impurity problem for the ilat-th site
          lanc_nstates_total = neigen_total_ineq(ilat)
          call ed_set_neigen_sector( neigen_sector_ineq(ilat,:) )
          !
          if(allocated(Hloc_anomalous_ineq))then
            call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:),Hloc_anomalous_ineq(ilat,:,:,:,:))
          else
            call ed_set_Hloc(Hloc_ineq(ilat,:,:,:,:))
          endif
          call ed_solve(bath(ilat,:),flag_gf=flag_gf_,flag_mpi=mpi_lanc_)
          !
          neigen_total_ineq(ilat)     = lanc_nstates_total
          call ed_get_neigen_sector ( neigen_sector_ineq(ilat,:) )
          call ed_get_dens( dens_ineq(ilat,1:Norb) )
          call ed_get_docc( docc_ineq(ilat,1:Norb) ) 
          call ed_get_mag( mag_ineq(ilat,:,1:Norb) )
          call ed_get_phi( phisc_ineq(ilat,1:Norb,1:Norb) )
          call ed_get_argphi( argsc_ineq(ilat,1:Norb,1:Norb) )
          call ed_get_exct( exct_ineq(ilat,1:4,1:Norb,1:Norb) )
          call ed_get_eimp( e_ineq(ilat,:) )
          call ed_get_doubles( dd_ineq(ilat,:) )
          call ed_get_sp_dm( single_particle_density_matrix_ineq(ilat,:,:,:,:) )
          if(RDM_FLAG)call ed_get_impurity_rdm( impurity_density_matrix_ineq(ilat,:,:) )
       enddo
       if(MPI_MASTER)call stop_timer
       call ed_reset_suffix
    end select
    !
  end subroutine ed_solve_lattice

  subroutine ed_solve_lattice_nobath(Nlat,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii,flag_gf)
    integer                                       :: Nlat
    real(8)                                       :: bath_dummy(Nlat,1) !user bath input array
    logical,optional                              :: mpi_lanc  !parallelization strategy flag: if :code:`.false.` each core serially solves an inequivalent site, if :code:`.true.` all cores parallely solve each site in sequence. Default :code:`.false.` .
    real(8),optional,dimension(Nlat,Norb) :: Uloc_ii !site-dependent values for :f:var:`uloc` , overriding the ones in the input file. It has dimension [ :f:var:`nlat` , :f:var:`norb` ].
    real(8),optional,dimension(Nlat)      :: Ust_ii  !site-dependent values for :f:var:`ust` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(Nlat)      :: Jh_ii   !site-dependent values for :f:var:`jh` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(Nlat)      :: Jp_ii   !site-dependent values for :f:var:`jp` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    real(8),optional,dimension(Nlat)      :: Jx_ii   !site-dependent values for :f:var:`jx` , overriding the ones in the input file. It has dimension [ :f:var:`nlat`].
    logical,optional                      :: flag_gf   !flag to calculate ( :code:`.true.` ) or not ( :code:`.false.` ) Green's functions and susceptibilities. Default :code:`.true.` .   
    logical                               :: mpi_lanc_,flag_gf_

    mpi_lanc_=.false.;if(present(mpi_lanc))mpi_lanc_=mpi_lanc
    flag_gf_=.true.;if(present(flag_gf))flag_gf=flag_gf_    
    bath_dummy=zero
    call ed_solve_lattice(bath_dummy,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii,flag_gf)
  end subroutine ed_solve_lattice_nobath


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: deallocate and finalize one or multiple baths -+!
  subroutine ed_finalize_solver_lattice(Nineq)
    integer                              :: ilat
    integer                              :: Nineq !number of inequivalent impurity sites for real-space DMFT
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(phisc_ineq))deallocate(phisc_ineq)
    if(allocated(argsc_ineq))deallocate(argsc_ineq)
    if(allocated(exct_ineq))deallocate(exct_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    if(allocated(Hloc_anomalous_ineq))deallocate(Hloc_anomalous_ineq)

    if(allocated(single_particle_density_matrix_ineq))deallocate(single_particle_density_matrix_ineq)
    if(allocated(impurity_density_matrix_ineq))deallocate(impurity_density_matrix_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       call ed_finalize_solver()
    enddo
    call ed_reset_suffix
    !
  end subroutine ed_finalize_solver_lattice

end module E2I_MAIN




