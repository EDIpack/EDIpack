MODULE ED_INPUT_VARS
  !:synopsis: User-accessible input variables
  !Contains all global input variables which can be set by the user through the input file. A specific preocedure :f:func:`ed_read_input` should be called to read the input file using :f:func:`parse_input_variable` procedure from SciFortran. All variables are automatically set to a default, looked for and updated by reading into the file and, sequentially looked for and updated from command line (std.input) using the notation `variable_name=variable_value(s)` (case independent).  
  !
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str,free_unit,to_upper,to_lower
  USE ED_VERSION
  use iso_c_binding
  implicit none

  !======================================!
  ! These variables are bound to c types ! 
  !======================================!

  integer(c_int),bind(c, name="Nbath")                               :: Nbath             !
  !Number of bath sites:
  ! * :f:var:`bath_type` = :code:`normal` : number of bath sites per orbital
  ! * :f:var:`bath_type` = :code:`hybrid` : total number of bath sites
  ! * :f:var:`bath_type` = :code:`replica/general` : number of replicas 
  ! :Default Nbath:`6`
  !
  integer(c_int),bind(c, name="Norb")                                :: Norb              !
  !Number of impurity orbitals (max :code:`5` )
  ! :Default Norb:`1`
  !
  integer(c_int),bind(c, name="Nspin")                               :: Nspin             !
  !If :code:`1`, assume :math:`H_{\downarrow}` =
  ! :math:`H_{\uparrow}` . If :code:`2` , the Hamiltonian
  ! needs to explicitly include spin-up and spin-down blocks
  ! :Default Nspin:`1`
  !
  integer(c_int),bind(c, name="Nloop")                               :: Nloop             !
  !Maximum number of DMFT loops 
  ! :Default Nloop:`100`
  !
  integer(c_int),bind(c, name="Nph")                                 :: Nph               !
  !Max number of phonons allowed (cut off)
  ! :Default Nph:`0`
  !
  real(c_double),dimension(:),allocatable                            :: Uloc_ 
  real(c_double),dimension(5),bind(c, name="Uloc")                   :: Uloc              !
  !Values of the local interaction per orbital (max :code:`5` )
  ! :Default Uloc:`(/( 2d0,i=1,Norb )/)`
  !
  real(c_double),bind(c, name="Ust")                                 :: Ust               !
  !Value of the inter-orbital interaction term
  ! :Default Ust:`0d0`
  !
  real(c_double),bind(c, name="Jh")                                  :: Jh                !
  !Hund's coupling constant 
  ! :Default Jh:`0.d0`
  !
  real(c_double),bind(c, name="Jx")                                  :: Jx                !
  !Coupling constant for the spin-eXchange interaction term 
  ! :Default Jx:`0.d0`
  !
  real(c_double),bind(c, name="Jp")                                  :: Jp                !
  !Coupling constant for the Pair-hopping interaction term 
  ! :Default Jp:`0.d0`
  !
  real(c_double),bind(c, name="xmu")                                 :: xmu               !
  !Chemical potential. If :f:var:`HFMODE` = :code:`T`, :f:var:`xmu` = :code:`0.d0` 
  !indicates half-filling condition.
  ! :Default xmu:`0d0`
  !
  real(c_double),bind(c, name="beta")                                :: beta              !
  !Inverse temperature, at zero temperature it is used as a IR cut-off.
  ! :Default beta:`1000d0`
  !
  integer(c_int),bind(c, name="Nsuccess")                            :: Nsuccess          !
  !Number of repeated success to fall below convergence threshold 
  ! :Default Nsuccess:`1`
  !
  real(c_double),bind(c, name="dmft_error")                          :: dmft_error        !
  !Error threshold for DMFT convergence
  ! :Default dmft_error:`1d-5`
  !
  real(c_double),bind(c, name="eps")                                 :: eps               !
  !Broadening on the real frequency axis for Green's function and Susceptibility 
  !calculations. :Default eps:`1d-2`
  !
  real(c_double),bind(c, name="wini")                                :: wini              !
  !Minimum value of the real frequency range
  ! :Default wini:`-5d0`
  !
  real(c_double),bind(c, name="wfin")                                :: wfin              !
  !Maximum value of the real frequency range
  ! :Default wfin:`5d0`
  !
  real(c_double),bind(c, name="xmin")                                :: xmin              !
  !Minimum of the x-range for the local lattice 
  ! probability distribution function (phonons)
  ! :Default xmin:`-3d0`
  !
  real(c_double),bind(c, name="xmax")                                :: xmax              !
  !Maximum of the x-range for the local lattice 
  !probability distribution function (phonons)
  ! :Default xmax:`3d0`
  real(c_double),bind(c, name="sb_field")                            :: sb_field          !
  !Value of a symmetry breaking field for magnetic solutions
  ! :Default sb_field:`1d-1`  
  real(c_double),bind(c, name="nread")                               :: nread             !
  !Target occupation value for fixed-density calculations. If set to :code:`0.0` 
  !the calculation is assumed to be at fixed :f:var:`xmu`
  ! :Default nread:`0d0`

  !-----------------------------------------------------------------!
  ! These variable need an equivalent "internal" one because        !
  ! parse_input_variable is not capable of reading logical(c_bool)  !
  !-----------------------------------------------------------------!

  logical(c_bool),bind(c, name="ed_total_ud")                        :: ed_total_ud       !
  !Flag to select which type of quantum numbers have to be considered (if :f:var:`ed_mode` = :code:`normal`)
  ! * :code:`T` : blocks have different total :math:`N_{\uparrow}` and :math:`N_{\downarrow}`
  ! * :code:`F` : blocks have different total :math:`N^{\alpha}_{\uparrow}` and :math:`N^{\alpha}_{\downarrow}`
  !where :math:`\alpha` is the orbital index. Speeds up calculation in the case where orbitals are not hybridized 
  !:Default ed_total_ud:`T`
  !
  logical                                                            :: ed_total_ud_
  logical(c_bool),bind(c, name="ed_twin")                            :: ed_twin           !
  !Flag to reduce ( :code:`T` ) or not (:code:`F` ) the number of visited sector using twin symmetry
  ! :Default ed_twin:`F`
  !
  logical                                                            :: ed_twin_
  
  logical(c_bool),bind(c, name="ed_read_umatrix")                    :: ed_read_umatrix   !
  !Flag to enable ( :code:`T` ) or not (:code:`F` ) reading the two-body terms from an external file
  !defined by :f:var:`umatrix_file`
  ! :Default ed_read_umatrix:`F`
  !
  logical                                                            :: ed_read_umatrix_
  
  logical(c_bool),bind(c, name="ed_use_kanamori")                    :: ed_use_kanamori   !
  !Flag to enable ( :code:`T` ) or not (:code:`F` ) the use of the input variables :f:var:`ULOC`, 
  !:f:var:`UST`, :f:var:`JH`, :f:var:`JX`, :f:var:`JP` for Hubbard-Kanamori coefficients.
  ! :Default ed_use_kanamori:`T`
  !
  logical                                                            :: ed_use_kanamori_



  !==========================================!
  ! These variables are not bound to c types ! 
  !==========================================!

  logical              :: HFmode            !
  !Flag to set the form of the Hubbard-Kanamori interaction
  ! * :code:`T` : :math:`U(n_{\uparrow}-\frac{1}{2})(n_{\downarrow}-\frac{1}{2})`
  ! * :code:`F` : :math:`Un_{\uparrow}n_{\downarrow}` 
  ! :Default HFmode:`T`
  !
  real(8)              :: cutoff                        !
  !Spectrum cut-off, used to determine the number states to be retained
  ! :Default cutoff:`1d-9`
  !
  real(8)              :: gs_threshold                  !
  !Energy threshold for ground state degeneracy loop up
  ! :Default gs_threshold:`1d-9`
  !
  real(8)              :: deltasc                       !
  !Value of the SC symmetry breaking term (only used if :f:var:`ed_mode` = :code:`superc` )
  ! :Default deltasc:`2d-2`
  !
  integer              :: ph_type                       !
  !Shape of the e part of the e-ph interaction: 
  ! * :code:`1` = orbital occupation
  ! * :code:`2` = orbital hybridization
  ! :Default ph_type:`1`
  !
  real(8)              :: A_ph                          !
  !Phonon forcing field coupled to displacement operator (constant)
  ! :Default A_ph:`0d0`
  !
  complex(8),allocatable,dimension(:,:)  :: g_ph        !
  !Electron-phonon coupling constant all
  ! :Default g_ph:`zero`
  !
  real(8)                                :: w0_ph       !
  !Phonon frequency
  ! :Default w0_ph:`0d0`
  !
  real(8),allocatable,dimension(:)       :: g_ph_diag  
  !Diagonal electron-phonon density coupling constant
  ! :Default g_ph_diag:`zero`
  !
  real(8),allocatable,dimension(:)  :: spin_field_x     !
  !Magnetic field per orbital coupling to X-spin component 
  ! :Default spin_field_x:`zero`
  !
  real(8),allocatable,dimension(:)  :: spin_field_y     !
  !Magnetic field per orbital coupling to Y-spin component 
  ! :Default spin_field_y:`zero`
  !
  real(8),allocatable,dimension(:)  :: spin_field_z     !
  !Magnetic field per orbital coupling to Z-spin component 
  ! :Default spin_field_z:`zero`
  !
  real(8),allocatable,dimension(:)  :: pair_field       !
  !Pair field per orbital coupling to s-wave order parameter component 
  !which explicitly appears in the impurity Hamiltonian 
  ! :Default pair_field:`zero`
  !
  real(8),dimension(4) :: exc_field                     !
  !External field coupling to exciton order parameter 
  ! :Default exc_field:`zero`
  !
  logical              :: rdm_flag          !
  !Flag to activate Reduced Density Matrix evaluation 
  ! :Default rdm_flag:`F`
  !
  logical              :: chispin_flag      !
  !Flag to activate spin susceptibility evaluation 
  ! :Default chispin_flag:`F`
  !
  logical              :: chidens_flag       !
  !Flag to activate charge susceptibility evaluation 
  ! :Default chidens_flag:`F`
  !
  logical              :: chipair_flag       !
  !Flag to activate pairing susceptibility evaluation 
  ! :Default chipair_flag:`F`
  !
  logical              :: chiexct_flag       !
  !Flag to activate excitonic susceptibility evaluation 
  ! :Default chiexct_flag:`F`
  !
  character(len=7)     :: ed_mode            !
  !Flag to set the ED mode 
  ! * :code:`normal` : normal
  ! * :code:`superc` : s-wave superconductive (singlet pairing)
  ! * :code:`nonsu2` : broken SU(2) symmetry
  ! :Default ed_mode:`normal`
  !
  logical              :: ed_finite_temp     !
  !Flag to set whether the problem is at finite temperature.
  !If :code:`T` then :f:var:`lanc_nstates_total` must be greater than 1.
  ! :Default ed_finite_temp:`F`
  !
  logical              :: ed_sparse_H        !
  !Flag to select  storage of the Fock space Hamiltonian as a sparse matrix
  ! * :code:`T` : :math:`H` is stored. CPU intensive
  ! * :code:`F` : on-the-fly :math:`H \cdot v` product is stored. Memory intensive
  ! :Default ed_sparse_H:`T`
  !
  logical              :: ed_solve_offdiag_gf !
  !Flag to select the calculation of the off-diagonal
  ! impurity GF. Set to :code:`T` by default if :f:var:`bath_type` is not :code:`normal`
  ! :Default ed_solve_offdiag_gf:`F`
  !
  logical              :: ed_print_Sigma    !
  !Flag to print impurity Self-energies 
  ! :Default ed_print_Sigma:`T`
  !
  logical              :: ed_print_G        !
  !Flag to print impurity Green`s functions 
  ! :Default ed_print_G:`T`
  !
  logical              :: ed_print_G0       !
  !Flag to print impurity non-interacting Green`s functions 
  ! :Default ed_print_G0:`T`
  !
  logical              :: ed_print_chispin  !
  !Flag to print impurity spin susceptibility, if calculated 
  ! :Default ed_print_chispin:`T`
  !
  logical              :: ed_print_chidens  !
  !Flag to print impurity dens susceptibility, if calculated 
  ! :Default ed_print_chidens:`T`
  !
  logical              :: ed_print_chipair  !
  !Flag to print impurity pair susceptibility, if calculated 
  ! :Default ed_print_chipair:`T`
  !
  logical              :: ed_print_chiexct  !
  !Flag to print impurity exct susceptibility, if calculated 
  ! :Default ed_print_chiexct:`T`
  !
  logical              :: ed_all_G          !
  !Flag to evaluate all the components of the impurity Green`s functions irrespective of the symmetries
  ! :Default ed_all_G:`F`
  !
  logical              :: ed_sectors        !
  !Flag to reduce sector scan for the spectrum to specific sectors
  ! :Default ed_sectors:`F`
  !
  integer              :: ed_sectors_shift  !
  !Additional sectors to be evaluated if :f:var:`ed_sectors` is set. These are sectors 
  !with all the quantum numbers varying of at most by :f:var:`ed_sectors_shift` around the 
  !sectors listed in :f:var:`ed_sectors`.
  ! :Default ed_sectors_shift:`1`
  !
  integer              :: ed_verbose        !
  !Verbosity level:
  ! * :code:`0` : almost nothing 
  ! * ...                       
  ! * :code:`3` : most of the verbose output
  ! * ...
  ! * :code:`5` : everything. Really, everything
  ! :Default ed_verbose:`3`
  !
  real(8)              :: ed_offset_bath    !
  !Offset for the initialization of diagonal terms if :f:var:`bath_type` = :code:`replica, general` . 
  !The replicas will be equally spaced in the range :code:`[-offset,offset]`
  ! :Default ed_offset_bath:`1d-1`
  !
  real(8)              :: ed_hw_bath        !
  !Half-bandwidth for bath level initialization if :f:var:`bath_type` = :code:`normal, hybrid` . 
  !The levels will be equispaced in the range  :code:`[-hw,hw]`
  ! :Default ed_hw_bath:`2d0`
  !
  logical              :: ed_obs_all        !
  !Flag to print observables for every loop
  ! :Default ed_obs_all:`T`
  !
  character(len=12)    :: lanc_method       !
  !Flag to select the Lanczos method to be used in the 
  !determination of the spectrum. 
  ! * :code:`ARPACK` : uses the Arnoldi algorithm 
  ! * :code:`LANCZOS`: uses an in-house Lanczos algorithm (limited to zero temperature)
  ! :Default lanc_method:`ARPACK`
  !
  real(8)              :: lanc_tolerance    !
  !Tolerance for the Lanczos iterations as used in Arpack and plain Lanczos
  ! :Default lanc_tolerance:`1d-18`
  !
  integer              :: lanc_niter         !
  !Max number of Lanczos iterations
  ! :Default lanc_niter:`512`
  !
  integer              :: lanc_ngfiter       !
  !Number of Lanczos iteration in GF determination. Number of moments.
  ! :Default lanc_ngfiter:`200`
  !
  integer              :: lanc_ncv_factor    !
  !Size of the block used in Lanczos-Arpack by multiplying the required :code:`Neigen` according to 
  ! :math:`N_{cv}=\mathrm{lanc\_ncv\_factor} \cdot \mathrm{Neigen} + \mathrm{lanc\_ncv\_add}`
  ! :Default lanc_ncv_factor:`10`
  !
  integer              :: lanc_ncv_add       !
  !Offset to add to the size of the block to prevent it to become too small according to
  ! :math:`N_{cv}=\mathrm{lanc\_ncv\_factor} \cdot \mathrm{Neigen} + \mathrm{lanc\_ncv\_add}`
  ! :Default lanc_ncv_add:`0`
  !
  integer              :: lanc_nstates_sector !
  !Max number of required eigenvalues per sector
  ! :Default lanc_nstates_sector:`2`
  !
  integer              :: lanc_nstates_total  !
  !Max number of states contributing to the partition function for finite-temperature calculations.
  !It must be set to :code:`1` for zero-temperature calculations, and to a greater value for 
  !finite-temperature calculations.
  ! :Default lanc_nstates_total:`1`
  !
  integer              :: lanc_nstates_step   !
  !Number of states added at each step for finite-temperature calculations: if the latest state included in the
  !partition function has a Boltzmann weight higher than :f:var:`cutoff`, :f:var:`lanc_nstates_total` will be increased 
  !by this amount at the next DMFT iteration.
  ! :Default lanc_nstates_step:`2`
  !
  integer              :: lanc_dim_threshold  !
  !Minimal sector dimension for Lanczos diagonalization. Smaller sectors will be solved with Exact Diagonalization 
  ! provided by Lapack
  ! :Default lanc_dim_threshold:`1024`
  !
  character(len=6)     :: cg_Scheme         !
  !Which quantity to use in the bath fitting routine:
  ! * :code:`WEISS` : the lattice Weiss field Green's function :math:`\mathcal{G}_{0}(i\omega_{n})` is fitted
  ! * :code:`DELTA` : the lattice Hybridization function :math:`\Delta(i\omega_{n})` is fitted
  ! :Default cg_Scheme:`Weiss`
  !
  integer              :: cg_method         !
  !Conjugate-gradient fitting routine to be used:
  ! * :code:`0` : Numerical Recipes
  ! * :code:`1` :minimize (FORTRAN77 code)
  ! :Default cg_method:`0`
  !
  integer              :: cg_grad           !
  !Flag to set the type of gradient evaluation (if :f:var:`cg_method` = :code:`0` ) 
  ! * :code:`0` : analytic
  ! * :code:`1` : numeric
  ! :Default cg_grad:`0`
  !
  integer              :: cg_Niter          !
  !Maximum number of iterations in the bath fitting procedure
  ! :Default cg_Niter:`500`
  !
  real(8)              :: cg_Ftol           !
  !Tolerance in the conjugate-gradient fitting procedure
  ! :Default cg_Ftol:`1d-5`
  !
  integer              :: cg_stop           !
  !Conjugate-gradient stopping condition
  ! * :code:`0` : :code:`1 .and. 2`
  ! * :code:`1` : :math:`\vert F_{n-1} -F_{n} \vert < \mathrm{tol} \cdot (1+F_{n})`
  ! * :code:`2` : :math:`\vert\vert x_{n-1} -x_{n} \vert\vert <\mathrm{tol} \cdot (1+ \vert\vert x_{n} \vert\vert`)
  ! :Default cg_stop:`0`
  !
  integer              :: cg_Weight         !
  !Weight assigned to the imaginary frequency data points in the calculation of the :math:`\chi^{2}`
  ! * :code:`0` : :math:`1`
  ! * :code:`1` : :math:`\frac{1}{n}`
  ! * :code:`2` : :math:`\frac{1}{\omega_{n}}`
  ! :Default cg_Weight:`0`
  !
  integer              :: cg_pow            !
  !Power exponent in the :math:`\chi^{2}` , according to 
  ! :math:`\vert \mathcal{G}_{0}(i\omega_{n}) - G_{0}^{\mathrm{And}}(i\omega_{n}) \vert ^{\mathrm{cg\_pow}}` or 
  ! :math:`\vert \Delta(i\omega_{n})- \Delta^{\mathrm{And}}(i\omega_{n}) \vert ^{\mathrm{cg\_pow}}`
  ! :Default cg_pow:`2`
  !
  character(len=12)    :: cg_norm           !
  !Which norm to use in the evaluation of the :math:`\chi^{2}` for matrix quantities. 
  ! * :code:`ELEMENTAL` : :math:`\chi^{2}` is the sum of each component's :math:`\chi^{2}`
  ! * :code:`FROBENIUS` : :math:`\chi^{2}` is calculated on the Frobenius norm (Matrix distance)
  ! :Default cg_norm:`ELEMENTAL`
  !.. warning:: 
  !   The Frobenius norm is currently implemented only for :f:var:`ed_bath` = :code:`replica, general`. 
  !   Also, for :f:var:`cg_pow` :math:`\neq 2`  the Frobenius norm is ill-defined, at least with respect 
  !   to its usual mathematical meaning. The behavior of :code:`cg_norm=frobenius` might be changed or 
  !   removed in future versions of the code, breaking back-compatibility.
  !
  logical              :: cg_minimize_ver   !
  !If :f:var:`cg_grad` = :code:`1` , select which version of :code:`minimize.f` to use
  ! * :code:`T` : Lichtenstein (newer)
  ! * :code:`F` : Krauth (older)
  ! :Default cg_minimize_ver:`F`
  !
  real(8)              :: cg_minimize_hh    !
  !Unknown parameter used in the CG minimize procedure ( for :f:var:`cg_grad` = :code:`1` )
  ! :Default cg_minimize_hh:`1d-4`
  !
  logical              :: finiteT           !
  !Flag to set finite-temperature calculation
  !
  character(len=8)     :: bath_type         !
  !Flag to select the bath geometry 
  ! * :code:`normal`  : each impurity orbital has a set of bath levels in a star geometry
  ! * :code:`hybrid`  : all impurity orbitals communicate with the same set of bath levels in a star geometry
  ! * :code:`replica` : the impurity communicates with clusters of the same form via an hybridization term :math:`V\mathbb{I}`  
  ! * :code:`general` : extends :code:`replica` so that each orbital has a different hybridization strength
  ! :Default bath_type:`normal`
  !
  real(8)              :: nerr              !
  !Error threshold for fixed-density calculations
  ! :Default nerr:`1d-4`
  !
  real(8)              :: ndelta            !
  !Initial chemical potential variation for fixed-density calculations
  ! :Default ndelta:`1d-1`
  !
  real(8)              :: ncoeff            !
  !Multiplier for the :code:`ndelta` value if :f:var:`xmu` and its error are 
  !read from a file ( :math:`\mathrm{ndelta} \rightarrow \mathrm{ndelta} \cdot \mathrm{ncoeff}` )
  ! :Default ncoeff:`1d0`
  !
  integer              :: niter             !
  !
  logical              :: Jz_basis          !
  !Flag to enable the :math:`J_{z}` basis in SOC calculations 
  ! :Default Jz_basis:`F`
  !
  logical              :: Jz_max            !
  !Flag to enable a maximum value for :math:`J_{z}` 
  ! :Default Jz_max:`F`
  !
  real(8)              :: Jz_max_value      !
  !Maximum value for Jz 
  ! :Default Jz_max_value:`1000d0`

  !Some parameters for function dimensions:
  !-----------------------------------------

  integer(c_int),bind(c, name="Lmats")             :: Lmats !
  !Number of Matsubara frequencies 
  ! :Default Lmats:`4096`
  !
  integer(c_int),bind(c, name="Lreal")             :: Lreal !
  !Number of real-axis frequencies 
  !:Default Lreal:`5000`
  !
  integer(c_int),bind(c, name="Lfit")              :: Lfit  !
  !Number of frequencies for bath fitting 
  !:Default Lfit:`1000`
  !
  integer(c_int),bind(c, name="Ltau")              :: Ltau  !
  !Number of imaginary time points 
  ! :Default Ltau:`1024`
  !
  integer(c_int),bind(c, name="Lpos")              :: Lpos  !
  !Number of points in Probability Distribution Function lattice for phonons 
  !:Default Lpos:`100`

  !LOG AND Hamiltonian UNITS
  !--------------------------
  character(len=100)                                           :: Hfile      !
  !File where to retrieve/store the bath parameters 
  ! :Default Hfile:`hamiltonian[.used/restart]`
  !
  character(len=100)                                           :: Bfile      !
  !File where to retrieve/store the bath parameters 
  ! :Default Bfile:`hbasis[.used/restart]`
  !
  character(len=100)                                           :: HLOCfile   !
  !File to read the input local H from
  ! :Default HLOCfile:`inputHLOC.in`
  !
  character(len=100)                                           :: Umatrix_File !
  !File containing the list of two-body operators
  ! :Default umatrix_file:`umatrix[.used/restart]`
  !
  character(len=100)                                           :: SectorFile !
  !File where to retrieve/store the sectors contributing to the spectrum 
  ! :Default SectorFile:`sectors[.used/restart]`
  !
  character(len=100)                                           :: GPHfile    !
  !File of Phonon couplings. Set to NONE to use only density couplings.
  ! :Default GPHfile:`NONE`
  !
  integer(c_int),bind(c, name="LOGfile"),save                  :: LOGfile    !
  !Logfile unit  
  ! :Default LOGfile:`6`
  !
  logical                                                      :: print_input_vars   !
  !Flag to toggle the printing on the terminal output of a list of input variables and their 
  !values
  ! :Default print_input_vars:`T`
  !
  logical                                                      :: print_sector_eigenvalues   !
  !Flag to toggle the printing of the eigenvalues of each sectors to a file. 
  ! :Default print_sector_eigenvalues:`T`
  !
  !THIS IS JUST A RELOCATED GLOBAL VARIABLE
  character(len=200)                                 :: ed_input_file=""    !Name of input file

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit)
    !
    !This functions reads the input file provided by :code:`INPUTunit` and sets the global variables accordingly
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    logical          :: master=.true.,bool
    integer          :: i,iorb,rank=0
    integer          :: unit_xmu, unit_gph
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A,A)")"DEBUG ed_read_input: read input from",trim(INPUTunit)
#endif
#ifdef _MPI    
    if(check_MPI())then
       master=get_Master_MPI()
       rank  =get_Rank_MPI()
    endif
#endif
    !
    !
    !
    !Store the name of the input file:
    ed_input_file=str(INPUTunit)
    !
    !
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals (max 5).")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=6,comment="Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)(replica/general=>Nbath=Nreplica/Ngeneral)")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(Nph,"NPH",INPUTunit,default=0,comment="Max number of phonons allowed (cut off)")
    call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal',comment="flag to set bath type: normal (1bath/imp), hybrid(1bath), replica(1replica/imp), general(replica++)")
    !
    !
    !
    allocate(Uloc_(Norb)) !#TODO: put me back!
    call parse_input_variable(uloc_,"ULOC",INPUTunit,default=(/( 2d0,i=1,Norb )/),comment="Values of the local interaction per orbital (Norb)")
    Uloc=0.d0 
    Uloc(1:Norb)=Uloc_
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    !
    !
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0,comment="Error threshold for DMFT convergence")
    call parse_input_variable(sb_field,"SB_FIELD",INPUTunit,default=0.1d0,comment="Value of a symmetry breaking field for magnetic solutions.")
    call parse_input_variable(deltasc,"DELTASC",INPUTunit,default=0.02d0,comment="Value of the SC symmetry breaking term.")

    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")


    if(allocated(g_ph_diag))deallocate(g_ph_diag)
    allocate(g_ph_diag(Norb)) ! THIS the diagonal part of the g_ph matrix (see below)
    call parse_input_variable(g_ph_diag,"G_PH",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="Electron-phonon coupling density constant: diagonal part (Norb) if GPHfile=None")
    call parse_input_variable(w0_ph,"W0_PH",INPUTunit,default=0.d0,comment="Phonon frequency")
    call parse_input_variable(A_ph,"A_PH",INPUTunit,default=0.d0,comment="Forcing field coupled to phonons displacement operator")
    call parse_input_variable(GPHfile,"GPHfile",INPUTunit,default="NONE",comment="File of Phonon couplings. Put NONE to use only density couplings.")

    if(allocated(spin_field_x))deallocate(spin_field_x)
    if(allocated(spin_field_y))deallocate(spin_field_y)
    if(allocated(spin_field_z))deallocate(spin_field_z)
    if(allocated(pair_field))deallocate(pair_field)
    allocate(spin_field_x(Norb))
    allocate(spin_field_y(Norb))
    allocate(spin_field_z(Norb))
    allocate(pair_field(Norb))
    call parse_input_variable(spin_field_x,"SPIN_FIELD_X",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="magnetic field per orbital coupling to X-spin component (Norb)")
    call parse_input_variable(spin_field_y,"SPIN_FIELD_Y",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="magnetic field per orbital coupling to Y-spin component (Norb)")
    call parse_input_variable(spin_field_z,"SPIN_FIELD_Z",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="magnetic field per orbital coupling to Z-spin component (Norb)")
    call parse_input_variable(pair_field,"PAIR_FIELD",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="pair field per orbital coupling to s-wave order parameter component (4)")
    call parse_input_variable(exc_field,"EXC_FIELD",INPUTunit,default=[0d0,0d0,0d0,0d0],comment="external field coupling to exciton order parameters (Norb)")
    !
    !
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chidens_flag,"CHIDENS_FLAG",INPUTunit,default=.false.,comment="Flag to activate density susceptibility calculation.")
    call parse_input_variable(chipair_flag,"CHIPAIR_FLAG",INPUTunit,default=.false.,comment="Flag to activate pair susceptibility calculation.")
    call parse_input_variable(chiexct_flag,"CHIEXCT_FLAG",INPUTunit,default=.false.,comment="Flag to activate excitonis susceptibility calculation.")
    !
    !
    call parse_input_variable(ed_mode,"ED_MODE",INPUTunit,default='normal',comment="Flag to set ED type: normal=normal, superc=superconductive, nonsu2=broken SU(2)")
    call parse_input_variable(ed_finite_temp,"ED_FINITE_TEMP",INPUTunit,default=.false.,comment="flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1")
    call parse_input_variable(ed_sectors,"ED_SECTORS",INPUTunit,default=.false.,comment="flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift.")
    call parse_input_variable(ed_sectors_shift,"ED_SECTORS_SHIFT",INPUTunit,1,comment="shift to ed_sectors")
    call parse_input_variable(ed_sparse_H,"ED_SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    !
    !
    call parse_input_variable(ed_total_ud_,"ED_TOTAL_UD",INPUTunit,default=.true.,comment="flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw")
    ed_total_ud = ed_total_ud_
    call parse_input_variable(ed_twin_,"ED_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    ed_twin = ed_twin_
    call parse_input_variable(ed_read_umatrix_,"ED_READ_UMATRIX",INPUTunit,default=.false.,comment="flag to read (T) or not (F,default) the two-body operators from an external file.")
    ed_read_umatrix = ed_read_umatrix_
    call parse_input_variable(ed_use_kanamori_,"ED_USE_KANAMORI",INPUTunit,default=.true.,comment="flag to use (T,default) or not (F) input variables for Kanamori coefficients.")
    ed_use_kanamori = ed_use_kanamori_
    if(ed_read_umatrix .and. ed_use_kanamori) STOP "ED_READ_UMATRIX and ED_USE_KANAMORI cannot be both true."
    call parse_input_variable(ed_obs_all,"ED_OBS_ALL",INPUTunit,default=.true.,comment="flag to print observables for every loop.")
    !
    !
    call parse_input_variable(ed_solve_offdiag_gf,"ED_SOLVE_OFFDIAG_GF",INPUTunit,default=.false.,comment="flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal") 
    call parse_input_variable(ed_print_Sigma,"ED_PRINT_SIGMA",INPUTunit,default=.true.,comment="flag to print impurity Self-energies")
    call parse_input_variable(ed_print_G,"ED_PRINT_G",INPUTunit,default=.true.,comment="flag to print impurity Greens function")
    call parse_input_variable(ed_print_G0,"ED_PRINT_G0",INPUTunit,default=.true.,comment="flag to print non-interacting impurity Greens function")
    call parse_input_variable(ed_print_chispin,"ED_PRINT_CHISPIN",INPUTunit,default=.true.,comment="flag to print impurity spin susceptibility")
    call parse_input_variable(ed_print_chidens,"ED_PRINT_CHIDENS",INPUTunit,default=.true.,comment="flag to print impurity dens susceptibility")
    call parse_input_variable(ed_print_chipair,"ED_PRINT_CHIPAIR",INPUTunit,default=.true.,comment="flag to print impurity pair susceptibility")
    call parse_input_variable(ed_print_chiexct,"ED_PRINT_CHIEXCT",INPUTunit,default=.true.,comment="flag to print impurity exct susceptibility")
    call parse_input_variable(ed_all_G,"ED_ALL_G",INPUTunit,default=.true.,comment="flag to evaluate all the components of the impurity Greens functions irrespective of the symmetries")
    call parse_input_variable(ed_verbose,"ED_VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")
    call parse_input_variable(ed_hw_bath,"ed_hw_bath",INPUTunit,default=2d0,comment="half-bandwidth for the bath initialization: flat in -ed_hw_bath:ed_hw_bath")
    call parse_input_variable(ed_offset_bath,"ed_offset_bath",INPUTunit,default=1d-1,comment="offset for the initialization of diagonal terms in replica/general bath: -offset:offset")

    !
    !
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=4096,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1024,comment="Number of imaginary time points.")
    call parse_input_variable(Lfit,"LFIT",INPUTunit,default=1000,comment="Number of Matsubara frequencies used in the \Chi2 fit.")
    call parse_input_variable(Lpos,"LPOS",INPUTunit,default=100,comment="Number of points for the lattice PDF.")
    !
    !
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    call parse_input_variable(ncoeff,"NCOEFF",INPUTunit,default=1d0,comment="multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff).")
    !
    !
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(xmin,"XMIN",INPUTunit,default=-3.d0,comment="Smallest position for the lattice PDF")
    call parse_input_variable(xmax,"XMAX",INPUTunit,default=3.d0,comment="Largest position for the lattice PDF")
    call parse_input_variable(rdm_flag,"RDM_FLAG",INPUTunit,default=.false.,comment="Flag to activate RDM calculation.")
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chidens_flag,"CHIDENS_FLAG",INPUTunit,default=.false.,comment="Flag to activate density susceptibility calculation.")
    call parse_input_variable(chipair_flag,"CHIPAIR_FLAG",INPUTunit,default=.false.,comment="Flag to activate pair susceptibility calculation.")
    call parse_input_variable(chiexct_flag,"CHIEXCT_FLAG",INPUTunit,default=.false.,comment="Flag to activate excitonis susceptibility calculation.")
    !
    !
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    !
    !    
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only), DVDSON (no MPI)")
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=2,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=10,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=0,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=1d-18,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=1024,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    !
    !
    call parse_input_variable(cg_method,"CG_METHOD",INPUTunit,default=0,comment="Conjugate-Gradient method: 0=NumericalRecipes, 1=minimize.")
    call parse_input_variable(cg_grad,"CG_GRAD",INPUTunit,default=0,comment="Gradient evaluation method: 0=analytic (default), 1=numeric.")
    call parse_input_variable(cg_ftol,"CG_FTOL",INPUTunit,default=0.00001d0,comment="Conjugate-Gradient tolerance.")
    call parse_input_variable(cg_stop,"CG_STOP",INPUTunit,default=0,comment="Conjugate-Gradient stopping condition: 0-2, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).")
    call parse_input_variable(cg_niter,"CG_NITER",INPUTunit,default=500,comment="Max. number of Conjugate-Gradient iterations.")
    call parse_input_variable(cg_weight,"CG_WEIGHT",INPUTunit,default=1,comment="Conjugate-Gradient weight form: 1=1.0, 2=1/n , 3=1/w_n.")
    call parse_input_variable(cg_scheme,"CG_SCHEME",INPUTunit,default='weiss',comment="Conjugate-Gradient fit scheme: delta or weiss.")
    call parse_input_variable(cg_norm,"CG_NORM",INPUTunit,default='elemental',comment="Conjugate-Gradient norm definition: elemental (default) or frobenius.")
    call parse_input_variable(cg_pow,"CG_POW",INPUTunit,default=2,comment="Fit power for the calculation of the generalized distance as |G0 - G0and|**cg_pow")
    call parse_input_variable(cg_minimize_ver,"CG_MINIMIZE_VER",INPUTunit,default=.false.,comment="Flag to pick old/.false. (Krauth) or new/.true. (Lichtenstein) version of the minimize CG routine")
    call parse_input_variable(cg_minimize_hh,"CG_MINIMIZE_HH",INPUTunit,default=1d-4,comment="Unknown parameter used in the CG minimize procedure.")
    !
    if (cg_norm=='frobenius' .and. (bath_type/='replica' .or. bath_type/='general')) then
       print*, "WARNING: The Frobenius norm is currently not implemented for normal and hybrid bath types."
       print*, "         The elemental (usual element-wise chi^2) norm will be used instead."
    endif
    if(cg_norm=='frobenius' .and. cg_pow/=2)then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         We'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call parse_input_variable(Jz_basis,"JZ_BASIS",INPUTunit,default=.false.,comment="Flag to enable the Jz basis")
    call parse_input_variable(Jz_max,"JZ_MAX",INPUTunit,default=.false.,comment="Whether to cutoff Jz")
    call parse_input_variable(Jz_max_value,"JZ_MAX_VALUE",INPUTunit,default=1000.d0,comment="Maximum Jz")
    !
    !
    call parse_input_variable(SectorFile,"SectorFile",INPUTunit,default="sectors",comment="File where to retrieve/store the sectors contributing to the spectrum.")
    call parse_input_variable(Hfile,"Hfile",INPUTunit,default="hamiltonian",comment="File where to retrieve/store the bath parameters.")
    call parse_input_variable(Bfile,"Bfile",INPUTunit,default="hbasis",comment="File where to retrieve/store the H bath matrix basis.")
    call parse_input_variable(HLOCfile,"HLOCfile",INPUTunit,default="inputHLOC.in",comment="File read the input local H.")
    call parse_input_variable(umatrix_file,"umatrix_file",INPUTunit,default="umatrix",comment="File read the two-body operator list from.")
    call parse_input_variable(print_input_vars,"PRINT_INPUT_VARS",INPUTunit,default=.true.,comment="Flag to toggle console printing of input variables list")
    call parse_input_variable(print_sector_eigenvalues,"PRINT_SECTOR_EIGENVALUES",INPUTunit,default=.true.,comment="Flag to toggle printing to file  of sector eigenvalues list")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")

    !Parameter checks:
    select case (to_lower(trim(bath_type)))
      case("normal","hybrid","replica","general")
      case default ; STOP "ed_read_input: bath_type "//trim(bath_type)//" not valid"
    end select
    
    select case (to_lower(trim(ed_mode)))
      case("normal","superc","nonsu2")
      case default ; STOP "ed_read_input: ed_mode "//trim(ed_mode)//" not valid"
    end select
    
    select case (to_lower(trim(lanc_method)))
      case("arpack","lanczos","dvdson")
      case default ; STOP "ed_read_input: lanc_method "//trim(lanc_method)//" not valid"
    end select

    select case (to_lower(trim(cg_scheme)))
      case("weiss","delta")
      case default ; STOP "ed_read_input: cg_scheme "//trim(cg_scheme)//" not valid"
    end select
    
    select case (to_lower(trim(cg_norm)))
      case("elemental","frobenius")
      case default ; STOP "ed_read_input: cg_norm "//trim(cg_norm)//" not valid"
    end select



    if(nph>0)then
       !
       if(allocated(g_ph))deallocate(g_ph)
       allocate(g_ph(Norb,Norb)) ! THIS SHOULD BE A MATRIX Norb*Norb
       !
       !Here the non-diagonal (non-density) phononic coupling are read
       g_ph=0.d0
       if(to_upper(str(GPHfile))=="NONE")then
          do iorb=1,Norb
             g_ph(iorb,iorb)=g_ph_diag(iorb)
          enddo
       else
          inquire(file=trim(GPHfile),EXIST=bool)
          if(bool)then
             open(free_unit(unit_gph),file=GPHfile)
             do iorb=1,Norb
                read(unit_gph,*) g_ph(iorb,:)
             enddo
             close(unit_gph)
             !
             !maybe an assert_hermitian would be globally useful
             if(any(g_ph /= transpose(conjg(g_ph))))then
                stop "ERROR: non hermitian phonon coupling matrix (g_ph) in input"
             end if
          else
             stop "GPHfile/=NONE but there is no GPHfile with the provided name"
          endif
       endif
       !
       !TO BE PUT SOMEWHERE ELSE
       open(free_unit(unit_gph),file="GPHinput.used")
       do iorb=1,Norb
          write(unit_gph,*) g_ph(iorb,:)
       enddo
       close(unit_gph)
    end if

#ifdef _MPI
    if(check_MPI())then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI()
             if(i==rank .and. ed_verbose>1)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    !
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
       if(ed_verbose>0) call print_logo()
       call scifor_version()
       call code_version(version,name="EDIpack")
       if(print_input_vars)call print_input()
       call save_input(INPUTunit)
    endif
    !
    !
    if(nread .ne. 0d0) then
       inquire(file="xmu.restart",EXIST=bool)
       if(bool)then
          open(free_unit(unit_xmu),file="xmu.restart")
          read(unit_xmu,*)xmu,ndelta
          ndelta=abs(ndelta)*ncoeff
          close(unit_xmu)
          write(*,"(A,F9.7,A)")"Adjusting XMU to ",xmu," as per provided xmu.restart "
       endif
    endif
    !
    !Act on the input variable only after printing.
    !
    !In the new parser variables are hard-linked into the list:
    !
    !any change to the variable is immediately copied into the list... (if you delete .ed it won't be printed out)
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".ed")
    call substring_delete(umatrix_file,".restart")
    call substring_delete(umatrix_file,".ed")
    !
    !
    !final cleanups
    deallocate(Uloc_)
  end subroutine ed_read_input

  subroutine ed_update_input(name,vals)
    !
    !This functions updates some variables in the input file, namely
    !
    !:f:var:`exc_field`, :f:var:`pair_field`, :f:var:`exc_field`,
    !
    !:f:var:`spin_field_x`, :f:var:`spin_field_y`, and :f:var:`spin_field_z`.
    character(len=*)      :: name !the name of the variable to update
    real(8),dimension(:)  :: vals !the new value of the variable
    select case (name)
    case default
       stop "WRONG NAME ON ED_UPDATE_INPUT"
    case ("EXC_FIELD")
       if(size(vals)/=4)stop "WRONG SIZE IN ED_UPDATE_EXC_FIELD"
       exc_field=vals
    case ("PAIR_FIELD")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_PAIR_FIELD"
       pair_field=vals
    case ("SPIN_FIELD_X")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_SPIN_FIELD_X"
       spin_field_x=vals
    case ("SPIN_FIELD_Y")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_SPIN_FIELD_Y"
       spin_field_y=vals
    case ("SPIN_FIELD_Z")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_SPIN_FIELD_Z"
       spin_field_z=vals
    end select

  end subroutine ed_update_input




  subroutine substring_delete (s,sub)
    !
    !! S_S_DELETE2 recursively removes a substring from a string.
    !
    !    The remainder is left justified and padded with blanks.
    !
    !    The substitution is recursive, so
    !
    !    that, for example, removing all occurrences of "ab" from
    !
    !    "aaaaabbbbbQ" results in "Q".
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !
    !  Example:
    !
    !    S = 'Fred is not a jerk!'
    !
    !    call s_chop ( S, 9, 12 )
    !
    !    S = 'Fred is a jerk!    '
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop
  
  
  subroutine print_logo()
      write(*,*)" "
      write(*,*)"            ╔════╗"
      write(*,*)"            ║    ║"
      write(*,*)"            ╚════╝"
      write(*,*)"           𜰟 𜱁 𜱁 𜰠"
      write(*,*)"          𜰟  𜱁 𜱁  ╭╴╮"
      write(*,*)"         𜰟   𜱁 𜱁  ╰╴╯"
      write(*,*)"        ╭╴╮ ╭╴╮╭╴╮"
      write(*,*)"        ╰╴╯ ╰╴╯╰╴╯"
      write(*,*)"╔═╗  ╔╦╗  ╦  ╔═╗  ╔═╗  ╔═╗  ╦╔═"
      write(*,*)"║╣    ║║  ║  ╠═╝  ╠═╣  ║    ╠╩╗"
      write(*,*)"╚═╝  ═╩╝  ╩  ╩    ╩ ╩  ╚═╝  ╩ ╩"
      write(*,*)"Version " // trim(GIT_TAG)
      write(*,*)" "
  end subroutine print_logo


END MODULE ED_INPUT_VARS
