MODULE E2I_IO
  !:synopsis: Retrieve observables and response functions, inequivalent sites version
  !Contains a set of routines that retrieve quantities such as Green's functions, self-energies (see :f:mod:`ed_greens_functions` ) and observables (from :f:mod:`ed_observables` ) and pass them to the user, as well ass routines to read and store Green's function and self-energies.
  USE EDIPACK
  !
  USE E2I_VARS_GLOBAL
  USE E2I_AUX_FUNX
  !
  USE SF_LINALG
  USE SF_SPIN
  USE SF_ARRAYS,  only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  USE SF_MISC,    only: assert_shape
  implicit none
  private

  interface ed_get_gimp
     !This subroutine gets from the EDIpack library the value of the impurity Green's function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity Green's function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
     !    :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`] 
     !  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_gimp_lattice_n3
     module procedure :: ed_get_gimp_lattice_n4
     module procedure :: ed_get_gimp_lattice_n6
  end interface ed_get_gimp

  interface ed_get_dimp
     !This subroutine gets from the EDIpack library the value of the impurity phonon's Green's function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity phonon's Green's function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat`,  :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_dimp_lattice_n2
  end interface ed_get_dimp

  interface ed_get_sigma
     !| This subrotine gets from the EDIpack library the value of the self-energy calculated 
     ! on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !| The self-energy is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
     !    :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`] 
     !  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_sigma_lattice_n3
     module procedure :: ed_get_sigma_lattice_n4
     module procedure :: ed_get_sigma_lattice_n6
  end interface ed_get_sigma



  interface ed_get_g0imp
     !| This subroutine gets from the EDIpack library the value of the impurity non-interacting Green's function calculated 
     ! on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !| It autonomously decides whether the system is single-impurity or real-space DMFT based on the :f:var:`bath` shape
     !
     !The impurity non-interacting Green's function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nlat` :math:`\cdot` :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nlat` :math:`\cdot` :f:var:`nspin` 
     !    :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`] 
     !  * [:f:var:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !  * [:f:var:`nlat`, :f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :f:var:`lmats` / :f:var:`lreal`]
     !
     !The bath is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat`, :f:var:`nb`] for real-space DMFT, with :f:var:`nlat` the number of inequivalent impurity sites
     !
     !Where :f:var:`nb` is the length of the :f:var:`bath` array.
     !
     module procedure :: ed_get_g0imp_lattice_n3
     module procedure :: ed_get_g0imp_lattice_n4
     module procedure :: ed_get_g0imp_lattice_n6
  end interface ed_get_g0imp


  interface ed_get_spinChi
     !This subroutine gets from the EDIpack library the value of the impurity spin susceptibility function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity spin susceptibility function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat`,  :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_spinChi_lattice_n3
  end interface ed_get_spinChi

  interface ed_get_densChi
     !This subroutine gets from the EDIpack library the value of the impurity dens susceptibility function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity dens susceptibility function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat`,  :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_densChi_lattice_n3
  end interface ed_get_densChi

  interface ed_get_pairChi
     !This subroutine gets from the EDIpack library the value of the impurity pair susceptibility function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity pair susceptibility function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat`,  :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_pairChi_lattice_n3
  end interface ed_get_pairChi

  interface ed_get_exctChi
     !This subroutine gets from the EDIpack library the value of the impurity exct susceptibility function calculated 
     !on the Matsubara or real-frequency axis, with number of frequencies :f:var:`lmats` or :f:var:`lreal` .
     !
     !The impurity exct susceptibility function is an array having the following possible dimensions:
     !
     !  * [:f:var:`nlat`,  :f:var:`lmats` / :f:var:`lreal`]
     !
     module procedure :: ed_get_exctChi_lattice_n3
  end interface ed_get_exctChi


  !Observables
  interface ed_get_dens
     !This subroutine gets from the EDIpack library the value of the charge density and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, density for that orbital
     !  * [:f:var:`nlat`]: if :f:var:`iorb` (default = 1) is provided for real-space DMFT with :f:var:`nlat` impurities, density for that orbital for all impurity sites
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, density for all impurity sites and orbitals
     !
     module procedure :: ed_get_dens_n1
     module procedure :: ed_get_dens_n2
  end interface ed_get_dens

  interface ed_get_mag
     !This subroutine gets from the EDIpack library the value of the magnetization and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`component` and :f:var:`iorb` are provided for single-impurity DMFT, given magnetization component for that orbital
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities, magnetization for that orbital for all impurity sites
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, one magnetization component for all orbitals and impurity sites
     !  * [:f:var:`nlat`, :code:`3`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, all magnetization components for all orbitals and sites
     !
     module procedure :: ed_get_mag_n1
     module procedure :: ed_get_mag_n2
     module procedure :: ed_get_mag_n3
  end interface ed_get_mag

  interface ed_get_docc
     !This subroutine gets from the EDIpack library the value of the double occupation and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * scalar: if :f:var:`iorb` is provided for single-impurity DMFT, dobule-occupation for that orbital
     !  * [:f:var:`nlat`]: if :f:var:`iorb` (default = 1) is provided for real-space DMFT with :f:var:`nlat` impurities, double-occupation for that orbital for all impurity sites
     !  * [:f:var:`nlat`, :f:var:`norb`]: if :f:var:`nlat` is provided for real-space DMFT, double-occupation for all impurity sites and orbitals
     !
     module procedure :: ed_get_docc_n1
     module procedure :: ed_get_docc_n2
  end interface ed_get_docc

  interface ed_get_phi
     !This subroutine gets from the EDIpack library the modulus value of the superconducting order parameter :math:`\phi` ( :f:var:`ed_mode` = :code:`superc` ) and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * [:f:var:`nlat`]:  :math:`\phi` for the specified :f:var:`iorb` orbital and all impurity sites
     !  * [:f:var:`nlat`, :f:var:`norb`]:  :math:`\phi` for all impurity sites and orbitals
     !
     module procedure :: ed_get_phisc_n1
     module procedure :: ed_get_phisc_n3
  end interface ed_get_phi


  interface ed_get_arg
     !This subroutine gets from the EDIpack library the argument value of the superconducting order parameter  :math:`\theta=tan^{-1}(Im\phi/Re\phi)` ( :f:var:`ed_mode` = :code:`superc` ) and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * [:f:var:`nlat`]:  :math:`\theta` for the specified :f:var:`iorb` orbital and all impurity sites
     !  * [:f:var:`nlat`, :f:var:`norb`]:  :math:`\theta` for all impurity sites and orbitals
     !
     module procedure :: ed_get_argsc_n1
     module procedure :: ed_get_argsc_n3
  end interface ed_get_arg


  interface ed_get_exct
     !This subroutine gets from the EDIpack library the value of the excitonic order parameters :math:`X^a` ( :f:var:`ed_mode` = :code:`normal`,:code:`nonsu2` ) and passes it to the user.
     !
     !The :f:var:`self` variable can have the following dimensions:
     ! 
     !  * [:f:var:`nlat`]: if :f:var:`iorb` is provided :math:`X` between orbital and orbital+1 for :f:var:`component` and and all impurity sites
     !  * [:f:var:`nlat` , :f:var:`norb` , :f:var:`norb`]:  :math:`X` for all orbitals, for  :f:var:`component`  and all impurity sites.
     !  * [:f:var:`nlat` , 4, :f:var:`norb` , :f:var:`norb`]: :math:`X` for all orbitals, :f:var:`component` and all impurity sites. 
     !
     module procedure :: ed_get_exct_n1
     module procedure :: ed_get_exct_n3
     module procedure :: ed_get_exct_n4
  end interface ed_get_exct


  !Get Energies
  interface ed_get_eimp
     !This subroutine gets from the EDIpack library and passes to the user the array [ :f:var:`ed_epot` , :f:var:`ed_eint` , :f:var:`ed_ehartree` , :f:var:`ed_eknot` ].
     !These are the expectation values various contribution to the internal energy
     !
     !  * :f:var:`ed_epot` = energy contribution from the interaction terms, **including** the Hartree term
     !  * :f:var:`ed_eint` = energy contribution from the interaction terms, **excluding** the Hartree term
     !  * :f:var:`ed_ehartree` = :math:`-\frac{U}{2} \sum_{i} \langle n_{i\uparrow} + n_{i\downarrow} \rangle 
     !    -\frac{2U^{'}-J_{H}}{2} \sum_{i < j} \langle n_{i\uparrow}+n_{i\downarrow} + n_{i\downarrow}+n_{j\downarrow} \rangle
     !    +\frac{U}{4} + \frac{2U^{'}-J_{H}}{2}` for :math:`i,j` orbitals
     !  * :f:var:`ed_eknot` = kinetic term from the **local** 1-body Hamiltonian
     !
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`, :code:`4`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_eimp_n2
  end interface ed_get_eimp

  interface ed_get_epot
     !This subroutine gets from the EDIpack library and passes to the user the value of 
     !:f:var:`ed_epot`, the energy contribution from the interaction terms, **including** the Hartree term.
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_epot_n1
  end interface ed_get_epot

  interface ed_get_eint
     !This subroutine gets from the EDIpack library and passes to the user the value of 
     !:f:var:`ed_int`, the energy contribution from the interaction terms, **excluding** the Hartree term.
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_eint_n1
  end interface ed_get_eint

  interface ed_get_ehartree
     !This subroutine gets from the EDIpack library and passes to the user the value of the Hartree potential 
     !:f:var:`ed_ehartree`. The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_ehartree_n1
  end interface ed_get_ehartree

  interface ed_get_eknot
     !This subroutine gets from the EDIpack library and passes to the user the value
     !:f:var:`ed_eknot`, the kinetic term from the **local** 1-body Hamiltonian
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_eknot_n1
  end interface ed_get_eknot


  interface ed_get_doubles
     !This subroutine gets from the EDIpack library and passes to the user the array [ :f:var:`ed_dust` , :f:var:`ed_dund` , :f:var:`ed_dse` , :f:var:`ed_dph` ].
     !These are the expectation values of the two-body operators associated with the density-density inter-orbital interaction (with opposite and parallel spins), 
     !spin-exchange and pair-hopping.
     !
     !  * :f:var:`ed_dust` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\downarrow} + n_{i\downarrow}n_{j\uparrow}` for :math:`i,j` orbitals
     !  * :f:var:`ed_dund` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\uparrow}  + n_{i\downarrow}n_{j\downarrow}` for :math:`i,j` orbitals
     !  * :f:var:`ed_dse` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{j\uparrow}c_{i\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !  * :f:var:`ed_dph` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{i\downarrow}c_{j\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`, :code:`4`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_doubles_n2
  end interface ed_get_doubles

  interface ed_get_dust
     !This subroutine gets from the EDIpack library and passes to the user the value of 
     !:f:var:`ed_dust` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\downarrow} + n_{i\downarrow}n_{j\uparrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_dust_n1
  end interface ed_get_dust

  interface ed_get_dund
     !This subroutine gets from the EDIpack library and passes to the user the value of 
     !:f:var:`ed_dund` = :math:`\sum_{i < j} n_{i\uparrow}n_{j\uparrow}  + n_{i\downarrow}n_{j\downarrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_dund_n1
  end interface ed_get_dund

  interface ed_get_dse
     !This subroutine gets from the EDIpack library and passes to the user the value of 
     !:f:var:`ed_dse` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{j\uparrow}c_{i\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_dse_n1
  end interface ed_get_dse

  interface ed_get_dph
     !This subroutine gets from the EDIpack library and passes to the user the value of 
     !:f:var:`ed_dph` = :math:`\sum_{i < j} c^{\dagger}_{i\uparrow}c^{\dagger}_{i\downarrow}c_{j\downarrow}c_{j\uparrow}` for :math:`i,j` orbitals
     !The returned array can have the following dimensions:
     !
     !  * [:f:var:`nlat`]: for real-space DMFT with :f:var:`nlat` impurities
     !
     module procedure :: ed_get_dph_n1
  end interface ed_get_dph



  ! interface ed_get_sp_dm
  !    !This subroutine returns to the user the impurity density matrix.
  !    !The density matrix is an array having the following possible dimensions:
  !    ! 
  !    !  * [:code:`nlat`, :f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`] for real-space DMFT
  !    !
  !    module procedure :: get_density_matrix_lattice
  ! end interface ed_get_sp_dm





  public :: ed_get_gimp
  public :: ed_get_dimp
  public :: ed_get_sigma
  public :: ed_get_g0imp
  !
  public :: ed_get_spinChi
  public :: ed_get_densChi
  public :: ed_get_pairChi
  public :: ed_get_exctChi

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_phi
  public :: ed_get_arg
  public :: ed_get_exct
  public :: ed_get_eimp
  public :: ed_get_epot
  public :: ed_get_eint 
  public :: ed_get_ehartree
  public :: ed_get_eknot
  public :: ed_get_doubles
  public :: ed_get_dust
  public :: ed_get_dund
  public :: ed_get_dse
  public :: ed_get_dph
  ! public :: ed_get_sp_dm
  ! public :: ed_get_neigen_total
  ! public :: ed_get_quantum_SOC_operators



  !****************************************************************************************!
  !****************************************************************************************!


  character(len=64)                           :: suffix


  integer                                     :: ilat,jlat
  integer                                     :: iorb,jorb
  integer                                     :: ispin,jspin
  integer                                     :: is,js
  integer                                     :: i,j
  integer                                     :: L
  complex(8),dimension(:,:,:,:,:),allocatable :: F




contains


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!

#include "get_gimp.f90"
#include "get_dimp.f90"


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity self-energy 
  !+-----------------------------------------------------------------------------+!

#include "get_sigma.f90"


  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve non-interacting green's functions 
  !+--------------------------------------------------------------------------+!

#include "get_g0imp.f90"



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve spin.dens.pair.exct susceptibilties
  !+--------------------------------------------------------------------------+!

#include "get_chi.f90"



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!

#include "get_dens.f90"
#include "get_mag.f90"
#include "get_docc.f90"
#include "get_phi.f90"
#include "get_exct.f90"
#include "get_energy.f90"
#include "get_doubles.f90"


  ! !+--------------------------------------------------------------------------+!
  ! ! PURPOSE: Single Particle DM
  ! !+--------------------------------------------------------------------------+!
  ! include "get_sp_dm.f90"


  !   !+--------------------------------------------------------------------------+!
  !   ! PURPOSE: Single Particle DM
  !   !+--------------------------------------------------------------------------+!
  ! #if __INTEL_COMPILER
  !   #include "get_sp_dm.f90"
  ! #else
  !   include "get_sp_dm.f90"
  ! #endif













  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : get SOC operators
  ! !+-------------------------------------------------------------------+
  ! subroutine ed_get_quantum_soc_operators()
  !   !This subroutine gets and prints the values of the components :math:`\overrightarrow{L}`, :math:`\overrightarrow{S}`, :math:`\overrightarrow{J}`
  !   !in the chosen basis depending on :f:var:`jz_basis`, and prints them on the files :code:`"L_imp_"//reg(str(ndx))//".dat"` , 
  !   !:code:`"S_imp_"//reg(str(ndx))//".dat"` and :code:`"J_imp_"//reg(str(ndx))//".dat"` , where :code:`ndx` is the inequivalent
  !   !impurity site for real-space DMFT (if that is the case). The ordering of the results in the output files is described by comments
  !   !in the files themselves
  !   !
  !   if (allocated(single_particle_density_matrix)) then
  !      call ed_get_quantum_soc_operators_lattice()
  !   else
  !      call ed_get_quantum_SOC_operators_single()
  !   endif
  ! end subroutine ed_get_quantum_soc_operators




  ! !+----------------------------------------------------------------------------+
  ! !PURPOSE  : get number of spectrum eigenstates for inequivalent impurity sites
  ! !+----------------------------------------------------------------------------+ 
  ! subroutine ed_get_neigen_total(nlii,Nlat) 
  !   !In the case of inequivalent impurity sites, this function returns the number of eigenstates per impurity
  !   !site in the ED spectrum.
  !   integer                      :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !   integer,dimension(Nlat)      :: nlii !array containing the number of eigenstates per inequivalent impurity site
  !   nlii=0d0
  !   if(allocated(neigen_total_ineq))then
  !      if(Nlat>size(neigen_total_ineq)) stop "ed_get_neigen_total error: required N_sites > evaluated N_sites"
  !      nlii=neigen_total_ineq
  !   endif
  ! end subroutine ed_get_neigen_total





  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(vr))allocate(vr(Lreal))
    if(.not.allocated(tau))allocate(tau(Ltau))
    wm  = pi/beta*(2*arange(1,Lmats)-1)
    vm  = pi/beta*(2*arange(1,Lmats)-2)
    wr  = linspace(wini,wfin,Lreal)
    vr  = linspace(0d0,wfin,Lreal)
    tau = linspace(0d0,beta,Ltau)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
    if(allocated(vr))deallocate(vr)
  end subroutine deallocate_grids





END MODULE E2I_IO







