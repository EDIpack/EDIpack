MODULE ED_OBSERVABLES_NONSU2
!:synopsis: Routines for observables calculation, :code:`NONSU2` case
  !This module calculates a series of observables, and stores them in aptly named plain-text files. :f:var:`ed_mode` = :code:`nonsu2`
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_BATH
  USE ED_HAMILTONIAN_NONSU2
  !
  !
  implicit none
  private
  !
  public :: observables_nonsu2
  public :: local_energy_nonsu2


  real(8),dimension(:),allocatable      :: dens ! orbital-resolved charge density
  real(8),dimension(:),allocatable      :: dens_up ! orbital-resolved spin-:math:`\uparrow` electron density
  real(8),dimension(:),allocatable      :: dens_dw ! orbital-resolved spin-:math:`\downarrow` electron density
  real(8),dimension(:),allocatable      :: docc ! orbital-resolved double occupation
  real(8),dimension(:),allocatable      :: magx ! orbital-resolved magnetization ( :code:`x` component )
  real(8),dimension(:),allocatable      :: magy ! orbital-resolved magnetization ( :code:`y` component )
  real(8),dimension(:),allocatable      :: magz ! orbital-resolved magnetization ( :code:`z` component )
  real(8),dimension(:,:),allocatable    :: n2 ! :math:`\langle n_{i} n_{j} \rangle` for i,j orbitals
  real(8),dimension(:,:),allocatable    :: sz2! :math:`\langle S^{z}_{i} S^{z}_{j} \rangle` for i,j orbitals
  real(8),dimension(:,:),allocatable    :: exct_s0 ! excitonic order parameter :math:`\langle c^{\dagger}_{is}\sigma^{0}c_{js^{'}} \rangle`
  real(8),dimension(:,:),allocatable    :: exct_tx ! excitonic order parameter :math:`\langle c^{\dagger}_{is}\sigma^{x}c_{js^{'}} \rangle`
  real(8),dimension(:,:),allocatable    :: exct_ty ! excitonic order parameter :math:`\langle c^{\dagger}_{is}\sigma^{y}c_{js^{'}} \rangle`
  real(8),dimension(:,:),allocatable    :: exct_tz ! excitonic order parameter :math:`\langle c^{\dagger}_{is}\sigma^{z}c_{js^{'}} \rangle`
  real(8)                               :: s2tot ! :math:`\langle S_{z}^{2} \rangle`
  real(8)                               :: Egs ! Ground-state energy
  real(8)                               :: Ei
  !
  integer                               :: iorb,jorb,istate
  integer                               :: ispin,jspin
  integer                               :: isite,jsite
  integer                               :: ibath
  integer                               :: r,m,k,k1,k2,k3,k4
  integer                               :: iup,idw
  integer                               :: jup,jdw
  integer                               :: mup,mdw
  integer                               :: iph,i_el,isz
  real(8)                               :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                               :: gs_weight
  !
  real(8)                               :: peso
  real(8)                               :: norm
  !
  integer                               :: i,j,ii
  integer                               :: isector,jsector
  !
  complex(8),dimension(:),allocatable   :: vvinit
  complex(8),dimension(:),allocatable   :: v_state
  logical                               :: Jcondition
  !
  type(sector)                          :: sectorI,sectorJ



contains 



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_nonsu2()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Calculate the values of the local observables
    integer,dimension(2*Ns)      :: ib
    integer,dimension(2,Ns)      :: Nud
    integer,dimension(Ns)        :: IbUp,IbDw
    real(8),dimension(Norb)      :: nup,ndw,Sz,nt
    real(8),dimension(Norb,Norb) :: theta_upup,theta_dwdw
    real(8),dimension(Norb,Norb) :: theta_updw,theta_dwup
    real(8),dimension(Norb,Norb) :: omega_updw,omega_dwup
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG observables_nonsu2"
#endif
    !
    !LOCAL OBSERVABLES:
    ! density, 
    ! double occupancy, 
    ! magnetization, 
    ! orbital//spin correlations  
    ! superconducting order parameter, etc..
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magZ(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(magX(Norb),magY(Norb))
    allocate(exct_S0(Norb,Norb),exct_Tz(Norb,Norb))
    allocate(exct_Tx(Norb,Norb),exct_Ty(Norb,Norb))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    magx    = 0.d0
    magy    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    exct_s0 = 0d0
    exct_tz = 0d0
    exct_tx = 0d0
    exct_ty = 0d0
    theta_upup = 0d0
    theta_dwdw = 0d0
    theta_updw = 0d0
    theta_dwup = 0d0
    omega_updw = 0d0
    omega_dwup = 0d0
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_nonsu2: get local observables"
#endif    
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_nonsu2: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(Mpimaster)then
          call build_sector(isector,sectorI)
          do i = 1,sectorI%Dim
             gs_weight=peso*abs(v_state(i))**2
             !
             m  = sectorI%H(1)%map(i)
             ib = bdecomp(m,2*Ns)
             do iorb=1,Norb
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
             enddo
             sz = (nup-ndw)/2d0
             nt =  nup+ndw
             !
             !Evaluate averages of observables:
             do iorb=1,Norb
                dens(iorb)     = dens(iorb)      +  nt(iorb)*gs_weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*gs_weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*gs_weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*gs_weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*gs_weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight
                do jorb=iorb+1,Norb
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
          enddo
       endif
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
    !
    !
    !EVALUATE <SX> AND <SY>
    do iorb=1,Norb
       !
#ifdef _DEBUG
       if(ed_verbose>2)write(Logfile,"(A)")&
            "DEBUG observables_nonsu2: eval in-plane magnetization <Sx>, <Sy>, a:"//str(iorb)
#endif
       do istate=1,state_list%size
          isector = es_return_sector(state_list,istate)
          Ei      = es_return_energy(state_list,istate)
          v_state    =  es_return_cvec(state_list,istate)
          !
#ifdef _DEBUG
          if(ed_verbose>3)write(Logfile,"(A)")&
               "DEBUG observables_nonsu2: get contribution from state:"//str(istate)
#endif
          !
          peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
          peso = peso/zeta_function
          !
          !GET <(CDG_UP + CDG_DW)(C_UP + C_DW)> = 
          !<CDG_UP*C_UP> + <CDG_DW*C_DW> + <CDG_UP*C_DW + CDG_DW*C_UP> = 
          !<N_UP> + <N_DW> + <Sx> 
          !since <Sx> = <CDG_UP*C_DW + CDG_DW*C_UP>
          !
          !APPLY (c_{iorb,up} + c_{iorb,dw})|gs> =  [1,1].[C_{-1},C_{-1}].[iorb,iorb].[1,2]
          jsector = getCsector(1,1,isector)
          if(jsector/=0)then
             vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,iorb],[1,2],isector,jsector)                
             magx(iorb) = magx(iorb) + dot_product(vvinit,vvinit)*peso
             deallocate(vvinit)
          endif
          !
          !GET <(-i*CDG_UP + CDG_DW)(i*C_UP + C_DW)> = 
          !<CDG_UP*C_UP> + <CDG_DW*C_DW> - i<CDG_UP*C_DW - CDG_DW*C_UP> = 
          !<N_UP> + <N_DW> + <Sy>
          !since <Sy> = -i<CDG_UP*C_DW - CDG_DW*C_UP>
          !
          !APPLY (C_DW + i*C_UP)|gs> =  [1,xi].[C_{-1},C_{-1}].[iorb,iorb].[2,1]
          jsector = getCsector(1,1,isector)
          if(jsector/=0)then
             vvinit = apply_Cops(v_state,[one,xi],[-1,-1],[iorb,iorb],[2,1],isector,jsector)                
             if(Mpimaster)magy(iorb) = magy(iorb) + dot_product(vvinit,vvinit)*peso
             deallocate(vvinit)
          endif
          !
          if(allocated(v_state))deallocate(v_state)
       enddo
    enddo
    !So we have:
    if(MpiMaster)then
       !<Sx> = <(CDG_UP + CDG_DW)(C_UP + C_DW)> - <N_UP> - <N_DW>
       magx = magx - dens_up - dens_dw
       !<Sy> = <(-i*CDG_UP + CDG_DW)(i*C_UP + C_DW)> - <N_UP> - <N_DW 
       magy = magy - dens_up - dens_dw
    endif

#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif



    !EVALUATE EXCITON OP <S_ab> AND <T^x,y,z_ab>
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_nonsu2: eval excitoninc OP <S_av>, <T^{x,y,z}_ab>"
#endif
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_nonsu2: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(Mpimaster)call build_sector(isector,sectorI)
       !
       !<S_ab>  :=   <C^+_{a,up}C_{b,up} + C^+_{a,dw}C_{b,dw}>
       !<T^z_ab>:=   <C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw}>
       ! O_uu  = a_up + b_up
       ! O_dd  = a_dw + b_dw
       !|v_up> = O_uu|v>
       !|v_dw> = O_dd|v> 
       ! Theta_uu = <v_up|v_up>
       ! Theta_dd = <v_dw|v_dw>
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !
             !|v_up> = (C_aup + C_bup)|>
             jsector = getCsector(1,1,isector)
             if(jsector/=0)then
                vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[1,1],isector,jsector)
                theta_upup(iorb,jorb) = theta_upup(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                deallocate(vvinit)
             endif
             !
             !|v_dw> = (C_adw + C_bdw)|>
             jsector = getCsector(1,2,isector)
             if(jsector/=0)then
                vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[2,2],isector,jsector)
                theta_dwdw(iorb,jorb) = theta_dwdw(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                deallocate(vvinit)
             endif
          enddo
       enddo
       !
       !<T^x_ab>:=   <C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}>
       ! O_ud  = a_up + b_dw
       ! O_du  = a_dw + b_up
       !|v_ud> = O_ud|v>
       !|v_du> = O_du|v> 
       ! Theta_ud = <v_ud|v_ud>
       ! Theta_du = <v_du|v_du>       
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !
             !|v_ud> = (C_aup + C_bdw)|>
             jsector = getCsector(1,1,isector)
             if(jsector/=0)then
                vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[1,2],isector,jsector)
                theta_updw(iorb,jorb) = theta_updw(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                deallocate(vvinit)
             endif
             !
             !|v_du> = (C_adw + C_bup)|>
             jsector = getCsector(1,1,isector)
             if(jsector/=0)then
                vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[2,1],isector,jsector)
                theta_dwup(iorb,jorb) = theta_dwup(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                deallocate(vvinit)
             endif
          enddo
       enddo
       !
       !<T^y_ab>:= -i<C^+_{a,up}C_{b,dw} - C^+_{a,dw}C_{b,up}>
       ! K_ud  = a_up - i*b_dw
       ! K_du  = a_dw - i*b_up
       !|w_ud> = K_ud|v>
       !|w_du> = K_du|v> 
       ! Omega_ud = <v_ud|v_ud>
       ! Omega_du = <v_du|v_du>    
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !|w_ud> = (C_aup - xi*C_bdw)|>
             jsector = getCsector(1,1,isector)
             if(jsector/=0)then
                vvinit = apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[1,2],isector,jsector)
                omega_updw(iorb,jorb) = omega_updw(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                if(allocated(vvinit))deallocate(vvinit)
             endif
             !
             !|w_du> = (C_adw - i*C_bup)|>
             jsector = getCsector(1,1,isector)
             if(jsector/=0)then
                vvinit = apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[2,1],isector,jsector)
                omega_dwup(iorb,jorb) = omega_dwup(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                if(allocated(vvinit))deallocate(vvinit)
             endif
             !
          enddo
       enddo
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    ! <S_ab>  = Theta_uu + Theta_dd - n_a - n_b
    ! <T^z_ab>= Theta_uu - Theta_dd - m_a - m_b
    ! <T^x_ab>= Theta_ud + Theta_du - n_a - n_b
    ! <T^y_ab>= Omega_ud - Omega_du - m_a - m_b
    if(MpiMaster)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             exct_s0(iorb,jorb) = theta_upup(iorb,jorb) + theta_dwdw(iorb,jorb) - dens(iorb) - dens(jorb)
             exct_tz(iorb,jorb) = theta_upup(iorb,jorb) - theta_dwdw(iorb,jorb) - magZ(iorb) - magZ(jorb)
             exct_tx(iorb,jorb) = theta_updw(iorb,jorb) + theta_dwup(iorb,jorb) - dens(iorb) - dens(jorb)
             exct_ty(iorb,jorb) = omega_updw(iorb,jorb) - omega_dwup(iorb,jorb) - magZ(iorb) + magZ(jorb)
          enddo
       enddo
    endif

    !SINGLE PARTICLE IMPURITY DENSITY MATRIX
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_nonsu2:  eval single particle density matrix <C^+_a C_b>"
#endif
    if(allocated(single_particle_density_matrix)) deallocate(single_particle_density_matrix)
    allocate(single_particle_density_matrix(Nspin,Nspin,Norb,Norb));single_particle_density_matrix=zero
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG observables_nonsu2: get contribution from state:"//str(istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(Mpimaster)then
          call build_sector(isector,sectorI)
          !Diagonal densities
          do ispin=1,Nspin
             do iorb=1,Norb
                isite=iorb + (ispin-1)*Norb
                do m=1,sectorI%Dim
                   i  = sectorI%H(1)%map(m)
                   ib = bdecomp(i,2*Ns)
                   single_particle_density_matrix(ispin,ispin,iorb,iorb) = &
                        single_particle_density_matrix(ispin,ispin,iorb,iorb) + &
                        peso*ib(isite)*conjg(v_state(m))*v_state(m)
                enddo
             enddo
          enddo
          !off-diagonal
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if((bath_type=="normal").and.(iorb/=jorb))cycle
                      ! if(bath_type=="replica".and.Jz_basis)then
                      !    if((.not.dmft_bath%mask(ispin,jspin,iorb,jorb,1)).AND.&
                      !         (.not.dmft_bath%mask(ispin,jspin,iorb,jorb,2)))cycle
                      ! endif
                      isite=iorb + (ispin-1)*Norb
                      jsite=jorb + (jspin-1)*Norb
                      do m=1,sectorI%Dim
                         i  = sectorI%H(1)%map(m)
                         ib = bdecomp(i,2*Ns)
                         if((ib(jsite)==1).and.(ib(isite)==0))then
                            call c(jsite,i,r,sgn1)
                            call cdg(isite,r,k,sgn2)
                            j=binary_search(sectorI%H(1)%map,k)
                            single_particle_density_matrix(ispin,jspin,iorb,jorb) = &
                                 single_particle_density_matrix(ispin,jspin,iorb,jorb) + &
                                 peso*sgn1*v_state(m)*sgn2*conjg(v_state(j))
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
          call delete_sector(sectorI)
       endif
       if(allocated(v_state))deallocate(v_state)
    enddo
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
    !
    write(LOGfile,"(A,10f18.12,f18.12)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    write(LOGfile,"(A,10f18.12)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12)")    "magX"//reg(ed_file_suffix)//"=",(magX(iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12)")    "magY"//reg(ed_file_suffix)//"=",(magY(iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12)")    "magZ"//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
    !
    write(LOGfile,"(A)",advance="no")"exS0"//reg(ed_file_suffix)//"="
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          write(LOGfile,"(20(F18.12,1X))")exct_s0(iorb,jorb)
       enddo
    enddo
    write(LOGfile,"(A)",advance="no")"exTX"//reg(ed_file_suffix)//"="
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          write(LOGfile,"(20(F18.12,1X))")exct_tx(iorb,jorb)
       enddo
    enddo
    write(LOGfile,"(A)",advance="no")"exTY"//reg(ed_file_suffix)//"="
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          write(LOGfile,"(20(F18.12,1X))")exct_ty(iorb,jorb)
       enddo
    enddo
    write(LOGfile,"(A)",advance="no")"exTZ"//reg(ed_file_suffix)//"="
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          write(LOGfile,"(20(F18.12,1X))")exct_tz(iorb,jorb)
       enddo
    enddo
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
       ed_mag(1,iorb)  =magX(iorb)
       ed_mag(2,iorb)  =magY(iorb)
       ed_mag(3,iorb)  =magZ(iorb)
       do jorb=1,Norb
          ed_exct(1,iorb,jorb) = exct_S0(iorb,jorb)
          ed_exct(2,iorb,jorb) = exct_Tx(iorb,jorb)
          ed_exct(3,iorb,jorb) = exct_Ty(iorb,jorb)
          ed_exct(4,iorb,jorb) = exct_Tz(iorb,jorb)
       enddo
    enddo
    !
    ed_imp_info=[s2tot,egs]
    !

#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_mag)
       call Bcast_MPI(MpiComm,ed_exct)
       call Bcast_MPI(MpiComm,ed_imp_info)
       if(allocated(single_particle_density_matrix))call Bcast_MPI(MpiComm,single_particle_density_matrix)
    endif
#endif
    !
    if(MPIMASTER)then
       call write_observables()
    endif
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
    deallocate(magX,magY)
    deallocate(exct_S0,exct_Tz)
    deallocate(exct_Tx,exct_Ty)
  end subroutine observables_nonsu2







  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
subroutine local_energy_nonsu2()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  !Calculate the value of the local energy components
  integer,dimension(2*Ns) :: ib
  integer,dimension(2,Ns) :: Nud
  integer,dimension(Ns)   :: IbUp,IbDw
  real(8),dimension(Norb)         :: nup,ndw
  real(8),dimension(Nspin,Norb)   :: eloc
  !
#ifdef _DEBUG
  write(Logfile,"(A)")"DEBUG local_energy_nonsu2"
#endif
  !
  Egs     = state_list%emin
  ed_Ehartree= 0.d0
  ed_Eknot   = 0.d0
  ed_Epot    = 0.d0
  ed_Dust    = 0.d0
  ed_Dund    = 0.d0
  ed_Dse     = 0.d0
  ed_Dph     = 0.d0
  !
  !Get diagonal part of Hloc
  do ispin=1,Nspin
     do iorb=1,Norb
        eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
     enddo
  enddo
  !
  do istate=1,state_list%size
     isector = es_return_sector(state_list,istate)
     Ei      = es_return_energy(state_list,istate)
     v_state    =  es_return_cvec(state_list,istate)
     !
#ifdef _DEBUG
     if(ed_verbose>3)write(Logfile,"(A)")&
          "DEBUG local_energy_nonsu2: get contribution from state:"//str(istate)
#endif
     !
     peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
     peso = peso/zeta_function
     !
     if(Mpimaster)then
        !
        call build_sector(isector,sectorI)
        do i=1,sectorI%Dim
           gs_weight=peso*abs(v_state(i))**2
           m  = sectorI%H(1)%map(i)
           ib = bdecomp(m,2*Ns)
           do iorb=1,Norb
              nup(iorb)= dble(ib(iorb))
              ndw(iorb)= dble(ib(iorb+Ns))
           enddo
           !
           !start evaluating the Tr(H_loc) to estimate potential energy
           !LOCAL ENERGY
           ed_Eknot = ed_Eknot + dot_product(eloc(1,:),nup)*gs_weight + dot_product(eloc(Nspin,:),ndw)*gs_weight
           !==> HYBRIDIZATION TERMS I: same or different orbitals, same spins.
           do iorb=1,Norb
              do jorb=1,Norb
                 !SPIN UP
                 if((ib(iorb)==0).AND.(ib(jorb)==1))then
                    call c(jorb,m,k1,sg1)
                    call cdg(iorb,k1,k2,sg2)
                    j=binary_search(sectorI%H(1)%map,k2)
                    if(Jz_basis.and.j==0)cycle
                    ed_Eknot = ed_Eknot + impHloc(1,1,iorb,jorb)*sg1*sg2*v_state(i)*conjg(v_state(j))*peso
                 endif
                 !SPIN DW
                 if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
                    call c(jorb+Ns,m,k1,sg1)
                    call cdg(iorb+Ns,k1,k2,sg2)
                    j=binary_search(sectorI%H(1)%map,k2)
                    if(Jz_basis.and.j==0)cycle
                    ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*v_state(i)*conjg(v_state(j))*peso
                 endif
              enddo
           enddo
           !==> HYBRIDIZATION TERMS II: same or different orbitals, opposite spins.
           do iorb=1,Norb
              do jorb=1,Norb
                 !UP-DW
                 if((impHloc(1,Nspin,iorb,jorb)/=zero).AND.(ib(iorb)==0).AND.(ib(jorb+Ns)==1))then
                    call c(jorb+Ns,m,k1,sg1)
                    call cdg(iorb,k1,k2,sg2)
                    j=binary_search(sectorI%H(1)%map,k2)
                    if(Jz_basis.and.j==0)cycle
                    ed_Eknot = ed_Eknot + impHloc(1,Nspin,iorb,jorb)*sg1*sg2*v_state(i)*conjg(v_state(j))*peso
                 endif
                 !DW-UP
                 if((impHloc(Nspin,1,iorb,jorb)/=zero).AND.(ib(iorb+Ns)==0).AND.(ib(jorb)==1))then
                    call c(jorb,m,k1,sg1)
                    call cdg(iorb+Ns,k1,k2,sg2)
                    j=binary_search(sectorI%H(1)%map,k2)
                    if(Jz_basis.and.j==0)cycle
                    ed_Eknot = ed_Eknot + impHloc(Nspin,1,iorb,jorb)*sg1*sg2*v_state(i)*conjg(v_state(j))*peso
                 endif
              enddo
           enddo
           !
           !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
           !Euloc=\sum=i U_i*(n_u*n_d)_i
           !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
           do iorb=1,Norb
              ed_Epot = ed_Epot + Uloc_internal(iorb)*nup(iorb)*ndw(iorb)*gs_weight
           enddo
           !
           !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
           !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
           !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
           if(Norb>1)then
              do iorb=1,Norb
                 do jorb=iorb+1,Norb
                    ed_Epot = ed_Epot + Ust_internal(iorb,jorb)*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                    ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                 enddo
              enddo
           endif
           !
           !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
           !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
           !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
           !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
           if(Norb>1)then
              do iorb=1,Norb
                 do jorb=iorb+1,Norb
                    ed_Epot = ed_Epot + (Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                    ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                 enddo
              enddo
           endif
           !
           !SPIN-EXCHANGE (S-E) TERMS
           !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
           if(Norb>1.AND.(any((Jx_internal/=0d0)).OR.any((Jp_internal/=0d0))))then
              do iorb=1,Norb
                 do jorb=1,Norb
                    Jcondition=((iorb/=jorb).AND.&
                         (ib(jorb)==1)      .AND.&
                         (ib(iorb+Ns)==1)   .AND.&
                         (ib(jorb+Ns)==0)   .AND.&
                         (ib(iorb)==0))
                    if(Jcondition)then
                       call c(jorb,m,k1,sg1)
                       call c(iorb+Ns,k1,k2,sg2)
                       call cdg(jorb+Ns,k2,k3,sg3)
                       call cdg(iorb,k3,k4,sg4)
                       j=binary_search(sectorI%H(1)%map,k4)
                       if(Jz_basis.and.j==0)cycle
                       ed_Epot = ed_Epot + Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                       ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                    endif
                 enddo
              enddo
           endif
           !
           !
           !PAIR-HOPPING (P-H) TERMS
           !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
           !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
           if(Norb>1.AND.(any((Jx_internal/=0d0)).OR.any((Jp_internal/=0d0))))then
              do iorb=1,Norb
                 do jorb=1,Norb
                    Jcondition=((iorb/=jorb).AND.&
                         (ib(jorb)==1)      .AND.&
                         (ib(jorb+Ns)==1)   .AND.&
                         (ib(iorb+Ns)==0)   .AND.&
                         (ib(iorb)==0))
                    if(Jcondition)then
                       call c(jorb,m,k1,sg1)
                       call c(jorb+Ns,k1,k2,sg2)
                       call cdg(iorb+Ns,k2,k3,sg3)
                       call cdg(iorb,k3,k4,sg4)
                       j=binary_search(sectorI%H(1)%map,k4)
                       if(Jz_basis.and.j==0)cycle
                       ed_Epot = ed_Epot + Jp_internal(iorb,jorb)*sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                       ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*v_state(i)*conjg(v_state(j))*peso
                    endif
                 enddo
              enddo
           endif
           !
           !
           !HARTREE-TERMS CONTRIBUTION:
           if(hfmode)then
              !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
              do iorb=1,Norb
                 ed_Ehartree=ed_Ehartree - 0.5d0*Uloc_internal(iorb)*(nup(iorb)+ndw(iorb))*gs_weight + 0.25d0*Uloc_internal(iorb)*gs_weight
              enddo
              if(Norb>1)then
                 do iorb=1,Norb
                    do jorb=iorb+1,Norb
                       ed_Ehartree=ed_Ehartree - 0.5d0*Ust_internal(iorb,jorb)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.5d0*Ust_internal(iorb,jorb)*gs_weight
                       ed_Ehartree=ed_Ehartree - 0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*gs_weight
                    enddo
                 enddo
              endif
           endif
        enddo
        call delete_sector(sectorI)
     endif
     !
     if(allocated(v_state))deallocate(v_state)
     !
  enddo
  !
#ifdef _DEBUG
  write(Logfile,"(A)")""
#endif
#ifdef _MPI
  if(MpiStatus)then
     call Bcast_MPI(MpiComm,ed_Epot)
     call Bcast_MPI(MpiComm,ed_Eknot)
     call Bcast_MPI(MpiComm,ed_Ehartree)
     call Bcast_MPI(MpiComm,ed_Dust)
     call Bcast_MPI(MpiComm,ed_Dund)
     call Bcast_MPI(MpiComm,ed_Dse)
     call Bcast_MPI(MpiComm,ed_Dph)
  endif
#endif
  !
  ed_Eint = ed_Epot
  ed_Epot = ed_Epot + ed_Ehartree
  !
  if(ed_verbose>=3)then
     write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
     write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
     write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
     write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
     write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
     write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
     write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
     write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
  endif
  !
  if(MPIMASTER)then
     call write_energy()
  endif
  !
  !
end subroutine local_energy_nonsu2



  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################








  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    !Write a plain-text file called :code:`observables_info.ed` detailing the names and contents of the observable output files.
    !Write the observable output files. Filenames with suffix :code:`_all` contain values for all DMFT interations, those with suffix :code:`_last` 
    !only values for the last iteration
    call write_obs_info()
    call write_obs_last()
    if(ed_obs_all)call write_obs_all()
  end subroutine write_observables




  subroutine write_obs_info()
    integer :: unit,iorb,jorb,ispin
    !Parameters used:
    if(.not.ed_read_umatrix)then
      unit = free_unit()
      open(unit,file="parameters_info.ed")
      write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
           (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
           reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
      close(unit)
    endif
    !
    !Generic observables 
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,*(A10,6X))")"#",&
         (str(iorb)//"dens_"//str(iorb),iorb=1,Norb),&
         (str(Norb+iorb)//"docc_"//str(iorb),iorb=1,Norb),&
         (str(2*Norb+iorb)//"nup_"//str(iorb),iorb=1,Norb),&
         (str(3*Norb+iorb)//"ndw_"//str(iorb),iorb=1,Norb),&
         str(4*Norb+1)//"s2tot",str(5*Norb+2)//"egs"
    close(unit)
    !
    !Magnetization along the 3 axis
    unit = free_unit()
    open(unit,file="magXYZ_info.ed")
    write(unit,"(A1,*(A10,6X))") "# ",(str(iorb)//"mag_"//str(iorb),iorb=1,Norb)
    close(unit)
    !
    !Spin-Spin correlation
    unit = free_unit()
    open(unit,file="Sz2_info.ed")
    write(unit,"(A1,2A6,A15)")"#","a","b","Sz.Sz(a,b)"
    close(unit)
    !
    !Density-Density correlation
    unit = free_unit()
    open(unit,file="N2_info.ed")
    write(unit,"(A1,2A6,A15)")"#","a","b","N.N(a,b)"
    close(unit)
    !
    !Exciton Singlet-Triplets
    unit = free_unit()
    open(unit,file="Exct_info.ed")
    write(unit,"(A1,*(A10,6X))") "# ",(("Exct_"//str(iorb)//str(jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
  end subroutine write_obs_info









  subroutine write_obs_last()
    integer :: unit,iorb,jorb,ispin
    !Parameters used:
    if(.not.ed_read_umatrix)then
      unit = free_unit()
      open(unit,file="parameters.ed")
      write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
      close(unit)
    endif
    !
    !Generic observables 
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         s2tot,egs
    close(unit)
    !
    !Magnetization along the 3 axis
    unit = free_unit()
    open(unit,file="magX_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))") (magx(iorb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="magY_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))") (magy(iorb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="magZ_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))") (magz(iorb),iorb=1,Norb)
    close(unit)
    !
    !Spin-Spin correlation
    unit = free_unit()
    open(unit,file="Sz2_last"//reg(ed_file_suffix)//".ed")
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,sz2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Density-Density correlation
    unit = free_unit()
    open(unit,file="N2_last"//reg(ed_file_suffix)//".ed")
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,n2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Exciton Singlet-Triplets
    unit = free_unit()
    open(unit,file="ExctS0_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))") ((exct_s0(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="ExctTx_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))") ((exct_tx(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="ExctTy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))") ((exct_ty(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="ExctTz_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))") ((exct_tz(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
  end subroutine write_obs_last









  subroutine write_obs_all()
    integer :: unit,iorb,jorb,ispin
    !
    !Generic observables 
    unit = free_unit()
    open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(*(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         s2tot,egs
    close(unit)
    !
    !
    !Magnetization along the 3 axis
    unit = free_unit()
    open(unit,file="magX_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))") (magx(iorb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="magY_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))") (magy(iorb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="magZ_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))") (magz(iorb),iorb=1,Norb)
    close(unit)
    !
    !Spin-Spin correlation
    unit = free_unit()
    open(unit,file="Sz2_all"//reg(ed_file_suffix)//".ed",position='append')
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,sz2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Density-Density correlation
    unit = free_unit()
    open(unit,file="N2_all"//reg(ed_file_suffix)//".ed",position='append')
    do iorb=1,Norb
       do jorb=1,Norb
          write(unit,"(1X,2I6,F15.9)")iorb,jorb,n2(iorb,jorb)
       enddo
    enddo
    close(unit)
    !
    !Exciton Singlet-Triplets
    unit = free_unit()
    open(unit,file="ExctS0_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(*(F15.9,1X))") ((exct_s0(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="ExctTx_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(*(F15.9,1X))") ((exct_tx(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="ExctTy_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(*(F15.9,1X))") ((exct_ty(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="ExctTz_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(*(F15.9,1X))") ((exct_tz(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
    close(unit)
  end subroutine write_obs_all





  subroutine write_energy()
    !Write the latest iteration values of energy observables
    integer :: unit
    !
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>",&
         reg(txtfy(7))//"<Dse>",&
         reg(txtfy(8))//"<Dph>"
    close(unit)
    !
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



end MODULE ED_OBSERVABLES_NONSU2
