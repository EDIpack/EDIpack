MODULE ED_GF_NORMAL
  !:synopsis: Routines for Green's function calculation, :code:`NORMAL` case
  !Constructs the interacting impurity electronic Green's functions. On request it evaluated the phononic Green's functions and and different susceptibilities.  
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy,to_lower
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private


  public :: build_impG_normal
  public :: get_impG_normal
  public :: get_impD_normal
  public :: get_Sigma_normal



  integer                             :: istate
  integer                             :: isector,jsector
  integer                             :: idim,idimUP,idimDW
  !
  integer                             :: ialfa,jalfa
  integer                             :: ipos,jpos
  integer                             :: i,j,m
  integer                             :: iph,i_el
  real(8)                             :: sgn,norm2
  integer                             :: in,jn
  integer                             :: inam,jnam
  integer                             :: ilat,jlat
  integer                             :: iorb,jorb
  integer                             :: ispin,jspin
  integer                             :: is,js
  integer                             :: io,jo
  integer                             :: iup,idw
  integer                             :: jup,jdw  
  integer                             :: mup,mdw

  !
#ifdef _CMPLX_NORMAL
  complex(8),allocatable              :: vvinit(:)
  complex(8),dimension(:),allocatable :: v_state
#else
  real(8),allocatable                 :: vvinit(:)
  real(8),dimension(:),allocatable    :: v_state
#endif
  real(8),allocatable                 :: alfa_(:),beta_(:)
  real(8)                             :: e_state



contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+

  subroutine build_impG_normal()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    !Evaluates the impurity electrons Green's function :math:`G(z)` and the phonons one :math:`D(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ] and rank-1 array :f:var:`impdmats`, :f:var:`impdreal`.    
    !
    !The off-diagonal components of :math:`G_{ab}` with :math:`a \neq b` are obtained using algebraic manipulation, see `j.cpc.2021.108261`_. 
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    integer                                     :: iorb,jorb,ispin,jspin,i
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_gf NORMAL: build GFs"
#endif
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_normal_diag(iorb,ispin)
       enddo
    enddo
    !
    if(offdiag_gf_flag)then
       call PrintHmask()
       !
       do ispin=1,Nspin
          do iorb=1,Norb
#ifdef _CMPLX_NORMAL
             do jorb=1,Norb
#else
             do jorb=iorb+1,Norb
#endif
                if(iorb==jorb)cycle
                if(.not.Gbool(ispin,ispin,iorb,jorb))cycle
                call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),Nstate=state_list%size)
                call lanc_build_gf_normal_mix(iorb,jorb,ispin)
             enddo
          enddo
       enddo
       !
       !> PHONONS
       if(DimPh>1)then
          call allocate_GFmatrix(impDmatrix,Nstate=state_list%size)
          call lanc_build_gf_phonon_main()
       endif
       !
    end if
    if(MPIMASTER)call stop_timer
    !
  end subroutine build_impG_normal







  subroutine lanc_build_gf_normal_diag(iorb,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer,intent(in)          :: iorb,ispin
    !
    if(ed_verbose>1)write(LOGfile,*)"Get G_l"//str(iorb,3)//"_m"//str(iorb,3)//"_s"//str(ispin)
    !
    ialfa=1 ; if(.not.ed_total_ud)ialfa = iorb
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,Nchan=2)
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
#ifdef _CMPLX_NORMAL
       v_state    =  es_return_cvec(state_list,istate)
#else
       v_state    =  es_return_dvec(state_list,istate)
#endif
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then 
          vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,1,iorb,iorb,ispin,ichan=1,istate=istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          vvinit =  apply_op_C(v_state,iorb,ispin,isector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,-1,iorb,iorb,ispin,ichan=2,istate=istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_build_gf_normal_diag

  




  subroutine lanc_build_gf_normal_mix(iorb,jorb,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                     :: iorb,jorb,ispin
    !
    if(ed_verbose>1)write(LOGfile,*)"Get G_l"//str(iorb,3)//"_m"//str(jorb,3)//"_s"//str(ispin)
    !
    ialfa=1 ; if(.not.ed_total_ud)ialfa = iorb
    !
    do istate=1,state_list%size
       !
#ifdef _CMPLX_NORMAL
       call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,Nchan=4)
#else
       call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,Nchan=2)
#endif
       !
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
#ifdef _CMPLX_NORMAL
       v_state    =  es_return_cvec(state_list,istate)
#else
       v_state    =  es_return_dvec(state_list,istate)
#endif
       !
       !
       !
       !(c^+_iorb + c^+_jorb)|gs> = [1,1].[C_{+1},C_{+1}].[iorb,jorb].[ispin,ispin]
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then
#ifdef _CMPLX_NORMAL
          vvinit = apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[ispin,ispin],isector,jsector)
#else        
          vvinit = apply_Cops(v_state,[1d0,1d0],[1,1],[iorb,jorb],[ispin,ispin],isector,jsector)
#endif
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,1,iorb,jorb,ispin,ichan=1,istate=istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,ichan=1,Nexc=0)
       endif
       !
       !
       !
       !(c_iorb + c_jorb)|gs> = [1,1].[C_{-1},C_{-1}].[iorb,jorb].[ispin,ispin]
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
#ifdef _CMPLX_NORMAL
          vvinit =  apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[ispin,ispin],isector,jsector)
#else  
          vvinit =  apply_Cops(v_state,[1d0,1d0],[-1,-1],[iorb,jorb],[ispin,ispin],isector,jsector)
#endif
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ispin,ichan=2,istate=istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,ichan=2,Nexc=0)
       endif
       !
       !
       !
#ifdef _CMPLX_NORMAL
       !(c^+_iorb + c^+_jorb)|gs> = [1,1].[C_{+1},C_{+1}].[iorb,jorb].[ispin,ispin]
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then
          vvinit = apply_Cops(v_state,[one,xi],[1,1],[iorb,jorb],[ispin,ispin],isector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(-xi*norm2,e_state,alfa_,beta_,1,iorb,jorb,ispin,ichan=3,istate=istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,ichan=3,Nexc=0)
       endif
       !
       !(c_iorb + c_jorb)|gs> = [1,1].[C_{-1},C_{-1}].[iorb,jorb].[ispin,ispin]
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          vvinit =  apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[ispin,ispin],isector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(-xi*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ispin,ichan=4,istate=istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,ichan=4,Nexc=0)
       endif
#endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix



  subroutine lanc_build_gf_phonon_main()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")&
         "DEBUG lanc_build_gf_phonon: build phonon GF"
#endif
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impDmatrix,istate=istate,Nchan=1)
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
#ifdef _CMPLX_NORMAL
       v_state    =  es_return_cvec(state_list,istate)
#else
       v_state    =  es_return_dvec(state_list,istate)
#endif
       !
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       if(MpiMaster.AND.ed_verbose>=3)write(LOGfile,"(A20,I6,20I4)")&
            'From sector',isector,Nups,Ndws
       !
       idim = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I12)")'Apply x',isector
          !
          allocate(vvinit(idim));vvinit=0d0
          !
          do i=1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             !apply destruction operator
             if(iph>1) then
                j = i_el + ((iph-1)-1)*iDimUp*iDimDw
                vvinit(j) = vvinit(j) + one*sqrt(dble(iph-1))*v_state(i)
             endif
             !
             !apply creation operator
             if(iph<DimPh) then
                j = i_el + ((iph+1)-1)*iDimUp*iDimDw
                vvinit(j) = vvinit(j) + one*sqrt(dble(iph))*v_state(i)
             endif
          enddo
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector_normal(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_phonon(one*norm2,e_state,alfa_,beta_,istate)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_build_gf_phonon_main







  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,ichan,istate)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8)                                 :: vnorm2,pesoBZ,peso,pesoF
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin,ichan,istate
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF NORMAL: add-up to GF. istate:"//str(istate)
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    !
    pesoF  = vnorm2/zeta_function
    pesoBZ = one
    if(finiteT)then
       if(beta*(Ei-Egs) < 200)then
          pesoBZ = one*exp(-beta*(Ei-Egs))
       else
          pesoBZ = zero
       endif
    endif
    !
    !Only the nodes in Mpi_Comm_Group did get the alanc,blanc.
    !However after delete_sectorHv MpiComm returns to be the global one
    !so we can safely Bcast the alanc,blanc (known only to the operative group)
    !to every nodes. The master is in charge of this (as a
    !participant of the operative group)
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF NORMAL: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),istate,ichan,Nlanc)
    !
    do j=1,nlanc
       de   = diag(j)-Ei
       peso = pesoF*pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%weight(j) = peso
       impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles(j)  = isign*de
    enddo
  end subroutine add_to_lanczos_gf_normal






  subroutine add_to_lanczos_phonon(vnorm2,Ei,alanc,blanc,istate)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8)                                 :: vnorm2,pesoF,pesoAB,pesoBZ,de,peso
    real(8)                                    :: Ei,Ej,Egs
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,istate
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function
    pesoBZ = one
    if(finiteT)then
       if(beta*(Ei-Egs) < 200)then
          pesoBZ = one*exp(-beta*(Ei-Egs))
       else
          pesoBZ = zero
       endif
    endif
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(impDmatrix,istate,1,Nexc=Nlanc) !ichan=1
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       !
       impDmatrix%state(istate)%channel(1)%weight(j) = peso
       impDmatrix%state(istate)%channel(1)%poles(j)  = de
    enddo
  end subroutine add_to_lanczos_phonon




  !################################################################
  !################################################################
  !################################################################
  !################################################################




  function get_impG_normal(zeta,axis) result(Gf)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in)                     :: zeta   !array of frequencies
    character(len=*),optional                              :: axis   !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Gf
    integer                                                :: iorb,jorb,ispin,jspin,i
    character(len=1)                                       :: axis_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impG_normal"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(impGmatrix))stop "get_Gimp_normal_array ERROR: impGmatrix not allocated!"
    !
    Gf = zero
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          call get_normal_Gab(iorb,iorb,ispin)
       enddo
    enddo
    !
    if(offdiag_gf_flag)then
       call PrintHmask()
       !
       do ispin=1,Nspin
          do iorb=1,Norb
#ifdef _CMPLX_NORMAL
             do jorb=1,Norb
#else
             do jorb=iorb+1,Norb
#endif
                if(iorb==jorb)cycle
                if(.not.Gbool(ispin,ispin,iorb,jorb))cycle
                call get_normal_Gab(iorb,jorb,ispin)
             enddo
          enddo
       enddo
#ifdef _CMPLX_NORMAL
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                if(.not.Gbool(ispin,ispin,iorb,jorb))cycle
                Gf(ispin,ispin,iorb,jorb,:) = 0.5d0*(Gf(ispin,ispin,iorb,jorb,:) &
                        - (one-xi)*Gf(ispin,ispin,iorb,iorb,:) &
                        - (one-xi)*Gf(ispin,ispin,jorb,jorb,:))
             enddo
          enddo
       enddo
#else
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                if(.not.Gbool(ispin,ispin,iorb,jorb))cycle
                Gf(ispin,ispin,iorb,jorb,:) = 0.5d0*(Gf(ispin,ispin,iorb,jorb,:) &
                     - Gf(ispin,ispin,iorb,iorb,:) - Gf(ispin,ispin,jorb,jorb,:))
                Gf(ispin,ispin,jorb,iorb,:) = Gf(ispin,ispin,iorb,jorb,:)
             enddo
          enddo
       enddo
#endif
    end if
    !
  contains
    !
    subroutine get_normal_Gab(iorb,jorb,ispin)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in)                 :: iorb,jorb,ispin
      integer                            :: Nstates,istate
      integer                            :: Nchannels,ichan
      integer                            :: Nexcs,iexc
      real(8)                            :: de
#ifdef _CMPLX_NORMAL
      complex(8)                         :: peso
#else
      real(8)                            :: peso
#endif
      !
      Gf(ispin,ispin,iorb,jorb,:)=zero
      !
#ifdef _DEBUG
      if(ed_verbose>2)write(LOGfile,"(A)")"DEBUG Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
#endif
      if(.not.allocated(impGmatrix(ispin,ispin,iorb,jorb)%state)) return
      !
      Nstates = size(impGmatrix(ispin,ispin,iorb,jorb)%state)
      do istate=1,Nstates
         if(.not.allocated(impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel))cycle
         Nchannels = size(impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel)
         do ichan=1,Nchannels
            Nexcs  = size(impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso  = impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
               de    = impGmatrix(ispin,ispin,iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
               Gf(ispin,ispin,iorb,jorb,:) = Gf(ispin,ispin,iorb,jorb,:) + peso/(zeta-de)
            enddo
         enddo
      enddo
      return
    end subroutine get_normal_Gab
    !
  end function get_impG_normal



  function get_impD_normal(zeta,axis) result(G)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    ! Reconstructs the phonon Green's functions using :f:var:`impdmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in) :: zeta !array of frequencies
    character(len=*),optional          :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(size(zeta))   :: G
    character(len=1)                   :: axis_
    !
    integer                            :: Nstates,istate
    integer                            :: Nchannels,ichan,i
    integer                            :: Nexcs,iexc
    real(8)                            :: de
#ifdef _CMPLX_NORMAL
    complex(8)                         :: peso
#else
    real(8)                            :: peso
#endif

#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impD_normal"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(impDmatrix%state)) return
    !
    G= zero
    Nstates = size(impDmatrix%state)
    do istate=1,Nstates
       if(.not.allocated(impDmatrix%state(istate)%channel))cycle
       Nchannels = size(impDmatrix%state(istate)%channel)
       do ichan=1,Nchannels
          Nexcs  = size(impDmatrix%state(istate)%channel(ichan)%poles)
          if(Nexcs==0)cycle
          do iexc=1,Nexcs
             peso  = impDmatrix%state(istate)%channel(ichan)%weight(iexc)
             de    = impDmatrix%state(istate)%channel(ichan)%poles(iexc)
             if(abs(beta*de) < 1e-8) then
                select case(axis_)
                   case("m","M")
                      do i=1,size(zeta)
                         if(abs(zeta(i))<1e-10) then ! \nu=0
                            G(i) = G(i) - peso*beta
                         endif
                      enddo
                   case("r","R")
                      ! A zero energy pole contributes only to D(Re(z) = 0) and
                      ! its contribution is taken to be the same as for D(\nu=0),
                      ! regardless of the imaginary shift in z
                      do i=1,size(zeta)
                         if(abs(dreal(zeta(i)))<1e-10) then
                            G(i) = G(i) + peso*beta
                         endif
                      enddo
                end select
                ! Ignore the negative energy poles and enforce the spectral symmetry by adding
                ! contributions from de and -de for each positive de
             elseif(de > 0) then
                select case(axis_)
                   case("m","M")
                      do i=1,size(zeta)
                         G(i) = G(i) - &
                            peso*(1d0-exp(-beta*de))*2d0*de/(dimag(zeta(i))**2 + de**2)
                      enddo
                   case("r","R")
                      do i=1,size(zeta)
                         G(i) = G(i) + &
                            peso*(1d0-exp(-beta*de)) * (1d0/(zeta(i)-de) - 1d0/(zeta(i)+de))
                      enddo
                end select
             endif
          enddo
       enddo
    enddo
    !
  end function get_impD_normal







  function get_Sigma_normal(zeta,axis) result(Sigma)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta      !array of frequencies
    character(len=*),optional                              :: axis      !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Sigma,invG0,invG
    complex(8),dimension(Norb,Norb)                        :: iGzeta
    character(len=1)                                       :: axis_
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_Sigma_normal"
#endif
    axis_="m";if(present(axis))axis_=str(axis)
    !
    !Get G0^-1
    invG0 = invg0_bath_function(zeta,axis_)
    !
    !Get G^-1
    invG  = get_impG_normal(zeta)
    !
    Sigma = zero
    !
    !Get Sigma= G0^-1 - G^-1
    do ispin=1,Nspin
       do i=1,size(zeta)
          select case(bath_type)
          case ("normal")
             do iorb=1,Norb
                iGzeta(iorb,iorb) = one/invG(ispin,ispin,iorb,iorb,i)
                invG(ispin,ispin,iorb,iorb,i)=iGzeta(iorb,iorb)
             enddo
          case default !Diagonal in spin
             iGzeta(:,:) = invG(ispin,ispin,:,:,i)
             call inv(iGzeta)
             invG(ispin,ispin,:,:,i)=iGzeta
          end select
       enddo
       Sigma(ispin,ispin,:,:,:) = invG0(ispin,ispin,:,:,:) - invG(ispin,ispin,:,:,:)
    enddo
    !
  end function get_Sigma_normal





  !################################################################
  !################################################################
  !################################################################
  !################################################################





  subroutine PrintHmask()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                     :: iorb,jorb,ispin,jspin
    Hmask= .true.
    if(.not.ed_all_g)then
       if(bath_type=="replica")Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
       if(bath_type=="general")Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
    endif
    write(LOGfile,"(A)")"Get mask(G):"
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,*)((Hmask(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
  end subroutine PrintHmask



  function Gbool(ispin,jspin,iorb,jorb) result(bool)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                     :: iorb,jorb,ispin,jspin
    logical                                     :: bool
    Hmask= .true.
    if(.not.ed_all_g)then
       if(bath_type=="replica")Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
       if(bath_type=="general")Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
    endif
    bool=.true.   
    if(bath_type=="replica".or.bath_type=="general")bool=Hmask(ispin,ispin,iorb,jorb)
  end function Gbool



END MODULE ED_GF_NORMAL


















