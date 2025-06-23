MODULE ED_GF_SUPERC
  !:synopsis: Routines for Green's functions calculation, :code:`SUPERC` case
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER
  USE SF_IOTOOLS, only: str,reg,txtfy, save_array, to_lower
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private

  public :: build_impG_superc
  public :: get_impG_superc , get_impF_superc, get_impD_superc
  public :: get_Sigma_superc, get_Self_superc

  integer                               :: istate
  integer                               :: isector,jsector,ksector
  integer                               :: idim,jdim
  integer                               :: isz,jsz
  integer                               :: ialfa,jalfa
  integer                               :: i,j,m,r
  integer                               :: iph,i_el
  real(8)                               :: sgn,norm2
  !
  complex(8),allocatable                :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)

  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: v_state
  real(8)                               :: e_state

contains



  !+------------------------------------------------------------------+
  !                        SUPERC
  !+------------------------------------------------------------------+
  subroutine build_impG_superc()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    !
    !Evaluates the impurity electrons Green's functions :math:`G(z)` and :math:`F(z)` and the phonons one :math:`D(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` , :f:var:`impfmats` , :f:var:`impfreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ] and rank-1 array :f:var:`impdmats`, :f:var:`impdreal`.    
    !
    !The off-diagonal components of :math:`G_{ab}` with :math:`a \neq b` as well as the anomalous Green's functions :math:`F_{ab}(z)\, \forall a,b` are obtained using algebraic manipulation to ensure working with hermitian conjugate operators in the dynamical Lanczos procedure.  
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations, normal or anomalous) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    !
    integer    :: iorb,jorb,ispin,i,isign
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_gf SUPERC: build GFs"
#endif
    !
    ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
    !
    !1: Get G^{aa}_{upup}-> saved in impG(1,1,a,a) and \bar{G}^{aa}_{dwdw}-> saved in barG
    ! (ispin,ispin) to components:
    ! (1,1) -> diag normal (diag e offdiag)
    ! (1,2) -> anomalous (diag e offdiag)
    ! (2,2) -> bar       
    if(MPIMASTER)call start_timer(unit=LOGfile)    
    do iorb=1,Norb
       call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),Nstate=state_list%size)
       call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),Nstate=state_list%size)
       call lanc_build_gf_superc_Gdiag(iorb)
    enddo
    if(MPIMASTER)call stop_timer
    !
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    select case(bath_type)
    case default
       !get F^{aa}_{updw} = <adg_up . adg_dw>
       do iorb=1,Norb
          call allocate_GFmatrix(impGmatrix(1,2,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_superc_Fmix(iorb,iorb)
       enddo
       !
    case ('hybrid','replica','general')
       !Get G^{ab}_{upup} --> saved in impG(1,1,a,b)
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_Gmix(iorb,jorb)
          enddo
       enddo
       !
       !get F^{ab}_{updw} = <adg_up . bdg_dw>
       do iorb=1,Norb
          do jorb=1,Norb
             call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_superc_Fmix(iorb,jorb)
          enddo
       enddo
    end select
    !
    !PHONONS
    if(DimPh>1)then
       call allocate_GFmatrix(impDmatrix,Nstate=state_list%size)
       call lanc_build_gf_phonon_main()
    endif
    !
    if(MPIMASTER)call stop_timer
    !
  end subroutine build_impG_superc










  subroutine lanc_build_gf_superc_Gdiag(iorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer      :: iorb
    !
    write(LOGfile,"(A)")"Get G & barG_l"//str(iorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,Nchan=2) ![c+_up,c_up]
       call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),istate,Nchan=2) ![c_dw,c+_dw]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE c^+_{up,iorb}|v> --> Gaux(1) = G^>_{up,up;iorb,iorb}
       jsector = getCDGsector(1,1,isector)
       if(jsector/=0)then 
          vvinit = apply_op_CDG(v_state,iorb,1,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,1,iorb,iorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c_{up,iorb}|v> --> Gaux(1) = G^<_{up,up;iorb,iorb}
       jsector = getCsector(1,1,isector)
       if(jsector/=0)then
          vvinit = apply_op_C(v_state,iorb,1,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,1,iorb,iorb,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE c_{dw,iorb}|v> --> Gaux(2) = G^>_{dw,dw;iorb,iorb}
       jsector = getCsector(1,2,isector)
       if(jsector/=0)then
          vvinit = apply_op_C(v_state,iorb,2,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,2,iorb,iorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c^+_{dw,iorb}|v> --> Gaux(2) = G^<_{dw,dw;iorb,iorb}
       jsector = getCDGsector(1,2,isector)
       if(jsector/=0)then
          vvinit = apply_op_CDG(v_state,iorb,2,isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,2,iorb,iorb,2,istate)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       else
          call allocate_GFmatrix(impGmatrix(2,2,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_Gdiag










  subroutine lanc_build_gf_superc_Gmix(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer      :: iorb,jorb
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,Nchan=4) !2*[QdgQ,RdgR]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE (cdg_iorb + cdg_jorb)|gs> = Qdg|gs> += <Q.Qdg>
       ! =[1,1].[C_{+1},C_{+1}].[iorb,jorb].[up,up]
       jsector = getCDGsector(1,1,isector)
       if(jsector/=0)then
          vvinit = apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,3,iorb,jorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE (c_iorb + c_jorb)|gs> = Q|gs> += <Qdg.Q>
       !=[1,1].[C_{-1},C_{-1}].[iorb,jorb].[up,up]
       jsector = getCsector(1,1,isector)
       if(jsector/=0)then
          vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,3,iorb,jorb,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !
       !EVALUATE (cdg_iorb + xi cdg_jorb)|gs> = Rdg|gs> += -xi*<Rdg.R>
       !=[1,1j].[C_{+1},C_{+1}].[iorb,jorb].[up,up]
       jsector = getCDGsector(1,1,isector)
       if(jsector/=0)then
          vvinit =  apply_Cops(v_state,[one,xi],[1,1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,3,iorb,jorb,3,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE (c_iorb - xi c_jorb)|gs>  = R|gs> += -xi*<R.Rdg>
       !=[1,-1j].[C_{-1},C_{-1}].[iorb,jorb].[up,up]
       jsector = getCsector(1,1,isector)
       if(jsector/=0)then
          vvinit =  apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[1,1],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,3,iorb,jorb,4,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,1,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_superc_Gmix







  subroutine lanc_build_gf_superc_Fmix(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                :: iorb,jorb
    type(sector)           :: sectorI,sectorJ
    !
    write(LOGfile,"(A)")"Get F_l"//str(iorb)//"_m"//str(jorb)
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,Nchan=4) !2*[Odg.O,Pdg.P]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE [cdg_{up,iorb} + c_{dw,jorb}]|gs> = Odg|gs> += <O.Odg>
       ! =[1,1].[C_{+1},C_{-1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit = apply_Cops(v_state,[one,one],[1,-1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,1,4,iorb,jorb,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} + cdg_{dw,jorb}]|gs> = O|gs> += <Odg.O>
       ! =[1,1].[C_{-1},C_{+1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit = apply_Cops(v_state,[one,one],[-1,1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(one*norm2,e_state,alfa_,beta_,-1,4,iorb,jorb,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !EVALUATE [c^+_{up,iorb} + xi*c_{dw,jorb}]|gs> = Pdg|gs> += -xi*<P.Pdg>
       ! =[1,1j].[C_{+1},C_{-1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz<Ns)then
          jsector = getsector(isz+1,1)
          vvinit = apply_Cops(v_state,[one,xi],[1,-1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,1,4,iorb,jorb,3,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !EVALUATE [c_{up,iorb} - xi*c^+_{dw,jorb}]|gs> = P|gs> += -xi*<Pdg.P>
       ! =[1,-1j].[C_{-1},C_{+1}].[iorb,jorb].[up,dw]
       isz = getsz(isector)
       if(isz>-Ns)then
          jsector = getsector(isz-1,1)
          vvinit = apply_Cops(v_state,[one,-xi],[-1,+1],[iorb,jorb],[1,2],isector,jsector)
          call tridiag_Hv_sector_superc(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_superc(-xi*norm2,e_state,alfa_,beta_,-1,4,iorb,jorb,4,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(1,2,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !
    return
  end subroutine lanc_build_gf_superc_Fmix




  subroutine lanc_build_gf_phonon_main()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    type(sector)                :: sectorI
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG lanc_build_gf_phonon: build phonon GF"
#endif
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impDmatrix,istate=istate,Nchan=1)
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,I6)")&
               'From sector  :',isector,sectorI%Sz
       endif
       !
       if(MpiMaster)then
          if(ed_verbose>=3)write(LOGfile,"(A20,I12)")'Apply x',isector
          !
          allocate(vvinit(sectorI%Dim));vvinit=0d0
          !
          do i=1,sectorI%Dim
             iph = (i-1)/(sectorI%DimEl) + 1
             i_el = mod(i-1,sectorI%DimEl) + 1
             !
             !apply destruction operator
             if(iph>1) then
                j = i_el + ((iph-1)-1)*sectorI%DimEl
                vvinit(j) = vvinit(j) + sqrt(dble(iph-1))*v_state(i)
             endif
             !
             !apply creation operator
             if(iph<DimPh) then
                j = i_el + ((iph+1)-1)*sectorI%DimEl
                vvinit(j) = vvinit(j) + sqrt(dble(iph))*v_state(i)
             endif
          enddo
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector_superc(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_phonon(norm2,e_state,alfa_,beta_,istate)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       if(allocated(v_state))deallocate(v_state)
    end do
    return
  end subroutine lanc_build_gf_phonon_main





  !################################################################
  !################################################################
  !################################################################
  !################################################################





  subroutine add_to_lanczos_gf_superc(vnorm2,Ei,alanc,blanc,isign,ichan,iorb,jorb,ic,istate)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,ichan,iorb,jorb,istate,ic
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr,chanI,chanJ
    complex(8)                                 :: iw
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF SUPERC: add-up to GF. istate:"//str(istate)
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if((finiteT).and.(beta*(Ei-Egs).lt.200))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !
    ! pesoBZ = vnorm2/zeta_function
    ! if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
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
         "DEBUG add_to_lanczos_GF SUPERC: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    select case(ichan)
    case(1,3)
       chanI=1
       chanJ=1
    case(2)
       chanI=2
       chanJ=2
    case(4)
       chanI=1
       chanJ=2
    case default
       STOP "add_to_lanczos_gf_superc: wrong ichan"
    end select
    !
    call allocate_GFmatrix(impGmatrix(chanI,chanJ,iorb,jorb),istate,ic,Nlanc)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(chanI,chanJ,iorb,jorb)%state(istate)%channel(ic)%weight(j) = peso
       impGmatrix(chanI,chanJ,iorb,jorb)%state(istate)%channel(ic)%poles(j)  = isign*de
    enddo
  end subroutine add_to_lanczos_gf_superc





  subroutine add_to_lanczos_phonon(vnorm2,Ei,alanc,blanc,istate)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
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
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
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
    call allocate_GFmatrix(impDmatrix,istate,1,Nexc=Nlanc)
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







  function get_impG_superc(zeta,axis) result(Gf)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    ! Reconstructs the system impurity electrons normal Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Gf
    integer                                                :: iorb,jorb,ispin,jspin,i
    character(len=1)                                       :: axis_
    complex(8)                                             :: auxG(4,size(zeta))
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impG_superc"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(impGmatrix))stop "get_impG_superc ERROR: impGmatrix not allocated!"
    !
    auxG = zero
    Gf   = zero
    !
    ispin=1 ! in this channel Nspin=2 is forbidden. check in ED_SETUP.
    !
    do iorb=1,Norb
       call get_superc_Gdiag(iorb)             !-> auxG(1:,:)
       Gf(1,1,iorb,iorb,:) = auxG(1,:)         !this is G_{up,up;iorb,iorb}
    enddo
    !
    select case(bath_type)
    case ("normal")
    case default
       ! Get G^{ab}_{upup} = <adg_up . b_up>
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             call get_superc_Gmix(iorb,jorb)     !-> auxG(3,:)
             Gf(1,1,iorb,jorb,:) = 0.5d0*( auxG(3,:) - (one-xi)*(Gf(1,1,iorb,iorb,:)+Gf(1,1,jorb,jorb,:)) )
          enddo
       enddo
    end select

  contains
    !
    subroutine get_superc_Gdiag(iorb) !get auxG(1:2)
#if __INTEL_COMPILER
      use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in)                 :: iorb
      integer                            :: Nstates,istate
      integer                            :: Nchannels,ic,ichan
      integer                            :: Nexcs,iexc
      real(8)                            :: de
      complex(8)                         :: peso
      !
#ifdef _DEBUG
      if(ed_verbose>2)write(LOGfile,"(A)")"DEBUG Get G_l"//str(iorb)//"_m"//str(iorb)//"_axis: "//str(axis_)
#endif
      if(.not.allocated(impGmatrix(1,1,iorb,iorb)%state))return
      !
      auxG(1,:) = zero
      Nstates = size(impGmatrix(1,1,iorb,iorb)%state)
      do istate=1,Nstates
         Nchannels = size(impGmatrix(1,1,iorb,iorb)%state(istate)%channel)     
         do ic=1,Nchannels        
            Nexcs  = size(impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso  = impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%weight(iexc)
               de    = impGmatrix(1,1,iorb,iorb)%state(istate)%channel(ic)%poles(iexc)
               auxG(1,:)     = auxG(1,:) + peso/(zeta-de)
            enddo
         enddo
      enddo

      return
    end subroutine get_superc_Gdiag

    subroutine get_superc_Gmix(iorb,jorb) !get auxG(3,:)
#if __INTEL_COMPILER
      use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in) :: iorb,jorb
      integer            :: Nstates,istate
      integer            :: Nchannels,ic,ichan
      integer            :: Nexcs,iexc
      real(8)            :: de
      complex(8)         :: peso
      !
#ifdef _DEBUG
      if(ed_verbose>2)write(LOGfile,"(A)")"DEBUG Get G_l"//str(iorb)//"_m"//str(jorb)//"_axis: "//str(axis_)
#endif
      if(.not.allocated(impGmatrix(1,1,iorb,jorb)%state)) return
      !
      auxG(3,:) = zero
      Nstates = size(impGmatrix(1,1,iorb,jorb)%state)
      do istate=1,Nstates
         Nchannels = size(impGmatrix(1,1,iorb,jorb)%state(istate)%channel)     
         do ic=1,Nchannels
            Nexcs  = size(impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso = impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%weight(iexc)
               de   = impGmatrix(1,1,iorb,jorb)%state(istate)%channel(ic)%poles(iexc)
               auxG(3,:)    = auxG(3,:) + peso/(zeta-de)
            enddo
         enddo
      enddo
      return
    end subroutine get_superc_Gmix
  end function get_impG_superc









  function get_impF_superc(zeta,axis) result(Ff)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    ! Reconstructs the system impurity anomalous electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Ff
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Gf
    integer                                                :: iorb,jorb,ispin,jspin,i
    character(len=1)                                       :: axis_
    complex(8)                                             :: barG(Norb,size(zeta))
    complex(8)                                             :: auxG(4,size(zeta))
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impF_superc"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(impGmatrix))stop "get_impF_superc ERROR: impGmatrix not allocated!"
    !
    auxG = zero
    barG = zero
    Gf   = zero
    Ff   = zero
    !
    ispin=1 ! in this channel Nspin=2 is forbidden. check in ED_SETUP.
    !
    do iorb=1,Norb
       call get_superc_Gdiag(iorb)               !-> auxG(1:2,:)
       Gf(1,1,iorb,iorb,:) = auxG(1,:)         !this is G_{up,up;iorb,iorb}
       barG(       iorb,:) = auxG(2,:) !this is G_{dw,dw;iorb,iorb}
    enddo
    !
    select case(bath_type)
    case ("normal")
       ! Get F^{aa}_{updw} = <adg_up . adg_dw>
       do iorb=1,Norb
          call get_superc_Fmix(iorb,iorb) !-> auxG(4,:)
          Ff(1,1,iorb,iorb,:) = 0.5d0*(auxG(4,:)-&
               (one-xi)*(Gf(1,1,iorb,iorb,:)+barG(iorb,:)))
       enddo
       !
    case default
       !
       ! Get F^{ab}_{updw} = <adg_up . bdg_dw>
       do iorb=1,Norb
          do jorb=1,Norb
             call get_superc_Fmix(iorb,jorb)
             Ff(1,1,iorb,jorb,:) = 0.5d0*(auxG(4,:)-&
                  (one-xi)*(Gf(1,1,iorb,iorb,:)+barG(jorb,:)))
          enddo
       enddo
    end select

  contains
    !
    subroutine get_superc_Gdiag(iorb) !get auxG(1:2)
#if __INTEL_COMPILER
      use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in) :: iorb
      integer            :: Nstates,istate
      integer            :: Nchannels,ic,is
      integer            :: Nexcs,iexc
      real(8)            :: de
      complex(8)         :: peso
      !
      auxG(1:2,:)=zero
      !
#ifdef _DEBUG
      if(ed_verbose>2)write(LOGfile,"(A)")"DEBUG Get G_l"//str(iorb)//"_m"//str(iorb)//"_axis: "//str(axis_)
#endif
      if(.not.allocated(impGmatrix(1,1,iorb,iorb)%state))return
      if(.not.allocated(impGmatrix(2,2,iorb,iorb)%state))return
      !
      do is=1,2
         Nstates = size(impGmatrix(is,is,iorb,iorb)%state)
         do istate=1,Nstates
            Nchannels = size(impGmatrix(is,is,iorb,iorb)%state(istate)%channel)     
            do ic=1,Nchannels        
               Nexcs  = size(impGmatrix(is,is,iorb,iorb)%state(istate)%channel(ic)%poles)
               if(Nexcs==0)cycle
               do iexc=1,Nexcs
                  peso  = impGmatrix(is,is,iorb,iorb)%state(istate)%channel(ic)%weight(iexc)
                  de    = impGmatrix(is,is,iorb,iorb)%state(istate)%channel(ic)%poles(iexc)
                  auxG(is,:)     = auxG(is,:) + peso/(zeta-de)
               enddo
            enddo
         enddo
      enddo
      !
      return
    end subroutine get_superc_Gdiag
    !
    subroutine  get_superc_Fmix(iorb,jorb)
#if __INTEL_COMPILER
      use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in) :: iorb,jorb
      integer            :: Nstates,istate
      integer            :: Nchannels,ic,ichan
      integer            :: Nexcs,iexc
      real(8)            :: de
      complex(8)         :: peso
      !
      auxG(4,:)=zero
      !
#ifdef _DEBUG
      if(ed_verbose>2)write(LOGfile,"(A)")"DEBUG Get F_l"//str(iorb)//"_m"//str(jorb)//"_axis: "//str(axis_)
#endif
      if(.not.allocated(impGmatrix(1,1,iorb,jorb)%state)) return
      !
      Nstates = size(impGmatrix(1,2,iorb,jorb)%state)
      do istate=1,Nstates
         Nchannels = size(impGmatrix(1,2,iorb,jorb)%state(istate)%channel)     
         do ic=1,Nchannels        
            Nexcs  = size(impGmatrix(1,2,iorb,jorb)%state(istate)%channel(ic)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso  = impGmatrix(1,2,iorb,jorb)%state(istate)%channel(ic)%weight(iexc)
               de    = impGmatrix(1,2,iorb,jorb)%state(istate)%channel(ic)%poles(iexc)
               auxG(4,:)     = auxG(4,:) + peso/(zeta-de)
            enddo
         enddo
      enddo
      return
    end subroutine get_superc_Fmix
    !
  end function get_impF_superc







  function get_impD_superc(zeta,axis) result(G)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    ! Reconstructs the phonon Green's functions using :f:var:`impdmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in) :: zeta
    character(len=*),optional          :: axis
    complex(8),dimension(size(zeta))   :: G
    character(len=1)                   :: axis_
    !
    integer                            :: Nstates,istate
    integer                            :: Nchannels,ichan,i
    integer                            :: Nexcs,iexc
    real(8)                            :: de
    complex(8)                         :: peso
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impD_superc"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    G = zero
    !
    write(LOGfile,"(A)")"Get D"
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
             ! ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
             ! ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
             ! ! This ensures that the correct null contribution is obtained.
             ! ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
             ! ! we do not include the contribution (because we are in the situation described above).
             ! ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
             select case(axis_)
             case("m","M")                
                if(beta*dE > 1d-3)G(1)=G(1) - peso*2*(1d0-exp(-beta*dE))/dE 
                do i=2,size(zeta)
                   G(i)=G(i) - peso*(1d0-exp(-beta*dE))*2d0*dE/(dreal(zeta(i))**2+dE**2)
                enddo
             case("r","R")
                do i=1,size(zeta)
                   G(i)=G(i) + peso*(1d0-exp(-beta*dE))*(1d0/(zeta(i) - dE) - 1d0/(zeta(i) + dE))
                enddo
             end select
          enddo
       enddo
    enddo
    !
  end function get_impD_superc








  function get_Sigma_superc(zeta,axis) result(Sigma)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Sigma
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: G,F
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: invG,invF,invG0,invF0
    real(8)                                                :: Gdet(size(zeta))
    complex(8)                                             :: detG(size(zeta))
    complex(8),dimension(2*Norb,2*Norb)                    :: M
    character(len=1)                                       :: axis_
    integer                                                :: L,ispin,iorb
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_Sigma_superc"
#endif
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    L = size(zeta)
    !
    !Get G0^-1,F0^-1
    ispin=1
    !
    invG0 = invg0_bath_function(zeta,axis_)
    invF0 = invf0_bath_function(zeta,axis_)
    !
    !Get G, F
    G     = get_impG_superc(zeta,axis_)
    F     = get_impF_superc(zeta,axis_)
    !
    !get G^{-1},F^{-1} --> Sigma
    Sigma = zero
    select case(bath_type)
    case ("normal")
       do iorb=1,Norb
          !
          select case(axis_)
          case default;stop "get_Sigma_superc error: axis_ != mats,real"
          case("m")
             Gdet =  dreal(abs(G(ispin,ispin,iorb,iorb,:))**2 + F(ispin,ispin,iorb,iorb,:)**2)
             invG(ispin,ispin,iorb,iorb,:)  =  conjg(G(ispin,ispin,iorb,iorb,:))/Gdet
          case("r")
             detG = -G(ispin,ispin,iorb,iorb,:)*conjg(G(ispin,ispin,iorb,iorb,L:1:-1)) - F(ispin,ispin,iorb,iorb,:)**2
             invG(ispin,ispin,iorb,iorb,:)  = -conjg(G(ispin,ispin,iorb,iorb,L:1:-1))/detG
          end select
          !
          Sigma(ispin,ispin,iorb,iorb,:) = invG0(ispin,ispin,iorb,iorb,:) - invG(ispin,ispin,iorb,iorb,:)
          !
       end do
       !
       !
    case ("hybrid","replica","general")
       do i=1,L
          M = zero
          select case(axis_)
          case default;stop "get_Sigma_superc error: axis_ != mats,real"
          case("m")
             M(1     :Norb  ,     1:Norb  ) = G(ispin,ispin,:,:,i)
             M(1     :Norb  ,Norb+1:2*Norb) = F(ispin,ispin,:,:,i)
             M(Norb+1:2*Norb,     1:Norb  ) = conjg(F(ispin,ispin,:,:,i)) !this is real so conjg does none, but it shouldn't be there
             M(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(G(ispin,ispin,:,:,i))
          case("r")
             M(1     :Norb  ,     1:Norb  ) = G(ispin,ispin,:,:,i)
             M(1     :Norb  ,Norb+1:2*Norb) = F(ispin,ispin,:,:,i)
             M(Norb+1:2*Norb,     1:Norb  ) = F(ispin,ispin,:,:,i)
             M(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(G(ispin,ispin,:,:,L-i+1))
          end select
          !
          call inv(M)
          invG(ispin,ispin,:,:,i) = M(1:Norb,1:Norb)
       enddo
       !
       Sigma(ispin,ispin,:,:,:)  = invG0(ispin,ispin,:,:,:) - invG(ispin,ispin,:,:,:)       
       !
    end select
    !
  end function get_Sigma_superc






  function get_Self_superc(zeta,axis) result(Self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Self
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: G,F
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: invG,invF,invG0,invF0
    complex(8)                                             :: detG(size(zeta))
    complex(8),dimension(2*Norb,2*Norb)                    :: M
    character(len=1)                                       :: axis_
    integer                                                :: L,ispin,iorb
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_Self_superc"
#endif
    !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    L = size(zeta)
    !
    !Get G0^-1,F0^-1
    ispin=1
    !
    invG0 = invg0_bath_function(zeta,axis_)
    invF0 = invf0_bath_function(zeta,axis_)
    !
    !Get G, F
    G     = get_impG_superc(zeta)
    F     = get_impF_superc(zeta)
    !
    !get G^{-1},F^{-1} --> Self
    Self = zero
    select case(bath_type)
    case ("normal")
       do iorb=1,Norb
          !
          select case(axis_)
          case default;stop "get_Sigma_superc error: axis_ != mats,real"
          case("m")
             detG =  dreal(abs(G(ispin,ispin,iorb,iorb,:))**2 + F(ispin,ispin,iorb,iorb,:)**2)
             invF(ispin,ispin,iorb,iorb,:)  =  F(ispin,ispin,iorb,iorb,:)/detG
          case("r")
             detG = -G(ispin,ispin,iorb,iorb,:)*conjg(G(ispin,ispin,iorb,iorb,L:1:-1)) - F(ispin,ispin,iorb,iorb,:)**2
             invF(ispin,ispin,iorb,iorb,:)  = -F(ispin,ispin,iorb,iorb,:)/detG
          end select
          !
          Self(ispin,ispin,iorb,iorb,:) = invF0(ispin,ispin,iorb,iorb,:) - invF(ispin,ispin,iorb,iorb,:)
       end do
       !
       !
    case ("hybrid","replica","general")
       do i=1,L
          M = zero
          select case(axis_)
          case default;stop "get_Sigma_superc error: axis_ != mats,real"
          case("m")
             M(1     :Norb  ,     1:Norb  ) = G(ispin,ispin,:,:,i)
             M(1     :Norb  ,Norb+1:2*Norb) = F(ispin,ispin,:,:,i)
             M(Norb+1:2*Norb,     1:Norb  ) = conjg(F(ispin,ispin,:,:,i)) !this is real so conjg does none, but it shouldn't be there
             M(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(G(ispin,ispin,:,:,i))
          case("r")
             M(1     :Norb  ,     1:Norb  ) = G(ispin,ispin,:,:,i)
             M(1     :Norb  ,Norb+1:2*Norb) = F(ispin,ispin,:,:,i)
             M(Norb+1:2*Norb,     1:Norb  ) = F(ispin,ispin,:,:,i)
             M(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(G(ispin,ispin,:,:,L-i+1))
          end select
          !
          call inv(M)
          invF(ispin,ispin,:,:,i) = M(1:Norb,Norb+1:2*Norb)
       enddo
       !
       Self(ispin,ispin,:,:,:)  = invF0(ispin,ispin,:,:,:) - invF(ispin,ispin,:,:,:)
    end select
    !
  end function get_Self_superc




















END MODULE ED_GF_SUPERC













