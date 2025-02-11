MODULE ED_GF_NONSU2
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
  USE ED_HAMILTONIAN_NONSU2
  implicit none
  private


  public :: build_impG_nonsu2
  public :: get_impG_nonsu2
  public :: get_Sigma_nonsu2



  integer                               :: istate
  integer                               :: isector,jsector,ksector
  integer                               :: idim,jdim
  integer                               :: isz,jsz
  integer                               :: ialfa,jalfa
  integer                               :: i,j,m,r
  integer                               :: iph,i_el
  real(8)                               :: sgn,norm2
  complex(8)                            :: cnorm2

  complex(8),allocatable                :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)



  !Lanczos shared variables
  !=========================================================
  complex(8),dimension(:),allocatable   :: v_state
  real(8)                               :: e_state



contains




  !+------------------------------------------------------------------+
  !                            NONSU2
  !+------------------------------------------------------------------+
  subroutine build_impG_nonsu2()
    !
    !
    !Evaluates the impurity electrons Green's functions :math:`G(z)` using dynamical Lanczos method. The result is stored in rank-5 arrays :f:var:`impgmats`, :f:var:`impgreal` , :f:var:`impfmats` , :f:var:`impfreal` of dimensions [ |Nspin| , |Nspin| , |Norb| , |Norb| , :f:var:`Lmats` / :f:var:`Lreal` ]
    !
    !The off-diagonal components of :math:`G_{ab\sigma\sigma'}` with :math:`a \neq b` and :math:`\sigma,\sigma'=\uparrow, \downarrow` are obtained using algebraic manipulation to ensure working with hermitian conjugate operators in the dynamical Lanczos procedure.  
    !
    !The weights and the poles obtained in this procedure are saved in a hierarchical data structure (for every state, every channel (creation or annihilation of excitations, normal or anomalous) and every degree of freedom) :f:var:`impgmatrix` of type :f:var:`gfmatrix`. 
    !
    ! .. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261
    !
    integer                                     :: iorb,jorb,ispin,jspin,i,io,jo
    logical                                     :: MaskBool
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG build_gf NONSU2: build GFs"
#endif
    !    
    call PrintHmask
    !
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !same orbital, same spin GF: G_{aa}^{ss}(z)    
    do ispin=1,Nspin
       do iorb=1,Norb
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),Nstate=state_list%size)
          call lanc_build_gf_nonsu2_diagOrb_diagSpin(iorb,ispin)
       enddo
    enddo
    !
    !
    !same orbital, different spin GF: G_{aa}^{ss'}(z)
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             if(.not.Gbool(ispin,jspin,iorb,iorb))cycle
             call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,iorb),Nstate=state_list%size)
             call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,iorb,ispin,jspin)
          enddo
       enddo
    enddo
    !
    !
    if(bath_type=="normal")return
    !
    !
    !different orbital, same spin GF: G_{ab}^{ss}(z)
    do ispin=1,Nspin
       do iorb=1,Norb
          do jorb=1,Norb
             if(iorb==jorb)cycle
             if(.not.Gbool(ispin,ispin,iorb,jorb))cycle
             call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,jorb),Nstate=state_list%size)
             call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,ispin)
          enddo
       enddo
    enddo
    !
    !
    !different orbital, different spin GF: G_{ab}^{ss'}(z)
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                if(.not.Gbool(ispin,jspin,iorb,jorb))cycle
                call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),Nstate=state_list%size)
                call lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,jspin)
             enddo
          enddo
       enddo
    enddo
    !
  end subroutine build_impG_nonsu2













  !PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin(iorb,ispin)
    integer      :: iorb,ispin,isite
    type(sector) :: sectorI,sectorJ
    integer :: ib(2*Ns)
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,Nchan=2) !2*[c,cdg]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !EVALUATE c^+_{a,s}|v> --> G^>_{s,s;a,a}
       jsector = getCDGsector(1,ispin,isector)
       if(Jz_basis)jsector = getCDGsector_Jz(iorb,ispin,isector)
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          vvinit = apply_op_CDG(v_state,iorb,ispin,isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,1,iorb,iorb,ispin,ispin,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,1,Nexc=0)
       endif
       !
       !EVALUATE c_{a,s}|v> --> G^<_{s,s;a,a}
       jsector = getCsector(1,ispin,isector)
       if(Jz_basis)jsector = getCsector_Jz(iorb,ispin,isector)
       if(getN(isector)/=0.and.jsector>=0)then
          vvinit =  apply_op_C(v_state,iorb,ispin,isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,-1,iorb,iorb,ispin,ispin,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,ispin,iorb,iorb),istate,2,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin




  !PURPOSE: Evaluate the  different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin(iorb,jorb,ispin,jspin)
    integer      :: iorb,jorb,ispin,jspin
    type(sector) :: sectorI,sectorJ
    !
    write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(jspin)
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,Nchan=4) !2*[aux1,aux2]
       !
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_cvec(state_list,istate)
       !
       !
       !APPLY (c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs> = [1,1].[C_{+1},C_{+1}].[iorb,jorb].[ispin,jspin]
       jsector = getCDGsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
          ksector = getCDGsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,one],[1,1],[iorb,jorb],[ispin,jspin],isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,1,iorb,jorb,ispin,jspin,1,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !
       !APPLY (c_{jorb,jspin} + c_{iorb,ispin})|gs> =  [1,1].[C_{-1},C_{-1}].[iorb,jorb].[ispin,jspin]
       jsector = getCsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
          ksector = getCsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       if(getN(isector)/=0.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,one],[-1,-1],[iorb,jorb],[ispin,jspin],isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(one*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ispin,jspin,2,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !
       !APPLY (+i*c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs> =  [1,1j].[C_{+1},C_{+1}].[iorb,jorb].[ispin,jspin]
       jsector = getCDGsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
          ksector = getCDGsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       !
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,xi],[1,1],[iorb,jorb],[ispin,jspin],isector,jsector)          
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(-xi*norm2,e_state,alfa_,beta_,1,iorb,jorb,ispin,jspin,3,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !
       !APPLY (-xi*c_{jorb,jspin} + c_{iorb,ispin})|gs> =  [1,-1j].[C_{-1},C_{-1}].[iorb,jorb].[ispin,jspin]
       jsector = getCsector(1,ispin,isector)
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
          ksector = getCsector_Jz(jorb,jspin,isector)
          if(getDim(jsector)/=getDim(ksector))stop "lanczos builgf dimensional error"
       endif
       if(getN(isector)/=0.and.jsector>=0)then
          vvinit = apply_Cops(v_state,[one,-xi],[-1,-1],[iorb,jorb],[ispin,jspin],isector,jsector)
          call tridiag_Hv_sector_nonsu2(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_nonsu2(-xi*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ispin,jspin,4,istate)
          deallocate(alfa_,beta_,vvinit)
       else
          call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin






  subroutine add_to_lanczos_gf_nonsu2(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,jspin,ichan,istate)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin,jspin,ichan,istate
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_GF NONSU2: add-up to GF. istate:"//str(istate)
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
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
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
         "DEBUG add_to_lanczos_GF NONSU2: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(impGmatrix(ispin,jspin,iorb,jorb),istate,ichan,Nlanc)
    !    
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       !
       impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%weight(j) = peso
       impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%poles(j)  = isign*de
       !       
    enddo
  end subroutine add_to_lanczos_gf_nonsu2



  !################################################################
  !################################################################
  !################################################################
  !################################################################





  function get_impG_nonsu2(zeta,axis) result(Gf)
    !
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Gf
    integer                                                :: iorb,jorb,ispin,jspin,i
    character(len=1)                                       :: axis_
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG get_impG_nonsu2: Get GFs on a input array zeta"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    !
    if(.not.allocated(impGmatrix))stop "get_impG_nonsu2 ERROR: impGmatrix not allocated!"
    !
    call PrintHmask()
    !
    Gf = zero
    !
    !same orbital, same spin GF: G_{aa}^{ss}(z)
    do ispin=1,Nspin
       do iorb=1,Norb
          call get_nonsu2_Gab(iorb,iorb,ispin,ispin)
       enddo
    enddo
    !
    !same orbital, different spin GF: G_{aa}^{ss'}(z)
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             if(.not.Gbool(ispin,jspin,iorb,iorb))cycle
             call get_nonsu2_Gab(iorb,iorb,ispin,jspin)
          enddo
       enddo
    enddo
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          if(ispin==jspin)cycle
          do iorb=1,Norb
             Gf(ispin,jspin,iorb,iorb,:) = 0.5d0*(Gf(ispin,jspin,iorb,iorb,:) &
                  - (one-xi)*Gf(ispin,ispin,iorb,iorb,:) &
                  - (one-xi)*Gf(jspin,jspin,iorb,iorb,:))
          enddo
       enddo
    enddo
    !
    select case(bath_type)
    case ("normal");
    case default;
       !
       !different orbital, same spin GF: G_{ab}^{ss}(z)
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                if(.not.Gbool(ispin,ispin,iorb,jorb))cycle
                call get_nonsu2_Gab(iorb,jorb,ispin,ispin)
             enddo
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if(iorb==jorb)cycle
                Gf(ispin,ispin,iorb,jorb,:) = 0.5d0*(Gf(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*Gf(ispin,ispin,iorb,iorb,:) &
                     - (one-xi)*Gf(ispin,ispin,jorb,jorb,:))
             enddo
          enddo
       enddo
       !
       !different orbital, different spin GF: G_{ab}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   if(.not.Gbool(ispin,jspin,iorb,jorb))cycle
                   call get_nonsu2_Gab(iorb,jorb,ispin,jspin)
                enddo
             enddo
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             if(ispin==jspin)cycle
             do iorb=1,Norb
                do jorb=1,Norb
                   if(iorb==jorb)cycle
                   Gf(ispin,jspin,iorb,jorb,:) = 0.5d0*(Gf(ispin,jspin,iorb,jorb,:) &
                        - (one-xi)*Gf(ispin,ispin,iorb,iorb,:) &
                        - (one-xi)*Gf(jspin,jspin,jorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
    !
  contains
    !
    subroutine get_nonsu2_Gab(iorb,jorb,ispin,jspin)
      integer,intent(in)                 :: iorb,jorb,ispin,jspin
      integer                            :: Nstates,istate
      integer                            :: Nchannels,ichan
      integer                            :: Nexcs,iexc
      real(8)                            :: peso,de
      !
      Gf(ispin,jspin,iorb,jorb,:)=zero
      !
      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
      if(.not.allocated(impGmatrix(ispin,jspin,iorb,jorb)%state)) return
      !
      associate(G => Gf(ispin,jspin,iorb,jorb,:)) !just an alias 
        G= zero
        Nstates = size(impGmatrix(ispin,jspin,iorb,jorb)%state)
        do istate=1,Nstates
           if(.not.allocated(impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel))cycle
           Nchannels = size(impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel)
           do ichan=1,Nchannels
              Nexcs  = size(impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%poles)
              if(Nexcs==0)cycle
              do iexc=1,Nexcs
                 peso  = impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
                 de    = impGmatrix(ispin,jspin,iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
                 G     = G + peso/(zeta-de)
              enddo
           enddo
        enddo
      end associate
      return
    end subroutine get_nonsu2_Gab
    !
  end function get_impG_nonsu2





  function get_Sigma_nonsu2(zeta,axis) result(Sigma)
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    character(len=1)                                       :: axis_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: Sigma,invG0,invG
    complex(8),dimension(Nspin*Norb,Nspin*Norb)            :: iGzeta
    !
    axis_="m";if(present(axis))axis_=str(axis)
    !
    !Get G0^-1
    invG0 = invg0_bath_function(zeta,dmft_bath,axis_)
    !
    !Get G^-1
    invG  = get_impG_nonsu2(zeta)
    !
    !Get Sigma= G0^-1 - G^-1
    do i=1,size(zeta)     
       iGzeta =  nn2so_reshape(invG(:,:,:,:,i),Nspin,Norb)
       call inv(iGzeta)
       invG(:,:,:,:,i)=so2nn_reshape(iGzeta,Nspin,Norb)
       Sigma(:,:,:,:,1) = invG0(:,:,:,:,1) - invG(:,:,:,:,i)
    enddo
    !
  end function get_Sigma_nonsu2












  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine PrintHmask()
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                     :: iorb,jorb,ispin,jspin
    Hmask=.true.
    if(.not.ed_all_g)then
       select case(bath_type)
       case("replica");Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
       case("general");Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
       case default;stop "ERROR: PrintHmask=FALSE AND BATH_TYPE!=REPLICA/GENERAL"
       end select
    end if
    write(LOGfile,"(A)")"Get mask(G):"
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,*)((Hmask(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
  end subroutine PrintHmask


  function Gbool(ispin,jspin,iorb,jorb) result(bool)
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                     :: iorb,jorb,ispin,jspin
    logical                                     :: bool
    Hmask=.true.
    if(.not.ed_all_g)then
       select case(bath_type)
       case("replica");Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
       case("general");Hmask=Hgeneral_mask(wdiag=.true.,uplo=.false.)
       case default;stop "ERROR: PrintHmask=FALSE AND BATH_TYPE!=REPLICA/GENERAL"
       end select
    end if
    bool=.true.   
    if(bath_type=="replica".or.bath_type=="general")bool=Hmask(ispin,ispin,iorb,jorb)
  end function Gbool





END MODULE ED_GF_NONSU2





