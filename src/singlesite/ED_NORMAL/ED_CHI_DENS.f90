MODULE ED_CHI_DENS
  !:synopsis: Routines for susceptibility calculation, :code:`NORMAL` case
  !Evaluates the impurity density susceptibility.
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
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


  public :: build_densChi_normal
  public :: get_densChi_normal

  integer                          :: istate,iorb,jorb,ispin,jspin
  integer                          :: isector
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ialfa
  integer                          :: jalfa
  integer                          :: ipos,jpos
  integer                          :: i,j
  integer                          :: iph,i_el
  real(8)                          :: sgn,norm2
  real(8),dimension(:),allocatable :: v_state
  real(8)                          :: e_state




contains


  !+------------------------------------------------------------------+
  !                            DENS
  !PURPOSE  : Evaluate the Dens susceptibility \Chi_dens for a 
  ! \chi_ab = <n_a(\tau)n_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_densChi_normal()
    ! Evaluates the impurity density susceptibility :math:`\chi^{n}=\langle T_\tau n_a(\tau) n_b\rangle` in the Matsubara :math:`i\omega_n` and Real :math:`\omega` frequency axis as well as imaginary time :math:`\tau`. 
    ! As for the Green's function, the off-diagonal component of the the susceptibility is determined using an algebraic manipulation to ensure use of Hermitian operator in the dynamical Lanczos.
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    write(LOGfile,"(A)")"Get dens Chi:"
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !
    do iorb=1,Norb
       call allocate_GFmatrix(densChimatrix(iorb,iorb),Nstate=state_list%size)
       call lanc_ed_build_densChi_diag(iorb)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             call allocate_GFmatrix(densChimatrix(iorb,jorb),Nstate=state_list%size)
             call lanc_ed_build_densChi_mix(iorb,jorb)
          end do
       end do
    endif
    !
    if(MPIMASTER)call stop_timer
    !
  end subroutine build_densChi_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_ed_build_densChi_diag(iorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                     :: iorb
    !
    write(LOGfile,"(A)")"Get Chi_dens_l"//reg(txtfy(iorb))
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(densChimatrix(iorb,iorb),istate,Nchan=1)
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       vvinit = apply_op_N(v_state,iorb,isector)
       call tridiag_Hv_sector_normal(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_densChi(norm2,e_state,alfa_,beta_,iorb,iorb)
       deallocate(alfa_,beta_,vvinit)
       if(allocated(v_state))deallocate(v_state)
    enddo
    !
    return
  end subroutine lanc_ed_build_densChi_diag



  subroutine lanc_ed_build_densChi_mix(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                          :: iorb,jorb
    real(8),dimension(:),allocatable :: vI,vJ
    !
    write(LOGfile,"(A)")"Get Chi_dens_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(densChimatrix(iorb,jorb),istate,Nchan=1)
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       !EVALUATE (N_jorb + N_iorb)|gs> = N_jorb|gs> + N_iorb|gs>
       vI = apply_op_N(v_state,iorb,isector)
       vJ = apply_op_N(v_state,jorb,isector)
       call tridiag_Hv_sector_normal(isector,vI+vJ,alfa_,beta_,norm2)
       call add_to_lanczos_densChi(norm2,e_state,alfa_,beta_,iorb,jorb)
       deallocate(alfa_,beta_,vI,vJ)
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_densChi_mix








  subroutine add_to_lanczos_densChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                                    :: iorb,jorb
    real(8)                                    :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_densChi: add-up to GF istate "//str(istate)
#endif
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    if((finiteT).and.(beta*(Ei-Egs) < 200))then
       pesoBZ = exp(-beta*(Ei-Egs))
    elseif(.not.finiteT)then
       pesoBZ = 1d0
    else
       pesoBZ = 0d0
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
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A)")&
         "DEBUG add_to_lanczos_densChi: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(densChiMatrix(iorb,jorb),istate,1,Nlanc)
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       !
       densChiMatrix(iorb,jorb)%state(istate)%channel(1)%weight(j) = peso
       densChiMatrix(iorb,jorb)%state(istate)%channel(1)%poles(j)  = de
       !
    enddo
    !
  end subroutine add_to_lanczos_densChi






  !################################################################
  !################################################################
  !################################################################
  !################################################################




  
  function get_densChi_normal(zeta,axis) result(Chi)
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)         :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                  :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(Norb,Norb,size(zeta)) :: Chi  !Charge susceptibility matrix
    integer                                    :: iorb,jorb,i
    character(len=1)                           :: axis_
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG get_densChi_normal"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(densChimatrix))stop "get_densChi_normal ERROR: densChimatrix not allocated!"
    !
    Chi = zero
    !
    do iorb=1,Norb
       call get_Chiab(iorb,iorb)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             call get_Chiab(iorb,jorb)
             Chi(iorb,jorb,:) = 0.5d0*(Chi(iorb,jorb,:)-Chi(iorb,iorb,:)-Chi(jorb,jorb,:))
             Chi(jorb,iorb,:) = Chi(iorb,jorb,:)
          enddo
       enddo
    end if
    !
  contains
    !
    subroutine get_Chiab(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in) :: iorb,jorb
      integer            :: Nstates,istate
      integer            :: Nchannels,ichan
      integer            :: Nexcs,iexc
      real(8)            :: peso,de
      !
      write(LOGfile,"(A)")"Get Chi_dens_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
      if(.not.allocated(densChimatrix(iorb,jorb)%state)) return
      !
      Chi(iorb,jorb,:)= zero
      Nstates = size(densChimatrix(iorb,jorb)%state)
      do istate=1,Nstates
         if(.not.allocated(densChimatrix(iorb,jorb)%state(istate)%channel))cycle
         Nchannels = size(densChimatrix(iorb,jorb)%state(istate)%channel)
         do ichan=1,Nchannels
            Nexcs  = size(densChimatrix(iorb,jorb)%state(istate)%channel(ichan)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso = densChimatrix(iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
               de   = densChimatrix(iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
               select case(axis_)
               case("m","M")
                  do i=1,size(zeta)
                    if (abs(zeta(i))<1e-10)then
                      if(beta*dE > 1d-3)Chi(iorb,jorb,i)=Chi(iorb,jorb,i) + &
                        peso*2*(1d0-exp(-beta*dE))/dE 
                     else
                       Chi(iorb,jorb,i)=Chi(iorb,jorb,i) + &
                            peso*(1d0-exp(-beta*dE))*2d0*dE/(dimag(zeta(i))**2 + dE**2)
                    endif
                  enddo
               case("r","R")
                  do i=1,size(zeta)
                     Chi(iorb,jorb,i)=Chi(iorb,jorb,i) - &
                          peso*(1d0-exp(-beta*dE))*(1d0/(zeta(i) - dE) - 1d0/(zeta(i) + dE))
                  enddo
               case("t","T")
                  do i=1,size(zeta)
                     Chi(iorb,jorb,i)=Chi(iorb,jorb,i) + peso*exp(-zeta(i)*dE)
                  enddo
               end select
            enddo
         enddo
      enddo
      return
    end subroutine get_Chiab
    !
  end function get_densChi_normal






END MODULE ED_CHI_DENS
























