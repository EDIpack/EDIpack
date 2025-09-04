MODULE ED_CHI_PAIR
  !:synopsis: Routines for pair susceptibility calculation, :code:`NORMAL` case
  !Evaluates the impurity pair susceptibility.
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh
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


  public :: build_pairChi_normal
  public :: get_pairChi_normal

  integer                          :: istate,iorb,jorb,ispin,jspin
  integer                          :: isector,jsector,ksector
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ialfa
  integer                          :: jalfa
  integer                          :: i,j,k
  real(8)                          :: sgn
  real(8)                          :: norm2
  real(8),dimension(:),allocatable :: v_state
  real(8)                          :: e_state




contains


  !+------------------------------------------------------------------+
  !                            PAIR
  !PURPOSE  : Evaluate the pair susceptibility \Chi_pair for a
  ! \chi_ab = <Delta*_a(\tau)Delta_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_pairChi_normal()
    ! Evaluates the impurity Pair susceptibility :math:`\chi^{\Delta}=\langle T_\tau \Delta^\dagger_a(\tau) \Delta_b\rangle` in the Matsubara :math:`i\omega_n` and Real :math:`\omega` frequency axis as well as imaginary time :math:`\tau`.
    ! As for the Green's function, the off-diagonal component of the the susceptibility is determined using an algebraic manipulation to ensure use of Hermitian operator in the dynamical Lanczos.
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    write(LOGfile,"(A)")"Get pair Chi:"
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !
    do iorb=1,Norb
       call allocate_GFmatrix(pairChimatrix(iorb,iorb),Nstate=state_list%size)
       call lanc_ed_build_pairChi_diag(iorb)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             call allocate_GFmatrix(pairChimatrix(iorb,jorb),Nstate=state_list%size)
             call lanc_ed_build_pairChi_mix(iorb,jorb)
          end do
       end do
    endif
    !
    if(MPIMASTER)call stop_timer
    !
  end subroutine build_pairChi_normal




  !################################################################
  !################################################################
  !################################################################
  !################################################################




  ! \chi_aa = <Delta*_a(\tau)Delta_a(0)>
  !         = <[C^+_a(\tau)C^+_a(\tau)][C_a(0)C_a(0)]>
  subroutine lanc_ed_build_pairChi_diag(iorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                          :: iorb
    real(8),dimension(:),allocatable :: vtmp
    !
    write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
    !
    if(ed_total_ud)then
       ialfa = 1
    else
       ialfa = iorb
    endif
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(pairChimatrix(iorb,iorb),istate,Nchan=2)
       isector  =  es_return_sector(state_list,istate)
       e_state  =  es_return_energy(state_list,istate)
       v_state  =  es_return_dvec(state_list,istate)
       !
       ! Lesser
       !
       ksector = getCsector(ialfa,2,isector);jsector=0
       if(ksector/=0)jsector = getCsector(ialfa,1,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !C_dw|gs>  = |tmp>
          vtmp   = apply_op_C(v_state,iorb,2,isector,ksector)
          !C_up|tmp> = C_up[C_dw|gs>] = |vvinit>
          vvinit = apply_op_C(vtmp,iorb,1,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_pairChi(norm2,e_state,alfa_,beta_,+1,iorb,iorb,ichan=1,istate=istate)
          deallocate(alfa_,beta_,vtmp,vvinit)
       else
          call allocate_GFmatrix(pairChiMatrix(iorb,iorb),istate,1,Nexc=0)
       endif
       !
       ! Greater
       !
       ksector = getCDGsector(ialfa,1,isector);jsector=0
       if(ksector/=0)jsector = getCDGsector(ialfa,2,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !C^+_up|gs>  = |tmp>
          vtmp   = apply_op_CDG(v_state,iorb,1,isector,ksector)
          !C^+_dw|tmp> = C^+_dw[C^+_up|gs>] = |vvinit>
          vvinit = apply_op_CDG(vtmp,iorb,2,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_pairChi(norm2,e_state,alfa_,beta_,-1,iorb,iorb,ichan=2,istate=istate)
          deallocate(alfa_,beta_,vtmp,vvinit)
       else
          call allocate_GFmatrix(pairChiMatrix(iorb,iorb),istate,2,Nexc=0)
       endif
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_pairChi_diag





  ! \chi_ab = <Delta*_a(\tau)Delta_b(0)>
  !         = <[C^+_a(\tau)C^+_a(\tau)][C_b(0)C_b(0)]>
  !from aux: <[C^+_a C^+_a + C^+_b C^+_b][C_a C_a + C_b C_b]>
  subroutine lanc_ed_build_pairChi_mix(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                             :: iorb,jorb
    real(8),dimension(:),allocatable    :: va,vb,vtmp
    !
    write(LOGfile,"(A)")"Get Chi_pair_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    if(.not.ed_total_ud)then
       write(LOGfile,"(A)")"ED_CHI_PAIR warning: can not evaluate \Chi_pair_ab with ed_total_ud=F"
       return
    endif
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(pairChimatrix(iorb,jorb),istate,Nchan=2) !Nchan+4 for complex case
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state  =  es_return_dvec(state_list,istate)
       !
       !First lesser
       ! --> Apply [C_b C_b + C_a C_a]|state>
       ksector = getCsector(1,2,isector);jsector=0
       if(ksector/=0)jsector = getCsector(1,1,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !Apply C_a,up*C_a,dw:
          !C_a.dw|gs>  = |tmp>
          vtmp = apply_op_C(v_state,iorb,2,isector,ksector)
          !C_a.up|tmp> = C_a.up[C_a.dw|gs>] = |vvinit>
          va   = apply_op_C(vtmp,iorb,1,ksector,jsector)
          !Apply + C_b,up*C_b,dw
          !C_b.dw|gs>  = |tmp>
          vtmp = apply_op_C(v_state,jorb,2,isector,ksector)
          !C_b.up|tmp> = C_b.up[C_b.dw|gs>] = |vvinit>
          vb   = apply_op_C(vtmp,jorb,1,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,va+vb,alfa_,beta_,norm2)
          call add_to_lanczos_pairChi(norm2,e_state,alfa_,beta_,+1,iorb,jorb,ichan=1,istate=istate)
          deallocate(alfa_,beta_,vtmp,va,vb)
       else
          call allocate_GFmatrix(pairChiMatrix(iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !First greater
       ! --> Apply [C^+_b C^+_b + C^+_a C^+_a]|state>
       ksector = getCDGsector(1,1,isector);jsector=0
       if(ksector/=0)jsector = getCDGsector(1,2,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !Apply C^+_a,dw*C^+_a,up:
          !C^+_a.up|gs>  = |tmp>
          vtmp = apply_op_CDG(v_state,iorb,1,isector,ksector)
          !C^+_a.dw|tmp> = C^+_a.dw[C^+_a.up|gs>] = |vvinit>
          va   = apply_op_CDG(vtmp,iorb,2,ksector,jsector)
          !Apply + C^+_b,dw*C^+_b,up
          !C^+_b.up|gs>  = |tmp>
          vtmp = apply_op_CDG(v_state,jorb,1,isector,ksector)
          !C^+_b.dw|tmp> = C^+_b.dw[C^+_b.up|gs>] = |vvinit>
          vb   = apply_op_CDG(vtmp,jorb,2,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,va+vb,alfa_,beta_,norm2)
          call add_to_lanczos_pairChi(norm2,e_state,alfa_,beta_,-1,iorb,jorb,ichan=2,istate=istate)
          deallocate(alfa_,beta_,vtmp,va,vb)
       else
          call allocate_GFmatrix(pairChiMatrix(iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !Uncomment for complex case
       !!Second lesser
       !! --> Apply [C_a C_a - xi* C_b C_b]|state>
       !ksector = getCsector(1,2,isector);jsector=0
       !if(ksector/=0)jsector = getCsector(1,1,ksector)
       !if(jsector/=0.AND.ksector/=0)then
       !   !Apply C_a,up*C_a,dw:
       !   !C_a.dw|gs>  = |tmp>
       !   vtmp = apply_op_C(v_state,iorb,2,isector,ksector)
       !   !C_a.up|tmp> = C_a.up[C_a.dw|gs>] = |vvinit>
       !   va   = apply_op_C(vtmp,iorb,1,ksector,jsector)
       !   !Apply + C_b,up*C_b,dw
       !   !C_b.dw|gs>  = |tmp>
       !   vtmp = apply_op_C(v_state,jorb,2,isector,ksector)
       !   !C_b.up|tmp> = C_b.up[C_b.dw|gs>] = |vvinit>
       !   vb   = apply_op_C(vtmp,jorb,1,ksector,jsector)
       !   call tridiag_Hv_sector_normal_complex(jsector,va - xi*vb,alfa_,beta_,norm2)
       !   call add_to_lanczos_pairChi(-xi*norm2,e_state,alfa_,beta_,+1,iorb,jorb,ichan=3,istate=istate)
       !   deallocate(alfa_,beta_,vtmp,va,vb)
       !else
       !   call allocate_GFmatrix(pairChiMatrix(iorb,jorb),istate,3,Nexc=0)
       !endif
       !
       !!Second greater
       !! --> Apply [C^+_a C^+_a + xi* C^+_b C^+_b]|state>
       !ksector = getCDGsector(1,1,isector);jsector=0
       !if(ksector/=0)jsector = getCDGsector(1,2,ksector)
       !if(jsector/=0.AND.ksector/=0)then
       !   !Apply C^+_a,dw*C^+_a,up:
       !   !C^+_a.up|gs>  = |tmp>
       !   vtmp = apply_op_CDG(v_state,iorb,1,isector,ksector)
       !   !C^+_a.dw|tmp> = C^+_a.dw[C^+_a.up|gs>] = |vvinit>
       !   va   = apply_op_CDG(vtmp,iorb,2,ksector,jsector)
       !   !Apply + C^+_b,dw*C^+_b,up
       !   !C^+_b.up|gs>  = |tmp>
       !   vtmp = apply_op_CDG(v_state,jorb,1,isector,ksector)
       !   !C^+_b.dw|tmp> = C^+_b.dw[C^+_b.up|gs>] = |vvinit>
       !   vb   = apply_op_CDG(vtmp,jorb,2,ksector,jsector)
       !   call tridiag_Hv_sector_normal_complex(jsector,va + xi*vb,alfa_,beta_,norm2)
       !   call add_to_lanczos_pairChi(-xi*norm2,e_state,alfa_,beta_,-1,iorb,jorb,ichan=4,istate=istate)
       !   deallocate(alfa_,beta_,vtmp,va,vb)
       !else
       !   call allocate_GFmatrix(pairChiMatrix(iorb,jorb),istate,4,Nexc=0)
       !endif
       !!
       !if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_pairChi_mix







  subroutine add_to_lanczos_pairChi(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ichan,istate)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                                    :: iorb,jorb,ichan,isign,istate
    real(8)                                    :: pesoF,pesoAB,pesoBZ,peso
    real(8)                                    :: vnorm2
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
         "DEBUG add_to_lanczos_pairChi: add-up to GF istate "//str(istate)
#endif
    !
    if(vnorm2==0d0)then
       call allocate_GFmatrix(pairChiMatrix(iorb,jorb),istate,ichan,Nexc=0)
       return
    endif
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
         "DEBUG add_to_lanczos_pairChi: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(pairChiMatrix(iorb,jorb),istate,ichan,Nlanc)
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       !
       pairChiMatrix(iorb,jorb)%state(istate)%channel(ichan)%weight(j) = peso
       pairChiMatrix(iorb,jorb)%state(istate)%channel(ichan)%poles(j)  = isign*de
    enddo
    !
  end subroutine add_to_lanczos_pairChi



  !################################################################
  !################################################################
  !################################################################
  !################################################################




  function get_pairChi_normal(zeta,axis) result(Chi)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
    complex(8),dimension(:),intent(in)         :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                  :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(Norb,Norb,size(zeta)) :: Chi  !Pairing susceptibility matrix
    integer                                    :: iorb,jorb,i
    character(len=1)                           :: axis_
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG get_pairChi_normal"
#endif
    !
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(pairChimatrix))stop "get_pairChi_normal ERROR: pairChimatrix not allocated!"
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
             !if (iorb /= jorb) then
                call get_Chiab(iorb,jorb)
                !Uncomment for complex case
                !Chi(iorb,jorb,:) = 0.5d0*(Chi(iorb,jorb,:) - (1-xi)*Chi(iorb,iorb,:) - (1-xi)*Chi(jorb,jorb,:))
                Chi(iorb,jorb,:) = 0.5d0*(Chi(iorb,jorb,:) - Chi(iorb,iorb,:) - Chi(jorb,jorb,:))
                Chi(jorb,iorb,:) = Chi(iorb,jorb,:)
             !endif
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
      write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
      if(.not.allocated(pairChimatrix(iorb,jorb)%state)) return
      !
      Chi(iorb,jorb,:)= zero
      Nstates = size(pairChimatrix(iorb,jorb)%state)
      do istate=1,Nstates
         if(.not.allocated(pairChimatrix(iorb,jorb)%state(istate)%channel))cycle
         Nchannels = size(pairChimatrix(iorb,jorb)%state(istate)%channel)
         do ichan=1,Nchannels
            Nexcs  = size(pairChimatrix(iorb,jorb)%state(istate)%channel(ichan)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso = pairChimatrix(iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
               de   = pairChimatrix(iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
               ! Zero energy poles are to be taken into account only once
               ! The coefficient 0.5 accounts for the two lesser/greater channels
               if(abs(beta*de) < 1e-8) then
                 select case(axis_)
                    case("m","M")
                       do i=1,size(zeta)
                          if(abs(zeta(i))<1e-10) then ! \nu=0
                             Chi(iorb,jorb,i) = Chi(iorb,jorb,i) + 0.5*peso*beta
                             !Uncomment for complex case
                             !if(iorb /= jorb) then
                             !   Chi(jorb,iorb,i) = Chi(jorb,iorb,i) + 0.5*peso*beta
                             !endif
                          endif
                       enddo
                    case("r","R")
                       ! A zero energy pole contributes only to \chi(Re(z) = 0) and
                       ! its contribution is taken to be the same as for \chi(\nu=0),
                       ! regardless of the imaginary shift in z
                       do i=1,size(zeta)
                          if(abs(dreal(zeta(i)))<1e-10) then
                             Chi(iorb,jorb,i) = Chi(iorb,jorb,i) + 0.5*peso*beta
                             !Uncomment for complex case
                             !if(iorb /= jorb) then
                             !   Chi(jorb,iorb,i) = Chi(jorb,iorb,i) + 0.5*peso*beta
                             !endif
                          endif
                       enddo
                    case("t","T")
                       Chi(iorb,jorb,:) = Chi(iorb,jorb,:) + 0.5*peso
                 end select
               ! Nonzero energy poles
               elseif(merge(de, -de, mod(ichan,2)==1)>0) then
                  select case(axis_)
                     case("m","M","r","R")
                        do i=1,size(zeta)
                           if(mod(ichan,2) == 1) then ! Lesser
                              Chi(iorb,jorb,i) = Chi(iorb,jorb,i) - &
                                 peso*(1d0-exp(-beta*de)) / (zeta(i) - de)
                              !Uncomment for complex case
                              !if(iorb /= jorb) then
                              !   Chi(jorb,iorb,i) = Chi(jorb,iorb,i) + &
                              !   peso*(1d0-exp(-beta*de)) / (zeta(i) + de)
                              !endif
                           else ! Greater
                              Chi(iorb,jorb,i) = Chi(iorb,jorb,i) + &
                                 peso*(1d0-exp(beta*de)) / (zeta(i) - de)
                              !Uncomment for complex case
                              !if(iorb /= jorb) then
                              !   Chi(jorb,iorb,i) = Chi(jorb,iorb,i) - &
                              !   peso*(1d0-exp(beta*de)) / (zeta(i) + de)
                              !endif
                           endif
                        enddo
                     case("t","T")
                        do i=1,size(zeta)
                           if(mod(ichan,2) == 1) then ! Lesser
                              Chi(iorb,jorb,i) = Chi(iorb,jorb,i) + peso*exp(-zeta(i)*de)
                              !Uncomment for complex case
                              !if(iorb /= jorb) then
                              !   Chi(jorb,iorb,i) = Chi(jorb,iorb,i) + peso*exp(-(beta-zeta(i))*de)
                              !endif
                           else ! Greater
                              Chi(iorb,jorb,i) = Chi(iorb,jorb,i) + peso*exp((beta-zeta(i))*de)
                              !Uncomment for complex case
                              !if(iorb /= jorb) then
                              !   Chi(jorb,iorb,i) = Chi(jorb,iorb,i) + peso*exp(zeta(i)*de)
                              !endif
                           endif
                        enddo
                  end select
               endif
            enddo
         enddo
      enddo
      return
    end subroutine get_Chiab
    !
  end function get_pairChi_normal



END MODULE ED_CHI_PAIR
