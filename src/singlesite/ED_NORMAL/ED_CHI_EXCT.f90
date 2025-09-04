MODULE ED_CHI_EXCT
  !:synopsis: Routines for excitonic susceptibility calculation, :code:`NORMAL` case
  !Evaluates the impurity excitonc susceptibility.
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  USE ED_AUX_FUNX

  implicit none
  private


  public :: build_exctChi_normal
  public :: get_exctChi_normal

  integer                          :: istate,iorb,jorb,ispin,jspin
  integer                          :: isector,jsector,ksector,lsector
  real(8),allocatable              :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: ipos,jpos
  integer                          :: i,j,k
  real(8)                          :: sgn,norm2
  real(8),dimension(:),allocatable :: v_state
  real(8)                          :: e_state



contains


  !+------------------------------------------------------------------+
  !                            EXCITON
  !PURPOSE  : Evaluate the Exciton susceptibility \Chi_exct for a
  ! \chi_ab = <O*_a(\tau)O_b(0)>
  ! a/=b
  ! Singlet: \sum_\sigma <C^+_{a\sigma}C_{b\sigma}
  ! Triplet: \sum_{\sigma\rho} C^+_{a\sigma} \tau_{\sigma\rho} C_{b\rho}
  !+------------------------------------------------------------------+
  subroutine build_exctChi_normal()
    ! Evaluates the impurity exciton-exciton susceptibility :math:`\chi^{X}_{ab}=\langle T_\tau X^\dagger_{ab}(\tau) X_{ab}\rangle` in the Matsubara :math:`i\omega_n` and Real :math:`\omega` frequency axis, the imaginary time :math:`\tau` as well as the singlet and triplet components of the operator.
    ! As for the Green's function, the off-diagonal component of the the susceptibility is determined using an algebraic manipulation to ensure use of Hermitian operator in the dynamical Lanczos.
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    write(LOGfile,"(A)")"Get exciton Chi:"
    if(MPIMASTER)call start_timer(unit=LOGfile)
    !
    do iorb=1,Norb
       do jorb=iorb,Norb
          call allocate_GFmatrix(exctChimatrix(1,iorb,jorb),Nstate=state_list%size)
          call lanc_ed_build_exctChi_Singlet(iorb,jorb)
          call allocate_GFmatrix(exctChimatrix(2,iorb,jorb),Nstate=state_list%size)
          call lanc_ed_build_exctChi_tripletXY(iorb,jorb)
          call allocate_GFmatrix(exctChimatrix(3,iorb,jorb),Nstate=state_list%size)
          call lanc_ed_build_exctChi_tripletZ(iorb,jorb)
       end do
    end do
    !
    if(MPIMASTER)call stop_timer
    !
  end subroutine build_exctChi_normal


  ! \chi_ab  = <Delta*_ab(\tau)Delta_ab(0)>
  !\Delta_ab = \sum_\sigma C^+_{a\sigma}C_{b\sigma}
  subroutine lanc_ed_build_exctChi_singlet(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                          :: iorb,jorb
    type(sector)                     :: sectorI
    !
    write(LOGfile,"(A)")"Get singlet Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    do istate=1,state_list%size
       !
       call allocate_GFmatrix(exctChimatrix(1,iorb,jorb),istate,Nchan=2)
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       ksector = getCsector(1,2,isector)
       lsector = getCsector(1,1,isector)
       !
       ! Lesser
       call auxiliary_vvinit_manipulation(ksector, &
                                          lsector, &
                                          iorb,    &
                                          jorb,    &  
                                          1,       & !indx: singlet component        
                                          1,       & !ichan: lesser
                                          istate)
       ! Greater
       call auxiliary_vvinit_manipulation(ksector, & 
                                          lsector, &
                                          iorb,    &
                                          jorb,    &       
                                          1,       & !indx: singlet component
                                          2,       & !ichan: greater
                                          istate)
       ! Cleanup
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_exctChi_singlet


  ! \chi_ab  = <O_ab(\tau)O_ab(0)>
  ! O_ab = \sum_sp C^+_{as}.tau^o_{sp}.C_{bp} with o=X,Y
  ! O_ab|0> X:=   [C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}]|0>
  !         Y:= -i[C^+_{a,up}C_{b,dw} - C^+_{a,dw}C_{b,up}]|0>
  !         X:=   [P_{up,dw} +  P_{dw,up}]|0> = |v> + |w>
  !         X:= -i[P_{up,dw} -  P_{dw,up}]|0> = |v> + |w>
  ! If |0>\in\SS(N_up,N_dw) => |v>\in\SS(N_up+1,N_dw-1), |w>\in\SS(N_up-1,N_dw+1)
  ! so that the sum |v> + |w> can not be accumulated onto a single vector |vvinit>
  ! Yet, we can recast the \Chi_ab expression in:
  ! \chi_ab = (<w| + <v|)[z-H]^{-1}(|v> + |w>)
  ! the direct terms: <v|[z-H]^{-1}|v> and <w|[z-H]^{-1}|w> are evaluated as usual.
  ! the mixed terms: <v|[z-H]^{-1}|w> and <w|[z-H]^{-1}|v> are indeed null.
  ! Proof:
  ! |v> and |w> belong to different sectors. H has a sector-block structure and so
  ! does its inverse H^{-1} made of the inverse of each block.
  ! The expected values <v|H^{-1}|w> are taken across different sectors, but unless
  ! spin-conservation is broken these sectors are not connected and, as such, these
  ! numbers have to be zero.
  ! Note that while H can have a sparse nature, its inverse in general does not.
  ! For this reason the same argument as above DOES NOT apply to the Z case as
  ! in that case |v> and |w> belong to the same sector (the same as |0>) and the
  ! mixed term is in general non null.
  subroutine lanc_ed_build_exctChi_tripletXY(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                          :: iorb,jorb
    real(8),dimension(:),allocatable :: vtmp
    !
    write(LOGfile,"(A)")"Get triplet XY Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(exctChimatrix(2,iorb,jorb),istate,Nchan=4)
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       !X - Component == Y -Component
       !X_{ab}= C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}
       !
       !C^+_{a,dw}C_{b,up}: First Lesser component
       ksector = getCsector(1,1,isector);jsector=0
       if(ksector/=0)jsector = getCDGsector(1,2,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !C_{b,up}|gs>   =|tmp>
          vtmp   = apply_op_C(v_state,jorb,1,isector,ksector)
          !C^+_{a,dw}|tmp>=|vvinit>
          vvinit = apply_op_CDG(vtmp,iorb,2,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,+1,iorb,jorb,2,1,istate)
          deallocate(alfa_,beta_,vtmp,vvinit)
       else
          call allocate_GFmatrix(exctChiMatrix(2,iorb,jorb),istate,1,Nexc=0)
       endif
       !
       !C^+_{b,up}C_{a,dw}: First Greater component
       ksector = getCsector(1,2,isector);jsector=0
       if(ksector/=0)jsector = getCDGsector(1,1,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !C_{a,dw}|gs>   =|tmp>
          vtmp   = apply_op_C(v_state,iorb,2,isector,ksector)
          !C^+_{b,up}|tmp>=|vvinit>
          vvinit = apply_op_CDG(vtmp,jorb,1,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,-1,iorb,jorb,2,2,istate)
          deallocate(alfa_,beta_,vtmp,vvinit)
       else
          call allocate_GFmatrix(exctChiMatrix(2,iorb,jorb),istate,2,Nexc=0)
       endif
       !
       !C^+_{a,up}C_{b,dw}: Second Lesser component
       ksector = getCsector(1,2,isector);jsector=0
       if(ksector/=0)jsector = getCDGsector(1,1,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !C_{b,dw}|gs>   =|tmp>
          vtmp   = apply_op_C(v_state,jorb,2,isector,ksector)
          !C^+_{a,up}|tmp>=|vvinit>
          vvinit = apply_op_CDG(vtmp,iorb,1,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,+1,iorb,jorb,2,3,istate)
          deallocate(alfa_,beta_,vtmp,vvinit)
       else
          call allocate_GFmatrix(exctChiMatrix(2,iorb,jorb),istate,3,Nexc=0)
       endif
       !
       !C^+_{b,dw}C_{a,up}: Second Greater component
       ksector = getCsector(1,1,isector);jsector=0
       if(ksector/=0)jsector = getCDGsector(1,2,ksector)
       if(jsector/=0.AND.ksector/=0)then
          !C_{a,up}|gs>   =|tmp>
          vtmp   = apply_op_C(v_state,iorb,1,isector,ksector)
          !C^+_{b,dw}|tmp>=|vvinit>
          vvinit = apply_op_CDG(vtmp,jorb,2,ksector,jsector)
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_exctChi(norm2,e_state,alfa_,beta_,-1,iorb,jorb,2,4,istate)
          deallocate(alfa_,beta_,vtmp,vvinit)
       else
          call allocate_GFmatrix(exctChiMatrix(2,iorb,jorb),istate,4,Nexc=0)
       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    return
  end subroutine lanc_ed_build_exctChi_tripletXY



  ! \chi_ab  = <Z_ab(\tau)Z_ab(0)>
  !Z_ab = \sum_sp C^+_{as}.tau^z_{sp}.C_{bp}
  subroutine lanc_ed_build_exctChi_tripletZ(iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                          :: iorb,jorb
    !
    write(LOGfile,"(A)")"Get triplet Z Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !
    !
    do istate=1,state_list%size
       call allocate_GFmatrix(exctChimatrix(3,iorb,jorb),istate,Nchan=2)
       isector    =  es_return_sector(state_list,istate)
       e_state    =  es_return_energy(state_list,istate)
       v_state    =  es_return_dvec(state_list,istate)
       !
       !Z - Component:
       !Z_{ab}= C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw}
       ksector = getCsector(1,2,isector)
       lsector = getCsector(1,1,isector)
       !
       ! Lesser
       call auxiliary_vvinit_manipulation(ksector,   &
                                          lsector,   &
                                          iorb,      &
                                          jorb,      &  
                                          3,         & !indx: singlet component        
                                          1,         & !ichan: lesser
                                          istate)
       ! Greater
       call auxiliary_vvinit_manipulation(ksector,   &
                                          lsector,   &
                                          iorb,      &
                                          jorb,      &  
                                          3,         & !indx: singlet component        
                                          2,         & !ichan: lesser
                                          istate)
       ! Cleanup                                   
       if(allocated(v_state))deallocate(v_state)
    enddo
    return
  end subroutine lanc_ed_build_exctChi_tripletZ



  subroutine add_to_lanczos_exctChi(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,indx,ichan,istate)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                                    :: iorb,jorb,ichan,indx,isign,istate
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
         "DEBUG add_to_lanczos_exctChi: add-up to GF istate "//str(istate)
#endif
    !
    if(vnorm2==0d0)then
       call allocate_GFmatrix(exctChiMatrix(indx,iorb,jorb),istate,ichan,Nexc=0)
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
         "DEBUG add_to_lanczos_exctChi: LApack tridiagonalization"
#endif
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    call allocate_GFmatrix(exctChiMatrix(indx,iorb,jorb),istate,ichan,Nlanc)
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       !
       exctChiMatrix(indx,iorb,jorb)%state(istate)%channel(ichan)%weight(j) = peso
       exctChiMatrix(indx,iorb,jorb)%state(istate)%channel(ichan)%poles(j)  = isign*de
    enddo
    !
  end subroutine add_to_lanczos_exctChi




  !################################################################
  !################################################################
  !################################################################
  !################################################################




  function get_exctChi_normal(zeta,axis) result(Chi)
    ! Reconstructs the system impurity electrons Green's functions using :f:var:`impgmatrix` to retrieve weights and poles.
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)             :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                      :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(3,Norb,Norb,size(zeta))   :: Chi  !Excitonic susceptibility matrix. dim-1, 1= :math:`S`, 2= :math:`T^x`, 3= :math:`T^z`
    integer                                        :: iorb,jorb,i,indx
    character(len=1)                               :: axis_
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG get_exctChi_normal"
#endif
    !
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1) !only for self-consistency, not used here
    !
    if(.not.allocated(exctChimatrix))stop "get_exctChi_normal ERROR: exctChimatrix not allocated!"
    !
    Chi = zero
    !
    do indx=1,3
       do iorb=1,Norb
          do jorb=iorb,Norb
             call get_Chiab(indx,iorb,jorb)
          enddo
       enddo
    enddo
    !
  contains
    !
    subroutine get_Chiab(indx,iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
      integer,intent(in) :: indx,iorb,jorb
      integer            :: Nstates,istate
      integer            :: Nchannels,ichan
      integer            :: Nexcs,iexc
      real(8)            :: peso,de
      !
      select case(indx)
      case(1); write(LOGfile,"(A)")"Get Chi_exct_S_l"//str(iorb)//str(jorb)
      case(2); write(LOGfile,"(A)")"Get Chi_exct_XY_l"//str(iorb)//str(jorb)
      case(3); write(LOGfile,"(A)")"Get Chi_exct_Z_l"//str(iorb)//str(jorb)
      end select
      if(.not.allocated(exctChimatrix(indx,iorb,jorb)%state)) return
      !
      Chi(indx,iorb,jorb,:) = zero
      if(iorb /= jorb) then
         Chi(indx,jorb,iorb,:) = zero
      endif
      Nstates = size(exctChimatrix(indx,iorb,jorb)%state)
      do istate=1,Nstates
         if(.not.allocated(exctChimatrix(indx,iorb,jorb)%state(istate)%channel))cycle
         Nchannels = size(exctChimatrix(indx,iorb,jorb)%state(istate)%channel)
         do ichan=1,Nchannels
            Nexcs  = size(exctChimatrix(indx,iorb,jorb)%state(istate)%channel(ichan)%poles)
            if(Nexcs==0)cycle
            do iexc=1,Nexcs
               peso = exctChimatrix(indx,iorb,jorb)%state(istate)%channel(ichan)%weight(iexc)
               de   = exctChimatrix(indx,iorb,jorb)%state(istate)%channel(ichan)%poles(iexc)
               ! Zero energy poles are to be taken into account only once
               ! The coefficient 0.5 accounts for the two lesser/greater channels
               if(abs(beta*de) < 1e-8) then
                 select case(axis_)
                    case("m","M")
                       do i=1,size(zeta)
                          if(abs(zeta(i))<1e-10) then ! \nu=0
                             Chi(indx,iorb,jorb,i) = Chi(indx,iorb,jorb,i) + 0.5*peso*beta
                             if(iorb /= jorb) then
                                Chi(indx,jorb,iorb,i) = Chi(indx,jorb,iorb,i) + 0.5*peso*beta
                             endif
                          endif
                       enddo
                    case("r","R")
                       ! A zero energy pole contributes only to \chi(Re(z) = 0) and
                       ! its contribution is taken to be the same as for \chi(\nu=0),
                       ! regardless of the imaginary shift in z
                       do i=1,size(zeta)
                          if(abs(dreal(zeta(i)))<1e-10) then
                             Chi(indx,iorb,jorb,i) = Chi(indx,iorb,jorb,i) + 0.5*peso*beta
                             if(iorb /= jorb) then
                                Chi(indx,jorb,iorb,i) = Chi(indx,jorb,iorb,i) + 0.5*peso*beta
                             endif
                          endif
                       enddo
                    case("t","T")
                       Chi(indx,iorb,jorb,:) = Chi(indx,iorb,jorb,:) + 0.5*peso
                       if(iorb /= jorb) then
                          Chi(indx,jorb,iorb,:) = Chi(indx,jorb,iorb,:) + 0.5*peso
                       endif
                 end select
               ! Nonzero energy poles
               elseif(merge(de, -de, mod(ichan,2)==1)>0) then
                  select case(axis_)
                     case("m","M","r","R")
                        do i=1,size(zeta)
                           if(mod(ichan,2) == 1) then ! Lesser
                              Chi(indx,iorb,jorb,i) = Chi(indx,iorb,jorb,i) - &
                                 peso*(1d0-exp(-beta*de)) / (zeta(i) - de)
                              if(iorb /= jorb) then
                                 Chi(indx,jorb,iorb,i) = Chi(indx,jorb,iorb,i) + &
                                 peso*(1d0-exp(-beta*de)) / (zeta(i) + de)
                              endif
                           else ! Greater
                              Chi(indx,iorb,jorb,i) = Chi(indx,iorb,jorb,i) + &
                                 peso*(1d0-exp(beta*de)) / (zeta(i) - de)
                              if(iorb /= jorb) then
                                 Chi(indx,jorb,iorb,i) = Chi(indx,jorb,iorb,i) - &
                                 peso*(1d0-exp(beta*de)) / (zeta(i) + de)
                              endif
                           endif
                        enddo
                     case("t","T")
                        do i=1,size(zeta)
                           if(mod(ichan,2) == 1) then ! Lesser
                              Chi(indx,iorb,jorb,i) = Chi(indx,iorb,jorb,i) + peso*exp(-zeta(i)*de)
                              if(iorb /= jorb) then
                                 Chi(indx,jorb,iorb,i) = Chi(indx,jorb,iorb,i) + peso*exp(-(beta-zeta(i))*de)
                              endif
                           else ! Greater
                              Chi(indx,iorb,jorb,i) = Chi(indx,iorb,jorb,i) + peso*exp((beta-zeta(i))*de)
                              if(iorb /= jorb) then
                                 Chi(indx,jorb,iorb,i) = Chi(indx,jorb,iorb,i) + peso*exp(zeta(i)*de)
                              endif
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
  end function get_exctChi_normal


  subroutine auxiliary_vvinit_manipulation(ksector,lsector,iorb,jorb,indx,ichan,istate)
    !this subroutine applies, when possible, the operators for the singlet and triplet_Z components
    !of the excitonic susceptibility
    integer                          :: iorb,jorb,ksector,lsector,isign,indx,ichan,istate
    real(8)                          :: comb_sign
    real(8),dimension(:),allocatable :: vup,vdw,vtmp
    !
    !See if we have to do anythin: if both c operators cannot be applied, exit.
    if(ksector==0 .AND. lsector==0)then
       call allocate_GFmatrix(exctChimatrix(indx,iorb,jorb),istate,1,Nexc=0)
       call allocate_GFmatrix(exctChimatrix(indx,iorb,jorb),istate,2,Nexc=0)
    else
      !startup
      if(allocated(alfa_))deallocate(alfa_)
      if(allocated(beta_))deallocate(beta_)
      if(allocated(vtmp))deallocate(vtmp)
      if(allocated(vup))deallocate(vup)
      if(allocated(vdw))deallocate(vdw)
      !
      !Set sign of vup+vdw linear combination depending if singlet (+) or triplet_z (-)
      if (indx==1)then
        comb_sign = 1.0
      elseif(indx==3)then
        comb_sign = -1.0
      else
        STOP "Wrong value of indx: only 1 or 3 are acceptable"
      endif
      !
      !Set isign depending if ichan=1 (lesser -> +) or ichan = 2 (greater-> -)
      if (ichan==1)then
        isign = +1
      elseif(ichan==2)then
        isign = -1
      else
        STOP "Wrong value of ichan: only 1 or 2 are acceptable"
      endif  
      !
      !try apply operators for down spin
      if(ksector/=0)then
        !C_b,dw|gs>=|tmp>
        vtmp = apply_op_C(v_state,jorb,2,isector,ksector)
        !C^+_a,dw|tmp>=|vvinit>
        vdw  = apply_op_CDG(vtmp,iorb,2,ksector,isector)
      endif
      !
      !try apply operators for up spin
      if(lsector/=0)then
        !C_b,up|gs>=|tmp>
        vtmp = apply_op_C(v_state,jorb,1,isector,lsector)
        !C^+_a,up|tmp>=|vvinit>
        vup  = apply_op_CDG(vtmp,iorb,1,lsector,isector)
      endif
      !
      !tridiagonalize
      if (allocated(vup) .and. allocated(vdw)) then
        call tridiag_Hv_sector_normal(isector, vup + comb_sign * vdw ,alfa_, beta_, norm2)
      elseif (allocated(vup)) then
        call tridiag_Hv_sector_normal(isector, vup, alfa_, beta_, norm2)
      elseif (allocated(vdw)) then
        call tridiag_Hv_sector_normal(isector, vdw, alfa_, beta_, norm2)
      else
        STOP "neither vup nor vdw are allocated"
      endif
      !
      !save weights and poles
      call add_to_lanczos_exctChi(norm2, e_state, alfa_, beta_, isign, iorb, jorb, indx, ichan, istate)
      !
      !cleanup
      deallocate(alfa_,beta_,vtmp)
      if(allocated(vup))deallocate(vup)
      if(allocated(vdw))deallocate(vdw)
    endif
  end subroutine auxiliary_vvinit_manipulation



END MODULE ED_CHI_EXCT
