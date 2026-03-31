MODULE ED_GREENS_FUNCTIONS
  !:synopsis: Routines for Green's functions calculation
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER
  USE SF_ARRAYS, only:arange
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,to_lower,splot
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_GF_NORMAL
  USE ED_GF_SUPERC
  USE ED_GF_NONSU2
  !
  implicit none
  private 

  public :: buildGf_impurity

  public :: get_impG
  public :: get_impF
  public :: get_impD
  public :: get_Sigma
  public :: get_Self
  !


  real(8),dimension(:,:),allocatable :: zimp !quasiparticle weight
  real(8),dimension(:,:),allocatable :: simp !scattering rate
  real(8)                            :: w_ph ! Renormalized phonon frequency



contains



  subroutine buildGF_impurity()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    ! 
    ! Build the quantum impurity electrons Green's functions :math:`\hat{G}` and the phonons Green's function :math:`\hat{D}` , calling the correct procedure according to the value of :f:var:`ed_mode` .
    ! Write the results on file according to the values of input variables :f:var:`ed_print_g` , :f:var:`ed_print_sigma`, :f:var:`ed_print_impD`  and :f:var:`ed_print_g0` .
    !
    ! * :code:`normal` : :f:func:`build_impg_normal`
    ! * :code:`superc` : :f:func:`build_impg_superc`
    ! * :code:`nonsu2` : :f:func:`build_impg_nonsu2`
    !
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG build_GF: build GFs"
#endif
    !
    if(.not.allocated(impGmatrix))then
       if(ed_mode=="superc")then
          allocate(impGmatrix(2*Nspin,2*Nspin,Norb,Norb))
       else
          allocate(impGmatrix(Nspin,Nspin,Norb,Norb))
       endif
    endif
    !
    call deallocate_GFmatrix(impGmatrix)
    call deallocate_GFmatrix(impDmatrix)
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default  ;STOP "buildGF_impurity error: ed_mode "//trim(ed_mode)//" not valid"
    case("normal");call build_impG_normal()
    case("superc");call build_impG_superc()
    case("nonsu2");call build_impG_nonsu2()
    end select
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG build_GF: writing results"
#endif
    if(MPIMASTER)then
       call print_impGmatrix()
       call print_impDmatrix()
       if(ed_print_G)          call print_impG()
       if(ed_print_Sigma)      call print_Sigma()
       if(ed_print_G0)         call print_impG0()
       if(ed_print_G.AND.Nph>0)call print_impD()
       !
       call get_szr
       call write_szr()
    endif
    !
  end subroutine buildGF_impurity




  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+

  function get_impG(zeta,axis) result(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impG"
#endif
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_impG error: not a valid ed_mode"
    case("normal");self = get_impG_normal(zeta,axis_)
    case("superc");self = get_impG_superc(zeta,axis_)
    case("nonsu2");self = get_impG_nonsu2(zeta,axis_)
    end select
  end function get_impG
  !
  function get_impF(zeta,axis) result(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impF"
#endif
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_impF error: not a valid ed_mode"
    case("superc");self = get_impF_superc(zeta,axis_)
    end select
  end function get_impF
  !
  function get_impD(zeta,axis) result(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in) :: zeta
    character(len=*),optional          :: axis
    complex(8),dimension(size(zeta))   :: self
    character(len=1)                   :: axis_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_impD"
#endif
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_impD error: not a valid ed_mode"
    case("normal");self = get_impD_normal(zeta,axis_)
    case("superc");self = get_impD_superc(zeta,axis_)
    case("nonsu2");self = get_impD_nonsu2(zeta,axis_)
    end select
  end function get_impD



  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+


  function get_Sigma(zeta,axis) result(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_Sigma"
#endif
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_Sigma error: not a valid ed_mode"
    case("normal");self = get_Sigma_normal(zeta,axis_)
    case("superc");self = get_Sigma_superc(zeta,axis_)
    case("nonsu2");self = get_Sigma_nonsu2(zeta,axis_)
    end select
  end function get_Sigma
  !
  function get_Self(zeta,axis) result(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                     :: zeta
    character(len=*),optional                              :: axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(zeta)) :: self
    character(len=1)                                       :: axis_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_Self"
#endif
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_Self error: not a valid ed_mode"
    case("superc");self = get_Self_superc(zeta,axis_)
    end select
  end function get_Self



  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Calculate the values of the scattering rate and quasiparticle weight
    integer                                       :: ispin,iorb
    real(8)                                       :: wm(2),vm(1)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,2) :: Smats
    complex(8),dimension(1)                       :: Dmats
    real(8)                                       :: iS1,iS2
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_Szr"
#endif
    wm    = pi/beta*(2*arange(1,2)-1)
    vm    = 0d0
    Smats = get_Sigma(xi*wm,'m')
    if(allocated(simp))deallocate(simp)
    if(allocated(zimp))deallocate(zimp)
    allocate(simp(Norb,Nspin))
    allocate(zimp(Norb,Nspin))
    w_ph  = w0_ph
    if(Nph>0)then
       Dmats = get_impD(xi*vm,'m')
       w_ph  = sqrt(-2.d0*w0_ph/Dmats(1)) !renormalized phonon frequency
    endif
    do ispin=1,Nspin
       do iorb=1,Norb
          iS1 = dimag(Smats(ispin,ispin,iorb,iorb,1))
          iS2 = dimag(Smats(ispin,ispin,iorb,iorb,2))
          simp(iorb,ispin) = iS1 - wm(1)*(iS2-iS1)/(wm(2)-wm(1))
          zimp(iorb,ispin) = 1d0/( 1d0 + abs( iS1/wm(1) ))
       enddo
    enddo
  end subroutine get_szr












  !+------------------------------------------------------------------+
  !                    PRINT FUNCTIONS
  !+------------------------------------------------------------------+
  subroutine print_impGmatrix(file)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !This subroutine prints weights and poles of the impurity Green's function by calling :f:func:`write_GFmatrix`. These are stored
    !in a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    if(.not.allocated(impGmatrix))stop "ED_PRINT_IMPGMATRIX ERROR: impGmatrix not allocated!"
    file_="gfmatrix";if(present(file))file_=str(file)
    call write_GFmatrix(impGmatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_impGmatrix

  subroutine print_impDmatrix(file)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !This subroutine prints weights and poles of the  phonon Green's function by calling :f:func:`write_GFmatrix`. These are stored
    !in a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    file_="dfmatrix";if(present(file))file_=str(file)
    call write_GFmatrix(impDmatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_impDmatrix






  subroutine print_impG
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !This subroutine print the impurity Green's function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}G,\mathrm{Re}G]` .
    !One file per Green'sfunction component, with the name
    !
    !  * :code:`"impG_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal G, Matsubara axis
    !  * :code:`"impG_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal G, real frequency axis
    !  * :code:`"impF_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous G, Matsubara axis
    !  * :code:`"impF_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous G, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impFmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impGreal
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impFreal
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG print_impG"
#endif
    call allocate_grids
    select case(ed_mode)
    case ("normal");
       impGmats  = get_impG(dcmplx(0d0,wm(:)),axis='m')
       impGreal  = get_impG(dcmplx(wr(:),eps),axis='r')
       call Gprint_normal("impG",impGmats,'m')
       call Gprint_normal("impG",impGreal,'r')
    case ("superc");
       impGmats = get_impG(dcmplx(0d0,wm(:)),axis='m')
       impFmats = get_impF(dcmplx(0d0,wm(:)),axis='m')
       impGreal = get_impG(dcmplx(wr(:),eps),axis='r')
       impFreal = get_impF(dcmplx(wr(:),eps),axis='r')
       call Gprint_superc("impG",impGmats,'m')
       call Gprint_superc("impF",impFmats,'m')
       call Gprint_superc("impG",impGreal,'r')
       call Gprint_superc("impF",impFreal,'r')
    case ("nonsu2");
       impGmats  = get_impG(dcmplx(0d0,wm(:)),axis='m')
       impGreal  = get_impG(dcmplx(wr(:),eps),axis='r')
       call Gprint_nonsu2("impG",impGmats,'m')
       call Gprint_nonsu2("impG",impGreal,'r')
    case default;stop "ed_print_impG error: ed_mode not valid"
    end select
    call deallocate_grids
  end subroutine print_impG




  

  subroutine print_impD
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Lmats,Lreal
#endif
    !This subroutine print the impurity phonon self-energy on the files
    !  * :code:`"impDph_iw.ed"`  matsubara axis
    !  * :code:`impDph_realw.ed"` real frequency axis
    !
    complex(8),dimension(Lmats) :: impDmats
    complex(8),dimension(Lreal) :: impDreal
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG print_impD"
#endif
    call allocate_grids()
    !Print the impurity functions:
    select case(ed_mode)
    case ("normal","superc","nonsu2")
       impDmats  = get_impD(dcmplx(0d0,wm(:)),axis='m')
       impDreal  = get_impD(dcmplx(wr(:),eps),axis='r')
       call splot("impDph_iw.ed"   ,vm,impDmats(:))
       call splot("impDph_realw.ed",vr,impDreal(:))
    end select
    call deallocate_grids()
  end subroutine Print_ImpD




  

  subroutine print_Sigma
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !This subroutine print the impurity self-energy on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}\Sigma,\mathrm{Re}\Sigma]` .
    !One file per self-energy component, with the name
    !
    !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal self-energy, Matsubara axis
    !  * :code:`"impSigma_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal self-energy, real frequency axis
    !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, Matsubara axis
    !  * :code:`"impSelf_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous self-energy, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impSmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impSAmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impSreal
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impSAreal
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG print_Sigma"
#endif
    call allocate_grids
    select case(ed_mode)
    case ("normal");
       impSmats  = get_Sigma(dcmplx(0d0,wm(:)),axis='m')
       impSreal  = get_Sigma(dcmplx(wr(:),eps),axis='r')
       call Gprint_normal("impSigma",impSmats,'m')
       call Gprint_normal("impSigma",impSreal,'r')
    case ("superc");
       impSmats  = get_Sigma(dcmplx(0d0,wm(:)),axis='m')
       impSreal  = get_Sigma(dcmplx(wr(:),eps),axis='r')
       impSAmats = get_Self(dcmplx(0d0,wm(:)),axis='m')
       impSAreal = get_Self(dcmplx(wr(:),eps),axis='r')
       call Gprint_superc("impSigma",impSmats,'m')
       call Gprint_superc("impSelf",impSAmats,'m')
       call Gprint_superc("impSigma",impSreal,'r')
       call Gprint_superc("impSelf",impSAreal,'r')
    case ("nonsu2");
       impSmats  = get_Sigma(dcmplx(0d0,wm(:)),axis='m')
       impSreal  = get_Sigma(dcmplx(wr(:),eps),axis='r')
       call Gprint_nonsu2("impSigma",impSmats,'m')
       call Gprint_nonsu2("impSigma",impSreal,'r')
    case default;stop "print_Sigma error: ed_mode not valid"
    end select
    call deallocate_grids
  end subroutine print_Sigma










  subroutine print_impG0
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !This subroutine print the non-interacting impurity Green's function on plain text files in the execution folder.
    !The files are formatted like :math:`[\omega,\mathrm{Im}G_{0},\mathrm{Re}G_{0}]` .
    !One file per Green's function component, with the name
    !
    !  * :code:`"impG0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_iw"//reg(ed_file_suffix)//".ed"` normal G, Matsubara axis
    !  * :code:`"impG0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)[str(jspin)]"_realw"//reg(ed_file_suffix)//".ed"` normal G, real frequency axis
    !  * :code:`"impF0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_iw"//reg(ed_file_suffix)//".ed"` anomalous G, Matsubara axis
    !  * :code:`"impF0_l"//str(iorb)[str(jorb)]//_s"//str(ispin)"_realw"//reg(ed_file_suffix)//".ed"` anomalous G, real frequency axis
    !
    !The variable :f:var:`ed_file_suffix` is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation.
    !
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impG0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impF0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impG0real
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: impF0real
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG print_impG0"
#endif
    call allocate_grids()
    select case(ed_mode)
    case ("normal");
       impG0mats  = g0and_bath_function(dcmplx(0d0,wm(:)),axis='mats')
       impG0real  = g0and_bath_function(dcmplx(wr(:),eps),axis='real')
       call Gprint_normal("impG0",impG0mats,'m')
       call Gprint_normal("impG0",impG0real,'r')
    case ("superc");
       impG0mats  = g0and_bath_function(dcmplx(0d0,wm(:)),axis='mats')
       impF0mats  = f0and_bath_function(dcmplx(0d0,wm(:)),axis='mats')
       impG0real  = g0and_bath_function(dcmplx(wr(:),eps),axis='real')
       impF0real  = f0and_bath_function(dcmplx(wr(:),eps),axis='real')
       call Gprint_superc("impG0",impG0mats,'m')
       call Gprint_superc("impF0",impF0mats,'m')
       call Gprint_superc("impG0",impG0real,'r')
       call Gprint_superc("impF0",impF0real,'r')
    case ("nonsu2");
       impG0mats  = g0and_bath_function(dcmplx(0d0,wm(:)),axis='mats')
       impG0real  = g0and_bath_function(dcmplx(wr(:),eps),axis='real')
       call Gprint_nonsu2("impG0",impG0mats,'m')
       call Gprint_nonsu2("impG0",impG0real,'r')
    case default;stop "ed_print_impG0 error: ed_mode not valid"
    end select
    call deallocate_grids()
  end subroutine Print_ImpG0






  subroutine write_szr()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer :: unit,iorb,jorb,ispin
    !Z renormalization constant
    unit = free_unit()
    open(unit,file="Z_info.ed")
    write(unit,"(A1,*(A10,6X))") "# ",&
         ((str(iorb+(ispin-1)*Norb)//"z_"//str(iorb)//"s"//str(ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    !\Sigma scattering rate
    unit = free_unit()
    open(unit,file="Sig_info.ed")
    write(unit,"(A1,*(A10,6X))") "#",&
         ((str(iorb+(ispin-1)*Norb)//"sig_"//str(iorb)//"s"//str(ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    !Z renormalization constant
    unit = free_unit()
    open(unit,file="Z_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))") ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    !\Sigma scattering rate
    unit = free_unit()
    open(unit,file="Sig_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(*(F15.9,1X))") ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    if(Nph>0)then
       unit = free_unit()
       open(unit,file="wph_last"//reg(ed_file_suffix)//".ed")
       write(unit,"(F15.9)") w_ph
       close(unit)
    endif
    !
    if(ed_obs_all)then
       !Z renormalization constant
       unit = free_unit()
       open(unit,file="Z_all"//reg(ed_file_suffix)//".ed",position='append')
       write(unit,"(90(F15.9,1X))") ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       !
       !\Sigma scattering rate
       unit = free_unit()
       open(unit,file="Sig_all"//reg(ed_file_suffix)//".ed",position='append')
       write(unit,"(*(F15.9,1X))") ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       !
       if(Nph>0)then
          unit = free_unit()
          open(unit,file="wph_all"//reg(ed_file_suffix)//".ed",position='append')
          write(unit,"(F15.9)") w_ph
          close(unit)
       end if
    endif
    !
  end subroutine write_szr







  !################################################################
  !################################################################
  !################################################################
  !################################################################







  subroutine Gprint_Normal(file,Self,axis)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*)                 :: file
    complex(8),dimension(:,:,:,:,:)  :: Self
    character(len=1)                 :: axis
    integer                          :: L,i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                :: suffix
    integer,dimension(:),allocatable :: getIorb,getJorb
    integer                          :: totNorb
    L = size(self,5)
    call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],"Gprint_Normal","Self")
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       totNorb=Norb
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
    case ("hybrid","replica","general")             !Diagonal in spin only. Full Orbital structure
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_normal error counting the orbitals"
    !!
    !Print the impurity functions:
    call allocate_grids
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          select case(to_lower(axis))
          case default;stop "Gprint_Normal error: axis not supported"
          case("m");call splot(reg(file)//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,Self(ispin,ispin,iorb,jorb,:))
          case("r");call splot(reg(file)//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Self(ispin,ispin,iorb,jorb,:))
          end select
       enddo
    enddo
    call deallocate_grids
    !
  end subroutine Gprint_Normal





  subroutine Gprint_superc(file,Self,axis)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*)                 :: file
    complex(8),dimension(:,:,:,:,:)  :: Self
    character(len=1)                 :: axis
    integer                          :: i,ispin,unit(4),iorb,jorb,isign
    character(len=20)                :: suffix
    integer,dimension(:),allocatable :: getIorb,getJorb
    integer                          :: totNorb,L
    L = size(self,5)
    call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],"Gprint_Superc","Self")
    !
    select case(bath_type)
    case default
       totNorb=Norb
       allocate(getIorb(Norb),getJorb(Norb))
       l=0
       do iorb=1,Norb
          l=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
    case ("hybrid","replica","general")
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_superc error counting the orbitals"
    !!
    !!PRINT OUT GF:
    call allocate_grids
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          select case(to_lower(axis))
          case default;stop "Gprint_Normal error: axis not supported"
          case("m");call splot(reg(file)//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,Self(ispin,ispin,iorb,jorb,:))
          case("r");call splot(reg(file)//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Self(ispin,ispin,iorb,jorb,:))
          end select
       enddo
    enddo
    call deallocate_grids
    !
  end subroutine Gprint_superc


  subroutine Gprint_nonsu2(file,Self,axis)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*)                 :: file
    complex(8),dimension(:,:,:,:,:)  :: Self
    character(len=1)                 :: axis
    integer                          :: i,isign,unit(2),iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: getIorb,getJorb,getIspin,getJspin
    integer                          :: totNso,totNorb,totNspin,l,io,jo
    character(len=20)                :: suffix
    L = size(self,5)
    call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],"Gprint_nonSU2","Self")
    !
    select case(bath_type)
    case default
       totNorb =Norb
       totNspin=Nspin*Nspin
       totNso  =totNorb*totNspin
       allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
       l=0
       do iorb=1,Norb
          do ispin=1,Nspin
             do jspin=1,Nspin
                l=l+1
                getIorb(l)=iorb
                getIspin(l)=ispin
                getJorb(l)=iorb
                getJspin(l)=jspin
             enddo
          enddo
       enddo
    case ("hybrid","replica","general")
       totNso  = (Norb*Nspin)**2
       allocate(getIorb(totNso),getJorb(totNso),getIspin(totNso),getJspin(totNso))
       l=0
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,Nspin
                do jspin=1,Nspin
                   l=l+1
                   getIorb(l)=iorb
                   getIspin(l)=ispin
                   getJorb(l)=jorb
                   getJspin(l)=jspin
                enddo
             enddo
          enddo
       enddo
    end select
    if(l/=totNso)stop "print_gf_nonsu2 error counting the spin-orbitals"
    !!
    !!PRINT OUT GF:
    call allocate_grids
    do l=1,totNso
       iorb=getIorb(l)
       ispin=getIspin(l)
       jorb=getJorb(l)
       jspin=getJspin(l)
       !
       suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
       select case(to_lower(axis))
       case default;stop "Gprint_Normal error: axis not supported"
       case("m");call splot(reg(file)//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,Self(ispin,jspin,iorb,jorb,:))
       case("r");call splot(reg(file)//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,Self(ispin,jspin,iorb,jorb,:))
       end select
    enddo
    call deallocate_grids
    !
  end subroutine Gprint_nonsu2


end MODULE ED_GREENS_FUNCTIONS










