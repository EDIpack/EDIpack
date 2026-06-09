MODULE ED_CHI_FUNCTIONS
  !:synopsis: Routines for susceptibilities calculation
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,to_lower,splot
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_SETUP
  USE ED_AUX_FUNX
  !
  USE ED_CHI_SPIN
  USE ED_CHI_DENS
  USE ED_CHI_PAIR
  USE ED_CHI_EXCT
  !
  implicit none
  private 

  public :: buildChi_impurity

  public :: get_spinChi
  public :: get_densChi
  public :: get_pairChi
  public :: get_exctChi


  

contains


  subroutine buildChi_impurity()
    ! 
    ! Build the quantum impurity electrons Susceptibilities :math:`\hat{\chi}` calling the correct procedure according to the value of :f:var:`ed_mode`.
    ! Write the results on file.
    !
    ! * :code:`normal` : :f:func:`build_spinChi_normal`, :f:func:`build_densChi_normal`, :f:func:`build_pairChi_normal`, :f:func:`build_exctChi_normal`
    ! * :code:`superc` : unavailable
    ! * :code:`nonsu2` : unavailable
    !
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
#ifdef _DEBUG
    if(any([chispin_flag,chidens_flag,chipair_flag,chiexct_flag]))&
         write(Logfile,"(A)")"DEBUG buildChi_impurity: build susceptibilities Chi"
#ifdef _CMPLX_NORMAL
    write(Logfile,"(A)")"WARNING: the calculation of susceptibilities with CMPLX_NORMAL activate is experimental and not tested!"
#endif
#endif
    !
    call deallocate_GFmatrix(spinChimatrix)
    call deallocate_GFmatrix(densChimatrix)
    call deallocate_GFmatrix(pairChimatrix)
    call deallocate_GFmatrix(exctChimatrix)
    !
    select case(ed_mode)
    case default;return
    case("normal")
       !BUILD SPIN SUSCEPTIBILITY
       if(chispin_flag)then
          call build_spinChi_normal()
          call print_spinChiMatrix()
          if(ed_print_chispin)call print_spinChi()
       endif
       !BUILD CHARGE SUSCEPTIBILITY
       if(chidens_flag)then
          call build_densChi_normal()
          call print_densChiMatrix()
          if(ed_print_chidens)call print_densChi()
       endif
       !BUILD PAIR SUSCEPTIBILITY
       if(chipair_flag)then
          call build_pairChi_normal()
          call print_pairChiMatrix()
          if(ed_print_chipair)call print_pairChi()
       endif
       !BUILD EXCITON SUSCEPTIBILITY
       if(chiexct_flag)then
          call build_exctChi_normal()
          call print_exctChiMatrix()
          if(ed_print_chiexct)call print_exctChi()
       endif
    end select
  end subroutine buildChi_impurity




  !+------------------------------------------------------------------+
  !+------------------------------------------------------------------+

  function get_spinChi(zeta,axis) result(self)
  !Returns the spin susceptibility for a provided values array
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)         :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                  :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(Norb,Norb,size(zeta)) :: self !The susceptibility matrix
    character(len=1)                           :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_spinChi error: not a valid ed_mode"
    case("normal");self = get_spinChi_normal(zeta,axis_)
    end select
  end function get_spinChi


  function get_densChi(zeta,axis) result(self)
  !Returns the charge susceptibility for a provided values array
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)         :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                  :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(Norb,Norb,size(zeta)) :: self !The susceptibility matrix
    character(len=1)                           :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_densChi error: not a valid ed_mode"
    case("normal");self = get_densChi_normal(zeta,axis_)
    end select
  end function get_densChi


  function get_pairChi(zeta,axis) result(self)
  !Returns the pair susceptibility for a provided values array
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)         :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                  :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(Norb,Norb,size(zeta)) :: self !The susceptibility matrix
    character(len=1)                           :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_pairChi error: not a valid ed_mode"
    case("normal");self = get_pairChi_normal(zeta,axis_)
    end select
  end function get_pairChi


  function get_exctChi(zeta,axis) result(self)
  !Returns the excitonic susceptibility for a provided values array
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)           :: zeta !Array of frequencies or imaginary times
    character(len=*),optional                    :: axis !Axis: can be :code:`m` for Matsubara, :code:`r` for real, :code:`t` for imaginary time
    complex(8),dimension(3,Norb,Norb,size(zeta)) :: self !The susceptibility matrix
    character(len=1)                             :: axis_
    axis_ = 'm' ; if(present(axis))axis_ = axis(1:1)
    select case(ed_mode)
    case default  ;stop "get_exctChi error: not a valid ed_mode"
    case("normal");self = get_exctChi_normal(zeta,axis_)
    end select
  end function get_exctChi








  !+------------------------------------------------------------------+
  !                    PRINT FUNCTIONS
  !+------------------------------------------------------------------+
  subroutine print_spinChimatrix(file)
    !This subroutine prints weights and poles of the impurity spin susceptibility function by calling :f:func:`write_GFmatrix`. These are stored
    !in a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    if(.not.allocated(spinChiMatrix))stop "ED_PRINT_SPINCHIMATRIX ERROR: spinChimatrix not allocated!"
    file_="spinchimatrix";if(present(file))file_=str(file)
    call write_GFmatrix(spinChiMatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_spinChimatrix


  subroutine print_densChimatrix(file)
    !This subroutine prints weights and poles of the impurity charge density susceptibility function by calling :f:func:`write_GFmatrix`. These are stored
    !in a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    if(.not.allocated(densChiMatrix))stop "ED_PRINT_DENSCHIMATRIX ERROR: densChimatrix not allocated!"
    file_="denschimatrix";if(present(file))file_=str(file)
    call write_GFmatrix(densChiMatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_densChimatrix


  subroutine print_pairChimatrix(file)
    !This subroutine prints weights and poles of the impurity pair susceptibility function by calling :f:func:`write_GFmatrix`. These are stored
    !in a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    if(.not.allocated(pairChiMatrix))stop "ED_PRINT_PAIRCHIMATRIX ERROR: pairChimatrix not allocated!"
    file_="pairchimatrix";if(present(file))file_=str(file)
    call write_GFmatrix(pairChiMatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_pairChimatrix


  subroutine print_exctChimatrix(file)
    !This subroutine prints weights and poles of the impurity exciton susceptibilities function by calling :f:func:`write_GFmatrix`. These are stored
    !in a file named :code:`"file"//str(ed_file_suffix)//.restart"` taking into account the value of the global variable :f:var:`ed_file_suffix` ,
    !which is :code:`"_ineq_Nineq"` padded with 4 zeros in the case of inequivalent sites, as per documentation
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    character(len=*),optional :: file !filename prefix (default :code:`gfmatrix`)
    character(len=256)        :: file_
    if(.not.allocated(exctChiMatrix))stop "ED_PRINT_EXCTCHIMATRIX ERROR: exctChimatrix not allocated!"
    file_="exctchimatrix";if(present(file))file_=str(file)
    call write_GFmatrix(exctChiMatrix,str(file_)//str(ed_file_suffix)//".restart")
  end subroutine print_exctChimatrix



















  ! SPIN-SPIN
  subroutine print_spinChi()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(Norb,Norb,Lmats) :: Cmats
    complex(8),dimension(Norb,Norb,Lreal) :: Creal
    complex(8),dimension(Norb,Norb,Ltau)  :: Ctau
    character(len=1)                      :: axis
    integer                               :: L,i,j,iorb,jorb,ispin,isign
    character(len=20)                     :: suffix
    !
    call allocate_grids
    !
    Cmats = get_spinChi(dcmplx(0d0,vm),axis='m')
    Creal = get_spinChi(dcmplx(vr,eps),axis='r')
    Ctau  = get_spinChi(dcmplx(tau,0d0),axis='t')
    !
    do iorb=1,Norb
       do jorb=1,Norb
          suffix="_l"//str(iorb)//str(jorb)
          call splot("spinChi"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Cmats(iorb,jorb,:))
          call splot("spinChi"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Creal(iorb,jorb,:))
          call splot("spinChi"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,dreal(Ctau(iorb,jorb,:)))
       enddo
    enddo
    call deallocate_grids
  end subroutine print_spinChi



  ! DENSITY-DENSITY
  subroutine print_densChi()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(Norb,Norb,Lmats) :: Cmats
    complex(8),dimension(Norb,Norb,Lreal) :: Creal
    complex(8),dimension(Norb,Norb,Ltau)  :: Ctau
    character(len=1)                      :: axis
    integer                               :: L,i,j,iorb,jorb,idens,isign
    character(len=20)                     :: suffix
    call allocate_grids
    !
    Cmats = get_densChi(dcmplx(0d0,vm),axis='m')
    Creal = get_densChi(dcmplx(vr,eps),axis='r')
    Ctau  = get_densChi(dcmplx(tau,0d0),axis='t')
    !
    do iorb=1,Norb
       do jorb=1,Norb
          suffix="_l"//str(iorb)//str(jorb)
          call splot("densChi"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Cmats(iorb,jorb,:))
          call splot("densChi"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Creal(iorb,jorb,:))
          call splot("densChi"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,dreal(Ctau(iorb,jorb,:)))
       enddo
    enddo
    call deallocate_grids
  end subroutine print_densChi


  ! PAIR-PAIR
  subroutine print_pairChi()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(Norb,Norb,Lmats) :: Cmats
    complex(8),dimension(Norb,Norb,Lreal) :: Creal
    complex(8),dimension(Norb,Norb,Ltau)  :: Ctau
    character(len=1)                      :: axis
    integer                               :: L,i,j,iorb,jorb,ipair,isign
    character(len=20)                     :: suffix
    call allocate_grids
    !
    Cmats = get_pairChi(dcmplx(0d0,vm),axis='m')
    Creal = get_pairChi(dcmplx(vr,eps),axis='r')
    Ctau  = get_pairChi(dcmplx(tau,0d0),axis='t')
    !
    do iorb=1,Norb
       do jorb=1,Norb
          suffix="_l"//str(iorb)//str(jorb)
          call splot("pairChi"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Cmats(iorb,jorb,:))
          call splot("pairChi"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Creal(iorb,jorb,:))
          call splot("pairChi"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,dreal(Ctau(iorb,jorb,:)))
       enddo
    enddo
    call deallocate_grids
  end subroutine print_pairChi



  ! EXCITON
  subroutine print_exctChi()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(3,Norb,Norb,Lmats) :: Cmats
    complex(8),dimension(3,Norb,Norb,Lreal) :: Creal
    complex(8),dimension(3,Norb,Norb,Ltau)  :: Ctau
    character(len=1)                      :: axis
    integer                               :: L,i,j,iorb,jorb,iexct,isign
    character(len=20)                     :: suffix
    call allocate_grids
    !
    Cmats = get_exctChi(dcmplx(0d0,vm),axis='m')
    Creal = get_exctChi(dcmplx(vr,eps),axis='r')
    Ctau  = get_exctChi(dcmplx(tau,0d0),axis='t')
    !
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          suffix="_l"//str(iorb)//str(jorb)
          call splot("exctChi_singlet"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Cmats(1,iorb,jorb,:))
          call splot("exctChi_tripletXY"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Cmats(2,iorb,jorb,:))
          call splot("exctChi_tripletZ"//str(suffix)//"_iv"//reg(ed_file_suffix)//".ed",vm,Cmats(3,iorb,jorb,:))
          !
          call splot("exctChi_singlet"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Creal(1,iorb,jorb,:))
          call splot("exctChi_tripletXY"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Creal(2,iorb,jorb,:))
          call splot("exctChi_tripletZ"//str(suffix)//"_realw"//reg(ed_file_suffix)//".ed",vr,Creal(3,iorb,jorb,:))
          !
          call splot("exctChi_singlet"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,dreal(Ctau(1,iorb,jorb,:)))
          call splot("exctChi_tripletXY"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,dreal(Ctau(2,iorb,jorb,:)))
          call splot("exctChi_tripletZ"//str(suffix)//"_tau"//reg(ed_file_suffix)//".ed",tau,dreal(Ctau(3,iorb,jorb,:)))
       enddo
    enddo
    call deallocate_grids
  end subroutine print_exctChi






end MODULE ED_CHI_FUNCTIONS
