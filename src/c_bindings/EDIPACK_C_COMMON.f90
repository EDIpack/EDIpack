MODULE EDIPACK_C_COMMON
!:synopsis: Fortran-C bindings: common routines and variables
  USE EDIPACK
  USE SCIFOR
  USE ISO_C_BINDING
  implicit none

#ifdef _WINEQ
  integer(c_int),bind(c, name="has_ineq") :: has_ineq=1 
#else
  integer(c_int),bind(c, name="has_ineq") :: has_ineq=0
#endif


!Dummy interfaces if the real-space DMFT routines are not compiled
#ifndef _WINEQ
  interface ed_finalize_solver
     module procedure :: ed_finalize_solver_dummy
  end interface ed_finalize_solver
  
  interface ed_get_spinChi
     module procedure :: ed_get_chi_dummy_n3
  end interface ed_get_spinChi

  interface ed_get_densChi
     module procedure :: ed_get_chi_dummy_n3
  end interface ed_get_densChi

  interface ed_get_pairChi
     module procedure :: ed_get_chi_dummy_n3
  end interface ed_get_pairChi

  interface ed_get_exctChi
     module procedure :: ed_get_chi_dummy2_n3
  end interface ed_get_exctChi  
  
  public :: ed_finalize_solver
  public :: ed_get_spinChi
  public :: ed_get_densChi
  public :: ed_get_pairChi
  public :: ed_get_exctChi
#endif

contains

  !integer to logical
  function i2l(var_integer) result (var_logical)
    use, intrinsic :: iso_c_binding
    integer        :: var_integer
    logical        :: var_logical   

    if (var_integer == 1) then
       var_logical = .true.
    else
       var_logical = .false.
    endif
  end function i2l

  !logical to integer
  function l2i(var_logical) result (var_integer)
    use, intrinsic :: iso_c_binding
    integer    :: var_integer
    logical    :: var_logical   

    if (var_logical) then
       var_integer = 1
    else
       var_integer = 0
    endif
  end function l2i

  !c string to fortran string
  subroutine c2f(c_str)
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*),intent(IN) :: c_str
    character(len=120), allocatable                 :: f_str
    integer                                         :: length
    integer                                         :: i

    length=0
    f_str=" "
    do
       if (c_str(length+1) == C_NULL_CHAR) exit
       length = length + 1
    end do
    do i = 1, length
       f_str(i:i) = c_str(i)
    enddo
    f_str=trim(f_str)
  end subroutine c2f

  !Get Bath Type
  integer(c_int) function get_bath_type_c() result(bt) bind(c, name='get_bath_type')
    use, intrinsic :: iso_c_binding
    select case(bath_type)
    case("normal");  bt = 1
    case("hybrid");  bt = 2
    case("replica"); bt = 3
    case("general"); bt = 4
    end select
  end function get_bath_type_c

  !Get ED mode
  integer(c_int) function get_ed_mode_c() result(edm) bind(c, name='get_ed_mode')
    use, intrinsic :: iso_c_binding
    select case(ed_mode)
    case("normal");  edm = 1
    case("superc");  edm = 2
    case("nonsu2");  edm = 3
    end select
  end function get_ed_mode_c



!Dummy routines in the case real-space DMFT is not compiled
#ifndef _WINEQ
  subroutine ed_finalize_solver_dummy(Nineq)
    integer                              :: Nineq 
    STOP "Cannot finalize solver for more than one site without R-DMFT support"
  end subroutine ed_finalize_solver_dummy
  
  subroutine ed_get_chi_dummy_n3(self,nlat,axis,z)
    complex(8),dimension(:,:,:,:),intent(inout)   :: self
    integer,intent(in)                            :: nlat
    character(len=*),optional                     :: axis
    complex(8),dimension(:),optional              :: z
    STOP "ed_get_chi: Only single-site available without r-DMFT module"
  end subroutine ed_get_chi_dummy_n3

  subroutine ed_get_chi_dummy2_n3(self,nlat,axis,z)
    complex(8),dimension(:,:,:,:,:),intent(inout) :: self
    integer,intent(in)                            :: nlat
    character(len=*),optional                     :: axis
    complex(8),dimension(:),optional              :: z
    STOP "ed_get_chi: Only single-site available without r-DMFT module"
  end subroutine ed_get_chi_dummy2_n3   
#endif
  

END MODULE EDIPACK_C_COMMON
