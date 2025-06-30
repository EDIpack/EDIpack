MODULE ED_BATH_FUNCTIONS
  !:synopsis: Routines for impurity Weiss field and hybridization functions
  !A comprehensive set of procedures to evaluate the non-interacting impurity Green's functions :math:`\hat{G}^{\rm And}` and hybridizations :math:`\hat{F}^{\rm And}` in the complex frequency domain given the :f:var:`effective_bath` instance.
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy,str,to_lower
  USE SF_LINALG, only: eye,inv,diag,zeye,inv_her,kron,det
  USE SF_SPIN, only: pauli_sigma_z,pauli_sigma_x,pauli_tau_z,pauli_tau_x
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH_AUX
  USE ED_BATH_REPLICA
  USE ED_BATH_DIM
  USE ED_BATH_USER
  USE ED_BATH_DMFT
  implicit none

  private


  !
  !\DELTA, THE HYBRIDIZATION FUNCTION
  interface delta_bath_function
     !
     ! Evaluates the normal hybridization function :math:`\Delta(x)`.
     !
     ! Output:
     !   * :f:var:`delta` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure delta_bath_array
  end interface delta_bath_function

  interface Fdelta_bath_function
     !
     ! Evaluates the anomalous hybridization function :math:`\Theta(x)`.
     !
     ! Output:
     !   * :f:var:`fdelta` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure Fdelta_bath_array
  end interface Fdelta_bath_function



  !NON-INTERACTING GREEN'S FUNCTION 
  interface g0and_bath_function
     !
     ! Evaluates the normal non-interacting Green's function :math:`G_0(x)`.
     !
     ! Output:
     !   * :f:var:`g0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure g0and_bath_array
  end interface g0and_bath_function

  interface f0and_bath_function
     !
     ! Evaluates the anomalous non-interacting Green's function :math:`F_0(x)`.
     !
     ! Output:
     !   * :f:var:`f0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure f0and_bath_array
  end interface f0and_bath_function



  !INVERSE NON-INTERACTING GREEN'S FUNCTION 
  interface invg0_bath_function
     !
     ! Evaluates the inverse of the normal non-interacting Green's function :math:`G^{-1}_0(x)`.
     !
     ! Output:
     !   * :f:var:`g0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure invg0_bath_array
  end interface invg0_bath_function

  interface invf0_bath_function
     !
     ! Evaluates the inverse of the anomalous non-interacting Green's function :math:`F^{-1}_0(x)`.
     !
     ! Output:
     !   * :f:var:`f0and` : complex rank-5 array with dimension [ |Nspin| , |Nspin| , |Norb| , |Norb| , :code:`size(x)` ]
     !
     module procedure invf0_bath_array
  end interface invf0_bath_function


  !Build Gand/Delta from a user bath
  interface ed_get_g0and
     !| This subroutine returns to the user the normal non-interacting Green's function :math:`G_0(x)` and
     ! the anomalous non-interacting Green's function :math:`F_0(x)` on a given set of frequencies. It does so
     ! by calling :f:func:`g0and_bath_function` and :f:func:`g0and_bath_function`.
     !
     !The non-interacting Green's function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(x)`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(x)`]
     !
     module procedure :: ed_get_g0and_n2
     module procedure :: ed_get_g0and_n4
  end interface ed_get_g0and

  interface ed_get_delta
     !| This subroutine returns to the user the normal hybridization function :math:`\Delta(x)` and
     ! the anomalous hybridization function :math:`\Theta(x)` on a given set of frequencies. It does so
     ! by calling :f:func:`delta_bath_function` and :f:func:`fdelta_bath_function`.
     !
     !The hybridization function is an array having the following possible dimensions:
     ! 
     !  * [:f:var:`nspin` :math:`\cdot` :f:var:`norb`, :f:var:`nspin`:math:`\cdot`:f:var:`norb`, :code:`size(x)`]  
     !  * [:f:var:`nspin`, :f:var:`nspin`, :f:var:`norb`, :f:var:`norb`, :code:`size(x)`]
     !
     module procedure :: ed_get_delta_n2
     module procedure :: ed_get_delta_n4
  end interface ed_get_delta





  public :: delta_bath_function
  public :: fdelta_bath_function
  public :: g0and_bath_function
  public :: f0and_bath_function
  public :: invg0_bath_function
  public :: invf0_bath_function
  !
  public :: ed_get_g0and
  public :: ed_get_delta





contains



  !##################################################################
  !             DELTA/ FDELTA
  !##################################################################
#if __INTEL_COMPILER
#include "delta_functions/delta_normal.f90"
#include "delta_functions/delta_hybrid.f90"
#include "delta_functions/delta_replica.f90"
#include "delta_functions/delta_general.f90"
#else
  include "delta_functions/delta_normal.f90"
  include "delta_functions/delta_hybrid.f90"
  include "delta_functions/delta_replica.f90"
  include "delta_functions/delta_general.f90"
#endif

  function delta_bath_array(x,axis) result(Delta)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                                :: x          !complex  array for the frequency
    character(len=*),optional                                         :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Delta
    character(len=1)                                                  :: axis_
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    Delta=zero
    !
    if(Nbath.eq.0)return
    !
    if(Nbath>0 .and. .not.dmft_bath%status)stop "delta_bath_array error: dmft_bath not allocated"
    !
    select case(bath_type)
    case default;stop "delta_bath_array error: bath_type not supported"
    case("normal")  ;Delta = delta_bath_array_normal(x,axis_)
    case("hybrid")  ;Delta = delta_bath_array_hybrid(x,axis_)
    case("replica") ;Delta = delta_bath_array_replica(x,axis_)
    case("general") ;Delta = delta_bath_array_general(x,axis_)
    end select
    !
  end function delta_bath_array


#if __INTEL_COMPILER
#include "delta_functions/fdelta_normal.f90"
#include "delta_functions/fdelta_hybrid.f90"
#include "delta_functions/fdelta_replica.f90"
#include "delta_functions/fdelta_general.f90"
#else
  include "delta_functions/fdelta_normal.f90"
  include "delta_functions/fdelta_hybrid.f90"
  include "delta_functions/fdelta_replica.f90"
  include "delta_functions/fdelta_general.f90"
#endif

  function fdelta_bath_array(x,axis) result(Fdelta)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                                :: x          !complex  array for the frequency
    character(len=*),optional                                         :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Fdelta
    character(len=1)                                                  :: axis_
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    Fdelta=zero
    !
    if(Nbath.eq.0)return
    !
    if(Nbath>0 .and. .not.dmft_bath%status)stop "fdelta_bath_array error: dmft_bath not allocated"
    !
    if(ed_mode/="superc")then
       print*,"fdelta_bath_array WARNING: called with ed_mode != superc. Return zero"
       call sleep(2000)
       return
    endif
    !
    select case(bath_type)
    case default;stop "fdelta_bath_array error: bath_type not supported"
    case("normal")  ;FDelta = fdelta_bath_array_normal(x,axis_)
    case("hybrid")  ;FDelta = fdelta_bath_array_hybrid(x,axis_)
    case("replica") ;FDelta = fdelta_bath_array_replica(x,axis_)
    case("general") ;FDelta = fdelta_bath_array_general(x,axis_)
    end select
    !
  end function fdelta_bath_array




  !##################################################################
  !##################################################################
  !##################################################################





  !##################################################################
  !             G0and / F0and
  !##################################################################
#if __INTEL_COMPILER
  #include "g0and_functions/g0and_normal.f90"
  #include "g0and_functions/g0and_hyrege.f90"
#else
  include "g0and_functions/g0and_normal.f90"
  include "g0and_functions/g0and_hyrege.f90"
#endif

  
  function g0and_bath_array(x,axis) result(G0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
    character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    character(len=1)                                    :: axis_
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    G0and = zero
    !
    !
    if(Nbath>0 .and. .not.dmft_bath%status)stop "g0and_bath_array error: dmft_bath not allocated"
    !
    select case(bath_type)
    case default;stop "G0and_bath_array error: bath_type not supported"
    case("normal")
       G0and = G0and_bath_array_normal(x,axis_)
    case("hybrid","replica","general")
       G0and = G0and_bath_array_hyrege(x,axis_)
    end select
    !
  end function g0and_bath_array


#if __INTEL_COMPILER
  #include "g0and_functions/f0and_normal.f90"
  #include "g0and_functions/f0and_hyrege.f90"
#else
  include "g0and_functions/f0and_normal.f90"
  include "g0and_functions/f0and_hyrege.f90"
#endif

  
  function f0and_bath_array(x,axis) result(F0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x !complex  array for the frequency
    character(len=*),optional                           :: axis!string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    character(len=1)                                    :: axis_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    F0and=zero
    !
    !
    if(Nbath>0 .and. .not.dmft_bath%status)stop "f0and_bath_array error: dmft_bath not allocated"
    !
    if(ed_mode/="superc")then
       print*,"f0and_bath_array WARNING: called with ed_mode != superc. Return zero"
       call sleep(2000)
       return
    endif
    !
    select case(bath_type)
    case default;stop "f0and_bath_array error: bath_type not supported"
    case("normal")
       F0and = f0and_bath_array_normal(x,axis_)
    case("hybrid","replica","general")
       F0and = f0and_bath_array_hyrege(x,axis_)
    end select
    !
  end function f0and_bath_array



  !##################################################################
  !##################################################################
  !##################################################################





  !##################################################################
  !              G0^-1 / F0^-1
  !##################################################################  
#if __INTEL_COMPILER
  #include "invg0_functions/invg0_normal.f90"
  #include "invg0_functions/invg0_hyrege.f90"
#else
  include "invg0_functions/invg0_normal.f90"
  include "invg0_functions/invg0_hyrege.f90"
#endif

  function invg0_bath_array(x,axis) result(G0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
    character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and    
    character(len=1)                                    :: axis_
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    G0and = zero
    !
    !
    if(Nbath>0 .and. .not.dmft_bath%status)stop "invg0_bath_array error: dmft_bath not allocated"
    !
    select case(bath_type)
    case default;stop "invg0_bath_array error: bath_type not supported"
    case("normal")
       G0and = invg0_bath_array_normal(x,axis_)
    case("hybrid","replica","general")
       G0and = invg0_bath_array_hyrege(x,axis_)
    end select
  end function invg0_bath_array


#if __INTEL_COMPILER
  #include "invg0_functions/invf0_normal.f90"
  #include "invg0_functions/invf0_hyrege.f90"
#else
  include "invg0_functions/invf0_normal.f90"
  include "invg0_functions/invf0_hyrege.f90"
#endif
  
  function invf0_bath_array(x,axis) result(F0and)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
    character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
    character(len=1)                                    :: axis_
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    !
    F0and=zero
    !
    !
    if(Nbath>0 .and. .not.dmft_bath%status)stop "invf0_bath_array error: dmft_bath not allocated"
    !
    if(ed_mode/="superc")then
       print*,"f0and_bath_array WARNING: called with ed_mode != superc. Return zero"
       call sleep(2000)
       return
    endif
    !
    select case(bath_type)
    case default;stop "f0and_bath_array error: bath_type not supported"
    case("normal")
       F0and = invf0_bath_array_normal(x,axis_)
    case("hybrid","replica","general")
       F0and = invf0_bath_array_hyrege(x,axis_)
    end select
    !
  end function invf0_bath_array



  !##################################################################
  !##################################################################
  !##################################################################







  !##################################################################
  !                bath --> G0and/F0and
  !##################################################################
  subroutine ed_get_g0and_n2(x,bath_,G0and,axis,type)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                  :: x !complex array of frequencies
    real(8),dimension(:)                                :: bath_ !user-accessible bath array
    complex(8),dimension(:,:,:)                         :: G0and !non-interacting Green's function
    character(len=*),optional                           :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis 
    character(len=*),optional                           :: type !string indicating the desired function, :code:`'n'` for normal (default), :code:`'a'` for anomalous
    !
    logical                                             :: check
    character(len=1)                                    :: axis_,type_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: g0
    integer                                             :: L
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    type_='n';if(present(type))type_=trim(type)
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    select case(type_)
    case default;stop "ed_get_g0and ERROR: type is wrong: either Normal or Anomalous"
    case ('n','N')
       g0 = g0and_bath_function(x)
    case('a','A')
       g0 = f0and_bath_function(x)
    end select
    call deallocate_dmft_bath()
    !
    L=size(x)
    call assert_shape(g0and,[Nspin*Norb,Nspin*Norb,L],'ed_get_g0and','g0and')
    g0and = nn2so_reshape(g0,Nspin,Norb,L)
  end subroutine ed_get_g0and_n2

  subroutine ed_get_g0and_n4(x,bath_,G0and,axis,type)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                  :: x
    real(8),dimension(:)                                :: bath_
    complex(8),dimension(:,:,:,:,:)                     :: G0and
    character(len=*),optional                           :: axis
    character(len=*),optional                           :: type
    !
    logical                                             :: check
    character(len=1)                                    :: axis_,type_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: g0
    integer                                             :: L
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    type_='n';if(present(type))type_=trim(type)
    check= check_bath_dimension(bath_)
    if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    select case(type_)
    case default;stop "ed_get_g0and ERROR: type is wrong: either Normal or Anomalous"
    case ('n','N')
       g0 = g0and_bath_function(x)
    case('a','A')
       g0 = f0and_bath_function(x)
    end select
    call deallocate_dmft_bath()
    !
    L=size(x)
    call assert_shape(g0and,[Nspin,Nspin,Norb,Norb,L],'ed_get_g0and','g0and')
    g0and = g0
  end subroutine ed_get_g0and_n4



  !##################################################################
  !##################################################################
  !##################################################################



  !##################################################################
  !                bath --> Delta / FDelta
  !##################################################################  
  subroutine ed_get_delta_n2(x,bath_,delta,axis,type)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                  :: x !complex array of frequencies
    real(8),dimension(:)                                :: bath_ !user-accessible bath array
    complex(8),dimension(:,:,:)                            :: delta !hybridization function
    character(len=*),optional                           :: axis !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis 
    character(len=*),optional                           :: type !string indicating the desired function, :code:`'n'` for normal (default), :code:`'a'` for anomalous
    !
    logical                                             :: check
    character(len=1)                                    :: axis_,type_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: d0
    integer                                             :: L
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    type_='n';if(present(type))type_=trim(type)
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    select case(type_)
    case default;stop "ed_get_delta ERROR: type is wrong: either Normal or Anomalous"
    case ('n','N')
       d0 = delta_bath_function(x,axis_)
    case('a','A')
       d0 = fdelta_bath_function(x,axis_)
    end select
    call deallocate_dmft_bath()
    !
    L=size(x)
    call assert_shape(delta,[Nspin*Norb,Nspin*Norb,L],'ed_get_delta','delta')
    delta = nn2so_reshape(d0,Nspin,Norb,L)
  end subroutine ed_get_delta_n2

  subroutine ed_get_delta_n4(x,bath_,delta,axis,type)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:),intent(in)                  :: x
    real(8),dimension(:)                                :: bath_
    complex(8),dimension(:,:,:,:,:)                            :: delta
    character(len=*),optional                           :: axis
    character(len=*),optional                           :: type
    !
    logical                                             :: check
    character(len=1)                                    :: axis_,type_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: d0
    integer                                             :: L
    !
    axis_="m";if(present(axis))axis_=str(to_lower(axis))
    type_='n';if(present(type))type_=trim(type)
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    select case(type_)
    case default;stop "ed_get_delta ERROR: type is wrong: either Normal or Anomalous"
    case ('n','N')
       d0 = delta_bath_function(x,axis_)
    case('a','A')
       d0 = fdelta_bath_function(x,axis_)
    end select
    call deallocate_dmft_bath()
    !
    L=size(x)
    call assert_shape(delta,[Nspin,Nspin,Norb,Norb,L],'ed_get_delta','delta')
    delta = d0
  end subroutine ed_get_delta_n4







  !##################################################################
  !##################################################################
  !##################################################################

  

  function zeta_superc(x,mu,axis) result(zeta)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    complex(8),dimension(:)                        :: x
    real(8)                                        :: mu
    character(len=*)                               :: axis
    complex(8),dimension(2*Nspin*Norb,size(x))     :: zeta
    integer                                        :: iorb,N,L
    N = Nspin*Norb
    L = size(x)
    do concurrent(iorb=1:N)
       select case(axis)
       case default
          zeta(iorb ,1:L)   = x(1:L) + mu
          zeta(N+iorb,1:L)  = x(1:L) - mu
       case ('real')
          zeta(iorb ,1:L)   = x(1:L) + mu
          zeta(N+iorb,1:L) = -conjg(x(L:1:-1) + mu)
       end select
    enddo
  end function zeta_superc


END MODULE ED_BATH_FUNCTIONS
