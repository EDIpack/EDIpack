MODULE E2I_AUX_FUNX
  !:synopsis: Assortment of auxiliary procedures, inequivalent sites version
  !Hosts a number of auxiliary procedures required in different parts of the code.
  !Specifically, it implements: creation/annihilation fermionic operators, binary decomposition of integer representation of Fock states and setup the local impurity Hamiltonian
  !
  USE EDIPACK
  USE E2I_VARS_GLOBAL
  !
  USE SF_TIMER
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS,  only: arange,linspace  
  implicit none
  private

  !> ED SET HLOC
  interface ed_set_Hloc
     !This subroutine sets the local Hamiltonian of the impurity problem. 
     !The local hamiltonian can have different shapes:
     !
     !   * [|Nlso|, |Nlso| ]: real-space DMFT case, rank-2 array.
     !   * [|Nlat| , |Nso| , |Nso| ]: real-space RDMFT case, rank-3 array.
     !   * [|Nlat| , |Nspin| , |Nspin| , |Norb| , |Norb| ]: real-space RDMFT case, rank-5 array.
     !
     !In the case of real-space DMFT, the number of impurities |Nlat| must be provided.
     !In case :f:var:`ED_MODE`=:code:`SUPERC` a matrix of the same shape containing the
     !anomalous terms can optionally be passed as second arsgument.
     !
     module procedure :: ed_set_Hloc_lattice_N2
     module procedure :: ed_set_Hloc_lattice_superc_N2
     module procedure :: ed_set_Hloc_lattice_N3
     module procedure :: ed_set_Hloc_lattice_superc_N3
     module procedure :: ed_set_Hloc_lattice_N5
     module procedure :: ed_set_Hloc_lattice_superc_N5
  end interface ed_set_Hloc



  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
     module procedure d_nso2nn_l
     module procedure c_nso2nn_l
  end interface so2nn_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
     module procedure d_nn2nso_l
     module procedure c_nn2nso_l
  end interface nn2so_reshape



  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
     module procedure d_nlso2nnn_l
     module procedure c_nlso2nnn_l
  end interface lso2nnn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
     module procedure d_nnn2nlso_l
     module procedure c_nnn2nlso_l
  end interface nnn2lso_reshape


  !SET local Impurity Hamiltonian (excluded the interaction part)
  public :: ed_set_Hloc

  !Internal use
  public :: set_impHloc

  !AUX RESHAPE FUNCTIONS (internal use)
  public :: so2nn_reshape
  public :: nn2so_reshape
  public :: lso2nnn_reshape
  public :: nnn2lso_reshape


contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Himpurity, the local part of the non-interacting Hamiltonian
  !+------------------------------------------------------------------+
  subroutine ed_set_Hloc_lattice_N2(Hloc,Nlat)
    complex(8),dimension(:,:),intent(in) :: Hloc
    integer                              :: Nlat !Number of impurities for real-space DMFT
    integer                              :: ilat
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_Hloc: set impHloc"
#endif
    !
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    allocate(Hloc_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
    !
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],'ed_set_Hloc','Hloc')
    Hloc_ineq  = lso2nnn_reshape(Hloc(1:Nlat*Nspin*Norb,1:Nlat*Nspin*Norb),Nlat,Nspin,Norb)
  end subroutine ed_set_Hloc_lattice_N2
  
  
  subroutine ed_set_Hloc_lattice_superc_N2(Hloc,Hloc_anomalous,Nlat)
    complex(8),dimension(:,:),intent(in) :: Hloc,Hloc_anomalous
    integer                              :: Nlat !Number of impurities for real-space DMFT
    integer                              :: ilat
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_Hloc: set impHloc"
#endif
    !
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    if(allocated(Hloc_anomalous_ineq))deallocate(Hloc_anomalous_ineq)
    allocate(Hloc_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
    allocate(Hloc_anomalous_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_anomalous_ineq=zero
    !
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],'ed_set_Hloc','Hloc')
    call assert_shape(Hloc_anomalous,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],'ed_set_Hloc','Hloc_anomalous')
    !
    Hloc_ineq  = lso2nnn_reshape(Hloc(1:Nlat*Nspin*Norb,1:Nlat*Nspin*Norb),Nlat,Nspin,Norb)
    Hloc_anomalous_ineq  = lso2nnn_reshape(Hloc_anomalous(1:Nlat*Nspin*Norb,1:Nlat*Nspin*Norb),Nlat,Nspin,Norb)
  end subroutine ed_set_Hloc_lattice_superc_N2


  subroutine ed_set_Hloc_lattice_N3(Hloc,Nlat)
    complex(8),dimension(:,:,:),intent(in) :: Hloc
    integer                                :: Nlat,ilat
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_Hloc: set impHloc"
#endif
    !
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    allocate(Hloc_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
    !
    call assert_shape(Hloc,[Nlat,Nspin*Norb,Nspin*Norb],'ed_set_Hloc','Hloc')
    do ilat=1,Nlat
       Hloc_ineq(ilat,:,:,:,:)  = so2nn_reshape(Hloc(ilat,1:Nspin*Norb,1:Nspin*Norb),Nspin,Norb)
    enddo
    !
  end subroutine ed_set_Hloc_lattice_N3
  
  subroutine ed_set_Hloc_lattice_superc_N3(Hloc,Hloc_anomalous,Nlat)
    complex(8),dimension(:,:,:),intent(in) :: Hloc,Hloc_anomalous
    integer                                :: Nlat,ilat
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_Hloc: set impHloc"
#endif
    !
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    if(allocated(Hloc_anomalous_ineq))deallocate(Hloc_anomalous_ineq)
    allocate(Hloc_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
    allocate(Hloc_anomalous_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_anomalous_ineq=zero
    !
    call assert_shape(Hloc,[Nlat,Nspin*Norb,Nspin*Norb],'ed_set_Hloc','Hloc')
    call assert_shape(Hloc_anomalous,[Nlat,Nspin*Norb,Nspin*Norb],'ed_set_Hloc','Hloc_anomalous')
    !
    do ilat=1,Nlat
       Hloc_ineq(ilat,:,:,:,:)  = so2nn_reshape(Hloc(ilat,1:Nspin*Norb,1:Nspin*Norb),Nspin,Norb)
       Hloc_anomalous_ineq(ilat,:,:,:,:)  = so2nn_reshape(Hloc_anomalous(ilat,1:Nspin*Norb,1:Nspin*Norb),Nspin,Norb)
    enddo
    !
  end subroutine ed_set_Hloc_lattice_superc_N3

  subroutine ed_set_Hloc_lattice_N5(Hloc,Nlat)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Hloc
    integer                                    :: Nlat,ilat
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_Hloc: set impHloc"
#endif
    !
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    allocate(Hloc_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
    !
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],'ed_set_Hloc','Hloc')
    Hloc_ineq  = Hloc
  end subroutine ed_set_Hloc_lattice_N5
  
  
  subroutine ed_set_Hloc_lattice_superc_N5(Hloc,Hloc_anomalous,Nlat)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Hloc,Hloc_anomalous
    integer                                    :: Nlat,ilat
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG ed_set_Hloc: set impHloc"
#endif
    !
    if(allocated(Hloc_ineq))deallocate(Hloc_ineq)
    if(allocated(Hloc_anomalous_ineq))deallocate(Hloc_anomalous_ineq)
    allocate(Hloc_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
    allocate(Hloc_anomalous_ineq(Nlat,Nspin,Nspin,Norb,Norb));Hloc_anomalous_ineq=zero
    !
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],'ed_set_Hloc','Hloc')
    call assert_shape(Hloc_anomalous,[Nlat,Nspin,Nspin,Norb,Norb],'ed_set_Hloc','Hloc_anomalous')
    Hloc_ineq  = Hloc
    Hloc_anomalous_ineq  = Hloc_anomalous
  end subroutine ed_set_Hloc_lattice_superc_N5


  subroutine set_impHloc(site)
    integer :: site
    ! if(allocated(impHloc))deallocate(impHloc)
    ! allocate(impHloc(Nspin,Nspin,Norb,Norb));impHloc=zero
    if(.not.allocated(Hloc_ineq))stop "set_impHloc error: called with Hloc_ineq not allocated"
    ! impHloc = Hloc_ineq(site,:,:,:,:)
    call ed_set_Hloc(Hloc_ineq(site,:,:,:,:))
  end subroutine set_impHloc










  !>> ALL THESE CAN BE EXPRESSED IN TERMS OF +DO CONCURRENT; ENDDO
  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix:
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn


  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso







  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix 
  ! _nso2nn   : from [Nso][Nso][L]   to [Nspin][Nspin][Norb][Norb][L]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb][L]   to [Nso][Nso][L]
  !+-----------------------------------------------------------------------------+!
  function d_nso2nn_l(Hso,Nspin,Norb,L) result(Hnn)
    integer                                    :: Nspin,Norb,L
    real(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hnn
    integer                                    :: iorb,ispin,is
    integer                                    :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb,:) = Hso(is,js,:)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn_l
  function c_nso2nn_l(Hso,Nspin,Norb,L) result(Hnn)
    integer                                       :: Nspin,Norb,L
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hso
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hnn
    integer                                       :: iorb,ispin,is
    integer                                       :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb,:) = Hso(is,js,:)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn_l


  function d_nn2nso_l(Hnn,Nspin,Norb,L) result(Hso)
    integer                                    :: Nspin,Norb,L
    real(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hso
    integer                                    :: iorb,ispin,is
    integer                                    :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js,:) = Hnn(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso_l

  function c_nn2nso_l(Hnn,Nspin,Norb,L) result(Hso)
    integer                                       :: Nspin,Norb,L
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hso
    integer                                       :: iorb,ispin,is
    integer                                       :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js,:) = Hnn(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso_l



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix:
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn


  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix:
  ! _nlso2nnn : from [Nlso][Nlso][L] to [Nlat][Nspin][Nspin][Norb][Norb][L]
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb][L] to [Nlso][Nlso][L]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn_l(Hlso,Nlat,Nspin,Norb,L) result(Hnnn)
    integer                                              :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hlso
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hnnn
    integer                                              :: iorb,ispin,ilat,is
    integer                                              :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb,:) = Hlso(is,js,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn_l
  function c_nlso2nnn_l(Hlso,Nlat,Nspin,Norb,L) result(Hnnn)
    integer                                                 :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hnnn
    integer                                                 :: iorb,ispin,ilat,is
    integer                                                 :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb,:) = Hlso(is,js,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn_l


  function d_nnn2nlso_l(Hnnn,Nlat,Nspin,Norb,L) result(Hlso)
    integer                                              :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hlso
    integer                                              :: iorb,ispin,ilat,is
    integer                                              :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js,:) = Hnnn(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso_l

  function c_nnn2nlso_l(Hnnn,Nlat,Nspin,Norb,L) result(Hlso)
    integer                                                 :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hlso
    integer                                                 :: iorb,ispin,ilat,is
    integer                                                 :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js,:) = Hnnn(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso_l

END MODULE E2I_AUX_FUNX
