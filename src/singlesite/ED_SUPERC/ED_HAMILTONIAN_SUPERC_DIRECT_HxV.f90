! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_SUPERC_DIRECT_HxV
  !:synopsis: Routines for direct matrix-vector product, :code:`SUPERC` case
  USE ED_HAMILTONIAN_SUPERC_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_superc_main
#ifdef _MPI
  public  :: directMatVec_MPI_superc_main
#endif




contains



  subroutine directMatVec_superc_main(Nloc,vin,Hv)
    !
    ! Serial version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm.
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to the input vector using serial algorithm.  
    !
    integer                                                         :: Nloc  !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc)                                      :: vin   !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                                      :: Hv    !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    integer                                                         :: i_el, j_el, iph, jj
    integer                                                         :: isector
    integer,dimension(Nlevels)                                      :: ib
    integer,dimension(Ns)                                           :: ibup,ibdw  
    real(8),dimension(Norb)                                         :: nup,ndw
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                                         :: first_state,last_state
    integer                                                         :: first_state_up,last_state_up
    integer                                                         :: first_state_dw,last_state_dw
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim   = Hsector%Dim
    DimEl = Hsector%DimEl
    !
    Hbath_tmp = zero
    !
    if(Nloc/=dim)stop "directMatVec_cc ERROR: Nloc /= dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0       
       do ibath=1,Nbath
          do ispin=1,Nspin             
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
             do iorb=1,Nfoo
                bath_diag(ispin,iorb,ibath)=dmft_bath%e(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    case ("replica")
       !H_p
       allocate(bath_diag(Nnambu*Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nnambu*Nspin
             do iorb=1,Norb
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
       !V_p
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
             enddo
          enddo
       enddo
    case ("general")
       !H_p
       allocate(bath_diag(Nnambu*Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nnambu*Nspin
             do iorb=1,Norb
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
       !V_p
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%vg(iorb+(ispin-1)*Norb)
             enddo
          enddo
       enddo
    end select

    !Check spin-singlet pairing
    if(bath_type=="replica" .or. bath_type=="general")then
      do ibath=1,Nbath
        do ispin=1,Nnambu*Nspin
          if(any(abs(Hbath_tmp(ispin,1+Nnambu*Nspin-ispin,:,:,ibath)-transpose(Hbath_tmp(ispin,1+Nnambu*Nspin-ispin,:,:,ibath)))>1d-10))then
            STOP "anomalous Hbath is not symmetric. Only spin-singlet s-wave pairing is allowed in SUPERC."
          endif
        enddo
      enddo
    endif
      
  
    !
    Hv=zero
    states: do i=1,Dim
       i_el = mod(i-1,DimEl) +1
       iph  = (i-1)/DimEl +1
       m    = Hsector%H(1)%map(i_el)
       ib   = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "direct/HxVimp_bath.f90"
       !
       if(DimPh>1)then
          !PHONON TERMS
          include "direct/HxV_ph.f90"
          !
          !ELECTRON-PHONON INTERACTION
          include "direct/HxV_eph.f90"
       endif
    enddo states
    !
  end subroutine directMatVec_superc_main



#ifdef _MPI
  subroutine directMatVec_MPI_superc_main(Nloc,v,Hv)
    !
    ! MPI parallel version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm.
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to a part of the vector own by the thread using parallel algorithm.  
    !
    integer                                                         :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)                                      :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                                      :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`

    integer                                                         :: i_el, j_el, iph, jj

    integer                                                         :: N
    complex(8),dimension(:),allocatable                             :: vin
    integer,allocatable,dimension(:)                                :: Counts,Offset
    integer                                                         :: isector
    integer,dimension(Nlevels)                                      :: ib
    integer,dimension(Ns)                                           :: ibup,ibdw  
    real(8),dimension(Norb)                                         :: nup,ndw
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                                         :: first_state,last_state
    integer                                                         :: first_state_up,last_state_up
    integer                                                         :: first_state_dw,last_state_dw
    integer                                                         :: mpiIerr
    !
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim   = Hsector%Dim
    DimEl = Hsector%DimEl
    !
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0       
       do ibath=1,Nbath
          do ispin=1,Nspin             
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
             do iorb=1,Nfoo
                bath_diag(ispin,iorb,ibath)=dmft_bath%e(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    case ("replica")
       !H_p
       allocate(bath_diag(Nnambu*Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nnambu*Nspin
             do iorb=1,Norb
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
       !V_p
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
             enddo
          enddo
       enddo
    case ("general")
       !H_p
       allocate(bath_diag(Nnambu*Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nnambu*Nspin
             do iorb=1,Norb
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
       !V_p
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%vg(iorb+(ispin-1)*Norb)
             enddo
          enddo
       enddo
    end select
    !
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)); vin  = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    Hv=zero
    !
    states: do i=MpiIstart,MpiIend
       i_el = mod(i-1,DimEl) +1
       iph  = (i-1)/DimEl +1
       m  = Hsector%H(1)%map(i_el)
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "direct/HxVimp_bath.f90"
       !
       if(DimPh>1)then
          !PHONON TERMS
          include "direct/HxV_ph.f90"
          !
          !ELECTRON-PHONON INTERACTION
          include "direct/HxV_eph.f90"
       endif
    enddo states
    !
  end subroutine directMatVec_MPI_superc_main
#endif

end MODULE ED_HAMILTONIAN_SUPERC_DIRECT_HXV






