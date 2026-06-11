MODULE ED_HAMILTONIAN_NORMAL_DIRECT_HxV
  !:synopsis: Routines for direct matrix-vector product, :code:`NORMAL` case
  !Constructs and applies on-the-fly each term of the sector Hamiltonian to the input
  !vector :math:`\vec{w} = H\times \vec{v}` in a Arpack/Lanczos framework.
  USE ED_HAMILTONIAN_NORMAL_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_normal_main
  public  :: directMatVec_normal_orbs
#ifdef _MPI
  public  :: directMatVec_MPI_normal_main
  public  :: directMatVec_MPI_normal_orbs
#endif



contains


  subroutine directMatVec_normal_main(Nloc,vin,Hv)
    !
    ! Serial version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True` 
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to the input vector using serial algorithm.  
    !
    integer                                           :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
#ifdef _CMPLX_NORMAL
    complex(8),dimension(Nloc)                        :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                        :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8),dimension(:),allocatable               :: vt,Hvt
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
#else
    real(8),dimension(Nloc)                           :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)                           :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    real(8),dimension(:),allocatable                  :: vt,Hvt
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)    :: Hbath_tmp
#endif    
    integer,dimension(2*Ns_Ud)                        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                             :: Nup,Ndw
    logical                                           :: nonloc_condition, sundry_condition, either_condition
    !
    nonloc_condition = (Norb>1 .AND. (any((Jx_internal/=0d0)) .OR. any((Jp_internal/=0d0))))
    sundry_condition = allocated(coulomb_sundry)
    either_condition = nonloc_condition .OR. sundry_condition
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc /= dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=zero      
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
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    case ("general")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%vg(iorb+Norb*(ispin-1))
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS
    include "direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "direct/HxV_dw.f90"
    !
    if(DimPh>1) then
       !PHONON TERMS
       include "direct/HxV_ph.f90"
       !ELECTRON-PHONON INTERACTION
       include "direct/HxV_eph.f90"
    end if
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(nonloc_condition)then
       include "direct/HxV_non_local.f90"
    endif
    !NON-LOCAL HAMILTONIAN TERMS
    if(sundry_condition)then
       include "direct/HxV_sundry.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_normal_main




  subroutine directMatVec_normal_orbs(Nloc,vin,Hv)
    !
    ! Serial version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm for :f:var:`ed_total_ud` = :code:`False` 
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to the input vector using serial algorithm.  
    !
    integer                                           :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
#ifdef _CMPLX_NORMAL
    complex(8),dimension(Nloc)                        :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                        :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8),dimension(:),allocatable               :: vt,Hvt
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
#else
    real(8),dimension(Nloc)                           :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)                           :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    real(8),dimension(:),allocatable                  :: vt,Hvt
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)    :: Hbath_tmp
#endif
    integer                                           :: isector
    integer,dimension(2*Ns_Ud)                        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                             :: Nup,Ndw
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc /= dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=zero      
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
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    case ("general")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%vg(iorb+Norb*(Nspin-1))!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/Orbs/HxV_local.f90" 
    !
    !UP HAMILTONIAN TERMS
    include "direct/Orbs/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS
    include "direct/Orbs/HxV_dw.f90"
    !
    if(DimPh>1)then
       !PHONON
       include "direct/Orbs/HxV_ph.f90"
       !ELECTRON-PHONON
       include "direct/Orbs/HxV_eph.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_normal_orbs







#ifdef _MPI
  subroutine directMatVec_MPI_normal_main(Nloc,vin,Hv)
    !
    ! MPI parallel version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True` 
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to a part of the vector own by the thread using parallel algorithm.  
    !
    integer                                           :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
#ifdef _CMPLX_NORMAL
    complex(8),dimension(Nloc)                        :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                        :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8),dimension(:),allocatable               :: vt,Hvt
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
#else
    real(8),dimension(Nloc)                           :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)                           :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    real(8),dimension(:),allocatable                  :: vt,Hvt
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)    :: Hbath_tmp
#endif
    integer                                           :: N
    !
    integer,dimension(2*Ns_Ud)                        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                             :: Nup,Ndw
    !
    integer                                           :: i_start,i_end
    logical                                           :: nonloc_condition, sundry_condition, either_condition  
    !
    nonloc_condition = (Norb>1 .AND. (any((Jx_internal/=0d0)) .OR. any((Jp_internal/=0d0))))
    sundry_condition = allocated(coulomb_sundry)
    either_condition = nonloc_condition .OR. sundry_condition
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=zero      
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
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    case ("general")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%vg(iorb+Norb*(ispin-1))
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/HxV_local.f90"
    !NON-LOCAL HAMILTONIAN TERMS
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    do iph=1,DimPh
       allocate(vt(mpiQup*DimDw))
       allocate(Hvt(mpiQup*DimDw))
       vt=zero
       Hvt=zero
       i_start = 1 + (iph-1)*DimUp*MpiQdw
       i_end = iph*DimUp*MpiQdw
       !
       call vector_transpose_MPI(DimUp,MpiQdw,Vin(i_start:i_end),DimDw,MpiQup,vt) !Vin^T --> Vt
       include "direct_mpi/HxV_dw.f90"
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=zero         !reallocate Vt
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
       Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
       !
       deallocate(vt,Hvt)
    end do
    !
    if(DimPh>1) then
       !PHONON TERMS
       include "direct_mpi/HxV_ph.f90"
       !ELECTRON-PHONON INTERACTION
       include "direct_mpi/HxV_eph.f90"
    end if
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(either_condition)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = zero
       call allgather_vector_MPI(MpiComm,vin,vt)
       !
       if(nonloc_condition)then
         include "direct_mpi/HxV_non_local.f90"
       endif
       if(sundry_condition)then
         include "direct_mpi/HxV_sundry.f90"
       endif
       !
       deallocate(vt)
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_normal_main



  subroutine directMatVec_MPI_normal_Orbs(Nloc,vin,Hv)
    !
    !
    ! MPI parallel version of the direct, on-the-fly matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm for :f:var:`ed_total_ud` = :code:`False` 
    ! This procedures evaluates the non-zero terms of any part of the global Hamiltonian and applies them to a part of the vector own by the thread using parallel algorithm.  
    !
    integer                                           :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
#ifdef _CMPLX_NORMAL
    complex(8),dimension(Nloc)                        :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)                        :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8),dimension(:),allocatable               :: vt,Hvt
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
#else
    real(8),dimension(Nloc)                           :: vin  !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)                           :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    real(8),dimension(:),allocatable                  :: vt,Hvt
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)    :: Hbath_tmp
#endif
    !
    integer,dimension(2*Ns_Ud)                        :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                   :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                             :: Nup,Ndw
    integer                                           :: i_start,i_end    
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=zero       
       do ibath=1,Nbath
          do ispin=1,Nspin             
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%v(ispin,iorb,ibath)
             enddo
             do iorb=1,Nfoo
                bath_diag(ispin,iorb,ibath)=one*dmft_bath%e(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    case ("replica")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    case ("general")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=zero
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=zero
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%vg(iorb+Norb*(ispin-1))
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    !    
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/Orbs/HxV_local.f90"
    !
    !NON-LOCAL TERMS:
    ! 
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/Orbs/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    do iph=1,DimPh
       allocate(vt(mpiQup*DimDw)) ;vt=zero
       allocate(Hvt(mpiQup*DimDw));Hvt=zero
       i_start = 1 + (iph-1)*DimUp*MpiQdw
       i_end = iph*DimUp*MpiQdw
       !
       call vector_transpose_MPI(DimUp,MpiQdw,Vin(i_start:i_end),DimDw,MpiQup,vt) !Vin^T --> Vt
       include "direct_mpi/Orbs/HxV_dw.f90"
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=zero         !reallocate Vt
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
       Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
       !
       deallocate(Hvt,vt)
    enddo
    !
    if(DimPh>1)then
       !PHONON TERMS
       include "direct_mpi/Orbs/HxV_ph.f90"
       !ELECTRON-PHONON INTERACTION
       include "direct_mpi/Orbs/HxV_eph.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_normal_Orbs
#endif


END MODULE ED_HAMILTONIAN_NORMAL_DIRECT_HXV
