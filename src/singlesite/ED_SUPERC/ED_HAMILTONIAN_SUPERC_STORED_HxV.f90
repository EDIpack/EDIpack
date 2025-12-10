! > BUILD STORED SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_SUPERC_STORED_HxV
  !:synopsis: Routines for sparse matrix-vector product, :code:`SUPERC` case
  USE ED_HAMILTONIAN_SUPERC_COMMON
  implicit none
  private

  !>Sparse Matric constructors
  public :: ed_buildH_superc_main


  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_superc_main
#ifdef _MPI
  public  :: spMatVec_MPI_superc_main
#endif





contains



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_superc_main(isector,Hmat)
    !
    !
    ! Builds the sector Hamiltonian :math:`H` and save each term in a suitable sparse matrix instance for :f:var:`ed_total_ed` = :code:`True`. If the dimension :f:var:`dim` of the sector are smaller than :f:var:`lanc_dim_threshold` the global matrix is dumped to the optional variable :f:var:`hmat`.
    !
    ! All the different electronic terms  are collected in the same sparse matrix, possibly using rows splitting and local / non-local blocks according to the :f:var:`MPI_Allgatherv` algorithm: 
    !  * :math:`H_{\rm int} \rightarrow` :f:var:`sph0` : interaction part of the electronic Hamiltonian
    !  * :math:`H_{\rm imp} \rightarrow` :f:var:`sph0` : impurity part of the eletronic Hamiltonian 
    !  * :math:`H_{\rm bath} \rightarrow` :f:var:`sph0`: bath levels part of the eletronic Hamiltonian
    !  * :math:`H_{\rm hyb} \rightarrow` :f:var:`sph0` : impurity - bath coupling part of the eletronic Hamiltonian
    !  * :math:`H_{ph} \rightarrow` :f:var:`sph0_ph` : phonon part of the of the global Hamiltonian
    !  * :math:`H_{e-eph} \rightarrow` :f:var:`sph0e_eph` : electron part of the electron-phonon term of the global Hamiltonian
    !  * :math:`H_{ph-eph} \rightarrow` :f:var:`sph0e_ph` : phonon part of the electron-phonon term of the global Hamiltonian
    !
    integer                                                         :: isector
    complex(8),dimension(:,:),optional                              :: Hmat  !optional dense matrix
    complex(8),dimension(:,:),allocatable                           :: Htmp_e,Htmp_ph,Htmp_eph_e,Htmp_eph_ph
    integer,dimension(Nlevels)                                      :: ib
    integer,dimension(Ns)                                           :: ibup,ibdw
    real(8),dimension(Norb)                                         :: nup,ndw
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    complex(8)                                                      :: htmp_el
    integer                                                         :: first_state,last_state
    integer                                                         :: first_state_up,last_state_up
    integer                                                         :: first_state_dw,last_state_dw
    integer                                                         :: jup, jdw, iph, i_el, j_el
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG ed_buildH_main SUPERC: build H"
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim   = Hsector%Dim
    DimEl = Hsector%DimEl
    !
    if(present(Hmat))call assert_shape(Hmat,[Dim,Dim],"ed_buildh_main","Hmat")
    !
    Hbath_tmp = zero
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
    do ibath=1,Nbath
      do ispin=1,Nnambu*Nspin
        if(any(abs(Hbath_tmp(ispin,1+Nnambu*Nspin-ispin,:,:,ibath)-transpose(Hbath_tmp(ispin,1+Nnambu*Nspin-ispin,:,:,ibath)))>1d-10))then
          STOP "anomalous Hbath is not symmetric. Only spin-singlet pairing is allowed in SUPERC."
        endif
      enddo
    enddo
    
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0,Dim)
       ! if(DimPh>1) then
       !    call sp_set_mpi_matrix(MpiComm,spH0e_eph,mpiIstart,mpiIend,mpiIshift)
       !    call sp_init_matrix(MpiComm,spH0e_eph,DimEl)
       ! endif
    else
       call sp_init_matrix(spH0,Dim)
       !if(DimPh>1) call sp_init_matrix(spH0e_eph,DimEl)
    endif
#else
    call sp_init_matrix(spH0,Dim)
    if(DimPh>1) then
       call sp_init_matrix(spH0e_eph,DimEl)
       call sp_init_matrix(spH0_ph,DimPh)
       call sp_init_matrix(spH0ph_eph,DimPh)
    end if
#endif
    !
    !-----------------------------------------------!
    !
#ifdef _MPI
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Himp"
#endif
    !IMPURITY  HAMILTONIAN
    include "stored_mpi/Himp.f90"
    !
    !LOCAL INTERACTION
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Hint"
#endif
    include "stored_mpi/Hint.f90"
    !
    !BATH HAMILTONIAN
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Hbath"
#endif
    include "stored_mpi/Hbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Himp_bath"
#endif
    include "stored_mpi/Himp_bath.f90"
    !
    if(DimPh>1) then
       !PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/H_ph"
#endif
       include "stored_mpi/H_ph.f90"
       !
       !ELECTRON-PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/H_e_ph"
#endif
       include "stored_mpi/H_e_ph.f90"
    endif


#else ! HERE NOT MPI
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Himp"
#endif
    !IMPURITY  HAMILTONIAN
    include "stored/Himp.f90"
    !
    !LOCAL INTERACTION
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Hint"
#endif
    include "stored/Hint.f90"
    !
    !BATH HAMILTONIAN
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Hbath"
#endif
    include "stored/Hbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/Himp_bath"
#endif
    include "stored/Himp_bath.f90"
    !
    if(DimPh>1) then
       !PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/H_ph"
#endif
       include "stored/H_ph.f90"
       !
       !ELECTRON-PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_SUPERC: stored/H_e_ph"
#endif
       include "stored/H_e_ph.f90"
    endif
#endif !IFdef _MPI
    !
    !-----------------------------------------------!
    !
    !
    if(present(Hmat))then
       Hmat=zero
#ifdef _MPI
       allocate(Htmp_e(Dim,Dim))
       Htmp_e=0.d0
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0,Htmp_e)
       else
          call sp_dump_matrix(spH0,Htmp_e)
       endif
       Hmat=Htmp_e
#else
       allocate(Htmp_e(DimEl,DimEl)); Htmp_e=0.d0
       call sp_dump_matrix(spH0,Htmp_e)
       if(DimPh>1) then
          allocate(Htmp_ph(DimPh,DimPh));Htmp_ph=0d0
          allocate(Htmp_eph_ph(DimPh,DimPh));Htmp_eph_ph=0d0
          allocate(Htmp_eph_e(DimEl,DimEl));Htmp_eph_e=0d0
          !
          call sp_dump_matrix(spH0_ph,Htmp_ph)
!#ifdef _MPI
          ! if(MpiStatus)then
          !    call sp_dump_matrix(MpiComm,spH0e_eph,Htmp_eph_e)
          ! else
          !    call sp_dump_matrix(spH0e_eph,Htmp_eph_e)
          ! endif
!#else
          call sp_dump_matrix(spH0e_eph,Htmp_eph_e)
!#endif
          call sp_dump_matrix(spH0ph_eph,Htmp_eph_ph)
          !
          Hmat = kronecker_product(zeye(DimPh),Htmp_e) + &
               kronecker_product(Htmp_ph,zeye(DimEl)) + &
               kronecker_product(Htmp_eph_ph,Htmp_eph_e)
          deallocate(Htmp_ph,Htmp_eph_e,Htmp_eph_ph)
       else
          Hmat = Htmp_e
       endif       
#endif
       !
       deallocate(Htmp_e)
    endif
    deallocate(diag_hybr,bath_diag)
    return
    !
  end subroutine ed_buildH_superc_main










  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
  subroutine spMatVec_superc_main(Nloc,v,Hv)
    !
    ! Serial version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm.
    ! This procedures applies one by one each term of the global Hamiltonian to an input vector using the stored sparse matrices.  
    !
    integer                         :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc)      :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)      :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    complex(8)                      :: val
    integer                         :: i,j,iph,jp, i_el,j_el
    !
    Hv=zero
    ! Electron Part
    do i=1,Nloc
       iph  = (i-1)/DimEl+1
       i_el = mod(i-1,DimEl) +1
       matmul: do j_el=1, spH0%row(i_el)%Size
          val = spH0%row(i_el)%cvals(j_el)
          j   = spH0%row(i_el)%cols(j_el) + (iph-1)*DimEl
          Hv(i) = Hv(i) + val*v(j)
       end do matmul
    end do
    !
    ! Phononic Part
    if(DimPh>1)then
       do iph=1,DimPh
          do i_el=1,DimEl
             i =i_el + (iph-1)*DimEl
             !
             !PHONON
             do jp=1,spH0_ph%row(iph)%size
                val = spH0_ph%row(iph)%cvals(jp)
                j = i_el + (spH0_ph%row(iph)%cols(jp)-1)*DimEl
                Hv(i) = Hv(i) +val*v(j)
             enddo
             !
             !ELECTRON-PHONON
             do j_el =1,spH0e_eph%row(i_el)%size
                do jp=1,spH0ph_eph%row(iph)%size
                   val = spH0e_eph%row(i_el)%cvals(j_el)*&
                        spH0ph_eph%row(iph)%cvals(jp)
                   j   = spH0e_eph%row(i_el)%cols(j_el)+&
                        (spH0ph_eph%row(iph)%cols(jp)-1)*DimEl
                   Hv(i) = Hv(i) + val*v(j)
                enddo
             enddo
          enddo
       enddo
    endif
    !
  end subroutine spMatVec_superc_main


#ifdef _MPI
  subroutine spMatVec_mpi_superc_main(Nloc,v,Hv)
    !
    !
    ! MPI parallel version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm.
    ! This procedures applies one by one each term of the global Hamiltonian to a part of the vector own by the thread using the stored sparse matrices.
    !
    integer                             :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)          :: v    !input vector part (passed by P-Arpack/P-Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)          :: Hv   !!output vector (required by P-Arpack/P-Lanczos) :math:`\vec{w}`
    complex(8)                          :: val
    integer                             :: i,j,mpiIerr,iph,jp, i_el,j_el, i_loc,j_loc,dimel_loc,j_max
    integer                             :: N,MpiShift
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: Counts,Offset
    !
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !if(DimPh>1) stop "ERROR: MPI superc phonon (Nph>0) is not supported"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    !
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    MpiShift = spH0%Ishift
    Hv=0d0
    do i_loc=1,Nloc
       !TO CHECK
       iph  = (MpiShift +i_loc-1)/DimEl +1 !GLOBAL/LOCAL phonon index
       i_el = mod(MpiShift +i_loc-1,DimEl) +1 !GLOBAL electron index
       local: do j_loc=1,spH0%loc(i_loc)%Size
          j = spH0%loc(i_loc)%cols(j_loc)
          Hv(i_loc) = Hv(i_loc) + spH0%loc(i_loc)%cvals(j_loc)*v(j-MpiShift)
       end do local
    end do
    !
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    ! Count(i_cpu) = number of elements in the i_cpu-th CPU
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)) ; vin = zero
    !GATHERVECTOR
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin      ,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    ! El,Ph and El-Ph Part
    do i_loc=1,Nloc                 !==spH0%Nrow
       iph  = (i_loc+MpiShift-1)/DimEl +1
       i_el = mod(i_loc+MpiShift-1,DimEl) +1
       matmul: do j_loc=1,spH0%row(i_loc)%Size
          j = spH0%row(i_loc)%cols(j_loc)
          Hv(i_loc) = Hv(i_loc) + spH0%row(i_loc)%cvals(j_loc)*vin(j)
       end do matmul
    end do
    !
  end subroutine spMatVec_mpi_superc_main
#endif











end MODULE ED_HAMILTONIAN_SUPERC_STORED_HXV
