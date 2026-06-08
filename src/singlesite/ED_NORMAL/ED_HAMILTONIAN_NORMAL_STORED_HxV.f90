MODULE ED_HAMILTONIAN_NORMAL_STORED_HxV
  !:synopsis: Routines for sparse matrix-vector product, :code:`NORMAL` case
  !Constructs each terms of the sector Hamiltonian storing them into different :f:var:`sparse_matrix` instances, implement the corresponding matrix-vector products using stored sparse matrices. 
  !
  USE ED_HAMILTONIAN_NORMAL_COMMON
  implicit none
  private


  !>Sparse Matric constructors
  public :: ed_buildh_normal_main
  public :: ed_buildh_normal_orbs

  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_normal_main
  public  :: spMatVec_normal_orbs
#ifdef _MPI
  public  :: spMatVec_MPI_normal_main
  public  :: spMatVec_MPI_normal_orbs
#endif


contains


  subroutine ed_buildh_normal_main(Hmat)
    !
    ! Builds the sector Hamiltonian :math:`H` and save each term in a suitable sparse matrix instance for :f:var:`ed_total_ed` = :code:`True`. If the dimension :f:var:`dim` of the sector are smaller than :f:var:`lanc_dim_threshold` the global matrix is dumped to the optional variable :f:var:`hmat`.
    !
    ! The sparse matrices are:
    !  * :math:`H_d \rightarrow` :f:var:`sph0d` : diagonal part of the electronic Hamiltonian
    !  * :math:`H_\uparrow \rightarrow` :f:var:`sph0ups` : :math:`\uparrow` spin terms of the eletronic Hamiltonian 
    !  * :math:`H_\downarrow \rightarrow` :f:var:`sph0dws` : :math:`\downarrow`  spin terms of the eletronic Hamiltonian
    !  * :math:`H_{nd} \rightarrow` :f:var:`sph0nd` : non-diagonal part of the eletronic Hamiltonian
    !  * :math:`H_{ph} \rightarrow` :f:var:`sph0_ph` : phonon part of the of the global Hamiltonian
    !  * :math:`H_{e-eph} \rightarrow` :f:var:`sph0e_eph` : electron part of the electron-phonon term of the global Hamiltonian
    !  * :math:`H_{ph_eph} \rightarrow` :f:var:`sph0e_ph` : phonon part of the electron-phonon term of the global Hamiltonian
    !
#ifdef _CMPLX_NORMAL
    complex(8),dimension(:,:),optional                :: Hmat !optional dense matrix
    complex(8),dimension(:,:),allocatable             :: Htmp_up,Htmp_dw,Hrdx,Hmat_tmp
    complex(8),dimension(:,:),allocatable             :: Htmp_ph,Htmp_eph_e,Htmp_eph_ph
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
#else    
    real(8),dimension(:,:),optional                   :: Hmat !optional dense matrix
    real(8),dimension(:,:),allocatable                :: Htmp_up,Htmp_dw,Hrdx,Hmat_tmp
    real(8),dimension(:,:),allocatable                :: Htmp_ph,Htmp_eph_e,Htmp_eph_ph
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)    :: Hbath_tmp
#endif
    integer                                           :: isector
    integer,dimension(2*Ns_Ud)                        :: Indices    ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                             :: Nup,Ndw    ![Ns]
    logical                                           :: nonloc_condition, sundry_condition, either_condition
    !
    nonloc_condition = (Norb>1 .AND. (any((Jx_internal/=0d0)) .OR. any((Jp_internal/=0d0))))
    sundry_condition = allocated(coulomb_sundry)
    either_condition = nonloc_condition .OR. sundry_condition
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG ed_buildH_main NORMAL: build H"
#endif
    !
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
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
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%v!(ispin)
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
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%vg(iorb+Norb*(ispin-1))!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,DimUp*DimDw)
       if(DimPh>1) then
          call sp_set_mpi_matrix(MpiComm,spH0e_eph,mpiIstart,mpiIend,mpiIshift)
          call sp_init_matrix(MpiComm,spH0e_eph,DimUp*DimDw)
       endif
       !
       if(either_condition)then
          call sp_set_mpi_matrix(MpiComm,spH0nd,mpiIstart,mpiIend,mpiIshift)
          call sp_init_matrix(MpiComm,spH0nd,DimUp*DimDw)
       endif
    else
       call sp_init_matrix(spH0d,DimUp*DimDw)
       if(DimPh>1) call sp_init_matrix(spH0e_eph,DimUp*DimDw)
       if(either_condition)call sp_init_matrix(spH0nd,DimUp*DimDw)
    endif
#else
    call sp_init_matrix(spH0d,DimUp*DimDw)
    if(DimPh>1) call sp_init_matrix(spH0e_eph,DimUp*DimDw)
    if(either_condition)call sp_init_matrix(spH0nd,DimUp*DimDw)
#endif
    call sp_init_matrix(spH0dws(1),DimDw)
    call sp_init_matrix(spH0ups(1),DimUp)
    if(DimPh>1) then
       call sp_init_matrix(spH0_ph,DimPh)
       call sp_init_matrix(spH0ph_eph,DimPh)
    end if
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/H_local"
#endif
    include "stored/H_local.f90"
    !
    !NON-LOCAL HAMILTONIAN TERMS
    if(nonloc_condition)then
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/H_non_local"
#endif
       include "stored/H_non_local.f90"
    endif
    !
    !NON-LOCAL HAMILTONIAN TERMS
    if(sundry_condition)then
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/user_defined_non_HK_terms"
#endif
       include "stored/H_sundry.f90"
    endif
    !
    !UP TERMS
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/H_up"
#endif
    include "stored/H_up.f90"
    !
    !DW TERMS
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/H_dw"
#endif
    include "stored/H_dw.f90"
    !
    if(DimPh>1) then
       !PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/H_ph"
#endif
       include "stored/H_ph.f90"
       !
       !ELECTRON-PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NORMAL: stored/H_e_ph"
#endif
       include "stored/H_e_ph.f90"
    endif
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = zero
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=zero
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=zero
       allocate(Hmat_tmp(DimUp*DimDw,DimUp*DimDw));Hmat_tmp=zero
       !
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat_tmp)
       else
          call sp_dump_matrix(spH0d,Hmat_tmp)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat_tmp)
#endif
       !
       if(either_condition)then
          allocate(Hrdx(DimUp*DimDw,DimUp*DimDw));Hrdx=zero
#ifdef _MPI
          if(MpiStatus)then
             call sp_dump_matrix(MpiComm,spH0nd,Hrdx)
          else
             call sp_dump_matrix(spH0nd,Hrdx)
          endif
#else
          call sp_dump_matrix(spH0nd,Hrdx)
#endif
          Hmat_tmp = Hmat_tmp + Hrdx
          deallocate(Hrdx)
       endif
       !
       call sp_dump_matrix(spH0ups(1),Htmp_up)
       call sp_dump_matrix(spH0dws(1),Htmp_dw)
       Hmat_tmp = Hmat_tmp + kronecker_product(Htmp_dw,eye(DimUp))
       Hmat_tmp = Hmat_tmp + kronecker_product(eye(DimDw),Htmp_up)
       !
       if(DimPh>1) then
          allocate(Htmp_ph(DimPh,DimPh));Htmp_ph=zero
          allocate(Htmp_eph_ph(DimPh,DimPh));Htmp_eph_ph=zero
          allocate(Htmp_eph_e(DimUp*DimDw,DimUp*DimDw));Htmp_eph_e=zero
          !
          call sp_dump_matrix(spH0_ph,Htmp_ph)
#ifdef _MPI
          if(MpiStatus)then
             call sp_dump_matrix(MpiComm,spH0e_eph,Htmp_eph_e)
          else
             call sp_dump_matrix(spH0e_eph,Htmp_eph_e)
          endif
#else
          call sp_dump_matrix(spH0e_eph,Htmp_eph_e)
#endif
          call sp_dump_matrix(spH0ph_eph,Htmp_eph_ph)
          !
          Hmat = kronecker_product(eye(Dimph),Hmat_tmp) +     &
               kronecker_product(Htmp_ph,eye(DimUp*DimDw)) +&
               kronecker_product(Htmp_eph_ph,Htmp_eph_e)
          !
          deallocate(Htmp_ph,Htmp_eph_e,Htmp_eph_ph)
       else
          Hmat = Hmat_tmp
       endif
       !
       deallocate(Htmp_up,Htmp_dw,Hmat_tmp)
    endif
    !
    deallocate(diag_hybr,bath_diag)
    return
    !
  end subroutine ed_buildh_normal_main





  subroutine ed_buildh_normal_orbs(Hmat)
    !
    ! Builds the sector Hamiltonian :math:`H` and save each term in a suitable sparse matrix instance for :f:var:`ed_total_ed` = :code:`False`. If the dimension :f:var:`dim` of the sector are smaller than :f:var:`lanc_dim_threshold` the global matrix is dumped to the optional variable :f:var:`hmat`. 
    !
    ! The sparse matrices are:
    !  * :math:`H_d \rightarrow` :f:var:`sph0d` : diagonal part of the electronic Hamiltonian
    !  * :math:`\vec{H}_\uparrow \rightarrow` :f:var:`sph0ups` : :math:`\uparrow` spin terms of the eletronic Hamiltonian 
    !  * :math:`\vec{H}_\downarrow \rightarrow` :f:var:`sph0dws` : :math:`\downarrow`  spin terms of the eletronic Hamiltonian
    !  * :math:`H_{nd} \rightarrow` :f:var:`sph0nd` : non-diagonal part of the eletronic Hamiltonian
    !  * :math:`H_{ph} \rightarrow` :f:var:`sph0_ph` : phonon part of the of the global Hamiltonian
    !  * :math:`H_{e-eph} \rightarrow` :f:var:`sph0e_eph` : electron part of the electron-phonon term of the global Hamiltonian
    !  * :math:`H_{ph_eph} \rightarrow` :f:var:`sph0e_ph` : phonon part of the electron-phonon term of the global Hamiltonian
    !
    ! where :math:`\vec{H}_\sigma = [H^1_\sigma,\dots,H^{Norb}_\sigma]` are the orbital and spin resolved Hamiltonian matrices.
    !
#ifdef _CMPLX_NORMAL
    complex(8),dimension(:,:),optional                :: Hmat !optional dense matrix
    complex(8),dimension(:,:),allocatable             :: Hmat_tmp,Htmp_ph,Htmp_eph_e,Htmp_eph_ph
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
#else
    real(8),dimension(:,:),optional                   :: Hmat !optional dense matrix
    real(8),dimension(:,:),allocatable                :: Hmat_tmp,Htmp_ph,Htmp_eph_e,Htmp_eph_ph
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)    :: Hbath_tmp
#endif
    integer                                           :: isector
    integer                                           :: mDimUp,mDimDw
    integer,dimension(2*Ns_Ud)                        :: Indices,Jndices
    integer,dimension(Ns_Ud,Ns_Orb)                   :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                             :: Nup,Ndw    ![Ns]
    integer                                           :: i,j,jj,iud
    integer                                           :: iup,idw,jup,jdw
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG ed_buildH_orbs NORMAL: build H"
#endif
    !
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT set"
    isector=Hsector%index
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
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
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%v!(ispin)
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
                diag_hybr(ispin,iorb,ibath)=one*dmft_bath%item(ibath)%vg(iorb+(ispin-1)*Nspin)!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,DimUp*DimDw)
       if(DimPh>1) then
          call sp_set_mpi_matrix(MpiComm,spH0e_eph,mpiIstart,mpiIend,mpiIshift)
          call sp_init_matrix(MpiComm,spH0e_eph,DimUp*DimDw)
       endif
    else
       call sp_init_matrix(spH0d,DimUp*DimDw)
       if(DimPh>1) call sp_init_matrix(spH0e_eph,DimUp*DimDw)
    endif
#else
    call sp_init_matrix(spH0d,DimUp*DimDw)
    if(DimPh>1) call sp_init_matrix(spH0e_eph,DimUp*DimDw)
#endif
    do iud=1,Ns_Ud
       call sp_init_matrix(spH0dws(iud),DimDws(iud))
       call sp_init_matrix(spH0ups(iud),DimUps(iud))
    enddo
    if(DimPh>1) then
       call sp_init_matrix(spH0_ph,DimPh)
       call sp_init_matrix(spH0ph_eph,DimPh)
    end if
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_orbs NORMAL: stored/Orbs/H_local"
#endif
    include "stored/Orbs/H_local.f90"
    !
    !UP TERMS
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_orbs NORMAL: stored/Orbs/H_up"
#endif
    include "stored/Orbs/H_up.f90"
    !
    !DW TERMS
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_orbs NORMAL: stored/Orbs/H_dw"
#endif
    include "stored/Orbs/H_dw.f90"
    !
    if(DimPh>1)then
       !PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_orbs NORMAL: stored/Orbs/H_ph"
#endif
       include "stored/Orbs/H_ph.f90"
       !
       !ELECTRON-PHONON TERMS
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_orbs NORMAL: stored/Orbs/H_e_ph"
#endif
       include "stored/Orbs/H_e_ph.f90"
    endif
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       Hmat = zero
       allocate(Hmat_tmp(DimUp*DimDw,DimUp*DimDw));Hmat_tmp=zero      
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat_tmp)
       else
          call sp_dump_matrix(spH0d,Hmat_tmp)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat_tmp)
#endif
       do i=1,DimUp*DimDw
          call state2indices(i,[DimUps,DimDws],Indices)
          do iud=1,Ns_Ud
             !UP:
             iup = Indices(iud)
             do jj=1,spH0ups(iud)%row(iup)%Size
                Jndices = Indices ; Jndices(iud) = spH0ups(iud)%row(iup)%cols(jj)
                call indices2state(Jndices,[DimUps,DimDws],j)
#ifdef _CMPLX_NORMAL
                Hmat_tmp(i,j) = Hmat_tmp(i,j) + spH0ups(iud)%row(iup)%cvals(jj)
#else
                Hmat_tmp(i,j) = Hmat_tmp(i,j) + spH0ups(iud)%row(iup)%dvals(jj)
#endif
             enddo
             !DW:
             idw = Indices(iud+Ns_Ud)
             do jj=1,spH0dws(iud)%row(idw)%Size
                Jndices = Indices ; Jndices(iud+Ns_Ud) = spH0dws(iud)%row(idw)%cols(jj)
                call indices2state(Jndices,[DimUps,DimDws],j)
#ifdef _CMPLX_NORMAL                
                Hmat_tmp(i,j) = Hmat_tmp(i,j) + spH0dws(iud)%row(idw)%dvals(jj)
#else                
                Hmat_tmp(i,j) = Hmat_tmp(i,j) + spH0dws(iud)%row(idw)%dvals(jj)
#endif                
             enddo
             !
          enddo
       enddo
       !
       if(DimPh>1) then
          allocate(Htmp_ph(DimPh,DimPh));Htmp_ph=zero
          allocate(Htmp_eph_ph(DimPh,DimPh));Htmp_eph_ph=zero
          allocate(Htmp_eph_e(DimUp*DimDw,DimUp*DimDw));Htmp_eph_e=zero
          !
          call sp_dump_matrix(spH0_ph,Htmp_ph)
#ifdef _MPI
          if(MpiStatus)then
             call sp_dump_matrix(MpiComm,spH0e_eph,Htmp_eph_e)
          else
             call sp_dump_matrix(spH0e_eph,Htmp_eph_e)
          endif
#else
          call sp_dump_matrix(spH0e_eph,Htmp_eph_e)
#endif
          call sp_dump_matrix(spH0ph_eph,Htmp_eph_ph)
          !
          Hmat = kronecker_product(eye(Dimph),Hmat_tmp) +     &
               kronecker_product(Htmp_ph,eye(DimUp*DimDw)) +&
               kronecker_product(Htmp_eph_ph,Htmp_eph_e)
          !
          deallocate(Htmp_ph,Htmp_eph_e,Htmp_eph_ph)
       else
          Hmat = Hmat_tmp
       endif
       !
       deallocate(Hmat_tmp)
    endif
    !
    deallocate(diag_hybr,bath_diag)
    return
    !
  end subroutine ed_buildh_normal_orbs












  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !+------------------------------------------------------------------+
  subroutine spMatVec_normal_main(Nloc,v,Hv)
    !
    ! Serial version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True` 
    ! This procedures applies one by one each term of the global Hamiltonian to an input vector using the stored sparse matrices.  
    !
    integer                 :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    real(8),dimension(Nloc) :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc) :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    real(8)                 :: val
    integer                 :: i,iup,idw,j,jup,jdw,jj,i_el,j_el
    logical                 :: nonloc_condition, sundry_condition, either_condition
    !
    nonloc_condition = (Norb>1 .AND. (any((Jx_internal/=0d0)) .OR. any((Jp_internal/=0d0))))
    sundry_condition = allocated(coulomb_sundry)
    either_condition = nonloc_condition .OR. sundry_condition
    !
    Hv=0d0
    !
    !Local:
    do i = 1,Nloc
       iph = (i-1)/(DimUp*DimDw) + 1   !phonon index [1:DimPh]
       i_el = mod(i-1,DimUp*DimDw) + 1 !electron index [1:DimUp*DimDw]
       do j_el=1,spH0d%row(i_el)%Size
          val = spH0d%row(i_el)%dvals(j_el)
          j = spH0d%row(i_el)%cols(j_el) + (iph-1)*DimUp*DimDw
          Hv(i) = Hv(i) + val*v(j)
       enddo
    enddo
    !
    do iph = 1,DimPh
       !DW:
       do iup=1,DimUp
          !
          do idw=1,DimDw
             i = iup + (idw-1)*DimUp + (iph-1)*DimUp*DimDw
             do jj=1,spH0dws(1)%row(idw)%Size
                jup = iup
                jdw = spH0dws(1)%row(idw)%cols(jj)
                val = spH0dws(1)%row(idw)%dvals(jj)
                j     = jup +  (jdw-1)*DimUp + (iph-1)*DimUp*DimDw
                Hv(i) = Hv(i) + val*V(j)
             enddo
          enddo
          !
       enddo
       !
       !UP:
       do idw=1,DimDw
          !
          do iup=1,DimUp
             i = iup + (idw-1)*DimUp + (iph-1)*DimUp*DimDw
             do jj=1,spH0ups(1)%row(iup)%Size
                jup = spH0ups(1)%row(iup)%cols(jj)
                jdw = idw
                val = spH0ups(1)%row(iup)%dvals(jj)
                j =  jup + (jdw-1)*DimUp + (iph-1)*DimUp*DimDw
                Hv(i) = Hv(i) + val*V(j)
             enddo
          enddo
          !
       enddo
       !
       if(DimPh>1) then
          do i_el = 1,DimUp*DimDw
             i = i_el + (iph-1)*DimUp*DimDw
             !
             !PHONON
             do jj = 1,spH0_ph%row(iph)%Size
                val = spH0_ph%row(iph)%dvals(jj)
                j = i_el + (spH0_ph%row(iph)%cols(jj)-1)*DimUp*DimDw
                Hv(i) = Hv(i) + val*v(j)
             enddo
             !
             !ELECTRON-PHONON
             do j_el = 1,spH0e_eph%row(i_el)%Size
                do jj = 1,spH0ph_eph%row(iph)%Size
                   val = spH0e_eph%row(i_el)%dvals(j_el)*&
                        spH0ph_eph%row(iph)%dvals(jj)
                   j = spH0e_eph%row(i_el)%cols(j_el) +&
                        (spH0ph_eph%row(iph)%cols(jj)-1)*DimUp*DimDw
                   Hv(i) = Hv(i) + val*v(j)
                enddo
             enddo
             !
          enddo
       endif
       !
    enddo
    !
    !Non-Local:
    if(either_condition)then
       do i = 1,Nloc
          iph = (i-1)/(DimUp*DimDw) + 1
          i_el = mod(i-1,DimUp*DimDw) + 1
          do j_el=1,spH0nd%row(i_el)%Size
             val = spH0nd%row(i_el)%dvals(j_el)
             j = spH0nd%row(i_el)%cols(j_el) + (iph-1)*DimUp*DimDw
             Hv(i) = Hv(i) + val*v(j)
          enddo
       enddo
    endif
    !
  end subroutine spMatVec_normal_main

  subroutine spMatVec_normal_orbs(Nloc,v,Hv)
    !
    ! Serial version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm for :f:var:`ed_total_ud` = :code:`False` 
    ! This procedures applies one by one each term of the global Hamiltonian to an input vector using the stored sparse matrices.
    !
    integer                    :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    real(8),dimension(Nloc)    :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)    :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    real(8)                    :: val
    integer                    :: i,iup,idw,j,jup,jdw,jj,i_el,j_el
    integer                    :: iud
    integer,dimension(2*Ns_Ud) :: Indices,Jndices
    !
    !
    Hv=0d0
    !

    do i = 1,Nloc
       i_el = mod(i-1,DimUp*DimDw) + 1
       !
       do j=1,spH0d%row(i_el)%Size
          Hv(i) = Hv(i) + spH0d%row(i_el)%dvals(j)*v(i)
       enddo
    enddo
    !
    !
    do i=1,Nloc
       i_el = mod(i-1,DimUp*DimDw) + 1
       iph = (i-1)/(DimUp*DimDw) + 1
       !
       call state2indices(i_el,[DimUps,DimDws],Indices)
       do iud=1,Ns_Ud
          !
          !UP:
          iup = Indices(iud)
          do jj=1,spH0ups(iud)%row(iup)%Size
             Jndices = Indices ; Jndices(iud) = spH0ups(iud)%row(iup)%cols(jj)
             call indices2state(Jndices,[DimUps,DimDws],j)
             !
             j = j + (iph-1)*DimUp*DimDw
             Hv(i) = Hv(i) + spH0ups(iud)%row(iup)%dvals(jj)*V(j)
          enddo
          !
          !DW:
          idw = Indices(iud+Ns_Ud)
          do jj=1,spH0dws(iud)%row(idw)%Size
             Jndices = Indices ; Jndices(iud+Ns_Ud) = spH0dws(iud)%row(idw)%cols(jj)
             call indices2state(Jndices,[DimUps,DimDws],j)
             !
             j = j + (iph-1)*DimUp*DimDw
             Hv(i) = Hv(i) + spH0dws(iud)%row(idw)%dvals(jj)*V(j)
          enddo
          !
       enddo
    enddo
    !
    if(DimPh>1)then
       do i=1,Nloc
          i_el = mod(i-1,DimUp*DimDw) + 1
          iph = (i-1)/(DimUp*DimDw) + 1
          !
          !PHONON
          do jj = 1,spH0_ph%row(iph)%Size
             val = spH0_ph%row(iph)%dvals(jj)
             j = i_el + (spH0_ph%row(iph)%cols(jj)-1)*DimUp*DimDw
             Hv(i) = Hv(i) + val*v(j)
          enddo
          !
          !ELECTRON-PHONON
          do j_el = 1,spH0e_eph%row(i_el)%Size
             do jj = 1,spH0ph_eph%row(iph)%Size
                val = spH0e_eph%row(i_el)%dvals(j_el)*&
                     spH0ph_eph%row(iph)%dvals(jj)
                j = spH0e_eph%row(i_el)%cols(j_el) +&
                     (spH0ph_eph%row(iph)%cols(jj)-1)*DimUp*DimDw
                Hv(i) = Hv(i) + val*v(j)
             enddo
          enddo
          !
       enddo
    endif
    !
  end subroutine spMatVec_normal_orbs


#ifdef _MPI
  subroutine spMatVec_mpi_normal_main(Nloc,v,Hv)
    !
    ! MPI parallel version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm for :f:var:`ed_total_ud` = :code:`True`. 
    ! This procedures applies one by one each term of the global Hamiltonian to a part of the vector own by the thread using the stored sparse matrices.
    !
    integer                          :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    real(8),dimension(Nloc)          :: v    !input vector part (passed by P-Arpack/P-Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)          :: Hv   !output vector (required by P-Arpack/P-Lanczos) :math:`\vec{w}`
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    real(8),dimension(:),allocatable :: vin
    real(8)                          :: val
    integer                          :: i,iup,idw,j,jup,jdw,jj
    integer                          :: i_el,j_el,i_start,i_end
    !local MPI
    integer                          :: irank
    logical                          :: nonloc_condition, sundry_condition, either_condition
    !
    nonloc_condition = (Norb>1 .AND. (any((Jx_internal/=0d0)) .OR. any((Jp_internal/=0d0))))
    sundry_condition = allocated(coulomb_sundry)
    either_condition = nonloc_condition .OR. sundry_condition
    !
    ! if(MpiComm==Mpi_Comm_Null)return
    ! if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_normal_main ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    do i=1,Nloc                
       i_el = mod(i-1,DimUp*MpiQdw) + 1
       ! do j_el=1,spH0d%row(i_el)%Size
       val = spH0d%row(i_el)%dvals(1)!(j_el)
       Hv(i) = Hv(i) + val*v(i)
       ! enddo
    end do
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do iph=1,DimPh
       do idw=1,MpiQdw
          do iup=1,DimUp
             i = iup + (idw-1)*DimUp + (iph-1)*DimUp*MpiQdw
             hxv_up: do jj=1,spH0ups(1)%row(iup)%Size
                jup = spH0ups(1)%row(iup)%cols(jj)
                jdw = idw
                val = spH0ups(1)%row(iup)%dvals(jj)
                j   = jup + (idw-1)*DimUp + (iph-1)*DimUp*MpiQdw
                Hv(i) = Hv(i) + val*v(j)
             end do hxv_up
          enddo
       end do
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    do iph=1,DimPh
       allocate(vt(mpiQup*DimDw))
       allocate(Hvt(mpiQup*DimDw))
       vt=0d0
       Hvt=0d0
       i_start = 1 + (iph-1)*DimUp*MpiQdw
       i_end = iph*DimUp*MpiQdw
       call vector_transpose_MPI(DimUp,MpiQdw,v(i_start:i_end),DimDw,MpiQup,vt)
       do idw=1,MpiQup             !<= Transposed order:  column-wise DW <--> UP  
          do iup=1,DimDw           !<= Transposed order:  column-wise DW <--> UP
             i = iup + (idw-1)*DimDw
             hxv_dw: do jj=1,spH0dws(1)%row(iup)%Size
                jup = spH0dws(1)%row(iup)%cols(jj)
                jdw = idw             
                j   = jup + (jdw-1)*DimDw
                val = spH0dws(1)%row(iup)%dvals(jj)
                Hvt(i) = Hvt(i) + val*vt(j)
             end do hxv_dw
          enddo
       end do
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=0d0
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
       Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
       deallocate(vt,Hvt)
    end do
    !
    if(DimPh>1)then
       do iph=1,DimPh
          do i_el = 1,DimUp*MpiQdw
             i = i_el + (iph-1)*DimUp*MpiQdw
             !
             !PHONON
             do jj = 1,spH0_ph%row(iph)%Size
                val = spH0_ph%row(iph)%dvals(jj)
                j = i_el + (spH0_ph%row(iph)%cols(jj)-1)*DimUp*MpiQdw
                Hv(i) = Hv(i) + val*v(j)
             enddo
             !
             !ELECTRON-PHONON
             do j_el = 1,spH0e_eph%row(i_el)%Size
                do jj = 1,spH0ph_eph%row(iph)%Size
                   val = spH0e_eph%row(i_el)%dvals(j_el)*&
                        spH0ph_eph%row(iph)%dvals(jj)
                   !interaction is diag from the electron point of view (coupling to the density)
                   j = spH0e_eph%row(i_el)%cols(j_el) + (spH0ph_eph%row(iph)%cols(jj)-1)*DimUp*MpiQdw
                   Hv(i) = Hv(i) + val*v(j)
                enddo
             enddo
             !
          enddo
       enddo
    end if
    !
    !Non-Local:
    if(either_condition)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       ! 
       allocate(vt(N)) ; vt = 0d0
       call allgather_vector_MPI(MpiComm,v,vt)
       !
       do i=1,Nloc
          iph  = (i-1)/(DimUp*MpiQdw)  + 1
          i_el = mod(i-1,DimUp*MpiQdw) + 1
          matmul: do j_el=1,spH0nd%row(i_el)%Size
             val = spH0nd%row(i_el)%dvals(j_el)
             j = spH0nd%row(i_el)%cols(j_el) + (iph-1)*DimUp*DimDw
             Hv(i) = Hv(i) + val*Vt(j)
          enddo matmul
       enddo
       deallocate(Vt)
    endif
    !
  end subroutine spMatVec_mpi_normal_main


  subroutine spMatVec_mpi_normal_orbs(Nloc,v,Hv)
    !
    ! MPI parallel version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm for :f:var:`ed_total_ud` = :code:`False`. 
    ! This procedures applies one by one each term of the global Hamiltonian to a part of the vector own by the thread using the stored sparse matrices.
    !
    integer                          :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    real(8),dimension(Nloc)          :: v    !input vector part (passed by P-Arpack/P-Lanczos) :math:`\vec{v}`
    real(8),dimension(Nloc)          :: Hv   !output vector (required by P-Arpack/P-Lanczos) :math:`\vec{w}`
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    real(8)                          :: val
    integer                          :: i,iup,idw,j,jup,jdw,jj
    integer                          :: iiup,iidw
    integer                          :: i_el,j_el,i_start,i_end
    integer                          :: iud
    integer,dimension(2*Ns_Ud)       :: Indices,Jndices
    !local MPI
    integer                          :: irank
    !
    if(MpiComm==Mpi_Comm_Null)return
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    do i=1,Nloc                 
       i_el = mod(i-1,DimUp*MpiQdw) + 1
       !
       do j=1,spH0d%row(i_el)%Size
          Hv(i) = Hv(i) + spH0d%row(i_el)%dvals(j)*v(i)
       end do
    end do
    !
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    do iph=1,DimPh
       do iidw=1,MpiQdw
          do iiup=1,DimUp
             i = iiup + (iidw-1)*DimUp
             call state2indices(i,[DimUps,DimDws],Indices)
             i = i + (iph-1)*DimUp*MpiQdw
             do iud=1,Ns_Ud
                !
                iup = Indices(iud)
                hxv_up: do jj=1,spH0ups(iud)%row(iup)%Size
                   Jndices      = Indices
                   Jndices(iud) = spH0ups(iud)%row(iup)%cols(jj)
                   call indices2state(Jndices,[DimUps,DimDws],j)
                   !
                   j = j + (iph-1)*DimUp*MpiQdw
                   Hv(i) = Hv(i) + spH0ups(iud)%row(iup)%dvals(jj)*v(j)
                end do hxv_up
                !
             enddo
          enddo
       enddo
    end do
    !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    !
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    do iph=1,DimPh
       allocate(vt(mpiQup*DimDw)) ;vt=0d0
       allocate(Hvt(mpiQup*DimDw));Hvt=0d0
       i_start = 1 + (iph-1)*DimUp*MpiQdw
       i_end = iph*DimUp*MpiQdw
       !
       call vector_transpose_MPI(DimUp,MpiQdw,v(i_start:i_end),DimDw,MpiQup,vt)
       Hvt=0d0    
       do iidw=1,MpiQup            !<= Transposed order:  column-wise DW <--> UP  
          do iiup=1,DimDw          !<= Transposed order:  column-wise DW <--> UP  
             i = iiup + (iidw-1)*DimDw
             call state2indices(i,[DimDws,DimUps],Indices)
             do iud=1,Ns_Ud
                !
                iup = Indices(iud)
                hxv_dw: do jj=1,spH0dws(iud)%row(iup)%Size
                   Jndices      = Indices
                   Jndices(iud) = spH0dws(iud)%row(iup)%cols(jj)
                   call indices2state(Jndices,[DimDws,DimUps],j)
                   Hvt(i) = Hvt(i) + spH0dws(iud)%row(iup)%dvals(jj)*vt(j)
                end do hxv_dw
                !
             enddo
          enddo
       end do
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ; vt=0d0
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
       Hv(i_start:i_end) = Hv(i_start:i_end) + vt
       deallocate(vt,Hvt)
    enddo
    !
    if(DimPh>1)then
       do iph=1,DimPh
          do i_el = 1,DimUp*MpiQdw
             i = i_el + (iph-1)*DimUp*MpiQdw
             !
             !PHONON
             do jj = 1,spH0_ph%row(iph)%Size
                val = spH0_ph%row(iph)%dvals(jj)
                j = i_el + (spH0_ph%row(iph)%cols(jj)-1)*DimUp*MpiQdw
                Hv(i) = Hv(i) + val*v(j)
             enddo
             !
             !ELECTRON-PHONON
             do j_el = 1,spH0e_eph%row(i_el)%Size
                do jj = 1,spH0ph_eph%row(iph)%Size
                   val = spH0e_eph%row(i_el)%dvals(j_el)*&
                        spH0ph_eph%row(iph)%dvals(jj)
                   !interaction is diag from the electron point of view (coupling to the density)
                   j = i_el + (spH0ph_eph%row(iph)%cols(jj)-1)*DimUp*MpiQdw
                   Hv(i) = Hv(i) + val*v(j)
                enddo
             enddo
             !
          enddo
       enddo
    end if
  end subroutine spMatVec_mpi_normal_orbs
#endif



end MODULE ED_HAMILTONIAN_NORMAL_STORED_HXV







