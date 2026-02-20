subroutine ed_get_g0imp_lattice_n3(self,bath,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout)     :: self !! [Nlso,Nlso,:]
  real(8),dimension(:,:)                        :: bath
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat,Nlat
  complex(8),dimension(:,:,:,:,:,:),allocatable :: gf
  integer                                       :: MPI_ID=0
  integer                                       :: MPI_SIZE=1
  integer                                       :: mpi_err 
#ifdef _MPI    
  if(check_MPI())then
     MPI_ID     = get_Rank_MPI()
     MPI_SIZE   = get_Size_MPI()
  endif
#endif
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        allocate(z_(lreal))
        z_=dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  Nlat = size(bath,1)
  call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,L],'ed_get_g0imp','self')
  !
  allocate(gf(Nlat,Nspin,Nspin,Norb,Norb,L))
  gf = zero
  !
  do ilat = 1 + MPI_ID, Nlat, MPI_SIZE
     call ed_set_suffix(ilat)
     call set_impHloc(ilat)
     call ed_get_g0and(z_,bath(ilat,:),gf(ilat,:,:,:,:,:),axis=axis_,type=type_)
  enddo
  !
#ifdef _MPI
  if(check_MPI())call MPI_AllReduce(MPI_IN_PLACE, gf, size(gf), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ERR)
#endif
  !
  self = nnn2lso_reshape(gf,Nlat,Nspin,Norb,L)
  !
  call ed_reset_suffix()
  call deallocate_grids()
  deallocate(gf)
  !
end subroutine ed_get_g0imp_lattice_n3

subroutine ed_get_g0imp_lattice_n4(self,bath,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:),intent(inout) :: self !! [Nlat,Nso,Nso,:]
  real(8),dimension(:,:)                      :: bath
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  character(len=1)                            :: type_
  complex(8),dimension(:),allocatable         :: z_
  integer                                     :: ilat,Nlat
  complex(8),dimension(:,:,:,:,:),allocatable :: gf
  integer                                       :: MPI_ID=0
  integer                                       :: MPI_SIZE=1
  integer                                       :: mpi_err 
#ifdef _MPI    
  if(check_MPI())then
     MPI_ID     = get_Rank_MPI()
     MPI_SIZE   = get_Size_MPI()
  endif
#endif
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        allocate(z_(lreal))
        z_=dcmplx(wr,eps)
     end select
  endif
  !
  L    = size(z_)
  Nlat = size(bath,1)
  !
  call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,L],'ed_get_g0imp','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  gf = zero
  do ilat = 1 + MPI_ID, Nlat, MPI_SIZE
     call ed_set_suffix(ilat)
     call set_impHloc(ilat)
     call ed_get_g0and(z_,bath(ilat,:),gf,axis=axis_,type=type_)
     !
     self(ilat,:,:,:) = nn2so_reshape(gf,Nspin,Norb,L)
     !
  enddo
  !
#ifdef _MPI
  if(check_MPI())call MPI_AllReduce(MPI_IN_PLACE, self, size(self), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ERR)
#endif
  !
  call ed_reset_suffix()
  call deallocate_grids()
  deallocate(gf)
  !
end subroutine ed_get_g0imp_lattice_n4

subroutine ed_get_g0imp_lattice_n6(self,bath,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: self
  real(8),dimension(:,:)                          :: bath
  character(len=*),optional                       :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                       :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional                :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                                :: axis_
  character(len=1)                                :: type_
  complex(8),dimension(:),allocatable             :: z_
  integer                                         :: ilat,Nlat
  integer                                       :: MPI_ID=0
  integer                                       :: MPI_SIZE=1
  integer                                       :: mpi_err 
#ifdef _MPI    
  if(check_MPI())then
     MPI_ID     = get_Rank_MPI()
     MPI_SIZE   = get_Size_MPI()
  endif
#endif
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        allocate(z_(lreal))
        z_=dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  Nlat=size(bath,1)
  !
  call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,L],'ed_get_g0imp','self')
  !
  do ilat = 1 + MPI_ID, Nlat, MPI_SIZE
     call ed_set_suffix(ilat)
     call set_impHloc(ilat)
     call ed_get_g0and(z_,bath(ilat,:),self(ilat,:,:,:,:,:),axis=axis_,type=type_)
  enddo
  !
#ifdef _MPI
  if(check_MPI())call MPI_AllReduce(MPI_IN_PLACE, self, size(self), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ERR)
#endif
  !
  call ed_reset_suffix()
  call deallocate_grids()
  !
end subroutine ed_get_g0imp_lattice_n6



