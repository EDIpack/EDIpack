!##################################################################
!LATTICE EXTENSION: need to read to the 
!##################################################################

subroutine ed_get_gimp_lattice_n3(self,nlat,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout)     :: self !! [Nlso,Nlso,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
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
        allocate(z_, source=dcmplx(0d0,wm))
     case ('r','R')
        allocate(z_, source=dcmplx(wr,eps))
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,L],'ed_get_gimp','self')
  !
  allocate(gf(Nlat,Nspin,Nspin,Norb,Norb,L))
  gf = zero
  !
  do ilat = 1 + MPI_ID, Nlat, MPI_SIZE
     call ed_set_suffix(ilat)
     call set_impHloc(ilat)
     call ed_read_impGmatrix()
     select case(type_)
     case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
     case ('n','N');gf(ilat,:,:,:,:,:) = ed_build_impG(z_,axis_)
     case ('a','A');gf(ilat,:,:,:,:,:) = ed_build_impF(z_,axis_)
     end select
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
  ! if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
  ! if(allocated(impGmatrix))deallocate(impGmatrix)
  !
end subroutine ed_get_gimp_lattice_n3

subroutine ed_get_gimp_lattice_n4(self,nlat,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:),intent(inout) :: self !! [Nlat,Nso,Nso,:]
  integer,intent(in)                          :: nlat
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  character(len=1)                            :: type_
  complex(8),dimension(:),allocatable         :: z_
  integer                                     :: ilat
  complex(8),dimension(:,:,:,:,:),allocatable :: gf
  integer                                     :: MPI_ID=0
  integer                                     :: MPI_SIZE=1
  integer                                     :: mpi_err 
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
        allocate(z_, source=dcmplx(0d0,wm))
     case ('r','R')
        allocate(z_, source=dcmplx(wr,eps))
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,Nspin*Norb,Nspin*Norb,L],'ed_get_gimp','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  gf = zero
  do ilat = 1 + MPI_ID, Nlat, MPI_SIZE
     call ed_set_suffix(ilat)
     call set_impHloc(ilat)
     call ed_read_impGmatrix()
     !
     select case(type_)
     case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
     case ('n','N');gf = ed_build_impG(z_,axis_)
     case ('a','A');gf = ed_build_impF(z_,axis_)
     end select
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
  ! if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
  ! if(allocated(impGmatrix))deallocate(impGmatrix)
  !
end subroutine ed_get_gimp_lattice_n4

subroutine ed_get_gimp_lattice_n6(self,nlat,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: self
  integer,intent(in)                              :: nlat
  character(len=*),optional                       :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                       :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional                :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                                :: axis_
  character(len=1)                                :: type_
  complex(8),dimension(:),allocatable             :: z_
  integer                                         :: ilat
  integer                                         :: MPI_ID=0
  integer                                         :: MPI_SIZE=1
  integer                                         :: mpi_err 
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
        allocate(z_, source=dcmplx(0d0,wm))
     case ('r','R')
        allocate(z_, source=dcmplx(wr,eps))
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,Nspin,Nspin,Norb,Norb,L],'ed_get_gimp','self')
  self = zero
  !
  do ilat = 1 + MPI_ID, Nlat, MPI_SIZE
     call ed_set_suffix(ilat)
     call set_impHloc(ilat)
     call ed_read_impGmatrix()
     select case(type_)
     case default; stop "ed_get_gimp ERROR: type is neither Normal, nor Anomalous"
     case ('n','N');self(ilat,:,:,:,:,:) = ed_build_impG(z_,axis_)
     case ('a','A');self(ilat,:,:,:,:,:) = ed_build_impF(z_,axis_)
     end select
  enddo
  !
#ifdef _MPI
  if(check_MPI())call MPI_AllReduce(MPI_IN_PLACE, self, size(self), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ERR)
#endif
  !
  call ed_reset_suffix()
  call deallocate_grids()
  ! if(allocated(impGmatrix))call deallocate_GFmatrix(impGmatrix)
  ! if(allocated(impGmatrix))deallocate(impGmatrix)
  !
end subroutine ed_get_gimp_lattice_n6



