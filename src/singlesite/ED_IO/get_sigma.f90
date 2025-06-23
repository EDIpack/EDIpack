subroutine ed_get_sigma_site_n3(self,axis,type,z)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout)   :: self ! Green's function matrix
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                   :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  character(len=1)                            :: type_
  complex(8),dimension(:),allocatable         :: z_
  complex(8),dimension(:,:,:,:,:),allocatable :: gf
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  !
  call allocate_grids
  if(.not.dmft_bath%status)call read_dmft_bath()
  if(.not.allocated(impGmatrix))call read_impGmatrix()
  if(.not.allocated(impHloc))stop "ed_get_sigma error: Hloc not allocated. Either call ed_set_Hloc or ed_init_solver before ed_get_sigma"
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
  call assert_shape(self,[Nspin*Norb,Nspin*Norb,L],'ed_get_sigma','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  gf = zero
  !
  select case(type_)
  case default; stop "ed_get_sigma ERROR: type is neither Normal, nor Anomalous"
  case ('n','N');gf = get_Sigma(z_,axis_)
  case ('a','A');gf = get_Self(z_,axis_)
  end select
  !
  self = nn2so_reshape( gf, Nspin,Norb,L)
  !
  call deallocate_grids
  call deallocate_dmft_bath()
  !
end subroutine ed_get_sigma_site_n3

subroutine ed_get_sigma_site_n5(self,axis,type,z)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  type_='n';if(present(type))type_=trim(type)
  !
  !
  call allocate_grids
  if(.not.dmft_bath%status)call read_dmft_bath()
  if(.not.allocated(impGmatrix))call read_impGmatrix()
  if(.not.allocated(impHloc))stop "ed_get_sigma error: Hloc not allocated. Either call ed_set_Hloc or ed_init_solver before ed_get_sigma"
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
  call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],'ed_get_sigma','self')
  self = zero
  !  
  select case(type_)
  case default; stop "ed_get_sigma ERROR: type is neither Normal, nor Anomalous"
  case ('n','N')
     self = get_Sigma(z_,axis_)
  case('a','A')
     self = get_Self(z_,axis_)
  end select
  !
  call deallocate_grids
  call deallocate_dmft_bath()
  !
end subroutine ed_get_sigma_site_n5


