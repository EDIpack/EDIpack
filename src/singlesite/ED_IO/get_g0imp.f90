subroutine ed_get_g0imp_site_n3(self,bath,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout)   :: self ! Green's function matrix
  real(8),dimension(:)                        :: bath ! The bath vector
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
        allocate(z_(Lreal))
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin*Norb,Nspin*Norb,L],'ed_get_g0imp','self')
  !
  allocate(gf(Nspin,Nspin,Norb,Norb,L))
  !
  call ed_get_g0and(z_,bath,gf,axis=axis_,type=type_)
  !
  self = nn2so_reshape( gf, Nspin,Norb,L)
  !
  call deallocate_grids
  !
end subroutine ed_get_g0imp_site_n3


subroutine ed_get_g0imp_site_n5(self,bath,axis,type,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self
  real(8),dimension(:)                          :: bath ! The bath vector
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  character(len=*),optional                     :: type ! Can be :f:var:`"n"` for Normal (default), :f:var:`"a"` for anomalous
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  character(len=1)                              :: type_
  complex(8),dimension(:),allocatable           :: z_
  !
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
        allocate(z_(Lreal))
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nspin,Nspin,Norb,Norb,L],'ed_get_g0imp','self')
  !
  call ed_get_g0and(z_,bath,self,axis=axis_,type=type_)
  !
  call deallocate_grids
  !
end subroutine ed_get_g0imp_site_n5

