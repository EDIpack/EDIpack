subroutine ed_get_spinChi_site_n3(self,axis,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis, nor Time"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,vm)
     case('r','R')
        allocate(z_(lreal))
        z_ = dcmplx(vr,eps)
     case ('t','T')
        allocate(z_(ltau))
        z_ = dcmplx(tau,0d0)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_spinChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_spinChi_site_n3


subroutine ed_get_densChi_site_n3(self,axis,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis, nor Time"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,vm)
     case('r','R')
        allocate(z_(lreal))
        z_ = dcmplx(vr,eps)
     case ('t','T')
        allocate(z_(ltau))
        z_ = dcmplx(tau,0d0)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_densChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_densChi_site_n3


subroutine ed_get_pairChi_site_n3(self,axis,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                 :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional          :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                          :: axis_
  complex(8),dimension(:),allocatable       :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis, nor Time"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,vm)
     case('r','R')
        allocate(z_(lreal))
        z_ = dcmplx(vr,eps)
     case ('t','T')
        allocate(z_(ltau))
        z_ = dcmplx(tau,0d0)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_pairChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_pairChi_site_n3

subroutine ed_get_exctChi_site_n3(self,axis,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:,:,:),intent(inout) :: self ! spin susceptibility 
  character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                            :: axis_
  complex(8),dimension(:),allocatable         :: z_
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_sigma ERROR: axis is neither Matsubara, nor Realaxis, nor Time"
     case ('m','M')
        allocate(z_(lmats))
        z_ = dcmplx(0d0,vm)
     case('r','R')
        allocate(z_(lreal))
        z_ = dcmplx(vr,eps)
     case ('t','T')
        allocate(z_(ltau))
        z_ = dcmplx(tau,0d0)
     end select
  endif
  !
  L = size(z_)
  call assert_shape(self,[3,Norb,Norb,L],'ed_get_spinChi','self')
  !
  self = get_exctChi(z_,axis_)
  call deallocate_grids
  !
end subroutine ed_get_exctChi_site_n3



