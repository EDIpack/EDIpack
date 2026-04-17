subroutine ed_get_spinChi_lattice_n3(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout)       :: self !! [Nlat,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
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
  !
  call assert_shape(self,[Nlat,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_read_spinChimatrix()
     self(ilat,:,:,:) =  ed_build_spinChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  ! if(allocated(spinChimatrix))call deallocate_GFmatrix(spinChimatrix)
  ! if(allocated(spinChimatrix))deallocate(spinChimatrix)
  !
end subroutine ed_get_spinChi_lattice_n3





subroutine ed_get_densChi_lattice_n3(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout)       :: self !! [Nlat,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
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
  !
  call assert_shape(self,[Nlat,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_read_densChimatrix()
     self(ilat,:,:,:) =  ed_build_densChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  ! if(allocated(densChimatrix))call deallocate_GFmatrix(densChimatrix)
  ! if(allocated(densChimatrix))deallocate(densChimatrix)
  !
end subroutine ed_get_densChi_lattice_n3



subroutine ed_get_pairChi_lattice_n3(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:),intent(inout)       :: self !! [Nlat,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
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
  !
  call assert_shape(self,[Nlat,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_read_pairChimatrix()
     self(ilat,:,:,:) =  ed_build_pairChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  ! if(allocated(pairChimatrix))call deallocate_GFmatrix(pairChimatrix)
  ! if(allocated(pairChimatrix))deallocate(pairChimatrix)
  !
end subroutine ed_get_pairChi_lattice_n3




subroutine ed_get_exctChi_lattice_n3(self,nlat,axis,z)
  complex(8),dimension(:,:,:,:,:),intent(inout) :: self !! [Nlat,3,Norb,Norb,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
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
  !
  call assert_shape(self,[Nlat,3,Norb,Norb,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_read_exctChimatrix()
     self(ilat,:,:,:,:) =  ed_build_exctChi(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  ! if(allocated(exctChimatrix))call deallocate_GFmatrix(exctChimatrix)
  ! if(allocated(exctChimatrix))deallocate(exctChimatrix)
  !
end subroutine ed_get_exctChi_lattice_n3
