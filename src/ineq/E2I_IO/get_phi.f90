!phi_i(a_,a_)
subroutine ed_get_phisc_n1(self,Nlat,iorb)
  real(8),dimension(:) :: self
  integer              :: Nlat ! number of inequivalent impurity sites for real-space DMFT
  integer              :: iorb
  if(iorb>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  if(.not.allocated(phisc_ineq))stop "ed_get_phisc error: phisc_ineq not allocated"
  if(Nlat>size(phisc_ineq,1))stop "ed_get_phisc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_phisc','phisc')
  self = phisc_ineq(:,iorb,iorb)
end subroutine ed_get_phisc_n1


!phi_i,ab
subroutine ed_get_phisc_n3(self,Nlat)
  real(8),dimension(:,:,:) :: self
  integer                  :: Nlat
  if(.not.allocated(phisc_ineq))stop "ed_get_phisc error: phisc_ineq not allocated"
  if(Nlat>size(phisc_ineq,1))stop "ed_get_phisc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb,Norb],'ed_get_phisc','phisc')
  self = phisc_ineq
end subroutine ed_get_phisc_n3











!arg_i(a_,a_)
subroutine ed_get_argsc_n1(self,Nlat,iorb)
  real(8),dimension(:) :: self
  integer              :: Nlat ! number of inequivalent impurity sites for real-space DMFT
  integer              :: iorb
  if(iorb>Norb)stop "ed_get_argsc error: orbital index > N_orbital"
  if(.not.allocated(argsc_ineq))stop "ed_get_argsc error: argsc_ineq not allocated"
  if(Nlat>size(argsc_ineq,1))stop "ed_get_argsc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_argsc','argsc')
  self = argsc_ineq(:,iorb,iorb)
end subroutine ed_get_argsc_n1


!arg_i,ab
subroutine ed_get_argsc_n3(self,Nlat)
  real(8),dimension(:,:,:) :: self
  integer                  :: Nlat
  if(.not.allocated(argsc_ineq))stop "ed_get_argsc error: argsc_ineq not allocated"
  if(Nlat>size(argsc_ineq,1))stop "ed_get_argsc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb,Norb],'ed_get_argsc','argsc')
  self = argsc_ineq
end subroutine ed_get_argsc_n3

