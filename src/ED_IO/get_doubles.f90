subroutine ed_get_doubles_n1(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self !array of two-body terms expectation values
  call assert_shape(self,[4],'ed_get_doubles','doubles')
  self = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
end subroutine ed_get_doubles_n1

subroutine ed_get_doubles_n2(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  integer                :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_doubles error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_doubles error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,4],'ed_get_doubles','doubles')
  self = dd_ineq
end subroutine ed_get_doubles_n2





subroutine ed_get_dust_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dust`
  self = ed_Dust
end subroutine ed_get_dust_n0

subroutine ed_get_dust_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dust error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dust error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dust','dust')
  self = dd_ineq(:,1)
end subroutine ed_get_dust_n1





subroutine ed_get_dund_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dund`
  self = ed_Dund
end subroutine ed_get_dund_n0

subroutine ed_get_dund_n1(self,Nlat)
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dund error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dund error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dund','dund')
  self = dd_ineq(:,2)
end subroutine ed_get_dund_n1



subroutine ed_get_dse_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dse`
  self = ed_Dse
end subroutine ed_get_dse_n0

subroutine ed_get_dse_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dse error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dse error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dse','dse')
  self = dd_ineq(:,3)
end subroutine ed_get_dse_n1




subroutine ed_get_dph_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dph`
  self = ed_Dph
end subroutine ed_get_dph_n0

subroutine ed_get_dph_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dph error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dph error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dph','dph')
  self = dd_ineq(:,4)
end subroutine ed_get_dph_n1
