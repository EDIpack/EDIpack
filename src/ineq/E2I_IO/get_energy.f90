subroutine ed_get_eimp_n2(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  integer                :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(e_ineq))stop "ed_get_eimp error: e_ineq not allocated"
  if(Nlat>size(e_ineq,1))stop "ed_get_eimp error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,4],'ed_get_eimp','eimp')
  self = e_ineq
end subroutine ed_get_eimp_n2


subroutine ed_get_epot_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DM
  !
  if(.not.allocated(e_ineq))stop "ed_get_epot error: e_ineq not allocated"
  if(Nlat>size(e_ineq,1))stop "ed_get_epot error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_epot','epot')
  self = e_ineq(:,1)
end subroutine ed_get_epot_n1

subroutine ed_get_eint_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DM
  !
  if(.not.allocated(e_ineq))stop "ed_get_eint error: e_ineq not allocated"
  if(Nlat>size(e_ineq,1))stop "ed_get_eint error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_eint','eint')
  self = e_ineq(:,2)
end subroutine ed_get_eint_n1

subroutine ed_get_ehartree_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DM
  !
  if(.not.allocated(e_ineq))stop "ed_get_ehartree error: e_ineq not allocated"
  if(Nlat>size(e_ineq,1))stop "ed_get_ehartree error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_ehartree','ehartree')
  self = e_ineq(:,3)
end subroutine ed_get_ehartree_n1


subroutine ed_get_eknot_n1(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DM
  !
  if(.not.allocated(e_ineq))stop "ed_get_eknot error: e_ineq not allocated"
  if(Nlat>size(e_ineq,1))stop "ed_get_eknot error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_eknot','eknot')
  self = e_ineq(:,4)
end subroutine ed_get_eknot_n1


