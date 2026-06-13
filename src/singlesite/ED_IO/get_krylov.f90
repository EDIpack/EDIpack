subroutine ed_get_Kcdg_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Kcdg))stop "ed_get_Kcdg error: ed_Kcdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Kcdg error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Kcdg','Kcdg')
  self = ed_Kcdg(ispin,iorb,:)
end subroutine ed_get_Kcdg_n1

subroutine ed_get_Kcdg_n3(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Kcdg))stop "ed_get_Kcdg error: ed_Kcdg not allocated"
  call assert_shape(self,shape(ed_Kcdg),'ed_get_Kcdg','Kcdg')
  self = ed_Kcdg
end subroutine ed_get_Kcdg_n3






subroutine ed_get_Scdg_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Scdg))stop "ed_get_Scdg error: ed_Scdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Scdg error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Scdg','Scdg')
  self = ed_Scdg(ispin,iorb,:)
end subroutine ed_get_Scdg_n1

subroutine ed_get_Scdg_n3(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Scdg))stop "ed_get_Scdg error: ed_Scdg not allocated"
  call assert_shape(self,shape(ed_Scdg),'ed_get_Scdg','Scdg')
  self = ed_Scdg
end subroutine ed_get_Scdg_n3







subroutine ed_get_Normcdg_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Normcdg))stop "ed_get_Normcdg error: ed_Normcdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Normcdg error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Normcdg','Normcdg')
  self = ed_Normcdg(ispin,iorb,:)
end subroutine ed_get_Normcdg_n1

subroutine ed_get_Normcdg_n3(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Normcdg))stop "ed_get_Normcdg error: ed_Normcdg not allocated"
  call assert_shape(self,shape(ed_Normcdg),'ed_get_Normcdg','Normcdg')
  self = ed_Normcdg
end subroutine ed_get_Normcdg_n3






subroutine ed_get_Pcdg_n2(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: ispin,iorb
  if(.not.allocated(ed_Pcdg))stop "ed_get_Pcdg error: ed_Pcdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Pcdg error: index out of bounds"
  call assert_shape(self,[size(ed_Pcdg,3),Ltimes],'ed_get_Pcdg','Pcdg')
  self = ed_Pcdg(ispin,iorb,:,:)
end subroutine ed_get_Pcdg_n2

subroutine ed_get_Pcdg_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_Pcdg))stop "ed_get_Pcdg error: ed_Pcdg not allocated"
  call assert_shape(self,shape(ed_Pcdg),'ed_get_Pcdg','Pcdg')
  self = ed_Pcdg
end subroutine ed_get_Pcdg_n4




!##################################################
!##################################################
!##################################################



subroutine ed_get_Kn_n1(self,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iorb
  if(.not.allocated(ed_Kn))stop "ed_get_Kn error: ed_Kn not allocated"
  if(iorb>Norb)stop "ed_get_Kn error: orbital index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Kn','Kn')
  self = ed_Kn(1,iorb,:)
end subroutine ed_get_Kn_n1

subroutine ed_get_Kn_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Kn))stop "ed_get_Kn error: ed_Kn not allocated"
  call assert_shape(self,shape(ed_Kn),'ed_get_Kn','Kn')
  self = ed_Kn
end subroutine ed_get_Kn_n3



subroutine ed_get_Sn_n1(self,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iorb
  if(.not.allocated(ed_Sn))stop "ed_get_Sn error: ed_Sn not allocated"
  if(iorb>Norb)stop "ed_get_Sn error: orbital index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Sn','Sn')
  self = ed_Sn(1,iorb,:)
end subroutine ed_get_Sn_n1

subroutine ed_get_Sn_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Sn))stop "ed_get_Sn error: ed_Sn not allocated"
  call assert_shape(self,shape(ed_Sn),'ed_get_Sn','Sn')
  self = ed_Sn
end subroutine ed_get_Sn_n3



subroutine ed_get_Normn_n1(self,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iorb
  if(.not.allocated(ed_Normn))stop "ed_get_Normn error: ed_Normn not allocated"
  if(iorb>Norb)stop "ed_get_Normn error: orbital index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Normn','Normn')
  self = ed_Normn(1,iorb,:)
end subroutine ed_get_Normn_n1

subroutine ed_get_Normn_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Normn))stop "ed_get_Normn error: ed_Normn not allocated"
  call assert_shape(self,shape(ed_Normn),'ed_get_Normn','Normn')
  self = ed_Normn
end subroutine ed_get_Normn_n3




subroutine ed_get_Pn_n2(self,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: iorb
  if(.not.allocated(ed_Pn))stop "ed_get_Pn error: ed_Pn not allocated"
  if(iorb>Norb)stop "ed_get_Pn error: orbital index out of bounds"
  call assert_shape(self,[size(ed_Pn,3),Ltimes],'ed_get_Pn','Pn')
  self = ed_Pn(1,iorb,:,:)
end subroutine ed_get_Pn_n2

subroutine ed_get_Pn_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_Pn))stop "ed_get_Pn error: ed_Pn not allocated"
  call assert_shape(self,shape(ed_Pn),'ed_get_Pn','Pn')
  self = ed_Pn
end subroutine ed_get_Pn_n4
