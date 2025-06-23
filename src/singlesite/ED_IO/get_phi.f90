subroutine ed_get_phisc_n0(self,arg,iorb,jorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8)          :: self ! :math:`\phi` value or array of values
  real(8),optional :: arg ! :math:`\theta` value or array of values
  integer,optional :: iorb ! first orbital index
  integer,optional :: jorb ! second orbital index
  integer          :: iorb_,jorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=1;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  self = ed_phisc(iorb_,jorb_)
  if(present(arg))arg=ed_argsc(iorb_,jorb_)
end subroutine ed_get_phisc_n0

!phi_aa
subroutine ed_get_phisc_n1(self,arg)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:)          :: self
  real(8),dimension(:),optional :: arg
  call assert_shape(self,[Norb],'ed_get_phisc','phisc')
  self = diagonal(ed_phisc)
  if(present(arg))arg = diagonal(ed_argsc)
end subroutine ed_get_phisc_n1


!phi_ab
subroutine ed_get_phisc_n2(self,arg)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:)          :: self
  real(8),dimension(:,:),optional :: arg
  call assert_shape(self,[Norb,Norb],'ed_get_phisc','phisc')
  self = ed_phisc
  if(present(arg))arg = ed_argsc
end subroutine ed_get_phisc_n2












subroutine ed_get_argsc_n0(self,iorb,jorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8)          :: self ! :math:`\theta` value or array of values
  integer,optional :: iorb ! first orbital index
  integer,optional :: jorb ! second orbital index
  integer          :: iorb_,jorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=1;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_argsc error: orbital index > N_orbital"
  self = ed_argsc(iorb_,jorb_)
end subroutine ed_get_argsc_n0

!arg_aa
subroutine ed_get_argsc_n1(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  call assert_shape(self,[Norb],'ed_get_argsc','argsc')
  self = diagonal(ed_argsc)
end subroutine ed_get_argsc_n1


!arg_ab
subroutine ed_get_argsc_n2(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  call assert_shape(self,[Norb,Norb],'ed_get_argsc','argsc')
  self = ed_argsc
end subroutine ed_get_argsc_n2

