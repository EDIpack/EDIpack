subroutine ed_get_phisc_n0(self,iorb,jorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8)          :: self ! :math:`\phi` value or array of values
  integer,optional :: iorb ! first orbital index
  integer,optional :: jorb ! second orbital index
  integer          :: iorb_,jorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=1;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_phisc error: orbital index > N_orbital"
  self = ed_phisc(iorb_,jorb_)
end subroutine ed_get_phisc_n0

!phi_aa
subroutine ed_get_phisc_n1(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  call assert_shape(self,[Norb],'ed_get_phisc','phisc')
  self = diagonal(ed_phisc)
end subroutine ed_get_phisc_n1


!phi_ab
subroutine ed_get_phisc_n2(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  call assert_shape(self,[Norb,Norb],'ed_get_phisc','phisc')
  self = ed_phisc
end subroutine ed_get_phisc_n2

