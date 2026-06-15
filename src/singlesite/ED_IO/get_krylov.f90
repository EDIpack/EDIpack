subroutine ed_get_KSC_n1(self,iop,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iop,ispin,iorb
  if(.not.allocated(ed_KSC))stop "ed_get_KSC error: ed_KSC not allocated"
  if(iop<1.OR.iop>size(ed_KSC,1).OR.ispin<1.OR.ispin>Nspin.OR.iorb<1.OR.iorb>Norb)stop "ed_get_KSC error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_KSC','KSC')
  self = ed_KSC(iop,ispin,iorb,:)
end subroutine ed_get_KSC_n1


subroutine ed_get_KSC_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_KSC))stop "ed_get_KSC error: ed_KSC not allocated"
  call assert_shape(self,shape(ed_KSC),'ed_get_KSC','KSC')
  self = ed_KSC
end subroutine ed_get_KSC_n4


subroutine ed_get_SSC_n1(self,iop,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iop,ispin,iorb
  if(.not.allocated(ed_SSC))stop "ed_get_SSC error: ed_SSC not allocated"
  if(iop<1.OR.iop>size(ed_SSC,1).OR.ispin<1.OR.ispin>Nspin.OR.iorb<1.OR.iorb>Norb)stop "ed_get_SSC error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_SSC','SSC')
  self = ed_SSC(iop,ispin,iorb,:)
end subroutine ed_get_SSC_n1


subroutine ed_get_SSC_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_SSC))stop "ed_get_SSC error: ed_SSC not allocated"
  call assert_shape(self,shape(ed_SSC),'ed_get_SSC','SSC')
  self = ed_SSC
end subroutine ed_get_SSC_n4


subroutine ed_get_PSC_n2(self,iop,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: iop,ispin,iorb
  if(.not.allocated(ed_PSC))stop "ed_get_PSC error: ed_PSC not allocated"
  if(iop<1.OR.iop>size(ed_PSC,1).OR.ispin<1.OR.ispin>Nspin.OR.iorb<1.OR.iorb>Norb)stop "ed_get_PSC error: index out of bounds"
  call assert_shape(self,[size(ed_PSC,4),Ltimes],'ed_get_PSC','PSC')
  self = ed_PSC(iop,ispin,iorb,:,:)
end subroutine ed_get_PSC_n2


subroutine ed_get_PSC_n5(self)
  real(8),dimension(:,:,:,:,:) :: self
  if(.not.allocated(ed_PSC))stop "ed_get_PSC error: ed_PSC not allocated"
  call assert_shape(self,shape(ed_PSC),'ed_get_PSC','PSC')
  self = ed_PSC
end subroutine ed_get_PSC_n5


subroutine ed_get_KOC_n1(self,iop,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iop,ispin,iorb
  if(.not.allocated(ed_KOC))stop "ed_get_KOC error: ed_KOC not allocated"
  if(iop<1.OR.iop>size(ed_KOC,1).OR.ispin<1.OR.ispin>Nspin.OR.iorb<1.OR.iorb>Norb)stop "ed_get_KOC error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_KOC','KOC')
  self = ed_KOC(iop,ispin,iorb,:)
end subroutine ed_get_KOC_n1


subroutine ed_get_KOC_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_KOC))stop "ed_get_KOC error: ed_KOC not allocated"
  call assert_shape(self,shape(ed_KOC),'ed_get_KOC','KOC')
  self = ed_KOC
end subroutine ed_get_KOC_n4


subroutine ed_get_SOC_n1(self,iop,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: iop,ispin,iorb
  if(.not.allocated(ed_SOC))stop "ed_get_SOC error: ed_SOC not allocated"
  if(iop<1.OR.iop>size(ed_SOC,1).OR.ispin<1.OR.ispin>Nspin.OR.iorb<1.OR.iorb>Norb)stop "ed_get_SOC error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_SOC','SOC')
  self = ed_SOC(iop,ispin,iorb,:)
end subroutine ed_get_SOC_n1


subroutine ed_get_SOC_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_SOC))stop "ed_get_SOC error: ed_SOC not allocated"
  call assert_shape(self,shape(ed_SOC),'ed_get_SOC','SOC')
  self = ed_SOC
end subroutine ed_get_SOC_n4


subroutine ed_get_POC_n2(self,iop,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: iop,ispin,iorb
  if(.not.allocated(ed_POC))stop "ed_get_POC error: ed_POC not allocated"
  if(iop<1.OR.iop>size(ed_POC,1).OR.ispin<1.OR.ispin>Nspin.OR.iorb<1.OR.iorb>Norb)stop "ed_get_POC error: index out of bounds"
  call assert_shape(self,[size(ed_POC,4),Ltimes],'ed_get_POC','POC')
  self = ed_POC(iop,ispin,iorb,:,:)
end subroutine ed_get_POC_n2


subroutine ed_get_POC_n5(self)
  real(8),dimension(:,:,:,:,:) :: self
  if(.not.allocated(ed_POC))stop "ed_get_POC error: ed_POC not allocated"
  call assert_shape(self,shape(ed_POC),'ed_get_POC','POC')
  self = ed_POC
end subroutine ed_get_POC_n5
