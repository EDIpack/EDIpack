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



subroutine ed_get_Kc_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Kc))stop "ed_get_Kc error: ed_Kc not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Kc error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Kc','Kc')
  self = ed_Kc(ispin,iorb,:)
end subroutine ed_get_Kc_n1

subroutine ed_get_Kc_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Kc))stop "ed_get_Kc error: ed_Kc not allocated"
  call assert_shape(self,shape(ed_Kc),'ed_get_Kc','Kc')
  self = ed_Kc
end subroutine ed_get_Kc_n3



subroutine ed_get_Sc_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Sc))stop "ed_get_Sc error: ed_Sc not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Sc error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Sc','Sc')
  self = ed_Sc(ispin,iorb,:)
end subroutine ed_get_Sc_n1

subroutine ed_get_Sc_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Sc))stop "ed_get_Sc error: ed_Sc not allocated"
  call assert_shape(self,shape(ed_Sc),'ed_get_Sc','Sc')
  self = ed_Sc
end subroutine ed_get_Sc_n3





subroutine ed_get_Pc_n2(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: ispin,iorb
  if(.not.allocated(ed_Pc))stop "ed_get_Pc error: ed_Pc not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Pc error: index out of bounds"
  call assert_shape(self,[size(ed_Pc,3),Ltimes],'ed_get_Pc','Pc')
  self = ed_Pc(ispin,iorb,:,:)
end subroutine ed_get_Pc_n2

subroutine ed_get_Pc_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_Pc))stop "ed_get_Pc error: ed_Pc not allocated"
  call assert_shape(self,shape(ed_Pc),'ed_get_Pc','Pc')
  self = ed_Pc
end subroutine ed_get_Pc_n4


subroutine ed_get_Kg1_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Kg1))stop "ed_get_Kg1 error: ed_Kg1 not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Kg1 error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Kg1','Kg1')
  self = ed_Kg1(ispin,iorb,:)
end subroutine ed_get_Kg1_n1

subroutine ed_get_Kg1_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Kg1))stop "ed_get_Kg1 error: ed_Kg1 not allocated"
  call assert_shape(self,shape(ed_Kg1),'ed_get_Kg1','Kg1')
  self = ed_Kg1
end subroutine ed_get_Kg1_n3


subroutine ed_get_Sg1_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_Sg1))stop "ed_get_Sg1 error: ed_Sg1 not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Sg1 error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_Sg1','Sg1')
  self = ed_Sg1(ispin,iorb,:)
end subroutine ed_get_Sg1_n1

subroutine ed_get_Sg1_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_Sg1))stop "ed_get_Sg1 error: ed_Sg1 not allocated"
  call assert_shape(self,shape(ed_Sg1),'ed_get_Sg1','Sg1')
  self = ed_Sg1
end subroutine ed_get_Sg1_n3



subroutine ed_get_Pg1_n2(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: ispin,iorb
  if(.not.allocated(ed_Pg1))stop "ed_get_Pg1 error: ed_Pg1 not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_Pg1 error: index out of bounds"
  call assert_shape(self,[size(ed_Pg1,3),Ltimes],'ed_get_Pg1','Pg1')
  self = ed_Pg1(ispin,iorb,:,:)
end subroutine ed_get_Pg1_n2

subroutine ed_get_Pg1_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_Pg1))stop "ed_get_Pg1 error: ed_Pg1 not allocated"
  call assert_shape(self,shape(ed_Pg1),'ed_get_Pg1','Pg1')
  self = ed_Pg1
end subroutine ed_get_Pg1_n4


subroutine ed_get_KOCcdg_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_KOCcdg))stop "ed_get_KOCcdg error: ed_KOCcdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_KOCcdg error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_KOCcdg','KOCcdg')
  self = ed_KOCcdg(ispin,iorb,:)
end subroutine ed_get_KOCcdg_n1

subroutine ed_get_KOCcdg_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_KOCcdg))stop "ed_get_KOCcdg error: ed_KOCcdg not allocated"
  call assert_shape(self,shape(ed_KOCcdg),'ed_get_KOCcdg','KOCcdg')
  self = ed_KOCcdg
end subroutine ed_get_KOCcdg_n3

subroutine ed_get_SOCcdg_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_SOCcdg))stop "ed_get_SOCcdg error: ed_SOCcdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_SOCcdg error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_SOCcdg','SOCcdg')
  self = ed_SOCcdg(ispin,iorb,:)
end subroutine ed_get_SOCcdg_n1

subroutine ed_get_SOCcdg_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_SOCcdg))stop "ed_get_SOCcdg error: ed_SOCcdg not allocated"
  call assert_shape(self,shape(ed_SOCcdg),'ed_get_SOCcdg','SOCcdg')
  self = ed_SOCcdg
end subroutine ed_get_SOCcdg_n3


subroutine ed_get_POCcdg_n2(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: ispin,iorb
  if(.not.allocated(ed_POCcdg))stop "ed_get_POCcdg error: ed_POCcdg not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_POCcdg error: index out of bounds"
  call assert_shape(self,[size(ed_POCcdg,3),Ltimes],'ed_get_POCcdg','POCcdg')
  self = ed_POCcdg(ispin,iorb,:,:)
end subroutine ed_get_POCcdg_n2

subroutine ed_get_POCcdg_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_POCcdg))stop "ed_get_POCcdg error: ed_POCcdg not allocated"
  call assert_shape(self,shape(ed_POCcdg),'ed_get_POCcdg','POCcdg')
  self = ed_POCcdg
end subroutine ed_get_POCcdg_n4

subroutine ed_get_KOCc_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_KOCc))stop "ed_get_KOCc error: ed_KOCc not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_KOCc error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_KOCc','KOCc')
  self = ed_KOCc(ispin,iorb,:)
end subroutine ed_get_KOCc_n1

subroutine ed_get_KOCc_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_KOCc))stop "ed_get_KOCc error: ed_KOCc not allocated"
  call assert_shape(self,shape(ed_KOCc),'ed_get_KOCc','KOCc')
  self = ed_KOCc
end subroutine ed_get_KOCc_n3

subroutine ed_get_SOCc_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_SOCc))stop "ed_get_SOCc error: ed_SOCc not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_SOCc error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_SOCc','SOCc')
  self = ed_SOCc(ispin,iorb,:)
end subroutine ed_get_SOCc_n1

subroutine ed_get_SOCc_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_SOCc))stop "ed_get_SOCc error: ed_SOCc not allocated"
  call assert_shape(self,shape(ed_SOCc),'ed_get_SOCc','SOCc')
  self = ed_SOCc
end subroutine ed_get_SOCc_n3


subroutine ed_get_POCc_n2(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: ispin,iorb
  if(.not.allocated(ed_POCc))stop "ed_get_POCc error: ed_POCc not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_POCc error: index out of bounds"
  call assert_shape(self,[size(ed_POCc,3),Ltimes],'ed_get_POCc','POCc')
  self = ed_POCc(ispin,iorb,:,:)
end subroutine ed_get_POCc_n2

subroutine ed_get_POCc_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_POCc))stop "ed_get_POCc error: ed_POCc not allocated"
  call assert_shape(self,shape(ed_POCc),'ed_get_POCc','POCc')
  self = ed_POCc
end subroutine ed_get_POCc_n4


subroutine ed_get_KOCg1_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_KOCg1))stop "ed_get_KOCg1 error: ed_KOCg1 not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_KOCg1 error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_KOCg1','KOCg1')
  self = ed_KOCg1(ispin,iorb,:)
end subroutine ed_get_KOCg1_n1

subroutine ed_get_KOCg1_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_KOCg1))stop "ed_get_KOCg1 error: ed_KOCg1 not allocated"
  call assert_shape(self,shape(ed_KOCg1),'ed_get_KOCg1','KOCg1')
  self = ed_KOCg1
end subroutine ed_get_KOCg1_n3

subroutine ed_get_SOCg1_n1(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:) :: self
  integer,intent(in)   :: ispin,iorb
  if(.not.allocated(ed_SOCg1))stop "ed_get_SOCg1 error: ed_SOCg1 not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_SOCg1 error: index out of bounds"
  call assert_shape(self,[Ltimes],'ed_get_SOCg1','SOCg1')
  self = ed_SOCg1(ispin,iorb,:)
end subroutine ed_get_SOCg1_n1

subroutine ed_get_SOCg1_n3(self)
  real(8),dimension(:,:,:) :: self
  if(.not.allocated(ed_SOCg1))stop "ed_get_SOCg1 error: ed_SOCg1 not allocated"
  call assert_shape(self,shape(ed_SOCg1),'ed_get_SOCg1','SOCg1')
  self = ed_SOCg1
end subroutine ed_get_SOCg1_n3


subroutine ed_get_POCg1_n2(self,ispin,iorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Ltimes
#endif
  real(8),dimension(:,:) :: self
  integer,intent(in)     :: ispin,iorb
  if(.not.allocated(ed_POCg1))stop "ed_get_POCg1 error: ed_POCg1 not allocated"
  if(ispin>Nspin.OR.iorb>Norb)stop "ed_get_POCg1 error: index out of bounds"
  call assert_shape(self,[size(ed_POCg1,3),Ltimes],'ed_get_POCg1','POCg1')
  self = ed_POCg1(ispin,iorb,:,:)
end subroutine ed_get_POCg1_n2

subroutine ed_get_POCg1_n4(self)
  real(8),dimension(:,:,:,:) :: self
  if(.not.allocated(ed_POCg1))stop "ed_get_POCg1 error: ed_POCg1 not allocated"
  call assert_shape(self,shape(ed_POCg1),'ed_get_POCg1','POCg1')
  self = ed_POCg1
end subroutine ed_get_POCg1_n4
