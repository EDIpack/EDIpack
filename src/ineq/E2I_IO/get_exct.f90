subroutine ed_get_exct_n1(self,Nlat,component,iorb,jorb)
  real(8),dimension(:) :: self ! :math:`[S_0,T_x,T_y,T_z]` value or array of values
  integer              :: Nlat !Number of inequivalent impurities for real-space DMFT
  character(len=1)     :: component
  integer,optional     :: iorb ! first orbital index
  integer,optional     :: jorb ! second orbital index
  integer              :: iorb_,jorb_,id
  self=0d0
  if(Norb<2)return
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=2;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_exct error: orbital index > N_orbital"
  if(.not.allocated(exct_ineq))stop "ed_get_exct error: exct_ineq not allocated"
  select case(component)     
  case default;stop "ed_get_exct error: wrong component, not in [0,x,y,z]"
  case('0','S','s');id=1
  case('X','x');id=2
  case('Y','y');id=3
  case('Z','z');id=4
  end select
  call assert_shape(self,[Nlat],'ed_get_exct','exct')
  self = exct_ineq(:,id,iorb_,jorb_)
end subroutine ed_get_exct_n1

!X_ab
subroutine ed_get_exct_n3(self,Nlat,component)
  real(8),dimension(:,:,:) :: self
  integer                  :: Nlat !Number of inequivalent impurities for real-space DMFT
  character(len=1)         :: component
  integer                  :: id
  if(.not.allocated(exct_ineq))stop "ed_get_exct error: exct_ineq not allocated"
  select case(component)
  case default;stop "ed_get_exct error: wrong component, not in [0,x,y,z]"
  case('0','S','s');id=1
  case('X','x');id=2
  case('Y','y');id=3
  case('Z','z');id=4
  end select
  call assert_shape(self,[Nlat,Norb,Norb],'ed_get_exct','exct')
  self = exct_ineq(:,id,:,:)
end subroutine ed_get_exct_n3


!\vecX_ab
subroutine ed_get_exct_n4(self,Nlat)
  real(8),dimension(:,:,:,:) :: self
  integer                    :: Nlat
  if(.not.allocated(exct_ineq))stop "ed_get_exct error: exct_ineq not allocated"
  call assert_shape(self,[Nlat,4,Norb,Norb],'ed_get_exct','exct')
  self = exct_ineq
end subroutine ed_get_exct_n4


