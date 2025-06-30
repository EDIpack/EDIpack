function fdelta_bath_array_hybrid(x,axis) result(Fdelta)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                        :: x !complex  array for the frequency
  character(len=*),optional                                 :: axis    !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))       :: Fdelta
  character(len=1)                                          :: axis_
  integer                                                   :: iorb,ispin,jorb,ibath
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: Vk
  real(8),dimension(Nbath)                                  :: eps,dps
  real(8),dimension(Norb,Nbath)                             :: vops
  integer                                                   :: i,L
  complex(8)                                                :: detH(Nbath),zz(2)
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  Fdelta=zero
  !
  L = size(x)
  !
  do ispin=1,Nspin
     eps  = dmft_bath%e(ispin,1     ,1:Nbath)
     dps  = dmft_bath%d(ispin,1     ,1:Nbath)
     vops = dmft_bath%v(ispin,1:Norb,1:Nbath)
     do iorb=1,Norb
        do jorb=1,Norb
           select case(axis_)
           case default ;stop "fdelta_bath_array_hybrid error: axis not supported"
           case ("m")
              do i=1,L
                 zz(1) = x(i)
                 zz(2) = x(i)
                 detH  = dps(:)**2 + (eps(:)+zz(2))*(eps(:)-zz(1))
                 Fdelta(ispin,ispin,iorb,jorb,i) = &
                      sum( vops(iorb,:)*vops(jorb,:)*dps(:)/detH)
                 ! Fdelta(ispin,ispin,iorb,jorb,i) = &
                 !      sum( dps(:)*vops(iorb,:)*vops(jorb,:)/(dimag(x(i))**2+eps(:)**2+dps(:)**2) )
              enddo
           case ("r")
              do i=1,L
                 zz(1) = x(i)
                 zz(2) = -conjg(x(L-i+1))
                 detH  = dps(:)**2 + (eps(:)+zz(2))*(eps(:)-zz(1))
                 Fdelta(ispin,ispin,iorb,jorb,i) = &
                      sum( vops(iorb,:)*vops(jorb,:)*dps(:)/detH)
                 ! Fdelta(ispin,ispin,iorb,jorb,i) = &
                 !      sum( dps(:)*vops(iorb,:)*vops(jorb,:)/(x(i)*(-x(i)) + eps(:)**2 + dps(:)**2) )
              enddo
           end select
        enddo
     enddo
  enddo
  !
end function fdelta_bath_array_hybrid
