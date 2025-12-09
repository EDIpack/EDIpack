function f0and_bath_array_normal(x,axis) result(F0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x !complex  array for the frequency
  character(len=*),optional                           :: axis!string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  character(len=1)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  !
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta,Fdelta
  integer                                             :: iorb,jorb,ispin,i,L
  real(8),dimension(size(x))                          :: ddet
  complex(8),dimension(size(x))                       :: cdet
  complex(8),dimension(size(x))                       :: fg,ff
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  F0and=zero
  !
  L = size(x)
  !
  Delta =  delta_bath_array(x,axis_)
  Fdelta= fdelta_bath_array(x,axis_)
  select case(axis_)
  case default;stop "f0and_bath_array_normal error: axis_ not support"
  case ("m")
     do ispin=1,Nspin
        do iorb=1,Norb
           fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
           ff(:) = -impHloc_anomalous(ispin,ispin,iorb,iorb)   - Fdelta(ispin,ispin,iorb,iorb,:)
           ddet(:)= abs(fg(:))**2 + ff(:)*ff(:)
           F0and(ispin,ispin,iorb,iorb,:) = ff(:)/ddet(:)
        enddo
     enddo
  case ("r")
     do ispin=1,Nspin
        do iorb=1,Norb
           fg(:)  =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
           ff(:)  =  -impHloc_anomalous(ispin,ispin,iorb,iorb)   - Fdelta(ispin,ispin,iorb,iorb,:)
           cdet(:) =  -fg(:)*conjg(fg(L:1:-1)) - ff(:)*ff(:)
           F0and(ispin,ispin,iorb,iorb,:) = -ff(:)/cdet(:)
        enddo
     enddo
  end select
  !
end function f0and_bath_array_normal
