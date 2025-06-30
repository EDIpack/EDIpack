function invf0_bath_array_hyrege(x,axis) result(F0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
  character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: F0and
  character(len=1)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Fdelta
  complex(8),dimension(2*Nspin*Norb,size(x))          :: z
  complex(8),dimension(Nspin*Norb,Nspin*Norb)         :: zeta
  integer                                             :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  F0and=zero
  Nso  = Nspin*Norb
  !
  L = size(x)
  !
  Fdelta= fdelta_bath_array(x,axis_)
  do i=1,L
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              F0and(ispin,ispin,iorb,jorb,i) = -Fdelta(ispin,ispin,iorb,jorb,i)
           enddo
        enddo
     enddo
  enddo
end function invf0_bath_array_hyrege
