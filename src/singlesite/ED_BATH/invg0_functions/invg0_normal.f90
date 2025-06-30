function invg0_bath_array_normal(x,axis) result(G0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
  character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and    
  character(len=1)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  complex(8),dimension(:,:),allocatable               :: fgorb,zeta
  integer                                             :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  G0and = zero
  Nso = Nspin*Norb
  !
  L=size(x)
  !
  select case(ed_mode)
  case default; stop "invg0_bath_array error: ed_mode not supported"
  case ("normal")
     Delta = delta_bath_array(x)
     do ispin=1,Nspin
        do iorb=1,Norb
           G0and(ispin,ispin,iorb,iorb,:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
        enddo
     enddo
     !
  case ("superc")
     Delta =  delta_bath_array(x,axis_)
     ! select case(axis_)
     ! case default ;stop "invg0_bath_array_normal error: axis not supported"         !mats
     ! case ("m")
     !    do ispin=1,Nspin
     !       do iorb=1,Norb
     !          G0and(ispin,ispin,iorb,iorb,:)  =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
     !       enddo
     !    enddo
     ! case("r")
     !    do ispin=1,Nspin
     !       do iorb=1,Norb
     !          G0and(ispin,ispin,iorb,iorb,:) = dreal(x(:)) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
     !       enddo
     !    enddo
     ! end select
     do ispin=1,Nspin
        do iorb=1,Norb
           G0and(ispin,ispin,iorb,iorb,:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
        enddo
     enddo
     !
  case ("nonsu2")
     Delta = delta_bath_array(x)
     allocate(zeta(Nspin,Nspin))
     do i=1,L
        zeta  = (x(i) + xmu)*zeye(Nspin)
        do iorb=1,Norb
           do ispin=1,Nspin
              do jspin=1,Nspin
                 G0and(ispin,jspin,iorb,iorb,i) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
              enddo
           enddo
        enddo
     enddo
     deallocate(zeta)
  end select
  !
end function invg0_bath_array_normal
