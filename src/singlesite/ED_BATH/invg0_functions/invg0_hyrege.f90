function invg0_bath_array_hyrege(x,axis) result(G0and)
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
  Nso   = Nspin*Norb
  !
  L=size(x)
  !
  select case(ed_mode)
  case default; stop "invg0_bath_array_hyrege error: ed_mode not supported"
  case ("normal")
     Delta = delta_bath_array(x)
     do ispin=1,Nspin
        do i=1,L
           G0and(ispin,ispin,:,:,i) = (x(i)+xmu)*zeye(Norb)-impHloc(ispin,ispin,:,:)-Delta(ispin,ispin,:,:,i)
        enddo
     enddo
     !
  case ("superc")
     allocate(zeta(Nso,Nso))
     Delta =  delta_bath_array(x,axis_)
     do ispin=1,Nspin
        do i=1,L
          zeta = (x(i)+xmu)*zeye(Nso)
          do iorb=1,Norb
             do jorb=1,Norb
                G0and(ispin,ispin,iorb,jorb,i) = zeta(iorb,jorb) - impHloc(ispin,ispin,iorb,jorb) - Delta(ispin,ispin,iorb,jorb,i)
             enddo
          enddo
       enddo
     enddo
     deallocate(zeta)
     !
  case ("nonsu2")
     Nso=Nspin*Norb
     allocate(zeta(Nso,Nso))
     Delta = delta_bath_array(x)
     do i=1,L
        zeta  = (x(i) + xmu)*zeye(Nso)
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    G0and(ispin,jspin,iorb,jorb,i) = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                 enddo
              enddo
           enddo
        enddo
     enddo
     deallocate(zeta)
  end select
  !
end function invg0_bath_array_hyrege
