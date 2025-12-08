function g0and_bath_array_normal(x,axis) result(G0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
  character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  character(len=1)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta,Fdelta
  integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
  real(8),dimension(size(x))                          :: ddet
  complex(8),dimension(size(x))                       :: cdet
  complex(8),dimension(size(x))                       :: fg,ff
  complex(8),dimension(:,:),allocatable               :: fgorb,zeta
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  G0and = zero
  Nso   = Nspin*Norb
  !
  L=size(x)
  !
  select case(ed_mode)
  case default;stop "g0and_bath_array_normal error: ed_mode not supported"
  case ("normal")
     Delta = delta_bath_array(x)
     do ispin=1,Nspin
        do iorb=1,Norb
           fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
           G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
        enddo
     enddo
  case ("superc")
     Delta =  delta_bath_array(x,axis_)
     Fdelta= fdelta_bath_array(x,axis_)
     select case(axis_)
     case default ;stop "g0and_bath_array_normal error: axis not supported"         !mats
     case ("m")
        do ispin=1,Nspin
           do iorb=1,Norb
              fg(:)  = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
              ff(:)  = -pair_field(iorb)                           - Fdelta(ispin,ispin,iorb,iorb,:)
              ddet(:) = abs(fg(:))**2 + ff(:)*ff(:)
              G0and(ispin,ispin,iorb,iorb,:) = conjg(fg(:))/ddet(:)
           enddo
        enddo
     case ("r")
        do ispin=1,Nspin
           do iorb=1,Norb
              fg(:)  =  x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(ispin,ispin,iorb,iorb,:)
              ff(:)  =  -pair_field(iorb)                           - Fdelta(ispin,ispin,iorb,iorb,:)
              cdet(:) =  fg(:)*conjg(fg(L:1:-1)) + ff(:)*ff(:)
              G0and(ispin,ispin,iorb,iorb,:) = conjg(fg(L:1:-1))/cdet(:)
           enddo
        enddo
     end select
  case ("nonsu2")
     allocate(fgorb(Nspin,Nspin),zeta(Nspin,Nspin))
     Delta = delta_bath_array(x)
     do i=1,L
        zeta  = (x(i) + xmu)*zeye(Nspin)
        fgorb = zero
        do iorb=1,Norb
           do ispin=1,Nspin
              do jspin=1,Nspin
                 fgorb(ispin,jspin) = zeta(ispin,jspin) - impHloc(ispin,jspin,iorb,iorb) - Delta(ispin,jspin,iorb,iorb,i)
              enddo
           enddo
           call inv(fgorb)
           do ispin=1,Nspin
              do jspin=1,Nspin
                 G0and(ispin,jspin,iorb,iorb,i) = fgorb(ispin,jspin)
              enddo
           enddo
        enddo
     enddo
     deallocate(fgorb,zeta)
  end select
  !
end function g0and_bath_array_normal
