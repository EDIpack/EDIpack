function g0and_bath_array_hyrege(x,dmft_bath_,axis) result(G0and)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                  :: x          !complex  array for the frequency
  type(effective_bath)                                :: dmft_bath_ !the current :f:var:`effective_bath` instance
  character(len=*),optional                           :: axis       !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis    
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  character(len=4)                                    :: axis_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta,Fdelta
  integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
  real(8),dimension(size(x))                          :: ddet
  complex(8),dimension(size(x))                       :: cdet
  complex(8),dimension(2*Nspin*Norb,size(x))          :: z
  complex(8),dimension(size(x))                       :: fg,ff
  complex(8),dimension(:,:),allocatable               :: fgorb,zeta
  !
  axis_="mats";if(present(axis))axis_=str(axis)
  !
  G0and = zero
  Nso   = Nspin*Norb
  !
  L=size(x)
  !
  select case(ed_mode)
     !
  case default;stop "g0and_bath_array error: ed_mode not supported"
     !
  case ("normal")
     allocate(fgorb(Norb,Norb),zeta(Norb,Norb))
     Delta = delta_bath_array(x,dmft_bath_)
     do ispin=1,Nspin
        do i=1,L
           fgorb = (x(i)+xmu)*zeye(Norb) - impHloc(ispin,ispin,:,:) - Delta(ispin,ispin,:,:,i)
           call inv(fgorb)
           G0and(ispin,ispin,:,:,i)=fgorb
        enddo
     enddo
     deallocate(fgorb,zeta)
     !
  case ("superc")
     allocate(fgorb(2*Norb,2*Norb),zeta(2*Norb,2*Norb)) !2==Nnambu
     Delta =  delta_bath_array(x,dmft_bath_,axis_)
     Fdelta= fdelta_bath_array(x,dmft_bath_,axis_)
     z     = zeta_superc(x,xmu,axis_)
     select case(axis_)
     case default;stop "g0and_bath_array error: axis_ not supported"
     case ("mats")
        do ispin=1,Nspin   !==1
           do i=1,L
              zeta = diag(z(:,i))
              fgorb= zero
              do iorb=1,Norb
                 do jorb=1,Norb
                    fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + conjg(impHloc(ispin,ispin,iorb,jorb)) + conjg( Delta(ispin,ispin,iorb,jorb,i) )
                 enddo
              enddo
              call inv(fgorb)
              G0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1:Norb)
           enddo
        enddo
     case("real")
        do ispin=1,Nspin   !==1
           do i=1,L
              zeta = diag(z(:,i))
              fgorb= zero
              do iorb=1,Norb
                 do jorb=1,Norb
                    fgorb(iorb,jorb)           = zeta(iorb,jorb)           - impHloc(ispin,ispin,iorb,jorb)  - Delta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb,jorb+Norb)      = zeta(iorb,jorb+Norb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb)      = zeta(iorb+Norb,jorb)                                        - Fdelta(ispin,ispin,iorb,jorb,i)
                    fgorb(iorb+Norb,jorb+Norb) = zeta(iorb+Norb,jorb+Norb) + conjg(impHloc(ispin,ispin,iorb,jorb))  + conjg( Delta(ispin,ispin,iorb,jorb,L-i+1) )
                 enddo
              enddo
              call inv(fgorb)
              G0and(ispin,ispin,:,:,i) = fgorb(1:Norb,1:Norb)
           enddo
        enddo
     end select
     deallocate(fgorb,zeta)
     !
  case ("nonsu2")
     Nso=Nspin*Norb
     allocate(fgorb(Nso,Nso),zeta(Nso,Nso));fgorb=zero
     Delta = delta_bath_array(x,dmft_bath_)
     do i=1,L
        zeta  = (x(i) + xmu)*zeye(Nso)
        fgorb = zero
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    fgorb(io,jo) = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                 enddo
              enddo
           enddo
        enddo
        call inv(fgorb)
        do ispin=1,Nspin
           do jspin=1,Nspin
              do iorb=1,Norb
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
                 enddo
              enddo
           enddo
        enddo
     enddo
     deallocate(fgorb,zeta)
     !
  end select
end function g0and_bath_array_hyrege
