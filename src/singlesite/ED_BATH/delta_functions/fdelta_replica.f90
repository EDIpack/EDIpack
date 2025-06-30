function fdelta_bath_array_replica(x,axis) result(Fdelta)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
  complex(8),dimension(:),intent(in)                                :: x     !complex  array for the frequency
  character(len=*),optional                                         :: axis  !string indicating the desired axis, :code:`'m'` for Matsubara (default), :code:`'r'` for Real-axis
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x))               :: Fdelta
  character(len=1)                                                  :: axis_
  integer                                                           :: ibath
  integer                                                           :: i,L
  real(8),dimension(Nspin*Norb)                                     :: V
  complex(8),dimension(Nnambu*Nspin*Norb,size(x))                   :: Z
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: Vk
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: Hk
  complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)         :: invH_k
  complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb)         :: invH_knn
  !
  axis_="m";if(present(axis))axis_=str(to_lower(axis))
  !
  Fdelta=zero
  !
  L = size(x)
  !
  Z = zeta_superc(x,0d0,axis_)
  !
  do ibath=1,Nbath
     V  = dmft_bath%item(ibath)%v
     Hk = nn2so_reshape(build_Hreplica(dmft_bath%item(ibath)%lambda),Nnambu*Nspin,Norb)
     Vk = kron( pauli_sigma_z, one*diag(v) )
     do i=1,L
        invH_k   = one*diag(Z(:,i)) - Hk
        call inv(invH_k)
        invH_k   = matmul(matmul(Vk,invH_k),Vk)
        invH_knn = so2nn_reshape(invH_k,Nnambu*Nspin,Norb)
        FDelta(1,1,:,:,i)=FDelta(1,1,:,:,i) + invH_knn(1,2,:,:)
     enddo
  enddo
end function fdelta_bath_array_replica
