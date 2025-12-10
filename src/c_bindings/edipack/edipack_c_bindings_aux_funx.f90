subroutine ed_set_Hloc_single_N2_c(Hloc,Hloc_anomalous,d) bind(c, name='ed_set_Hloc_single_N2')
  !
  !Sets the local Hamiltonian of the impurity problem.
  !The passed input is a rank-2 array with dimension :code:`[d(1),d(2)]`
  !
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                        :: d(2) !Array dimensions
  complex(c_double_complex),dimension(d(1),d(2)),intent(in) :: Hloc !local Hamiltonian rank-2
  complex(c_double_complex),dimension(d(1),d(2)),intent(in) :: Hloc_anomalous !local anomalous Hamiltonian rank-2
  if(any(abs(Hloc_anomalous)/=0.0))then
    call ed_set_Hloc(Hloc,Hloc_anomalous)
  else
    call ed_set_Hloc(Hloc)
  endif
end subroutine ed_set_Hloc_single_N2_c

subroutine ed_set_Hloc_single_N4_c(Hloc,Hloc_anomalous,d) bind(c, name='ed_set_Hloc_single_N4')
  !
  !Sets the local Hamiltonian of the impurity problem.
  !The passed input is a rank-4 array with dimension :code:`[d(1),...,d(4)]`
  !
  use, intrinsic :: iso_c_binding
  integer(c_int64_t)                                                   :: d(4) !Array dimensions
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4)),intent(in)  :: Hloc !local Hamiltonian rank-4
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4)),intent(in)  :: Hloc_anomalous !local anomalous Hamiltonian rank-4
  if(any(abs(Hloc_anomalous)/=0.0))then
    call ed_set_Hloc(Hloc,Hloc_anomalous)
  else
    call ed_set_Hloc(Hloc)
  endif
end subroutine ed_set_Hloc_single_N4_c


!SEARCH VARIABLE:
subroutine search_variable(var,ntmp,converged) bind(c, name='search_variable')
  use, intrinsic :: iso_c_binding
  real(c_double),dimension(1)         :: var(1)
  real(c_double),dimension(1)         :: ntmp(1)
  integer(c_int),dimension(1)         :: converged(1)
  logical                             :: bool
  converged(1)=0
  call ed_search_variable(var(1),ntmp(1),bool)
  if (bool) converged(1)=1
end subroutine search_variable






