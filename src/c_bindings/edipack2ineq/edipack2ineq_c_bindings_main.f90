!ED_MAIN:
subroutine init_solver_ineq_c(bath,dim_bath) bind(c, name='init_solver_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t),dimension(2),intent(in)                           :: dim_bath
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(inout)      :: bath
  call ed_init_solver(bath)
end subroutine init_solver_ineq_c

subroutine init_solver_ineq_nobath_c(Nlat) bind(c, name='init_solver_ineq_nobath')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                      :: Nlat
  call ed_init_solver(Nlat)
end subroutine init_solver_ineq_nobath_c

subroutine solve_ineq_c(bath,dim_bath,flag_gf,mpi_lanc) bind(c, name='solve_ineq')
  use, intrinsic :: iso_c_binding
  integer(c_int64_t),dimension(2),intent(in)                                                           :: dim_bath
  integer(c_int),value                                                                                 :: mpi_lanc,flag_gf
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(in)                                         :: bath
  integer                                                                                              :: Nineq
  Nineq = size(bath,1)
  call ed_solve(bath,mpi_lanc=i2l(mpi_lanc),flag_gf=i2l(flag_gf))
end subroutine solve_ineq_c

subroutine solve_ineq_nobath_c(Nlat,flag_gf,mpi_lanc) bind(c, name='solve_ineq_nobath')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                                                 :: Nlat,mpi_lanc,flag_gf
  call ed_solve(Nlat,mpi_lanc=i2l(mpi_lanc),flag_gf=i2l(flag_gf))
end subroutine solve_ineq_nobath_c

