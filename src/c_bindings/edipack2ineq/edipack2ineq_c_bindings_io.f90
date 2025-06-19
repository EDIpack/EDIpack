!ED_IO:

!OBSERVABLES

!density
subroutine ed_get_dens_n1_c(self) bind(c,name="ed_get_dens_n1")
  use, intrinsic :: iso_c_binding
  real(c_double)     :: self(Norb)
  call ed_get_dens(self)
end subroutine ed_get_dens_n1_c

subroutine ed_get_dens_n2_c(self,Nlat) bind(c,name="ed_get_dens_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(Nlat,Norb)
  integer(c_int),value     :: Nlat
  call ed_get_dens(self,Nlat)
end subroutine ed_get_dens_n2_c

!magnetization
subroutine ed_get_mag_n2_c(self) bind(c,name="ed_get_mag_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(3,Norb)
  integer(c_int)           :: iorb
  do iorb = 1,Norb
     call ed_get_mag(self(1,iorb),component="x",iorb=iorb)
     call ed_get_mag(self(2,iorb),component="y",iorb=iorb)
     call ed_get_mag(self(3,iorb),component="z",iorb=iorb)
  enddo
end subroutine ed_get_mag_n2_c

subroutine ed_get_mag_n3_c(self,Nlat) bind(c,name="ed_get_mag_n3")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(Nlat,3,Norb)
  integer(c_int),value     :: Nlat
  call ed_get_mag(self(:,1,:),Nlat,"x")
  call ed_get_mag(self(:,2,:),Nlat,"y")
  call ed_get_mag(self(:,3,:),Nlat,"z")
end subroutine ed_get_mag_n3_c

!double occupation
subroutine ed_get_docc_n1_c(self) bind(c,name="ed_get_docc_n1")
  use, intrinsic :: iso_c_binding
  real(c_double)     :: self(Norb)
  call ed_get_docc(self)
end subroutine ed_get_docc_n1_c

subroutine ed_get_docc_n2_c(self,Nlat) bind(c,name="ed_get_docc_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)           :: self(Nlat,Norb)
  integer(c_int),value     :: Nlat
  call ed_get_docc(self,Nlat)
end subroutine ed_get_docc_n2_c

!superconductive phi
subroutine ed_get_phisc_n2_c(self) bind(c,name="ed_get_phisc_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)     :: self(Norb,Norb)
  call ed_get_phi(self)
end subroutine ed_get_phisc_n2_c

subroutine ed_get_phisc_n3_c(self,Nlat) bind(c,name="ed_get_phisc_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value     :: Nlat
  real(c_double)           :: self(Nlat,Norb,Norb)
  call ed_get_phi(self,Nlat)
end subroutine ed_get_phisc_n3_c

!superconductive arg
subroutine ed_get_argsc_n2_c(self) bind(c,name="ed_get_argsc_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)     :: self(Norb,Norb)
  call ed_get_arg(self)
end subroutine ed_get_argsc_n2_c

subroutine ed_get_argsc_n3_c(self,Nlat) bind(c,name="ed_get_argsc_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value     :: Nlat
  real(c_double)           :: self(Nlat,Norb,Norb)
  call ed_get_arg(self,Nlat)
end subroutine ed_get_argsc_n3_c

!energy
subroutine ed_get_eimp_n1_c(self) bind(c,name="ed_get_eimp_n1")
  use, intrinsic :: iso_c_binding
  real(c_double) :: self(4)
  call ed_get_eimp(self)
end subroutine ed_get_eimp_n1_c

subroutine ed_get_eimp_n2_c(self,Nlat) bind(c,name="ed_get_eimp_n2")
  use, intrinsic :: iso_c_binding
  real(c_double)                   :: self(Nlat,4)
  integer(c_int),value             :: Nlat
  call ed_get_eimp(self,Nlat)
end subroutine ed_get_eimp_n2_c

!---------!
!GET SIGMA!
!---------!

!SITE
subroutine get_sigma_site_n3_c(self,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_site_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_sigma(self,axis_,typ_,zeta)
  else
     call ed_get_sigma(self,axis_,typ_)
  endif
end subroutine get_sigma_site_n3_c

subroutine get_sigma_site_n5_c(self,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_site_n5")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_sigma(self,axis_,typ_,zeta)
  else
     call ed_get_sigma(self,axis_,typ_)
  endif
end subroutine get_sigma_site_n5_c

!LATTICE
subroutine get_sigma_lattice_n3_c(self,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_lattice_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_sigma(self,Nineq,axis_,typ_,zeta)
  else
     call ed_get_sigma(self,Nineq,axis_,typ_)
  endif
end subroutine get_sigma_lattice_n3_c

subroutine get_sigma_lattice_n4_c(self,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_lattice_n4")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nineq,Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_sigma(self,Nineq,axis_,typ_,zeta)
  else
     call ed_get_sigma(self,Nineq,axis_,typ_)
  endif
end subroutine get_sigma_lattice_n4_c

subroutine get_sigma_lattice_n6_c(self,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_sigma_lattice_n6")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: self(Nineq,Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_sigma(self,Nineq,axis_,typ_,zeta)
  else
     call ed_get_sigma(self,Nineq,axis_,typ_)
  endif
end subroutine get_sigma_lattice_n6_c

!---------!
!GET GIMP !
!---------!

!SITE
subroutine get_gimp_site_n3_c(gimp,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_site_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_gimp(gimp,axis_,typ_,zeta)
  else
     call ed_get_gimp(gimp,axis_,typ_)
  endif
end subroutine get_gimp_site_n3_c

subroutine get_gimp_site_n5_c(gimp,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_site_n5")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_gimp(gimp,axis_,typ_,zeta)
  else
     call ed_get_gimp(gimp,axis_,typ_)
  endif
end subroutine get_gimp_site_n5_c

!LATTICE
subroutine get_gimp_lattice_n3_c(gimp,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_lattice_n3")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nineq*Nspin*Norb,Nineq*Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_gimp(gimp,Nineq,axis_,typ_,zeta)
  else
     call ed_get_gimp(gimp,Nineq,axis_,typ_)
  endif
end subroutine get_gimp_lattice_n3_c

subroutine get_gimp_lattice_n4_c(gimp,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_lattice_n4")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nineq,Nspin*Norb,Nspin*Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_gimp(gimp,Nineq,axis_,typ_,zeta)
  else
     call ed_get_gimp(gimp,Nineq,axis_,typ_)
  endif
end subroutine get_gimp_lattice_n4_c

subroutine get_gimp_lattice_n6_c(gimp,Nineq,axis,typ,zeta,dz,zflag) bind(c,name="get_gimp_lattice_n6")
  use, intrinsic :: iso_c_binding
  integer(c_int),value                          :: dz,axis,Nineq,typ,zflag
  character(len=1)                              :: axis_
  character(len=1)                              :: typ_
  complex(c_double_complex)                     :: zeta(dz)
  complex(c_double_complex)                     :: gimp(Nineq,Nspin,Nspin,Norb,Norb,dz)
  !
  axis_="m"
  if(axis==1)axis_="r"
  typ_="n"
  if(typ==1)typ_="a"
  if(zflag==1)then
     call ed_get_gimp(gimp,Nineq,axis_,typ_,zeta)
  else
     call ed_get_gimp(gimp,Nineq,axis_,typ_)
  endif
end subroutine get_gimp_lattice_n6_c

!Get G0and given a bath array
subroutine ed_get_g0and_d3_c(warray,dim_warray,bath,dim_bath,G0and,dim_g0and,axis,typ) bind(c, name='get_g0and_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                        :: dim_warray, dim_bath
  integer(c_int64_t)                                                          :: dim_G0and(3)
  complex(c_double_complex),dimension(dim_warray)                             :: warray
  real(c_double),dimension(dim_bath)                                          :: bath
  complex(c_double_complex),dimension(dim_g0and(1),dim_g0and(2),dim_g0and(3)) :: g0and
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)

  call ed_get_g0and(warray,bath,G0and,axis_,typ_)

end subroutine ed_get_g0and_d3_c

subroutine ed_get_g0and_d5_c(warray,dim_warray,bath,dim_bath,G0and,dim_g0and,axis,typ) bind(c, name='get_g0and_n5')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                        :: dim_warray, dim_bath
  integer(c_int64_t)                                                          :: dim_G0and(5)
  complex(c_double_complex),dimension(dim_warray)                             :: warray
  real(c_double),dimension(dim_bath)                                          :: bath
  complex(c_double_complex),dimension(dim_g0and(1),dim_g0and(2),dim_g0and(3),dim_g0and(4),dim_g0and(5)) :: g0and
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)

  call ed_get_g0and(warray,bath,G0and,axis_,typ_)

end subroutine ed_get_g0and_d5_c

!Get Delta given a bath array
subroutine ed_get_delta_d3_c(warray,dim_warray,bath,dim_bath,Delta,dim_Delta,axis,typ) bind(c, name='get_delta_n3')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                        :: dim_warray, dim_bath
  integer(c_int64_t)                                                          :: dim_Delta(3)
  complex(c_double_complex),dimension(dim_warray)                             :: warray
  real(c_double),dimension(dim_bath)                                          :: bath
  complex(c_double_complex),dimension(dim_Delta(1),dim_Delta(2),dim_Delta(3)) :: Delta
  character(kind=c_char), dimension(1),optional                               :: axis,typ
  character(len=1)                                                            :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)

  call ed_get_Delta(warray,bath,Delta,axis_,typ_)

end subroutine ed_get_delta_d3_c

subroutine ed_get_delta_d5_c(warray,dim_warray,bath,dim_bath,Delta,dim_Delta,axis,typ) bind(c, name='get_delta_n5')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                                                  :: dim_warray, dim_bath
  integer(c_int64_t)                                                                                    :: dim_Delta(5)
  complex(c_double_complex),dimension(dim_warray)                                                       :: warray
  real(c_double),dimension(dim_bath)                                                                    :: bath
  complex(c_double_complex),dimension(dim_Delta(1),dim_Delta(2),dim_Delta(3),dim_Delta(4),dim_Delta(5)) :: Delta
  character(kind=c_char), dimension(1),optional                                                         :: axis,typ
  character(len=1)                                                                                      :: axis_,typ_
  typ_(1:1)=typ(1)
  axis_(1:1)=axis(1)

  call ed_get_Delta(warray,bath,Delta,axis_,typ_)

end subroutine ed_get_delta_d5_c


!----------------!
!SUSCEPTIBILITIES!
!----------------!

subroutine get_spinChi_c(self,zeta,dim_zeta,zetaflag,axis,Nsites,latticeflag) bind(c, name='ed_get_spinchi')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                :: dim_zeta
  complex(c_double_complex),dimension(dim_zeta)                       :: zeta
  integer(c_int),value                                                :: axis  ! 0="m", 1="r", 2="t"
  integer(c_int),value                                                :: latticeflag, Nsites, zetaflag
  character(len=1)                                                    :: axis_
  complex(c_double_complex),dimension(Nsites,Norb,Norb,dim_zeta)      :: self
  !
  if(axis==0)axis_="m"
  if(axis==1)axis_="r"
  if(axis==2)axis_="t"
  !
  if(zetaflag==0)then
     if (latticeflag==0)then
        call ed_get_spinChi(self(1,:,:,:),axis_)
     else
        call ed_get_spinChi(self,Nsites,axis_)
     endif
  else
     if (latticeflag==0)then
        call ed_get_spinChi(self(1,:,:,:),axis_,zeta)
     else
        call ed_get_spinChi(self,Nsites,axis_,zeta)
     endif
  endif
end subroutine get_spinChi_c

subroutine get_densChi_c(self,zeta,dim_zeta,zetaflag,axis,Nsites,latticeflag) bind(c,name='ed_get_denschi')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                       :: dim_zeta
  complex(c_double_complex),dimension(dim_zeta)                              :: zeta
  integer(c_int),value                                                       :: axis  ! 0="m", 1="r", 2="t"
  integer(c_int),value                                                       :: latticeflag, Nsites, zetaflag
  character(len=1)                                                           :: axis_
  complex(c_double_complex),dimension(Nsites,Norb,Norb,dim_zeta)             :: self
  !
  if(axis==0)axis_="m"
  if(axis==1)axis_="r"
  if(axis==2)axis_="t"
  !
  if(zetaflag==0)then
     if (latticeflag==0)then
        call ed_get_densChi(self(1,:,:,:),axis_)
     else
        call ed_get_densChi(self,Nsites,axis_)
     endif
  else
     if (latticeflag==0)then
        call ed_get_densChi(self(1,:,:,:),axis_,zeta)
     else
        call ed_get_densChi(self,Nsites,axis_,zeta)
     endif
  endif
  !
end subroutine get_densChi_c

subroutine get_pairChi_c(self,zeta,dim_zeta,zetaflag,axis,Nsites,latticeflag) bind(c,name='ed_get_pairchi')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                :: dim_zeta
  complex(c_double_complex),dimension(dim_zeta)                       :: zeta
  integer(c_int),value                                                :: axis  ! 0="m", 1="r", 2="t"
  integer(c_int),value                                                :: latticeflag, Nsites, zetaflag
  character(len=1)                                                    :: axis_
  complex(c_double_complex),dimension(Nsites,Norb,Norb,dim_zeta)      :: self
  !
  if(axis==0)axis_="m"
  if(axis==1)axis_="r"
  if(axis==2)axis_="t"
  !
  if(zetaflag==0)then
     if (latticeflag==0)then
        call ed_get_pairChi(self(1,:,:,:),axis_)
     else
        call ed_get_pairChi(self,Nsites,axis_)
     endif
  else
     if (latticeflag==0)then
        call ed_get_pairChi(self(1,:,:,:),axis_,zeta)
     else
        call ed_get_pairChi(self,Nsites,axis_,zeta)
     endif
  endif
  !
end subroutine get_pairChi_c

subroutine get_exctChi_c(self,zeta,dim_zeta,zetaflag,axis,Nsites,latticeflag) bind(c,name='ed_get_exctchi')
  use, intrinsic :: iso_c_binding
  integer(c_int),value                                                :: dim_zeta
  complex(c_double_complex),dimension(dim_zeta)                       :: zeta
  integer(c_int),value                                                :: axis  ! 0="m", 1="r", 2="t"
  integer(c_int),value                                                :: latticeflag, Nsites, zetaflag
  character(len=1)                                                    :: axis_
  complex(c_double_complex),dimension(Nsites,3,Norb,Norb,dim_zeta)    :: self
  !
  if(axis==0)axis_="m"
  if(axis==1)axis_="r"
  if(axis==2)axis_="t"
  !
  if(zetaflag==0)then
     if (latticeflag==0)then
        call ed_get_exctChi(self(1,:,:,:,:),axis_)
     else
        call ed_get_exctChi(self,Nsites,axis_)
     endif
  else
     if (latticeflag==0)then
        call ed_get_exctChi(self(1,:,:,:,:),axis_,zeta)
     else
        call ed_get_exctChi(self,Nsites,axis_,zeta)
     endif
  endif
  !
end subroutine get_exctChi_c

!----------------!
!       RDM      !
!----------------!

subroutine ed_get_impurity_rdm_c(rdm,doprint) bind(c,name='ed_get_impurity_rdm')
  use, intrinsic :: iso_c_binding
  complex(8),dimension(4**Norb,4**norb),intent(inout) :: rdm
  integer(c_int),value                                :: doprint
  !
  if (doprint==0)then
     call ed_get_impurity_rdm(rdm)
  else
     call ed_get_impurity_rdm(rdm,.true.)
  endif
  !
end subroutine ed_get_impurity_rdm_c
