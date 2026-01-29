program ed_bhz
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  
  integer                                 :: iloop,Lk,Nso
  logical                                 :: converged
  !Bath:
  integer                                 :: Nb
  real(8),allocatable                     :: Bath(:)
  !The local hybridization function:
  !The local hybridization function:
  complex(8),dimension(:,:,:),allocatable :: Weiss, Weiss_
  complex(8),dimension(:,:,:),allocatable :: Smats,Sreal
  complex(8),dimension(:,:,:),allocatable :: Gmats,Greal
  !hamiltonian input:
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:),allocatable   :: Hloc
  complex(8),dimension(:,:),allocatable   :: sigmaBHZ
  real(8),dimension(:,:),allocatable      :: Zbhz
  !variables for the model:
  integer                                 :: Nx,Nkpath
  real(8)                                 :: mh,lambda,wmixing,z2
  character(len=16)                       :: finput
  !>Gamma matrices:
  complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma5,GammaN
  complex(8),dimension(4,4)               :: GammaE0,GammaEx,GammaEy,GammaEz
  !>Replica bath
  real(8),dimension(:,:),allocatable      :: lambdasym_vector
  complex(8),dimension(:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                 :: comm,rank
  logical                                 :: master,fhtop

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.conf')  
  call parse_input_variable(Nx,"NX",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  !
  call ed_read_input(trim(finput))
  !
  !
  !Add DMFT CTRL Variables used in DMFT_TOOLS:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
 call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb
  Lk=Nx**2

  !Allocate local functions:
  allocate(Weiss(Nso,Nso,Lmats))
  allocate(Weiss_(Nso,Nso,Lmats))
  allocate(Smats(Nso,Nso,Lmats))
  allocate(Gmats(Nso,Nso,Lmats))
  allocate(Sreal(Nso,Nso,Lreal))
  allocate(Greal(Nso,Nso,Lreal))
  allocate(SigmaBHZ(Nso,Nso))
  allocate(Zbhz(Nso,Nso))

  gamma1=kron( pauli_sigma_z, pauli_tau_x)
  gamma2=kron( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron( pauli_sigma_0, pauli_tau_z)
  gammaN=kron( pauli_sigma_0, pauli_tau_0)
  !
  gammaE0=kron( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron( pauli_sigma_z, pauli_tau_x )



  !> Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])  ! real-space lattice basis vectors
  call TB_set_bk([pi2,0d0],[0d0,pi2])   ! k-space lattice basis vectors
  !> Set the self-energy correction to zero
  call set_SigmaBHZ()
  !> Allocate and build the lattice Hamiltonian using model function above.
  allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
  call TB_build_model(Hk,hk_bhz,Nso,[Nx,Nx])
  !> Solve the band structure
  call solve_Hk()


  !> Get local Hamiltonian summing over k (one can do better)
  allocate(Hloc(Nso,Nso))
  Hloc = sum(Hk,dim=3)/Lk
  where(abs(dreal(Hloc))<1d-6)Hloc=zero
  !> Set H_{loc} in EDIpack
  call ed_set_hloc(Hloc)


  !> Setup the replica bath basis for the case E0EzEx(singlet,tripletZ,tripletX)
  allocate(lambdasym_vector(Nbath,4))
  allocate(Hsym_basis(Nso,Nso,4))
  Hsym_basis(:,:,1)=Gamma5  ;lambdasym_vector(:,1)= Mh
  Hsym_basis(:,:,2)=GammaE0 ;lambdasym_vector(:,2)= sb_field
  Hsym_basis(:,:,3)=GammaEz ;lambdasym_vector(:,3)= sb_field
  Hsym_basis(:,:,4)=GammaEx ;lambdasym_vector(:,4)=-sb_field
  !> Set the replica bath
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  !> Get bath dimension and allocate user bath to this size
  Nb=ed_get_bath_dimension(4)!(Hsym_basis)
  allocate(Bath(Nb))
  !
  !> Initialize the ED solver (bath is guessed or read from file) 
  call ed_init_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !> Solve the impurity problem, retrieve matsubara self-energy 
     call ed_solve(bath)
     call ed_get_sigma(Smats,axis='mats')

     !> Get Gloc (using DMFT_TOOLS)
     call get_gloc(Hk,Gmats,Smats,axis='m')
     call write_gf(Gmats,"Gloc",axis='mats',iprint=3)

     !> Update the Weiss field: (using DMFT_TOOLS). Linear mixing.
     call dmft_self_consistency(Gmats,Smats,Weiss)
     if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_

     !> Fit the new bath, starting from the current bath + the update Weiss field
     call ed_chi2_fitgf(Weiss,bath)

     !Check convergence (using DMFT_TOOLS)
     converged = check_convergence(Weiss(1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  !> retrieve real-axis self-energy and build local Green's function
  call ed_get_sigma(Sreal,axis='real')
  call get_gloc(Hk,Greal,Sreal,axis='r')
  call write_gf(Greal,"Gloc",axis='real',iprint=3)
  !
  !> Set the self-energy correction, build the topological Hamiltonian and get corresponding invariant
  call set_SigmaBHZ(Smats(:,:,1))
  call TB_build_model(Hk,hk_bhz,Nso,[Nx,Nx])  !this is now Htop = Z.(Hk + ReSigma)
  !> Solve the topological Hamiltonian Band structure (a proxy for GF poles).
  call solve_hk()


  call finalize_MPI()


  

contains






  function hk_bhz(kvec,N) result(hk)     
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky
    !
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx = kvec(1) ; ky=kvec(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
    !
    !> Include the self-energy correction, if previously defined
    if(fhtop)then
       Hk = Hk + dreal(SigmaBHZ)
       Hk = matmul(Zbhz,Hk)
    endif
  end function hk_bhz


  subroutine set_SigmaBHZ(sigma)    
    complex(8),dimension(Nso,Nso),optional :: sigma
    integer                                :: ii
    !
    sigmaBHZ = zero ; Zbhz=eye(Nso); fhtop=.false.
    !
    if(present(sigma))then
       sigmaBHZ=sigma
       !
       Zbhz=zero
       do ii=1,Nso
          Zbhz(ii,ii)  = 1.d0/abs( 1.d0 +  abs(dimag(sigmaBHZ(ii,ii))/pi*beta)) 
       end do
       !
       fhtop=.true.
       !
    endif
  end subroutine set_SigmaBHZ




  !--------------------------------------------------------------------!
  !PURPOSE: Solve the topological Hamiltonian
  !--------------------------------------------------------------------!
  subroutine solve_hk()
    integer                            :: Npts
    real(8),dimension(:,:),allocatable :: kpath
    !
    if(master)then
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_gamma
       kpath(2,:)=kpoint_x1
       kpath(3,:)=kpoint_m1
       kpath(4,:)=kpoint_gamma
       write(LOGfile,*)"Build H(k) BHZ along path \Gamma-X-M-\Gamma:"
       call TB_solve_model(hk_bhz,Nso,kpath,Nkpath,&
            colors_name=[red,blue,red,blue],&
            points_name=[character(len=20) :: "{\Symbol G}","X","M","{\Symbol G}"],&
            file="Eig_Hk_"//str(fhtop)//".ed")
    endif
  end subroutine solve_hk


end program ed_bhz
