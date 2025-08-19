program ed_bhz
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                 :: Nso,iloop,Lk
  logical                                 :: converged
  !Bath:
  integer                                 :: Nb
  real(8),allocatable                     :: Bath(:)
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
  !
  complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma5,GammaN
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


  !Allocate Weiss Field:
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




  !> Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])  ! real-space lattice basis vectors
  call TB_set_bk([pi2,0d0],[0d0,pi2])   ! k-space lattice basis vectors
  !> Set the self-energy correction to zero
  call set_SigmaBHZ()
  !> Allocate and build the lattice Hamiltonian using model function above.
  allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
  call TB_build_model(Hk,hk_bhz,Nso,[Nx,Nx])
  !> Get the topological invariant Z_2 = 1/2[C_up - C_dw] using discretized Berry flux in the BZ
  z2 = get_spinChern(Hk,[Nx,Nx])
  if(master)print*,"Z2 = ", z2
  !> Solve the band structure
  call solve_Hk()





  !> Get local Hamiltonian summing over k (one can do better)
  allocate(Hloc(Nso,Nso))
  Hloc = sum(Hk,dim=3)/Lk
  where(abs(dreal(Hloc))<1d-6)Hloc=zero
  !> Set H_{loc} in EDIpack
  call ed_set_hloc(Hloc)
  !> Get bath dimension and allocate user bath to this size
  Nb=ed_get_bath_dimension()
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
     call write_gf(Gmats,"Gloc",axis='mats',iprint=1)


     !> Update the Weiss field: (using DMFT_TOOLS). Linear mixing.
     call dmft_self_consistency(Gmats,Smats,Weiss)
     if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_


     !> Fit the new bath, starting from the current bath + the update Weiss field
     call ed_chi2_fitgf(Weiss,bath,ispin=1)
     call ed_spin_symmetrize_bath(bath,save=.true.)

     !Check convergence (using DMFT_TOOLS)
     converged = check_convergence(Weiss(1,1,:),dmft_error,nsuccess,nloop)
     Weiss_=Weiss

     call end_loop
  enddo


  !> retrieve real-axis self-energy and build local Green's function
  call ed_get_sigma(Sreal,axis='real')
  call get_gloc(Hk,Greal,Sreal,axis='r')
  call write_gf(Greal,"Gloc",axis='real',iprint=1)
  !
  !> Set the self-energy correction, build the topological Hamiltonian and get corresponding invariant
  call set_SigmaBHZ(Smats(:,:,1))
  call TB_build_model(Hk,hk_bhz,Nso,[Nx,Nx])  !this is now Htop = Z.(Hk + ReSigma)
  z2 = get_spinChern(Hk,[Nx,Nx])
  if(master)print*,"Z2* = ", z2
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



  !get spin Chern number of H(k) using discretized Berry flux in the BZ
  function get_spinChern(Hk,Nkvec,spin,berry) result(sp_chern)
    !get spin Chern number of H(k) using discretized Berry flux in the BZ
    complex(8),intent(in),dimension(:,:,:)      :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)             :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    integer,optional                            :: spin
    real(8),dimension(:,:),allocatable,optional :: berry
    real(8)                                     :: sp_chern
    !
    integer                                     :: spin_
    integer                                     :: Nocc
    integer                                     :: Nktot
    integer                                     :: Nkx,Nky
    integer                                     :: ispin
    integer                                     :: ikx,iky
    integer                                     :: ikxP,ikyP
    integer                                     :: ik,iocc
    complex(8),dimension(:,:),allocatable       :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable            :: Eigval ![Nlso]
    complex(8),dimension(:,:,:,:),allocatable   :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                     :: Ulink
    real(8)                                     :: berry_phase,chern(2)
    integer                                     :: unit
    !
    spin_ = 0;if(present(spin))spin_=spin
    Nocc  = Nso/2
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    Nktot = product(Nkvec)
    call assert_shape(Hk,[Nso,Nso,Nktot],"Get_Chern_NUmber","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nso,Nso))
    allocate(Eigval(Nso))
    allocate(BlochStates(Nkx,Nky,Nocc,Nso)) 
    if(present(berry))then
       if(allocated(berry))deallocate(berry)
       allocate(berry(Nkx,Nky))
    endif
    !
    ik=0
    do ikx=1,Nkx
       do iky=1,Nky
          ik=ik+1
          Eigvec = Hk(:,:,ik)
          call eigh(Eigvec,Eigval)
          do iocc=1,Nocc
             BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
          enddo
       enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
       ikxP = modulo(ikx,Nkx) + 1
       do iky= 1, Nky
          ikyP = modulo(iky,Nky) + 1
          !
          select case(spin_)
          case default;stop "hk_to_spin_Chern error: spin != {0,1,2}: 0=Z2/Both, 1=UP, 2=DW"
          case(0)
             do ispin=1,2
                Ulink(1)    = dot_product(BlochStates(ikx,iky,ispin,:)  , BlochStates(ikx,ikyP,ispin,:))
                Ulink(2)    = dot_product(BlochStates(ikx,ikyP,ispin,:) , BlochStates(ikxP,ikyP,ispin,:))
                Ulink(3)    = dot_product(BlochStates(ikxP,ikyP,ispin,:), BlochStates(ikxP,iky,ispin,:))
                Ulink(4)    = dot_product(BlochStates(ikxP,iky,ispin,:) , BlochStates(ikx,iky,ispin,:))
                berry_phase = -dimag(zlog( product(Ulink(:))  ))/pi2
                chern(ispin)= chern(ispin) + berry_phase
             enddo
             sp_chern = 0.5d0*(chern(1)-chern(2))
             if(master)then
                open(unit=free_unit(unit),file="Hk_to_z2.dat")
                write(unit,*)sp_chern
                close(unit)
             endif
          case(1,2)
             ispin=spin_
             Ulink(1)    = dot_product(BlochStates(ikx,iky,ispin,:)  , BlochStates(ikx,ikyP,ispin,:))
             Ulink(2)    = dot_product(BlochStates(ikx,ikyP,ispin,:) , BlochStates(ikxP,ikyP,ispin,:))
             Ulink(3)    = dot_product(BlochStates(ikxP,ikyP,ispin,:), BlochStates(ikxP,iky,ispin,:))
             Ulink(4)    = dot_product(BlochStates(ikxP,iky,ispin,:) , BlochStates(ikx,iky,ispin,:))
             berry_phase = dimag(zlog( product(Ulink(:))  ))/pi2
             chern(ispin)= chern(ispin) + berry_phase
             if(present(berry))&
                  berry(ikx,iky) = berry_phase!*one_over_area=Nkx/pi2*Nkx/pi2
             sp_chern = chern(ispin)
             if(master)then
                open(unit=free_unit(unit),file="Hk_to_spin_Chern_s"//str(spin_)//".dat")
                write(unit,*)sp_chern
                close(unit)
             endif
          end select
       enddo
    enddo
    !
    !
  end function get_spinChern

end program ed_bhz
