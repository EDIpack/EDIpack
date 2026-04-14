program ed_normal_superc
  USE EDIPACK
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  USE COMMON
  implicit none
  integer                :: i,js,Nso,Nmomenta
  integer                :: Nb,iorb,jorb,ispin,jspin,Ns_full
  integer, parameter     :: Nnambu=2
  complex(8),allocatable :: Smats(:,:,:,:,:,:)
  complex(8),allocatable :: Hloc(:,:,:,:)
  complex(8),allocatable :: denmat4(:,:,:,:)
  complex(8),allocatable :: denmat2(:,:)
  real(8)                :: Delta
  character(len=16)      :: finput
  logical                :: dsave
  !
  !
  ! MPI initialization
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !Parse additional variables && read Input
  call parse_cmd_variable(finput,"FINPUT",default="inputED.in")
  call parse_input_variable(delta,"DELTA",finput,default=0.d0)
  call parse_cmd_variable(dsave,"dsave",default=.false.)
  !
  !
  call ed_read_input(trim(finput))
  !
  if(bath_type/="normal")stop "Wrong setup from input file: non normal bath"
  if(ed_mode/="superc")stop "Wrong setup from input file: non superc ed_mode"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Nso=Nspin*Norb
  Nmomenta=4
  !
  !Allocate Weiss Field:
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Norb))
  allocate(docc(Norb))
  allocate(energy(4))
  allocate(doubles(4))
  allocate(phisc(Norb,Norb))
  allocate(imp(2))
  allocate(Smom(Norb,Nmomenta))
  allocate(ASmom(Norb,Nmomenta))
  !
  allocate(Wlist(Lmats))  
  Wlist = pi/beta*(2*arange(1,Lmats)-1)
  !
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = zero
  if(Norb==2)then
     do js=1,Nspin
        Hloc(js,js,:,:)= Delta*pauli_sigma_z
     end do
  endif
  !
  !
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  !
  !
  call run_test(sparse=.true.,umatrix=.false.,hk=.true.)  
  call run_test(sparse=.false.,umatrix=.false.,hk=.true.)
  call run_test(sparse=.true.,umatrix=.true.,hk=.false.)  
  call run_test(sparse=.false.,umatrix=.true.,hk=.false.)
  call run_test(sparse=.true.,umatrix=.false.,hk=.false.)  
  call run_test(sparse=.false.,umatrix=.false.,hk=.false.)
  !
  call finalize_MPI()



contains


  subroutine run_test(sparse,umatrix,hk)
    logical :: sparse,umatrix,hk
    ED_SPARSE_H    =sparse
    ED_READ_UMATRIX=umatrix
    ED_USE_KANAMORI=hk
    if(.not.umatrix.AND..not.hk)call set_twobody_hk()
    call print_status()
    call ed_init_solver(bath)
    call ed_set_Hloc(hloc)
    call ed_solve(bath)
    call ed_get_sigma(Smats(1,:,:,:,:,:),axis="m",type="n")
    call ed_get_sigma(Smats(2,:,:,:,:,:),axis="m",type="a")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_phi(phisc)
    call ed_get_eimp(energy)
    call ed_get_doubles(doubles)
    call ed_get_imp_info(imp)
    call ed_get_evals(evals)
    ! In SUPERC/NORMAL mode full_denmat covers all sites in Nambu space:
    ! shape [Nnambu,Nnambu,Ns,Ns] with Ns=(Nbath+1)*Norb for normal bath
    Ns_full = (Nbath+1)*Norb
    if(allocated(denmat4))deallocate(denmat4)
    if(allocated(denmat2))deallocate(denmat2)
    allocate(denmat4(Nnambu,Nnambu,Ns_full,Ns_full))
    allocate(denmat2(Nnambu*Ns_full,Nnambu*Ns_full))
    call ed_get_denmat(denmat4)
    call ed_get_denmat(denmat2)
    do i=1,Nmomenta
       do iorb=1,Norb
          call compute_momentum(Wlist,Smats(1,1,1,iorb,iorb,:),i,Smom(iorb,i))
          call compute_momentum(Wlist,Smats(2,1,1,iorb,iorb,:),i,ASmom(iorb,i))
       enddo
    enddo
    call ed_finalize_solver()
    if(dsave)then
       write(*,*)"Saving results to .check files and exit"
       call save_results()
       stop
    endif
    call test_results()
    call test_denmat_checks()
  end subroutine run_test



  subroutine read_results()
    integer :: L
    L = file_length("evals.check")
    if(allocated(evals))deallocate(evals)
    allocate(evals(L))
    call read_array("evals.check",evals)
    call read_array("dens.check",dens)
    call read_array("docc.check",docc)
    call read_array("phisc.check",phisc)
    call read_array("energy.check",energy)
    call read_array("doubles.check",doubles)
    call read_array("imp.check",imp)
    call read_array("Sigma_momenta.check",Smom)
    call read_array("Self_momenta.check",ASmom)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(phiscR))deallocate(phiscR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(impR))deallocate(impR)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(SmomR))deallocate(SmomR)
    if(allocated(ASmomR))deallocate(ASmomR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(phiscR, source=phisc)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(impR, source=imp)
    allocate(evalsR,source=evals)
    allocate(SmomR, source=Smom)
    allocate(ASmomR, source=ASmom)
  end subroutine read_results



  subroutine test_results()
    call read_results()
    write(*,*)
    write(*,"(A50)") "Summary RESULTS:"
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(phisc,phiscR,"phisc")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(imp,impR,"imp")
    call assert(evals,evalsR,"evals")
    call assert(Smom/SmomR,dble(ones(Norb,Nmomenta)),"Sigma_momenta 1:4",tol=1.0d-8)
    call assert(ASmom/ASmomR,dble(ones(Norb,Nmomenta)),"Self_momenta 1:4",tol=1.0d-8)
    call print_status()
    write(*,*)""
    write(*,*)""
    write(*,*)""
    call wait(1000)    
  end subroutine test_results



  subroutine save_results()
    call save_array("evals.check",evals)
    call save_array("dens.check",dens)
    call save_array("docc.check",docc)
    call save_array("phisc.check",phisc)
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("imp.check",imp)
    call save_array("Sigma_momenta.check",Smom)
    call save_array("Self_momenta.check",ASmom)
  end subroutine save_results




  subroutine test_denmat_checks()
    !
    ! Verify the one-body density matrix getters in SUPERC/Nambu mode.
    !
    ! full_denmat has shape [Nnambu=2, Nnambu=2, Ns_, Ns_] where Ns_=(Nbath+1)*Norb.
    ! The Nambu spinor is Psi = (c_up, c†_dn), so:
    !   dm4(1,1,io,io) = <c†_up c_up>  = <n_up(io)>         in [0,1]
    !   dm4(2,2,io,io) = <c_dn  c†_dn> = 1 - <n_dn(io)>     in [0,1]
    !   dm4(1,2,io,jo) = <c†_up c†_dn> (anomalous/pairing)
    !   dm4(2,1,io,jo) = <c_dn  c_up>  (conjugate)
    !
    ! Checks (no reference files needed):
    !   1. Both normal-sector diagonals are real
    !   2. The full Nambu dm4 is Hermitian: dm4(a,b,io,jo) = conj(dm4(b,a,jo,io))
    !   3. n2 variant matches the block-interleaved reshape of n4:
    !         dm2(io+(a-1)*Ns_, jo+(b-1)*Ns_) = dm4(a,b,io,jo)
    !   4. Density check for impurity orbitals (io=1..Norb are the impurity sites):
    !         dens(io) = <n_up> + <n_dn>
    !                  = dm4(1,1,io,io) + (1 - dm4(2,2,io,io))
    !                  = dm4(1,1,io,io) - dm4(2,2,io,io) + 1
    !      With SU2: <n_up>=<n_dn>, so dens = 2*dm4(1,1,io,io) = 2*(1-dm4(2,2,io,io)).
    !
    complex(8),allocatable :: dm4_hc(:,:,:,:)
    complex(8),allocatable :: dm2_from_dm4(:,:)
    integer :: ia,ib,io,jo
    !
    allocate(dm4_hc(Nnambu,Nnambu,Ns_full,Ns_full))
    allocate(dm2_from_dm4(Nnambu*Ns_full,Nnambu*Ns_full))
    !
    write(*,*)
    write(*,"(A50)") "Summary DENMAT CHECKS (SUPERC/Nambu):"
    !
    ! 1. Normal-sector diagonal imaginary parts vanish
    do io=1,Ns_full
       call assert(dimag(denmat4(1,1,io,io)),0.d0, &
            "denmat4: Im[dm(1,1,"//str(io)//","//str(io)//")]")
       call assert(dimag(denmat4(2,2,io,io)),0.d0, &
            "denmat4: Im[dm(2,2,"//str(io)//","//str(io)//")]")
    enddo
    !
    ! 2. Hermiticity: dm4(a,b,io,jo) = conj(dm4(b,a,jo,io))
    !    Holds even for the anomalous blocks because
    !    <c†_up c†_dn>† = <c_dn c_up>, i.e. dm4(1,2)† = dm4(2,1).
    do ia=1,Nnambu
       do ib=1,Nnambu
          do io=1,Ns_full
             do jo=1,Ns_full
                dm4_hc(ia,ib,io,jo) = conjg(denmat4(ib,ia,jo,io))
             enddo
          enddo
       enddo
    enddo
    call assert(denmat4,dm4_hc,"denmat4: Hermitian (Nambu)")
    !
    ! 3. n2 variant is block-interleaved reshape of n4:
    !    dm2(io+(a-1)*Ns_, jo+(b-1)*Ns_) = dm4(a,b,io,jo)
    !    Block layout of dm2:  [ particle-particle | particle-hole ]
    !                          [ hole-particle     | hole-hole     ]
    do ia=1,Nnambu
       do ib=1,Nnambu
          do io=1,Ns_full
             do jo=1,Ns_full
                dm2_from_dm4(io+(ia-1)*Ns_full,jo+(ib-1)*Ns_full) = denmat4(ia,ib,io,jo)
             enddo
          enddo
       enddo
    enddo
    call assert(denmat2,dm2_from_dm4,"denmat: n2 vs n4 consistency (Nambu)")
    !
    ! 4. Density check for impurity orbitals (io=1..Norb):
    !    dens(io) = dm4(1,1,io,io) + 1 - dm4(2,2,io,io)
    do io=1,Norb
       call assert(dreal(denmat4(1,1,io,io))+1.d0-dreal(denmat4(2,2,io,io)),dens(io), &
            "denmat4 vs dens: Re[dm_pp-dm_hh+1]("//str(io)//") == dens")
    enddo
    !
    deallocate(dm4_hc,dm2_from_dm4)
    !
  end subroutine test_denmat_checks


  subroutine set_twobody_hk()
    call ed_add_twobody_operator(1,"u",1,"d",1,"u",1,"d",-2.000000d0)
    call ed_add_twobody_operator(1,"d",1,"u",1,"d",1,"u",-2.000000d0)
    call ed_add_twobody_operator(2,"u",2,"d",2,"u",2,"d",-2.000000d0)
    call ed_add_twobody_operator(2,"d",2,"u",2,"d",2,"u",-2.000000d0)
    call ed_add_twobody_operator(1,"d",2,"u",1,"d",2,"u",-1.500000d0)
    call ed_add_twobody_operator(1,"u",2,"d",1,"u",2,"d",-1.500000d0)
    call ed_add_twobody_operator(2,"d",1,"u",2,"d",1,"u",-1.500000d0)
    call ed_add_twobody_operator(2,"u",1,"d",2,"u",1,"d",-1.500000d0)
    call ed_add_twobody_operator(1,"u",2,"u",1,"u",2,"u",-1.500000d0)
    call ed_add_twobody_operator(1,"d",2,"d",1,"d",2,"d",-1.500000d0)
    call ed_add_twobody_operator(2,"u",1,"u",2,"u",1,"u",-1.500000d0)
    call ed_add_twobody_operator(2,"d",1,"d",2,"d",1,"d",-1.500000d0)
    call ed_add_twobody_operator(1,"u",2,"u",2,"u",1,"u",0.250000d0)
    call ed_add_twobody_operator(1,"d",2,"d",2,"d",1,"d",0.250000d0)
    call ed_add_twobody_operator(2,"u",1,"u",1,"u",2,"u",0.250000d0)
    call ed_add_twobody_operator(2,"d",1,"d",1,"d",2,"d",0.250000d0)
    call ed_add_twobody_operator(1,"d",2,"u",2,"d",1,"u",0.250000d0)
    call ed_add_twobody_operator(1,"u",2,"d",2,"u",1,"d",0.250000d0)
    call ed_add_twobody_operator(2,"d",1,"u",1,"d",2,"u",0.250000d0)
    call ed_add_twobody_operator(2,"u",1,"d",1,"u",2,"d",0.250000d0)
    call ed_add_twobody_operator(1,"d",1,"u",2,"d",2,"u",0.250000d0)
    call ed_add_twobody_operator(1,"u",1,"d",2,"u",2,"d",0.250000d0)
    call ed_add_twobody_operator(2,"d",2,"u",1,"d",1,"u",0.250000d0)
    call ed_add_twobody_operator(2,"u",2,"d",1,"u",1,"d",0.250000d0)
  end subroutine set_twobody_hk




end program



