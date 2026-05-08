program ed_normal_normal
  USE EDIPACK
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  USE COMMON
  implicit none
  integer                :: i,js,Nso,Nmomenta
  integer                :: Nb,iorb,jorb,ispin,jspin,Ns_full
  complex(8),allocatable :: Smats(:,:,:,:,:)
  complex(8),allocatable :: Hloc(:,:,:,:)
  complex(8),allocatable :: denmat4(:,:,:,:)
  complex(8),allocatable :: denmat2(:,:)

  !variables for the model:
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
  if(ed_mode/="normal")stop "Wrong setup from input file: non normal ed_mode"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Nso=Nspin*Norb
  Nmomenta=4
  !
  !Allocate Weiss Field:
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Norb))
  allocate(docc(Norb))
  allocate(energy(4))
  allocate(doubles(4))
  allocate(imp(2))
  allocate(Smom(Norb,Nmomenta))
  ! allocate(rdm(4**Norb,4**Norb))
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
    call ed_get_sigma(Smats,axis="m",type="n")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_eimp(energy)
    call ed_get_doubles(doubles)
    call ed_get_imp_info(imp)
    call ed_get_evals(evals)
    ! In NORMAL mode full_denmat covers all (impurity+bath) sites:
    ! shape [Nspin,Nspin,Ns,Ns] with Ns=(Nbath+1)*Norb for normal bath
    Ns_full = (Nbath+1)*Norb
    if(allocated(denmat4))deallocate(denmat4)
    if(allocated(denmat2))deallocate(denmat2)
    allocate(denmat4(Nspin,Nspin,Ns_full,Ns_full))
    allocate(denmat2(Nspin*Ns_full,Nspin*Ns_full))
    call ed_get_denmat(denmat4)
    call ed_get_denmat(denmat2)
    ! call ed_get_impurity_rdm(rdm)
    do i=1,Nmomenta
       do iorb=1,Norb
          call compute_momentum(Wlist,Smats(1,1,iorb,iorb,:),i,Smom(iorb,i))
       enddo
    enddo
    call ed_finalize_solver()
    if(dsave)then
       write(*,*)"Saving results to .check files and exit"
       call save_results()
       stop "Results saved to *.check files"
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
    call read_array("energy.check",energy)
    call read_array("doubles.check",doubles)
    call read_array("imp.check",imp)
    call read_array("Sigma_momenta.check",Smom)
    ! call read_array("rdm.check",rdm)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(impR))deallocate(impR)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(SmomR))deallocate(SmomR)
    ! if(allocated(rdmR))deallocate(rdmR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(impR, source=imp)
    allocate(evalsR,source=evals)
    allocate(SmomR, source=Smom)
    ! allocate(rdmR, source=rdm)
  end subroutine read_results


  subroutine test_results()
    call read_results()
    write(*,*)
    write(*,"(A50)") "Summary RESULTS:"
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(imp,impR,"imp")
    call assert(evals,evalsR,"evals")
    call assert(Smom/SmomR,dble(ones(Norb,Nmomenta)),"Sigma_momenta 1:4",tol=1.0d-8)
    ! call assert(rdm-rdmR,zeros(4**Norb,4**Norb),"RDM",tol=1.0d-8)
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
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("imp.check",imp)
    call save_array("Sigma_momenta.check",Smom)
    ! call save_array("rdm.check",rdm)
  end subroutine save_results


  subroutine test_denmat_checks()
    !
    ! Verify the one-body density matrix getters in NORMAL mode (Nspin=1, SU2).
    ! full_denmat covers impurity orbitals only: shape [Nspin,Nspin,Norb,Norb].
    ! Four in-memory checks (no reference files needed):
    !   1. Diagonal elements are purely real
    !   2. Density check: dens(io) = <nup>+<ndw> = 2*dm4(1,1,io,io) (SU2 symmetry)
    !   3. dm4 is Hermitian: dm4(is,js,io,jo) == conj(dm4(js,is,jo,io))
    !   4. n2 getter matches block-interleaved reshape of n4:
    !         dm2(io+(is-1)*Norb, jo+(js-1)*Norb) == dm4(is,js,io,jo)
    !
    complex(8),dimension(Nspin,Nspin,Ns_full,Ns_full) :: dm4_hc
    complex(8),dimension(Nspin*Ns_full,Nspin*Ns_full) :: dm2_from_dm4
    integer :: is,js,io,jo
    !
    write(*,*)
    write(*,"(A50)") "Summary DENMAT CHECKS:"
    !
    ! 1. Diagonal imaginary parts vanish
    do io=1,Ns_full
       call assert(dimag(denmat4(1,1,io,io)),0.d0, &
            "denmat4: Im[dm(1,1,"//str(io)//","//str(io)//")]")
    enddo
    !
    ! 2. Density check: dens(io) = <nup>+<ndw> = 2*dm4(1,1,io,io) (by SU2 symmetry).
    !    Loose tolerance: dm4 and dens are accumulated via independent sums so
    !    floating-point rounding can cause O(1e-9) disagreement.
    do io=1,Norb
       call assert(2.d0*dreal(denmat4(1,1,io,io)),dens(io), &
            "denmat4 vs dens: 2*Re[dm(1,1,"//str(io)//","//str(io)//")]==dens",tol=1.d-8)
    enddo
    !
    ! 3. Hermiticity of dm4
    do is=1,Nspin
       do js=1,Nspin
          do io=1,Ns_full
             do jo=1,Ns_full
                dm4_hc(is,js,io,jo) = conjg(denmat4(js,is,jo,io))
             enddo
          enddo
       enddo
    enddo
    call assert(denmat4,dm4_hc,"denmat4: Hermitian")
    !
    ! 4. n2 variant is the block-interleaved reshape of n4:
    !    dm2(io+(is-1)*Norb, jo+(js-1)*Norb) = dm4(is,js,io,jo)
    do is=1,Nspin
       do js=1,Nspin
          do io=1,Ns_full
             do jo=1,Ns_full
                dm2_from_dm4(io+(is-1)*Norb,jo+(js-1)*Norb) = denmat4(is,js,io,jo)
             enddo
          enddo
       enddo
    enddo
    call assert(denmat2,dm2_from_dm4,"denmat: n2 vs n4 consistency")
    !
  end subroutine test_denmat_checks


  subroutine set_twobody_hk()
    call ed_add_twobody_operator(1,"u",1,"d",1,"u",1,"d",2.00000000d0)
    call ed_add_twobody_operator(1,"d",1,"u",1,"d",1,"u",2.00000000d0)
    call ed_add_twobody_operator(2,"u",2,"d",2,"u",2,"d",2.00000000d0)
    call ed_add_twobody_operator(2,"d",2,"u",2,"d",2,"u",2.00000000d0)
    call ed_add_twobody_operator(1,"d",2,"u",1,"d",2,"u",2.00000000d0)
    call ed_add_twobody_operator(1,"u",2,"d",1,"u",2,"d",2.00000000d0)
    call ed_add_twobody_operator(2,"d",1,"u",2,"d",1,"u",2.00000000d0)
    call ed_add_twobody_operator(2,"u",1,"d",2,"u",1,"d",2.00000000d0)
    call ed_add_twobody_operator(1,"u",2,"u",1,"u",2,"u",2.00000000d0)
    call ed_add_twobody_operator(1,"d",2,"d",1,"d",2,"d",2.00000000d0)
    call ed_add_twobody_operator(2,"u",1,"u",2,"u",1,"u",2.00000000d0)
    call ed_add_twobody_operator(2,"d",1,"d",2,"d",1,"d",2.00000000d0)
    call ed_add_twobody_operator(1,"u",2,"u",2,"u",1,"u",0.12500000d0)
    call ed_add_twobody_operator(1,"d",2,"d",2,"d",1,"d",0.12500000d0)
    call ed_add_twobody_operator(2,"u",1,"u",1,"u",2,"u",0.12500000d0)
    call ed_add_twobody_operator(2,"d",1,"d",1,"d",2,"d",0.12500000d0)
    call ed_add_twobody_operator(1,"d",2,"u",2,"d",1,"u",0.12500000d0)
    call ed_add_twobody_operator(1,"u",2,"d",2,"u",1,"d",0.12500000d0)
    call ed_add_twobody_operator(2,"d",1,"u",1,"d",2,"u",0.12500000d0)
    call ed_add_twobody_operator(2,"u",1,"d",1,"u",2,"d",0.12500000d0)
    call ed_add_twobody_operator(1,"d",1,"u",2,"d",2,"u",0.12500000d0)
    call ed_add_twobody_operator(1,"u",1,"d",2,"u",2,"d",0.12500000d0)
    call ed_add_twobody_operator(2,"d",2,"u",1,"d",1,"u",0.12500000d0)
    call ed_add_twobody_operator(2,"u",2,"d",1,"u",1,"d",0.12500000d0)
  end subroutine set_twobody_hk




end program ed_normal_normal
