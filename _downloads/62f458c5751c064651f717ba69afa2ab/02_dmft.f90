program ed_hm_bethe
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  !> Number of discretization points for the Bethe DOS
  integer,parameter                           :: Le=5000
  !> Bethe half-bandwidth = energy unit
  real(8),parameter                           :: D=1d0
  !> Bethe DOS and linear dispersion
  real(8),dimension(Le)                       :: DOS
  real(8),dimension(Le)                       :: ene
  !>Bath:
  real(8),allocatable                         :: Bath(:)
  !> local dynamical functions [Nspin,Nspin,Norb,Norb,L] (could use other format)
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Smats,Sreal,Gmats,Greal,Weiss_
  !> input file
  character(len=16)                           :: finput
  !Aux variables:
  integer                                     :: iloop,Nb,i
  logical                                     :: converged
  complex(8)                                  :: zeta
  real(8)                                     :: de,error
  real(8),dimension(:),allocatable            :: wfreq
  real(8)                                     :: wmixing

  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))


  !> Check dimensions of the problem
  if(Nspin*Norb>1)stop "Nso>1. This is test."


  !> Build Bethe dos
  Ene = linspace(-D,D,Le,mesh=de)
  DOS = dens_bethe(Ene,D)


  !Allocate local Functions:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats)) ;Weiss = zero
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats));Weiss_= zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats)) ;Smats = zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats)) ;Gmats = zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal)) ;Sreal = zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal)) ;Greal = zero



  !> query bath dimensions and allocate the user bath
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb))

  !> init the ED solver
  call ed_init_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !> Solve the effective impurity problem
     call ed_solve(bath)



     !> Impurity Self-energy on Matsubara axis
     call ed_get_sigma(Smats,'m')

     !> Build a local Green's function using the Impurity Self-energy
     ! if(allocated(wfreq))deallocate(wfreq);allocate(wfreq(Lmats))
     wfreq = pi/beta*(2*arange(1,Lmats)-1)   !automatic Fortran allocation
     open(100,file="Gloc_iw.dat")
     do i=1,Lmats
        zeta= xi*wfreq(i)+xmu - Smats(1,1,1,1,i)
        Gmats(1,1,1,1,i) = sum(DOS(:)/( zeta-Ene(:) ))*de
        write(100,*)wfreq(i),dimag(Gmats(1,1,1,1,i)),dreal(Gmats(1,1,1,1,i))
     enddo
     close(100)


     !> Self-consistency: get the new Weiss field:
     Weiss(1,1,1,1,:) = one/Gmats(1,1,1,1,:) + Smats(1,1,1,1,:)
     Weiss(1,1,1,1,:) = one/Weiss(1,1,1,1,:)

     !> Mix to avoid trapping:
     if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_


     !> Close the self-consistency fitting the new bath:
     call ed_chi2_fitgf(Weiss,bath,ispin=1)


     !>Check convergence
     error = ( sum(abs(Weiss(1,1,1,1,:)-Weiss_(1,1,1,1,:))) / sum(abs(Weiss(1,1,1,1,:))) ) 

     converged = error < dmft_error
     write(*,*) error, dmft_error, converged
     
     !> save Weiss for next error check
     Weiss_ = Weiss
     
     call end_loop
  enddo


  call ed_get_sigma(Sreal,'r')
  wfreq = linspace(wini,wfin,Lreal) !automatic Fortran allocation
  open(100,file="Gloc_realw.dat")
  do i=1,Lreal
     zeta= dcmplx(wfreq(i),eps)+xmu - Sreal(1,1,1,1,i)
     Greal(1,1,1,1,i) = sum(DOS(:)/( zeta-Ene(:) ))*de
     write(100,*)wfreq(i),dimag(Greal(1,1,1,1,i)),dreal(Greal(1,1,1,1,i))
  enddo
  close(100)



end program ed_hm_bethe
