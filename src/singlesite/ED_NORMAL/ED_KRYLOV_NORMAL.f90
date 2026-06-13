!KSC=Krylov-state complexity
!KOC=Krylov-operator complexity
MODULE ED_KRYLOV_NORMAL
  !:synopsis: Low-level Krylov routines, :code:`NORMAL` case
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private

  public :: KSC_Krylov_Basis_normal

contains

  subroutine KSC_Krylov_Basis_normal(isector,vvinit,alanc,blanc,norm2)
    integer,intent(in)                  :: isector
    real(8),dimension(:),intent(inout)    :: vvinit
    real(8),dimension(:),allocatable     :: alanc
    real(8),dimension(:),allocatable     :: blanc
    real(8)                             :: norm2

    call tridiag_Hv_sector_normal(isector,vvinit,alanc,blanc,norm2)
  end subroutine KSC_Krylov_Basis_normal

end MODULE ED_KRYLOV_NORMAL
