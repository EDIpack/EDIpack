MODULE ED_KRYLOV_SUPERC
  !:synopsis: Low-level Krylov routines, :code:`SUPERC` case
  USE ED_HAMILTONIAN_SUPERC
  implicit none
  private

  public :: KSC_Krylov_Basis_superc

contains

  subroutine KSC_Krylov_Basis_superc(isector,vvinit,alanc,blanc,norm2)
    integer,intent(in)                    :: isector
    complex(8),dimension(:),intent(inout) :: vvinit
    real(8),dimension(:),allocatable      :: alanc
    real(8),dimension(:),allocatable      :: blanc
    real(8)                               :: norm2

    call tridiag_Hv_sector_superc(isector,vvinit,alanc,blanc,norm2)
  end subroutine KSC_Krylov_Basis_superc

end MODULE ED_KRYLOV_SUPERC
