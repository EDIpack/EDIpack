MODULE ED_KRYLOV_NONSU2
  !:synopsis: Low-level Krylov routines, :code:`NONSU2` case
  USE ED_HAMILTONIAN_NONSU2
  implicit none
  private

  public :: KSC_Krylov_Basis_nonsu2

contains

  subroutine KSC_Krylov_Basis_nonsu2(isector,vvinit,alanc,blanc,norm2)
    integer,intent(in)                    :: isector
    complex(8),dimension(:),intent(inout) :: vvinit
    real(8),dimension(:),allocatable      :: alanc
    real(8),dimension(:),allocatable      :: blanc
    real(8)                               :: norm2

    call tridiag_Hv_sector_nonsu2(isector,vvinit,alanc,blanc,norm2)
  end subroutine KSC_Krylov_Basis_nonsu2

end MODULE ED_KRYLOV_NONSU2
