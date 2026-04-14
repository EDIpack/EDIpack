subroutine ed_get_Dimp_site(self,axis,z)
    complex(8),dimension(:),intent(inout)       :: self ! Green's function matrix
    character(len=*),optional                   :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
    complex(8),dimension(:),optional            :: z    ! User provided array of complex frequency where to evaluate Self
    character(len=1)                            :: axis_
    complex(8),dimension(:),allocatable         :: z_
    !TO BE COMPLETED
end subroutine  ed_get_Dimp_site