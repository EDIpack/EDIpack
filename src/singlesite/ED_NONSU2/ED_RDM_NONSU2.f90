MODULE ED_RDM_NONSU2
  !:synopsis: Routines for Reduced Density Matrix calculation, :code:`NONSU2` case
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_BATH
  USE ED_HAMILTONIAN_NONSU2
  implicit none
  private
  !
  public                           :: imp_rdm_nonsu2


  real(8)                             :: Egs ! Ground-state energy
  integer                             :: iorb,jorb,iorb1,jorb1
  integer                             :: r,m,k,k1,k2,k3,k4
  integer                             :: iup,idw
  integer                             :: jup,jdw
  integer                             :: mup,mdw
  integer                             :: iph,i_el,isectorDim
  real(8)                             :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                             :: gs_weight
  !
  real(8)                             :: peso
  real(8)                             :: norm
  !
  integer                             :: i,j,ii,jj,io,jo
  integer                             :: isector,jsector
  !
  real(8)                             :: e_state
  complex(8),dimension(:),allocatable :: v_state
  logical                             :: Jcondition
  !
  integer                             :: iImpUp,iImpDw,iImp
  integer                             :: jImpUp,jImpDw,jImp
  integer                             :: ib,iBath
  integer                             :: LenBath
  integer,allocatable                 :: Bath(:)
  type(sector)                        :: sectorI,sectorJ
  character(len=128)                  :: fmt

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Lanc method
  !+-------------------------------------------------------------------+
  subroutine imp_rdm_nonsu2()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    !Evaluates the RDM using the saved eigen-states in the :f:var:`state_list` and an efficient sparse algorithm.
    !For any given eigen-state :math:`|N\rangle` we proceed as follows.
    !Such state is a linear combination of the basis state in a given sector:
    !
    !:math:`|N\rangle = \sum_I c_I |I\rangle = \sum_{I}C_I |I_\uparrow\rangle|B_\uparrow\rangle|I_\downarrow\rangle|B_\downarrow\rangle`
    !
    !the goal is to build:
    !
    !:math:`\rho_{IJ} = \sum_{B_\sigma}|I_\uparrow\rangle|B_\uparrow\rangle|I_\downarrow\rangle|B_\downarrow\rangle \langle B_\downarrow|\langle J_\downarrow| \langle B_\uparrow| \langle J_\uparrow|`
    !
    !However, not all the bath configurations are admissibile in this sum. In fact if we store in a sparse map which bath configuration :math:`B_\sigma` (as integer)
    !corresponds to a given value of the impurity configuration :math:`I_\sigma`, then we can look for those value of :math:`B_\sigma` which *couples* to both
    !:math:`I_\sigma` and :math:`J_\sigma`. This reduces the sum to all and only the terms contributing to the RDM.
    !
    !In the  :code:`nonsu2` mode we need to enforce the link between different spin orientation imposed by the symmetry of the sectors.
    integer                 :: istate,val
    integer                 :: ibathUp,ibathDw
    integer                 :: IsgnUp,IsgnDw,JsgnUp,JsgnDw,nIdw,nJdw,nBup
    integer                 :: Isign,Jsign,signI,signJ,Sign
    real(8),dimension(Norb) :: nup,ndw
    !
    Egs     = state_list%emin
    !
    !IMPURITY DENSITY MATRIX
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_nonsu2: eval impurity density matrix \rho_IMP = Tr_BATH(\rho)"
#endif
    if(allocated(impurity_density_matrix))deallocate(impurity_density_matrix)
    allocate(impurity_density_matrix(4**Norb,4**Norb))
    impurity_density_matrix=zero

    do istate=1,state_list%size
       !
       isector = es_return_sector(state_list,istate)
       e_state =  es_return_energy(state_list,istate)
       v_state =  es_return_cvec(state_list,istate)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(e_state-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then

          select case(Nbath)
          case(0)
             !
             call build_sector(isector,sectorI)
             do IimpUp=0,2**Norb-1
                do IimpDw=0,2**Norb-1
                   !I: get the Fock state ii, search the corresponding sector i
                   ii= iImpUp + iImpDw*2**Ns
                   i = binary_search(sectorI%H(1)%map,ii)
                   if(i==0)cycle
                   do JimpUp=0,2**Norb-1
                      do JimpDw=0,2**Norb-1
                         !J: get the Fock state jj, search the corresponding sector j
                         jj= jImpUp + jImpDw*2**Ns
                         j = binary_search(sectorI%H(1)%map,jj)
                         if(j==0)cycle
                         !
                         !Build (i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                         io = (iImpUp+1) + 2**Norb*iImpDw
                         jo = (jImpUp+1) + 2**Norb*jImpDw
                         impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                              v_state(i)*conjg(v_state(j))*peso
                      enddo
                   enddo
                enddo
             enddo
             call delete_sector(sectorI)
             !
          case default
             !
             call build_sector(isector,sectorI,itrace=.true.)
             do IimpUp=0,2**Norb-1
                do IimpDw=0,2**Norb-1
                   do JimpUp=0,2**Norb-1
                      do JimpDw=0,2**Norb-1
                         !
                         !Build indices of the RDM in 1:4**Norb
                         iImp = iImpUp + iImpDw*2**Norb
                         jImp = jImpUp + jImpDw*2**Norb
                         call sp_return_intersection(sectorI%H(1)%sp,iImp,jImp,Bath,lenBath)
                         if(lenBATH==0)cycle
                         !
                         !=== >>> TRACE over bath states <<< =======
                         do ib=1,lenBath
                            iBath = Bath(ib)
                            !
                            !I: get the Fock state ii, search the corresponding sector i
                            ii= iImpUp + iImpDw*2**Ns + 2**Norb*iBath
                            i = binary_search(sectorI%H(1)%map,ii)
                            !
                            !J: get the Fock state jj, search the corresponding sector j
                            jj= jImpUp + jImpDw*2**Ns + 2**Norb*iBath
                            j = binary_search(sectorI%H(1)%map,jj)
                            !
                            !this is a safety measure which should never ever apply
                            if(i==0.OR.j==0)cycle
                            !
                            !Construct the sign of each components of RDM(io,jo)
                            nBup  = popcnt(Ibits(ii,Norb,Norb*Nbath))
                            nIdw  = popcnt(Ibits(ii,Ns,Norb))
                            nJdw  = popcnt(Ibits(jj,Ns,Norb))
                            signI = (-1)**(nIdw*nBup)
                            signJ = (-1)**(nJdw*nBup)
                            sign  = signI*signJ
                            !
                            !Build (i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                            io = (iImpUp+1) + 2**Norb*iImpDw
                            jo = (jImpUp+1) + 2**Norb*jImpDw
                            impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                                 v_state(i)*conjg(v_state(j))*peso*sign
                         enddo
                      enddo
                   enddo
                enddo
             enddo
             call delete_sector(sectorI)
             !
          end select

       endif
       !
       if(allocated(v_state))deallocate(v_state)
       !
    enddo
    !Useful prints to test:
    if(MPIMASTER .and. ed_verbose>2)then
       if(Norb==1)then
          write(LOGfile,*)"RDM:"
          do i=1,4**Norb
             write(LOGfile,"(*(F5.2,1x))")(dreal(impurity_density_matrix(i,j)),j=1,4**Norb)
          enddo
#ifdef _DEBUG
          ! Cfr Eq. 4 in Mod Phys Lett B 2013 27:05
          write(LOGfile,*)1-ed_dens_up(1)-ed_dens_dw(1)+ed_docc(1),abs(1-ed_dens_up(1)-ed_dens_dw(1)+ed_docc(1)-impurity_density_matrix(1,1))
          write(LOGfile,*)ed_dens_up(1)-ed_docc(1),abs(ed_dens_up(1)-ed_docc(1)-impurity_density_matrix(2,2))
          write(LOGfile,*)ed_dens_dw(1)-ed_docc(1),abs(ed_dens_dw(1)-ed_docc(1)-impurity_density_matrix(3,3))
          write(LOGfile,*)ed_docc(1),abs(ed_docc(1)-impurity_density_matrix(4,4))
#endif
       endif
    endif
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")""
#endif
  end subroutine imp_rdm_nonsu2






end MODULE ED_RDM_NONSU2

















