MODULE ED_RDM_NORMAL
  !:synopsis: Routines for Reduced Density Matrix, :code:`NORMAL` case
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
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private
  !
  public                           :: imp_rdm_normal
  
  real(8)                          :: Egs ! Ground-state energy
  integer                          :: iorb,jorb,iorb1,jorb1
  integer                          :: ispin,jspin
  integer                          :: isite,jsite
  integer                          :: ibath
  integer                          :: r,m,k,k1,k2,k3,k4
  integer                          :: iup,idw
  integer                          :: jup,jdw
  integer                          :: mup,mdw
  integer                          :: iph,i_el,isectorDim
  real(8)                          :: sgn,sgn1,sgn2,sg1,sg2,sg3,sg4
  real(8)                          :: gs_weight
  !
  real(8)                          :: peso
  real(8)                          :: norm
  !
  integer                          :: i,j,ii,io,jo
  integer                          :: isector,jsector
  !
  real(8)                          :: e_state
  real(8),dimension(:),allocatable :: v_state
  logical                          :: Jcondition
  !
  integer                          :: iImpUp,iImpDw
  integer                          :: jImpUp,jImpDw
  integer                          :: iBathUp,iBathDw
  integer                          :: iBup,iBdw

  integer                          :: LenBathUp,LenBathDw
  integer,allocatable              :: BathUp(:),BathDw(:)

  type(sector)                     :: sectorI,sectorJ
  character(len=128)               :: fmt

contains 


  subroutine imp_rdm_normal()
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
    !In the  :code:`normal` mode we operate each :math:`\sigma` independently as follows from the symmetry of the sectors.
    integer                         :: istate,val
    real(8),dimension(Norb)         :: nup,ndw
    !
    Egs     = state_list%emin
    !
    !IMPURITY DENSITY MATRIX
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")&
         "DEBUG observables_normal: eval impurity density matrix \rho_IMP = Tr_BATH(\rho)"
#endif
    if(allocated(impurity_density_matrix))deallocate(impurity_density_matrix)
    allocate(impurity_density_matrix(4**Norb,4**Norb))
    impurity_density_matrix=zero
    do istate=1,state_list%size
       !
       isector = es_return_sector(state_list,istate)
       e_state =  es_return_energy(state_list,istate)
       v_state =  es_return_dvec(state_list,istate)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(e_state-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then
          !
          select case(Nbath)
          case(0)
             !
             call build_sector(isector,sectorI) !no bath to trace
             !
             do IimpUp=0,2**Norb-1
                iUP = binary_search(sectorI%H(1)%map,IimpUp)
                if(iUP==0)cycle
                do JimpUp=0,2**Norb-1
                   jUP = binary_search(sectorI%H(1)%map,JimpUp)
                   if(jUP==0)cycle
                   do IimpDw=0,2**Norb-1
                      iDW = binary_search(sectorI%H(2)%map,IimpDw)
                      if(iDW==0)cycle
                      do JimpDw=0,2**Norb-1
                         jDW = binary_search(sectorI%H(2)%map,JimpDw)
                         if(jDW==0)cycle
                         !
                         i  = iUP + (iDW-1)*sectorI%DimUp
                         j  = jUP + (jDW-1)*sectorI%DimUp
                         !
                         io = (IimpUp + 2**Norb*IimpDw) + 1
                         jo = (JimpUp + 2**Norb*JimpDw) + 1
                         !
                         !(i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                         impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                              v_state(i)*v_state(j)*peso
                      enddo
                   enddo
                enddo
             enddo
             call delete_sector(sectorI)
             !
          case default
             !
             call build_sector(isector,sectorI,itrace=.true.)
             !
             do IimpUp=0,2**Norb-1
                do JimpUp=0,2**Norb-1
                   !
                   !Finding the unique bath states connecting IimpUp and JimpUp -> BATHup(:)
                   call sp_return_intersection(sectorI%H(1)%sp,IimpUp,JimpUp,BATHup,lenBATHup)
                   if(lenBATHup==0)cycle  !there are no bath states intersecting IimpUp,JimpUp
                   !
                   do IimpDw=0,2**Norb-1
                      do JimpDw=0,2**Norb-1
                         !
                         !Finding the unique bath states connecting IimpDw and JimpDw -> BATHdw(:)
                         call sp_return_intersection(sectorI%H(2)%sp,IimpDw,JimpDw,BATHdw,lenBATHdw)
                         if(lenBATHdw==0)cycle  !there are no bath states intersecting IimpDw,JimpDw
                         !
                         !=== >>> TRACE over bath states <<< 
                         do ibUP=1,lenBATHup
                            IbathUp = BATHup(ibUP)
                            do ibDW=1,lenBATHdw
                               IbathDw = BATHdw(ibDW)
                               !-----------------------------------------------------
                               !Allowed spin Fock space Istates:
                               !Iup = IimpUp +  2^Norb * IbathUp
                               !Idw = IimpDw +  2^Norb * IbathDw
                               !
                               !Corresponding sector indices per spin:
                               iUP= binary_search(sectorI%H(1)%map,IimpUp + 2**Norb*IbathUp)
                               iDW= binary_search(sectorI%H(2)%map,IimpDw + 2**Norb*IbathDw)
                               !
                               !Global sector index:
                               i  = iUP + (iDW-1)*sectorI%DimUp
                               !-----------------------------------------------------
                               !Allowed spin Fock space Jstates:
                               !Jup = JimpUp +  2^Norb * IbathUp
                               !Jdw = JimpDw +  2^Norb * IbathDw
                               !
                               !Corresponding sector jndices per spin:
                               jUP= binary_search(sectorI%H(1)%map,JimpUp + 2**Norb*IbathUp)
                               jDW= binary_search(sectorI%H(2)%map,JimpDw + 2**Norb*IbathDw)
                               !
                               !Global sector jndex:
                               j  = jUP + (jDW-1)*sectorI%DimUp
                               !-----------------------------------------------------
                               !Final [Up,Dw] composition for the impurity:
                               io = (IimpUp + 2**Norb*IimpDw) + 1
                               jo = (JimpUp + 2**Norb*JimpDw) + 1
                               !-----------------------------------------------------------------
                               !(i,j)_th contribution to the (io,jo)_th element of \rho_IMP
                               impurity_density_matrix(io,jo) = impurity_density_matrix(io,jo) + &
                                    v_state(i)*v_state(j)*peso
                               !-----------------------------------------------------------------
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
             call delete_sector(sectorI)
          end select
          !

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
  end subroutine imp_rdm_normal






end MODULE ED_RDM_NORMAL

















