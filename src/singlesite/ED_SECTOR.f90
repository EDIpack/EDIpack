MODULE ED_SECTOR
  !:synopsis: Routines for Fock space sectors creation and manipulation
  !Contains procedures to construct the symmetry sectors corresponding to a given set of quantum numbers :math:`\vec{Q}`, in particular it allocated and build the  :f:var:`sector_map` connecting the states of a given sector with the corresponding Fock ones. 
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  interface map_allocate
     !
     ! Allocate the map(s) connecting the states in the symmetry sector to thos in the Fock space. 
     !
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate

  interface flip_state
     module procedure :: flip_state_normal
     module procedure :: flip_state_other
  end interface flip_state

  interface get_Sector
     !
     ! Returns the index of :f:var:`isector` of the symmetry sector given the quantum numbers :f:var:`qn`
     !
     module procedure :: get_Sector_normal
     module procedure :: get_Sector_superc
     module procedure :: get_Sector_nonsu2
  end interface get_Sector

  interface get_QuantumNumbers
     !
     !Returns the quantum numbers :f:var:`qn` given the  index of :f:var:`isector` of the symmetry sector
     !
     module procedure :: get_QuantumNumbers_normal
     module procedure :: get_QuantumNumbers_other
  end interface get_QuantumNumbers



  interface apply_op_C
     !
     ! Apply the Fock operator :f:math:`C` to an input Fock state :f:var:`V` in a given sector :f:var:`sectorI`.
     ! Output: a Fock state :f:var:`OV` in the sector :f:var:`sectorJ` with one less particle as allocatable double real or double complex allocatable type.  
     !
     module procedure :: apply_op_C_d
     module procedure :: apply_op_C_c
  end interface apply_op_C



  interface apply_op_CDG
     !
     ! Apply the Fock operator :f:math:`C` to an input Fock state :f:var:`V` in a given sector :f:var:`sectorI`.
     ! Output: a Fock state :f:var:`OV` in the sector :f:var:`sectorJ` with one less particle as allocatable double real or double complex allocatable type.  
     !
     module procedure :: apply_op_CDG_d
     module procedure :: apply_op_CDG_c
  end interface apply_op_CDG



  interface apply_COps
     !
     ! Apply a linear combination of creation/destriction Fock operators :f:math:`\vec{C}=\sum_i a_i C^{o_i}_{p_i}` where :f:math:`a_i` are real or cmplx coeffients, :f:math:`o_i=\pm 1` indicate
     ! respectively creation or destruction operators and :f:math:`p_i` are the respective positions. 
     ! The input Fock state :f:var:`V` in a given sector :f:var:`sectorI`.
     ! Output: a Fock state :f:var:`OV` in the sector :f:var:`sectorJ` with one more particle as allocatable double real or double complex allocatable type.  
     !
     module procedure :: apply_COps_d
     module procedure :: apply_COps_c
  end interface apply_COps



  interface apply_op_Sz
     !
     ! Apply the Fock operator :f:math:`S_z` to an input Fock state :f:var:`V` in a given sector :f:var:`sectorI`.
     ! Output: a Fock state :f:var:`OV`.  
     !
     module procedure :: apply_op_Sz_d
     module procedure :: apply_op_Sz_c
  end interface apply_op_Sz



  interface apply_op_N
     !
     ! Apply the Fock operator :f:math:`N` to an input Fock state :f:var:`V` in a given sector :f:var:`sectorI`.
     ! Output: a Fock state :f:var:`OV`.
     !
     module procedure :: apply_op_N_d
     module procedure :: apply_op_N_c
  end interface apply_op_N


  public :: build_sector
  public :: delete_sector
  !
  public :: apply_Cops
  public :: apply_op_C
  public :: apply_op_CDG
  public :: apply_op_Sz
  public :: apply_op_N
  public :: build_op_Ns
  !
  public :: get_Sector
  public :: get_QuantumNumbers
  public :: get_Nup
  public :: get_Ndw
  public :: get_Sz
  public :: get_Ntot
  public :: get_DimUp
  public :: get_DimDw
  public :: get_Dim
  !
  public :: indices2state
  public :: state2indices
  public :: iup_index
  public :: idw_index


  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state






contains




  !##################################################################
  !##################################################################
  !BUILD SECTORS
  !
  !About sparse map:
  ! iImp = Iup + Bup*2**Norb + (Idw + Bdw*2**Norb)*2**Ns
  !      = Iup + Idw*2**Ns + (Bup+Bdw*2**Ns)*2**Norb
  !      = I_ + B_*2**Norb
  ! so that iImp should be: 
  ! iImp = ImpUp + ImpDw*2^Ns (1)
  ! and not
  ! iImp = ImpUp + ImpDw*2^Nimp (2)
  ! However (2) is used to have Imp states in 1:4**Norb, the size of RDM.
  !##################################################################
  !##################################################################
  subroutine build_sector(isector,self,iTrace)
    !
    ! This procedure build the sector :f:type:`sector` :math:`{\cal S}(\vec{Q})` given the index :f:var:`isector` which univocally determines the quantum numbers :math:`\vec{Q}`. All the components of the data structure :f:type:`sector` are determined here. Moreover the map connecting the states of the sector to thos in the Fock space is constructed.
    ! To avoid integer overflow the loop over Fock space states is decomposed in spin up and down integer so that :math:`I = I_\uparrow + 2^{N}I_\downarrow`.  
    !
    integer,intent(in) :: isector
    type(sector)       :: self
    logical,optional   :: itrace
    logical            :: itrace_
    integer            :: iup,idw
    integer            :: nup_,ndw_,sz_,nt_
    integer            :: twoSz_,twoLz_,twoJz
    integer            :: dim,iud,i,iorb,impDIM
    integer            :: ivec(Ns),jvec(Ns)
    integer            :: iImp,iImpUp,iImpDw
    integer            :: iBath,iBathUp,iBathDw
    !
    itrace_=.false. ; if(present(itrace))itrace_=itrace
    !
    if(self%status)call delete_sector(self)
    !
    self%index = isector
    !
    select case(ed_mode)
    case default
       !
       impDIM = 2**(Norb/Ns_ud)
       !
       allocate(self%H(2*Ns_Ud))
       allocate(self%DimUps(Ns_Ud))
       allocate(self%DimDws(Ns_Ud))
       allocate(self%Nups(Ns_Ud))
       allocate(self%Ndws(Ns_Ud))
       !
       call get_Nup(isector,self%Nups);self%Nup=sum(self%Nups)
       call get_Ndw(isector,self%Ndws);self%Ndw=sum(self%Ndws)
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A,2"//str(Ns_Ud)//"I3)")&
            "DEBUG build_sector: sector:"//str(isector)//"- Nups,Ndws",self%Nups,self%Ndws
#endif
       call get_DimUp(isector,self%DimUps);self%DimUp=product(self%DimUps)
       call get_DimDw(isector,self%DimDws);self%DimDw=product(self%DimDws)
       self%DimEl=self%DimUp*self%DimDw
       self%DimPh=Nph+1
       self%Dim=self%DimEl*self%DimPh
       !
       if(itrace_)then
          call map_allocate(self%H,[self%DimUps,self%DimDws],impDIM)
       else
          call map_allocate(self%H,[self%DimUps,self%DimDws])
       endif
       !
       do iud=1,Ns_Ud
          !UP    
          dim=0
          do iup=0,2**Ns_Orb-1
             nup_ = popcnt(iup)
             if(nup_ /= self%Nups(iud))cycle
             dim  = dim+1
             self%H(iud)%map(dim) = iup
             if(.not.itrace_)cycle
             iImp  = ibits(iup,0,Norb)
             iBath = ibits(iup,Norb,Norb*Nbath)
             call sp_insert_state(self%H(iud)%sp,iImp,iBath,dim)
          enddo
          !DW
          dim=0
          do idw=0,2**Ns_Orb-1
             ndw_= popcnt(idw)
             if(ndw_ /= self%Ndws(iud))cycle
             dim = dim+1
             self%H(iud+Ns_Ud)%map(dim) = idw
             if(.not.itrace_)cycle
             iIMP  = ibits(idw,0,Norb)
             iBATH = ibits(idw,Norb,Norb*Nbath)
             call sp_insert_state(self%H(iud+Ns_Ud)%sp,iImp,iBath,dim)
          enddo
       enddo
       !
    case ("superc")
       !
       impDIM = 4**Norb
       !
       allocate(self%H(1))
       self%Sz    = getSz(isector)
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A,I4)")&
            "DEBUG build_sector: sector:"//str(isector)//"- Sz",self%Sz
#endif
       self%DimPh = Nph+1
       self%DimEl = getDim(isector)/(Nph+1)
       self%Dim   = self%DimEl*self%DimPh
       if(itrace_)then
          call map_allocate(self%H,[self%DimEl],impDIM)
       else
          call map_allocate(self%H,[self%DimEl])
       endif
       dim=0
       do idw=0,2**Ns-1
          ndw_= popcnt(idw)
          do iup=0,2**Ns-1
             nup_ = popcnt(iup)
             sz_  = nup_ - ndw_
             if(sz_ == self%Sz)then
                dim=dim+1
                self%H(1)%map(dim) = iup + idw*2**Ns
                if(.not.itrace_)cycle
                iImpUp  = ibits(iup,0,Norb)
                iImpDw  = ibits(idw,0,Norb)
                iBathUp = ibits(iup,Norb,Norb*Nbath)
                iBathDw = ibits(idw,Norb,Norb*Nbath)
                iImp    = iImpUp  + iImpDw*(2**Norb)
                iBath   = iBathUp + iBathDw*(2**Ns)
                call sp_insert_state(self%H(1)%sp,iImp,iBath,dim)
             endif
          enddo
       enddo

    case ("nonsu2")
       !
       impDIM = 4**Norb
       !
       allocate(self%H(1))
       !
       if(Jz_basis)then
          self%Ntot  = getN(isector)
#ifdef _DEBUG
          if(ed_verbose>3)write(Logfile,"(A,I4)")&
               "DEBUG build_sector: sector:"//str(isector)//"- N",self%Ntot
#endif
          self%DimPh = Nph+1
          self%DimEl = getDim(isector)/(Nph+1)
          self%Dim   = self%DimEl*self%DimPh
          self%twoJz = gettwoJz(isector)
          if(itrace_)then
             call map_allocate(self%H,[self%DimEl],impDIM)
          else
             call map_allocate(self%H,[self%DimEl])
          endif
          dim=0
          do idw=0,2**Ns-1
             ndw_= popcnt(idw)
             jvec = bdecomp(idw,Ns)
             do iup=0,2**Ns-1
                nup_ = popcnt(iup)
                ivec = bdecomp(iup,Ns)
                nt_     = nup_ + ndw_
                twoSz_  = nup_ - ndw_
                twoLz_  = 0
                do ibath=0,Nbath
                   do iorb=1,Norb
                      twoLz_ = twoLz_ + 2 * Lzdiag(iorb) * ivec(iorb+Norb*ibath)  &
                           + 2 * Lzdiag(iorb) * jvec(iorb+Norb*ibath)
                   enddo
                enddo
                if(self%Ntot == nt_  .AND. self%twoJz==(twoSz_+twoLz_) )then
                   dim=dim+1
                   self%H(1)%map(dim) = iup + idw*2**Ns
                   if(.not.itrace_)cycle
                   iImpUp  = ibits(iup,0,Norb)
                   iImpDw  = ibits(idw,0,Norb)
                   iBathUp = ibits(iup,Norb,Norb*Nbath)
                   iBathDw = ibits(idw,Norb,Norb*Nbath)
                   iImp    = iImpUp  + iImpDw*(2**Norb)
                   iBath   = iBathUp + iBathDw*(2**Ns)
                   call sp_insert_state(self%H(1)%sp,iImp,iBath,dim)
                endif
             enddo
          enddo
          !
       else
          !
          self%Ntot  = getN(isector)
#ifdef _DEBUG
          if(ed_verbose>3)write(Logfile,"(A,I4)")&
               "DEBUG build_sector: sector:"//str(isector)//"- N",self%Ntot
#endif
          self%DimPh = Nph+1
          self%DimEl = getDim(isector)/(Nph+1)
          self%Dim   = self%DimEl*self%DimPh
          if(itrace_)then
             call map_allocate(self%H,[self%DimEl],impDIM)
          else
             call map_allocate(self%H,[self%DimEl])
          endif
          dim=0
          do idw=0,2**Ns-1
             ndw_= popcnt(idw)
             do iup=0,2**Ns-1
                nup_ = popcnt(iup)
                nt_  = nup_ + ndw_
                if(nt_ /= self%Ntot)cycle
                dim=dim+1
                self%H(1)%map(dim) = iup + idw*2**Ns
                if(.not.itrace_)cycle
                iImpUp  = ibits(iup,0,Norb)
                iImpDw  = ibits(idw,0,Norb)
                iBathUp = ibits(iup,Norb,Norb*Nbath)
                iBathDw = ibits(idw,Norb,Norb*Nbath)
                iImp    = iImpUp  + iImpDw*(2**Norb)
                iBath   = iBathUp + iBathDw*(2**Ns)
                call sp_insert_state(self%H(1)%sp,iImp,iBath,dim)
             enddo
          enddo
       endif
    end select
    self%Nlanc  = min(self%Dim,lanc_nGFiter)
    self%status = .true.
  end subroutine build_sector


  subroutine delete_sector(self)
    type(sector) :: self
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,I4)")"DEBUG delete_sector"
#endif
    call map_deallocate(self%H)
    if(allocated(self%H))deallocate(self%H)
    if(allocated(self%DimUps))deallocate(self%DimUps)
    if(allocated(self%DimDws))deallocate(self%DimDws)
    if(allocated(self%Nups))deallocate(self%Nups)
    if(allocated(self%Ndws))deallocate(self%Ndws)
    self%index=0
    self%DimUp=0
    self%DimDw=0
    self%Dim=0
    self%Nup=0
    self%Ndw=0
    self%Sz=-1000
    self%Ntot=-1
    self%twoJz=-1000
    self%Nlanc=0
    self%status=.false.
  end subroutine delete_sector







  subroutine map_allocate_scalar(H,N,Nsp)
    type(sector_map) :: H
    integer          :: N
    integer,optional :: Nsp
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    if(present(Nsp))call sp_init_map(H%sp,Nsp)
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N,Nsp)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer,optional                    :: Nsp
    integer                             :: i
    do i=1,size(H)
       if(present(Nsp))then
          call map_allocate_scalar(H(i),N(i),Nsp)
       else
          call map_allocate_scalar(H(i),N(i))
       endif
    enddo
  end subroutine map_allocate_vector



  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(LOGfile,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector






  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################






  function apply_op_C_d(V,iorb,ispin,isector,jsector) result(OV)
    real(8),dimension(:),intent(in)  :: V
    integer, intent(in)              :: iorb,ispin,isector,jsector
    type(sector)                     :: sectorI,sectorJ
    real(8),dimension(:),allocatable :: OV
    real(8)                          :: sgn
    integer                          :: ialfa,ibeta,ipos,isite
    integer                          :: i,j,r
    integer                          :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)       :: Indices
    integer,dimension(2*Ns_Ud)       :: Jndices
    integer,dimension(2,Ns_Orb)      :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)             :: Iud
    integer,dimension(2*Ns)          :: ib
    character(len=100)               :: fmt_print
    !
    if(MpiMaster)then
       !
       if(ed_total_ud)then
          ialfa = 1
          ipos  = iorb
       else
          ialfa = iorb
          ipos  = 1
       endif
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_C ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=0d0
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case ("normal")
             write(fmt_print, '(A,I0,A,I0,A)') '(A,I6,', 2*Ns_ud, 'I4,A,I6,', 2*Ns_ud, 'I4)' 
             write(LOGfile,fmt_print)&
                  'From:',sectorI%index,sectorI%Nups,sectorI%Ndws,&
                  ' -> apply C  :',sectorJ%index,sectorJ%Nups,sectorJ%Ndws
          case default;stop "apply_Op_C ERROR: called with ed_mode != normal"
          end select
       endif
       !
       do i=1,sectorI%Dim
          ibeta = ialfa + (ispin-1)*Ns_Ud
          iph   = (i-1)/(sectorI%DimEl) + 1
          i_el  = mod(i-1,sectorI%DimEl) + 1
          !
          call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
          iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
          iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
          nud(1,:) = Bdecomp(iud(1),Ns_Orb)
          nud(2,:) = Bdecomp(iud(2),Ns_Orb)
          if(Nud(ispin,ipos)/=1)cycle
          call c(ipos,iud(ispin),r,sgn)
          Jndices        = Indices
          Jndices(ibeta) = binary_search(sectorJ%H(ibeta)%map,r)
          call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
          j     = j + (iph-1)*sectorJ%DimEl
          OV(j) = sgn*V(i)
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=0d0
    end if
    !
  end function apply_op_C_d



  function apply_op_C_c(V,iorb,ispin,isector,jsector) result(OV)
    complex(8),dimension(:),intent(in)  :: V
    integer, intent(in)                 :: iorb,ispin,isector,jsector
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: OV
    real(8)                             :: sgn
    integer                             :: ialfa,ibeta,isite
    integer                             :: i,j,r
    integer                             :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(2*Ns_Ud)          :: Jndices
    integer,dimension(2,Ns_Orb)         :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                :: Iud
    integer,dimension(2*Ns)             :: ib
    !
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_C ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default;stop "apply_op_C ERROR: called with ed_mode!=superc/nonsu2"
          case ("superc")
             write(LOGfile,"(A,I6,I3,A,I6,I3)")&
                  'From:',sectorI%index,sectorI%Sz,&
                  ' -> apply C  :',sectorJ%index,sectorJ%Sz
          case ("nonsu2")
             if(Jz_basis)then
                write(LOGfile,"(A,I6,I3,A,I6,I3)")&
                     'From:',sectorI%index,sectorI%twoJz/2.,&
                     ' -> apply C  :',sectorJ%index,sectorJ%twoJz/2
             else
                write(LOGfile,"(A,I6,I3,A,I6,I3)")&
                     'From:',sectorI%index,sectorI%Ntot,&
                     ' -> apply C  :',sectorJ%index,sectorJ%Ntot
             endif
          end select
       endif
       !
       do i=1,sectorI%Dim
          isite= iorb + (ispin-1)*Ns
          iph  = (i-1)/(sectorI%DimEl)+1
          i_el = mod(i-1,sectorI%DimEl)+1
          ei   = sectorI%H(1)%map(i_el)
          ib   = bdecomp(ei,2*Ns)
          if(ib(isite)/=1)cycle
          call c(isite,ei,r,sgn)
          j_el = binary_search(sectorJ%H(1)%map,r)
          j    = j_el + (iph-1)*sectorJ%DimEl
          !
          OV(j) = sgn*V(i)
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
    !
  end function apply_op_C_c













  function apply_op_CDG_d(V,iorb,ispin,isector,jsector) result(OV)
    real(8),dimension(:),intent(in)  :: V
    integer, intent(in)              :: iorb,ispin,isector,jsector
    type(sector)                     :: sectorI,sectorJ
    real(8),dimension(:),allocatable :: OV
    real(8)                          :: sgn
    integer                          :: ialfa,ibeta,ipos,isite
    integer                          :: i,j,r
    integer                          :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)       :: Indices
    integer,dimension(2*Ns_Ud)       :: Jndices
    integer,dimension(2,Ns_Orb)      :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)             :: Iud
    integer,dimension(2*Ns)          :: ib
    character(len=100)               :: fmt_print
    !
    !
    if(MpiMaster)then
       !
       if(ed_total_ud)then
          ialfa = 1
          ipos  = iorb
       else
          ialfa = iorb
          ipos  = 1
       endif
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_CDG ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=0d0
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case ("normal")
             write(fmt_print, '(A,I0,A,I0,A)') '(A,I6,', 2*Ns_ud, 'I4,A,I6,', 2*Ns_ud, 'I4)'
             write(LOGfile,fmt_print)&
                  'From:',sectorI%index,sectorI%Nups,sectorI%Ndws,&
                  ' -> apply C^+:',sectorJ%index,sectorJ%Nups,sectorJ%Ndws
          case default;stop "apply_op_CDG ERROR: called with ed_mode/=normal"
          end select
       endif
       !
       do i=1,sectorI%Dim
          ibeta  = ialfa + (ispin-1)*Ns_Ud
          iph = (i-1)/(sectorI%DimEl) + 1
          i_el = mod(i-1,sectorI%DimEl) + 1
          !
          call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
          iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
          iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
          nud(1,:) = Bdecomp(iud(1),Ns_Orb)
          nud(2,:) = Bdecomp(iud(2),Ns_Orb)
          if(Nud(ispin,ipos)/=0)cycle
          call cdg(ipos,iud(ispin),r,sgn)
          Jndices        = Indices
          Jndices(ibeta) = binary_search(sectorJ%H(ibeta)%map,r)
          call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j_el)
          j     = j_el + (iph-1)*sectorJ%DimEl
          OV(j) = sgn*V(i)
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=0d0
    end if
    !
  end function apply_op_CDG_d


  function apply_op_CDG_c(V,iorb,ispin,isector,jsector) result(OV)
    complex(8),dimension(:),intent(in)  :: V
    integer, intent(in)                 :: iorb,ispin,isector,jsector
    type(sector)                        :: sectorI,sectorJ
    complex(8),dimension(:),allocatable :: OV
    real(8)                             :: sgn
    integer                             :: ialfa,ibeta,isite
    integer                             :: i,j,r
    integer                             :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(2*Ns_Ud)          :: Jndices
    integer,dimension(2,Ns_Orb)         :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                :: Iud
    integer,dimension(2*Ns)             :: ib
    !
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_CDG ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default;stop "apply_op_CDG ERROR: called with ed_mode!=superc/nonsu2"
          case ("superc")
             write(LOGfile,"(A,I6,I3,A,I6,I3)")&
                  'From:',sectorI%index,sectorI%Sz,&
                  ' -> apply C^+:',sectorJ%index,sectorJ%Sz
          case ("nonsu2")
             if(Jz_basis)then
                write(LOGfile,"(A,I6,I3,A,I6,I3)")&
                     'From:',sectorI%index,sectorI%twoJz/2.,&
                     ' -> apply C^+:',sectorJ%index,sectorJ%twoJz/2
             else
                write(LOGfile,"(A,I6,I3,A,I6,I3)")&
                     'From:',sectorI%index,sectorI%Ntot,&
                     ' -> apply C^+:',sectorJ%index,sectorJ%Ntot
             endif
          end select
       endif
       !
       do i=1,sectorI%Dim
          isite= iorb + (ispin-1)*Ns
          iph  = (i-1)/(sectorI%DimEl)+1
          i_el = mod(i-1,sectorI%DimEl) + 1
          ei   = sectorI%H(1)%map(i_el)
          ib   = bdecomp(ei,2*Ns)
          if(ib(isite)/=0)cycle
          call cdg(isite,ei,r,sgn)
          j_el = binary_search(sectorJ%H(1)%map,r)
          j    = j_el + (iph-1)*sectorJ%DimEl
          OV(j) = sgn*V(i)
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
    !
  end function apply_op_CDG_c











  function apply_COps_d(V,As,Os,Pos,Spin,isector,jsector) result(OV)
    real(8),dimension(:),intent(in)        :: V
    real(8),dimension(:),intent(in)        :: As
    integer,dimension(size(As)),intent(in) :: Os
    integer,dimension(size(As)),intent(in) :: Pos,Spin
    integer, intent(in)                    :: isector,jsector
    type(sector)                           :: sectorI,sectorJ
    real(8),dimension(:),allocatable       :: OV
    integer                                :: ipos,ispin,ios,isite
    integer                                :: i,j,Nos,is,N,iph,i_el,j_el,ei
    real(8)                                :: sgn
    integer                                :: r
    integer                                :: fi
    integer,dimension(2*Ns_Ud)             :: Indices
    integer,dimension(2*Ns_Ud)             :: Jndices
    integer,dimension(2,Ns_Orb)            :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                   :: Iud
    integer,dimension(2*Ns)                :: ib
    character(2),dimension(-1:1)           :: Cstr = ["C ","  ","C*"]
    character(:),allocatable               :: a,sg,Ostr
    !
    if(ed_mode == "normal" .and. .not.ed_total_ud)stop "apply_COps ERROR: called with ed_total_ud=F"
    !       
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_COps ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          Ostr = ""
          do is=1,size(As)
             ipos  = Pos(is)
             ispin = Spin(is)
             ios   = Os(is)
             sg = "";if(is>1)sg=" + "
             a  = str(As(is),2)
             if(As(is)==one)then
                a=""
                if(is>1)sg = " + "
             endif
             if(As(is)==-one)then
                a=""
                sg = " - "
             endif
             if(As(is)==xi)then
                a="i."
                if(is>1)sg = " + "
             endif
             if(As(is)==-xi)then
                a="i."
                sg = " - "
             endif
             if(As(is)==one+xi)then
                a="(1+i)."
                if(is>1)sg = " + "
             endif
             if(As(is)==-one+xi)then
                a="(-1+i)."
                sg = " + "
             endif
             if(As(is)==one-xi)then
                a="(1-i)."
                if(is>1)sg = " + "
             endif
             if(As(is)==-one-xi)then
                a="(1+i)."
                sg = " - "
             endif
             Ostr  = Ostr//sg//a//str(Cstr(ios))//"_l"//str(ipos)//"s"//str(ispin)
          enddo
          N = max(20,len(Ostr))
          !
          select case(ed_mode)
          case default;stop "apply_op_Cops ERROR: called with ed_mode/=normal"
          case ("normal")
             write(LOGfile,"(A,I6,2I4,A,A"//str(N)//",I6,2I4)")&
                  'From:',sectorI%index,sectorI%Nups,sectorI%Ndws,&
                  ' -> apply:',Ostr,sectorJ%index,sectorJ%Nups,sectorJ%Ndws
          end select
       endif
       !
       do is=1,size(As)
          ipos  = Pos(is)
          ispin = Spin(is)
          ios   = Os(is)
          do i=1,sectorI%Dim
             iph  = (i-1)/(sectorI%DimEl) + 1
             i_el = mod(i-1,sectorI%DimEl) + 1
             call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
             iud(1)   = sectorI%H(1)%map(Indices(1))
             iud(2)   = sectorI%H(2)%map(Indices(2))
             Nud(1,:) = Bdecomp(iud(1),Ns_Orb)
             Nud(2,:) = Bdecomp(iud(2),Ns_Orb)             
             select case(ios)
             case default;stop "apply_COps ERROR: ios sign not \in [-1,1]"
             case(-1)
                if(Nud(ispin,ipos)/=1)cycle
                call c(ipos,iud(ispin),r,sgn)
             case(1)
                if(Nud(ispin,ipos)/=0)cycle
                call cdg(ipos,iud(ispin),r,sgn)
             end select
             Jndices        = Indices
             Jndices(ispin) = binary_search(sectorJ%H(ispin)%map,r)
             call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j_el)
             j     = j_el + (iph-1)*sectorJ%DimEl
             OV(j) = OV(j) + sgn*V(i)*As(is)
          enddo
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=0d0
    end if
  end function apply_COps_d





  function apply_COps_c(V,As,Os,Pos,Spin,isector,jsector) result(OV)
    complex(8),dimension(:),intent(in)     :: V
    complex(8),dimension(:),intent(in)     :: As
    integer,dimension(size(As)),intent(in) :: Os
    integer,dimension(size(As)),intent(in) :: Pos,Spin
    integer, intent(in)                    :: isector,jsector
    type(sector)                           :: sectorI,sectorJ
    complex(8),dimension(:),allocatable    :: OV
    integer                                :: ipos,ispin,ios,isite
    integer                                :: i,j,Nos,is,N,iph,i_el,j_el,ei
    real(8)                                :: sgn
    integer                                :: r
    integer                                :: fi
    integer,dimension(2*Ns_Ud)             :: Indices
    integer,dimension(2*Ns_Ud)             :: Jndices
    integer,dimension(2,Ns_Orb)            :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                   :: Iud
    integer,dimension(2*Ns)                :: ib
    character(2),dimension(-1:1)           :: Cstr = ["C ","  ","C*"]
    character(:),allocatable               :: a,sg,Ostr
    !
    if(ed_mode == "normal" .and. .not.ed_total_ud)stop "apply_COps ERROR: called with ed_total_ud=F"
    !
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_COps ERROR: size(V) != sectorI.Dim"
       !
       call build_sector(jsector,sectorJ)
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorJ%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          Ostr = ""
          do is=1,size(As)
             ipos  = Pos(is)
             ispin = Spin(is)
             ios   = Os(is)
             sg = "";if(is>1)sg=" + "
             a  = str(As(is),2)
             if(As(is)==one)then
                a=""
                if(is>1)sg = " + "
             endif
             if(As(is)==-one)then
                a=""
                sg = " - "
             endif
             if(As(is)==xi)then
                a="i."
                if(is>1)sg = " + "
             endif
             if(As(is)==-xi)then
                a="i."
                sg = " - "
             endif
             if(As(is)==one+xi)then
                a="(1+i)."
                if(is>1)sg = " + "
             endif
             if(As(is)==-one+xi)then
                a="(-1+i)."
                sg = " + "
             endif
             if(As(is)==one-xi)then
                a="(1-i)."
                if(is>1)sg = " + "
             endif
             if(As(is)==-one-xi)then
                a="(1+i)."
                sg = " - "
             endif
             Ostr  = Ostr//sg//a//str(Cstr(ios))//"_l"//str(ipos)//"s"//str(ispin)
          enddo
          N = max(20,len(Ostr))
          !
          select case(ed_mode)
          case default;stop "apply_op_CDG ERROR: called with ed_mode!=superc/nonsu2"
          case ("superc")
             write(LOGfile,"(A,I6,I3,A,A"//str(N)//",I6,I3))")&
                  'From:',sectorI%index,sectorI%Sz,&
                  ' -> apply:',Ostr,sectorJ%index,sectorJ%Sz
          case ("nonsu2")
             if(Jz_basis)then                    
                write(LOGfile,"(A,I6,I3,A,A"//str(N)//",I6,I3))")&
                     'From:',sectorI%index,sectorI%twoJz/2.,&
                     ' -> apply:',Ostr,sectorJ%index,sectorJ%twoJz/2
             else
                write(LOGfile,"(A,I6,I3,A,A"//str(N)//",I6,I3))")&
                     'From:',sectorI%index,sectorI%Ntot,&
                     ' -> apply:',Ostr,sectorJ%index,sectorJ%Ntot
             endif
          end select
       endif
       !
       do is=1,size(As)
          ipos  = Pos(is)
          ispin = Spin(is)
          ios   = Os(is)
          do i=1,sectorI%Dim
             isite= ipos + (ispin-1)*Ns
             iph  = (i-1)/(sectorI%DimEl)+1
             i_el = mod(i-1,sectorI%DimEl) + 1
             fi   = sectorI%H(1)%map(i_el)
             ib   = bdecomp(fi,2*Ns)
             select case(ios)
             case default;stop "apply_COps ERROR: ios sign not \in [-1,1]"
             case(-1)
                if(ib(isite)/=1)cycle
                call c(isite,fi,r,sgn)
             case(1)
                if(ib(isite)/=0)cycle
                call cdg(isite,fi,r,sgn)
             end select
             j_el = binary_search(sectorJ%H(1)%map,r)
             j    = j_el + (iph-1)*sectorJ%DimEl
             OV(j) = OV(j) + sgn*V(i)*As(is)
          enddo
       enddo
       call delete_sector(sectorI)
       call delete_sector(sectorJ)
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
  end function apply_COps_c











  function apply_op_Sz_d(V,iorb,isector) result(OV)
    real(8),dimension(:),intent(in)  :: V
    integer, intent(in)              :: iorb,isector
    real(8),dimension(:),allocatable :: OV
    type(sector)                     :: sectorI
    real(8)                          :: sgn
    integer                          :: ialfa,ibeta,ipos,isite
    integer                          :: i,j,r
    integer                          :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)       :: Indices
    integer,dimension(2*Ns_Ud)       :: Jndices
    integer,dimension(2,Ns_Orb)      :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)             :: Iud
    integer,dimension(2*Ns)          :: ib
    !
    if(ed_total_ud)then
       ialfa = 1
       ipos  = iorb
    else
       ialfa = iorb
       ipos  = 1
    endif
    !
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_Sz ERROR: size(V) != sectorI.Dim"
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorI%Dim)) ; OV=0d0
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default;stop "apply_op_Sz ERROR: called with ed_mode/=normal"
          case ("normal")
             write(LOGfile,"(A,I6,2I4)")&
                  'apply Sz :',sectorI%index,sectorI%Nups,sectorI%Ndws
          end select
       endif
       !
       do i=1,sectorI%Dim
          iph = (i-1)/(sectorI%DimEl) + 1
          i_el = mod(i-1,sectorI%DimEl) + 1
          call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
          iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
          iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
          nud(1,:) = Bdecomp(iud(1),Ns_Orb)
          nud(2,:) = Bdecomp(iud(2),Ns_Orb)
          sgn      = dble(nud(1,ipos))-dble(nud(2,ipos))
          sgn      = sgn/2d0
          OV(i)    = sgn*V(i)
       enddo
       !
       call delete_sector(sectorI)
       !
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=0d0
    end if
    !
  end function apply_op_Sz_d




  function apply_op_Sz_c(V,iorb,isector) result(OV)
    complex(8),dimension(:),intent(in)  :: V
    integer, intent(in)                 :: iorb,isector
    type(sector)                        :: sectorI
    complex(8),dimension(:),allocatable :: OV
    real(8)                             :: sgn
    integer                             :: ialfa,ibeta,ipos,isite
    integer                             :: i,j,r
    integer                             :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)          :: Indices
    integer,dimension(2*Ns_Ud)          :: Jndices
    integer,dimension(2,Ns_Orb)         :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)                :: Iud
    integer,dimension(2*Ns)             :: ib
    !
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_Sz ERROR: size(V) != sectorI.Dim"
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorI%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default;stop "apply_op_Sz ERROR: called with ed_mode/=superc/nonsu2"
          case ("superc")
             write(LOGfile,"(A,I6,I3)")&
                  'apply Sz :',sectorI%index,sectorI%Sz
          case ("nonsu2")
             if(Jz_basis)then
                write(LOGfile,"(A,I6,I3)")&
                     'apply Sz :',sectorI%index,sectorI%twoJz/2.
             else
                write(LOGfile,"(A,I6,I3)")&
                     'apply Sz :',sectorI%index,sectorI%Ntot
             endif
          end select
       endif
       !
       do i=1,sectorI%Dim
          iph   = (i-1)/(sectorI%DimEl)+1
          i_el  = mod(i-1,sectorI%DimEl)+1
          ei    = sectorI%H(1)%map(i_el)
          ib    = bdecomp(ei,2*Ns)
          sgn   = dble(ib(ipos))-dble(ib(ipos+Ns))
          sgn   = sgn/2d0
          OV(i) = sgn*V(i)
       enddo
       !
       call delete_sector(sectorI)
       !
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
    !
  end function apply_op_Sz_c










  function apply_op_N_d(V,iorb,isector) result(OV)
    real(8),dimension(:),intent(in) :: V
    integer, intent(in)                :: iorb,isector
    type(sector)                  :: sectorI
    real(8),dimension(:),allocatable :: OV
    real(8)                            :: sgn
    integer                            :: ialfa,ibeta,ipos,isite
    integer                            :: i,j,r
    integer                            :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)         :: Indices
    integer,dimension(2*Ns_Ud)         :: Jndices
    integer,dimension(2,Ns_Orb)        :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)               :: Iud
    integer,dimension(2*Ns)            :: ib
    !
    if(MpiMaster)then
       !
       if(ed_total_ud)then
          ialfa = 1
          ipos  = iorb
       else
          ialfa = iorb
          ipos  = 1
       endif
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_N ERROR: size(V) != sectorI.Dim"
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorI%Dim)) ; OV=0d0
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default;stop "apply_op_N ERROR: called with ed_mode/=normal"
          case ("normal")
             write(LOGfile,"(A,I6,2I4)")&
                  'apply N  :',sectorI%index,sectorI%Nups,sectorI%Ndws
          end select
       endif
       !
       do i=1,sectorI%Dim
          iph = (i-1)/(sectorI%DimEl) + 1
          i_el = mod(i-1,sectorI%DimEl) + 1
          call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
          iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
          iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
          nud(1,:) = Bdecomp(iud(1),Ns_Orb)
          nud(2,:) = Bdecomp(iud(2),Ns_Orb)
          sgn      = dble(nud(1,ipos))+dble(nud(2,ipos))
          OV(i)    = sgn*V(i)
       enddo
       !
       call delete_sector(sectorI)       
       !
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=0d0
    end if
    !
  end function apply_op_N_d
  !



  function apply_op_N_c(V,iorb,isector) result(OV)
    complex(8),dimension(:),intent(in) :: V
    integer, intent(in)                :: iorb,isector
    type(sector)                        :: sectorI
    complex(8),dimension(:),allocatable :: OV
    real(8)                            :: sgn
    integer                            :: ialfa,ibeta,ipos,isite
    integer                            :: i,j,r
    integer                            :: iph,i_el,j_el,ei
    integer,dimension(2*Ns_Ud)         :: Indices
    integer,dimension(2*Ns_Ud)         :: Jndices
    integer,dimension(2,Ns_Orb)        :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)               :: Iud
    integer,dimension(2*Ns)            :: ib
    !
    if(MpiMaster)then
       !
       call build_sector(isector,sectorI)
       !
       if(size(V)/=sectorI%Dim)stop "apply_op_N ERROR: size(V) != sectorI.Dim"
       !
       if(allocated(OV))deallocate(OV)
       allocate(OV(sectorI%Dim)) ; OV=zero
       !
       if(ed_verbose>2)then
          select case(ed_mode)
          case default;stop "apply_op_N ERROR: called with ed_mode/=superc/nonsu2"
          case ("superc")
             write(LOGfile,"(A,I6,I3)")&
                  'apply N  :',sectorI%index,sectorI%Sz
          case ("nonsu2")
             if(Jz_basis)then
                write(LOGfile,"(A,I6,I3)")&
                     'apply N  :',sectorI%index,sectorI%twoJz/2.
             else
                write(LOGfile,"(A,I6,I3)")&
                     'apply N  :',sectorI%index,sectorI%Ntot
             endif
          end select
       endif
       !
       do i=1,sectorI%Dim
          iph   = (i-1)/(sectorI%DimEl)+1
          i_el  = mod(i-1,sectorI%DimEl)+1
          ei    = sectorI%H(1)%map(i_el)
          ib    = bdecomp(ei,2*Ns)
          sgn   = dble(ib(ipos))+dble(ib(ipos+Ns))
          OV(i) = sgn*V(i)
       enddo
       !
       call delete_sector(sectorI)
       !
    else
       if(allocated(OV))deallocate(OV)
       allocate(OV(1)) ; OV=zero
    end if
    !
  end function apply_op_N_c
  !





  ! subroutine apply_op_Sz(i,sgn,ipos,ialfa,sectorI) 
  !   integer, intent(in)               :: i,ipos,ialfa
  !   type(sector),intent(in)           :: sectorI
  !   real(8),intent(out)               :: sgn
  !   integer                           :: iph,i_el
  !   integer,dimension(2*Ns_Ud)        :: Indices
  !   integer,dimension(2*Ns_Ud)        :: Jndices
  !   integer,dimension(2,Ns_Orb)       :: Nud !Nbits(Ns_Orb)
  !   integer,dimension(2)              :: Iud
  !   !
  !   sgn=0d0
  !   !
  !   iph = (i-1)/(sectorI%DimEl) + 1
  !   i_el = mod(i-1,sectorI%DimEl) + 1
  !   !
  !   call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
  !   iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
  !   iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
  !   nud(1,:) = Bdecomp(iud(1),Ns_Orb)
  !   nud(2,:) = Bdecomp(iud(2),Ns_Orb)
  !   !
  !   sgn = dble(nud(1,ipos))-dble(nud(2,ipos))
  !   sgn = sgn/2d0
  ! end subroutine apply_op_Sz


  ! subroutine apply_op_N(i,sgn,ipos,ialfa,sectorI) 
  !   integer, intent(in)         :: i,ipos,ialfa
  !   type(sector),intent(in)     :: sectorI
  !   real(8),intent(out)         :: sgn
  !   integer                     :: iph,i_el
  !   integer,dimension(2*Ns_Ud)  :: Indices
  !   integer,dimension(2*Ns_Ud)  :: Jndices
  !   integer,dimension(2,Ns_Orb) :: Nud !Nbits(Ns_Orb)
  !   integer,dimension(2)        :: Iud
  !   !
  !   sgn=0d0
  !   !
  !   iph = (i-1)/(sectorI%DimEl) + 1
  !   i_el = mod(i-1,sectorI%DimEl) + 1
  !   !
  !   call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
  !   iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
  !   iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
  !   nud(1,:) = Bdecomp(iud(1),Ns_Orb)
  !   nud(2,:) = Bdecomp(iud(2),Ns_Orb)
  !   !
  !   sgn = dble(nud(1,ipos))+dble(nud(2,ipos))
  ! end subroutine apply_op_N








  subroutine build_op_Ns(i,Nup,Ndw,sectorI) 
    integer, intent(in)             :: i
    type(sector),intent(in)         :: sectorI
    integer,dimension(Ns)           :: Nup,Ndw  ![Ns]
    integer                         :: iph,i_el,ii,iorb
    integer,dimension(2*Ns_Ud)      :: Indices
    integer,dimension(Ns_Ud,Ns_Orb) :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(2*Ns)         :: Ib
    integer,dimension(2)            :: Iud
    !
    select case(ed_mode)
    case default
       iph = (i-1)/(sectorI%DimEl) + 1
       i_el = mod(i-1,sectorI%DimEl) + 1
       !
       call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
       do ii=1,Ns_Ud
          iud(1) = sectorI%H(ii)%map(Indices(ii))
          iud(2) = sectorI%H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
          Nups(ii,:) = Bdecomp(iud(1),Ns_Orb) ![Norb,1+Nbath]
          Ndws(ii,:) = Bdecomp(iud(2),Ns_Orb)
       enddo
       Nup = Breorder(Nups)
       Ndw = Breorder(Ndws)
       !
    case("superc","nonsu2")
       ii = sectorI%H(1)%map(i)
       Ib = bdecomp(ii,2*Ns)
       Nup = Ib(1:Ns)
       Ndw = Ib(Ns+1:)
       ! do ii=1,Ns
       !    Nup(ii)= ib(ii)
       !    Ndw(ii)= ib(ii+Ns)
       ! enddo
    end select
    !
  end subroutine build_op_Ns







  subroutine get_Sector_normal(QN,N,isector)
    integer,dimension(:) :: QN
    integer              :: N
    integer              :: isector
    integer              :: i,Nind,factor
    Nind = size(QN)
    Factor = N+1
    isector = 1
    do i=Nind,1,-1
       isector = isector + QN(i)*(Factor)**(Nind-i)
    enddo
  end subroutine get_Sector_normal
  !
  subroutine get_Sector_superc(QN,isector)
    integer :: QN
    integer :: isector
    isector=getSector(QN,1)
  end subroutine get_Sector_superc
  !
  subroutine get_Sector_nonsu2(QN,twoJz,isector)
    integer :: QN,twoJz
    integer :: isector
    if(Jz_basis)then
       isector=getSector(QN,twoJz)
    else
       isector=getSector(Qn,1)
    end if
  end subroutine get_Sector_nonsu2


  subroutine get_QuantumNumbers_normal(isector,N,QN)
    integer                          :: isector,N
    integer,dimension(:)             :: QN
    integer                          :: i,count,Dim
    integer,dimension(size(QN)) :: QN_
    !
    Dim = size(QN)
    if(mod(Dim,2)/=0)stop "get_QuantumNumbers error: Dim%2 != 0"
    count=isector-1
    do i=1,Dim
       QN_(i) = mod(count,N+1)
       count      = count/(N+1)
    enddo
    QN = QN_(Dim:1:-1)
  end subroutine get_QuantumNumbers_normal
  subroutine get_QuantumNumbers_other(isector,QN)
    integer                     :: isector
    integer                     :: QN
    select case(ed_mode)
    case ("superc")
       QN=getSz(isector)
    case("nonsu2")
       QN=getN(isector)
    case default
       stop "get_QuantumNumbers_other ERROR: invoked with ed_mode=normal"
    end select
  end subroutine get_QuantumNumbers_other


  subroutine get_Nup(isector,Nup)
    integer                   :: isector,Nup(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud)  :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Nup = indices_(2*Ns_Ud:Ns_Ud+1:-1)
  end subroutine get_Nup
  !
  subroutine get_Ndw(isector,Ndw)
    integer                   :: isector,Ndw(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud) :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Ndw = indices_(Ns_Ud:1:-1)
  end subroutine get_Ndw

  subroutine get_Sz(isector,Sz)
    integer :: isector,Sz
    Sz = getSz(isector)
  end subroutine get_Sz

  subroutine get_Ntot(isector,Ntot)
    integer :: isector,Ntot
    Ntot = getN(isector)
  end subroutine get_Ntot



  subroutine  get_DimUp(isector,DimUps)
    integer                :: isector,DimUps(Ns_Ud)
    integer                :: Nups(Ns_Ud),iud
    call get_Nup(isector,Nups)
    do iud=1,Ns_Ud
       DimUps(iud) = binomial(Ns_Orb,Nups(iud))
    enddo
  end subroutine get_DimUp

  subroutine get_DimDw(isector,DimDws)
    integer                :: isector,DimDws(Ns_Ud)
    integer                :: Ndws(Ns_Ud),iud
    call get_Ndw(isector,Ndws)
    do iud=1,Ns_Ud
       DimDws(iud) = binomial(Ns_Orb,Ndws(iud))
    enddo
  end subroutine get_DimDw

  subroutine  get_Dim(isector,Dim)
    !
    !Returns the dimension of the symmetry sector with index :f:var:`isector`
    !
    integer                :: isector,Dim
    Dim=getDim(isector)
  end subroutine get_Dim


  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index














  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !
  !- build the map from the A-sector to \HHH
  !- get the list of states in \HHH corresponding to sector B twin of A
  !- return the ordering of B-states in \HHH with respect to those of A
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                             :: isector
    integer,dimension(:)                :: order
    type(sector)                        :: sectorH
    type(sector_map),dimension(2*Ns_Ud) :: H
    integer,dimension(2*Ns_Ud)          :: Indices,Istates
    integer                             :: i,iud,iph,i_el
    !
    if(size(Order)/=GetDim(isector))stop "twin_sector_order error: wrong dimensions of *order* array"
    !
    call build_sector(isector,sectorH)
    select case(ed_mode)
    case ("normal")
       do i=1,sectorH%Dim
          iph = (i-1)/(sectorH%DimEl) + 1   !find number of phonons
          i_el = mod(i-1,sectorH%DimEl) + 1 !electronic index
          call state2indices(i_el,[sectorH%DimUps,sectorH%DimDws],Indices)
          forall(iud=1:2*Ns_Ud)Istates(iud) = sectorH%H(iud)%map(Indices(iud))        
          Order(i) = flip_state( Istates ) + (iph-1)*2**(2*Ns)          
       enddo
    case default
       do i=1,sectorH%dim
          iph = (i-1)/(sectorH%DimEl) + 1   !find number of phonons
          i_el = mod(i-1,sectorH%DimEl) + 1 !electronic index
          Order(i) = flip_state(sectorH%H(1)%map(i_el)) + (iph-1)*2**(2*Ns)     
       enddo
    end select
    call delete_sector(sectorH)
    call sort_array(Order) 
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  ! superc: j=|{dw}>|{up}>  , sz  --> -sz
  ! nonsu2: j=|{!up}>|{!dw}>, n   --> 2*Ns-n
  !+------------------------------------------------------------------+
  function flip_state_normal(istate) result(j)
    integer,dimension(2*Ns_Ud) :: istate
    integer                    :: j
    integer,dimension(Ns_Ud)   :: jups,jdws
    integer,dimension(2*Ns_Ud) :: dims
    jups = istate(Ns_Ud+1:2*Ns_Ud)
    jdws = istate(1:Ns_Ud)
    dims = 2**Ns_Orb
    call indices2state([jups,jdws],Dims,j)
  end function flip_state_normal
  function flip_state_other(istate) result(j)
    integer          :: istate
    integer          :: j
    integer          :: ivec(2*Ns),foo(2*Ns)
    !
    Ivec = bdecomp(istate,2*Ns)
    select case(ed_mode)
    case("superc")  !Invert the overall spin sign: |{up}> <---> |{dw}>
       foo(1:Ns)     = Ivec(Ns+1:2*Ns)
       foo(Ns+1:2*Ns)= Ivec(1:Ns)
    case ("nonsu2") !Exchange Occupied sites (1) with Empty sites (0)
       where(Ivec==1)foo=0
       where(Ivec==0)foo=1
    case default    !Exchange UP-config |{up}> with DW-config |{dw}>
       stop "flip_state_other error: called with ed_mode==normal"
    end select
    !
    j = bjoin(foo,2*Ns)
    !
  end function flip_state_other


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in)       :: isector
    integer                  :: jsector    
    integer,dimension(Ns_Ud) :: Iups,Idws
    integer                  :: Sz,Ntot
    select case(ed_mode)
    case default
       call get_Nup(isector,iups)
       call get_Ndw(isector,idws)
       call get_Sector([idws,iups],Ns_Orb,jsector)
    case ("superc")
       call get_Sz(isector,Sz)
       Sz = -Sz
       call get_Sector(Sz,jsector)
    case ("nonsu2")
       call get_Ntot(isector,Ntot)
       Ntot = Nlevels-Ntot
       call get_Sector(Ntot,jsector)
    end select
  end function get_twin_sector













  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



end MODULE ED_SECTOR












