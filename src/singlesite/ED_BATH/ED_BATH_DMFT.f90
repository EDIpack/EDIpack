MODULE ED_BATH_DMFT
  !:synopsis: Routines to allocate, read and write the bath
  !A class for the the :f:var:`effective_bath` data structure describing the effective bath in the code. 
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH_AUX
  USE ED_BATH_REPLICA
  USE ED_BATH_DIM 
  implicit none



contains



  subroutine allocate_dmft_bath()
    !
    ! Allocate the  :f:var:`effective_bath` input  :f:var:`dmft_bath_` according to the values of the global input parameters |Nspin| , |Norb|, |Nbath|, |ed_mode| and |bath_type|.
    !
    ! The correct components of the :f:var:`effective_bath` data structure are allocated.   
    !
    ! .. list-table:: allocated components of :f:var:`dmft_bath_`
    !    :widths: auto
    !    :header-rows: 1
    !    :stub-columns: 1
    !
    !    * - 
    !      - :code:`normal`
    !      - :code:`superc`
    !      - :code:`nonsu2`
    !
    !    * - :code:`normal`
    !      - :f:var:`e` , :f:var:`v`            
    !      - :f:var:`e` , :f:var:`v`, :f:var:`d`
    !      - :f:var:`e` , :f:var:`v`, :f:var:`u`
    !
    !    * - :code:`hybrid`
    !      - :f:var:`e` , :f:var:`v`            
    !      - :f:var:`e` , :f:var:`v`, :f:var:`d`
    !      - :f:var:`e` , :f:var:`v`, :f:var:`u`
    !
    !    * - :code:`replica`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda`
    !
    !    * - :code:`general`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda` , :f:var:`item` % :f:var:`vg`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda` , :f:var:`item` % :f:var:`vg`
    !      - :f:var:`nbasis` , :f:var:`item`  , :f:var:`item` % :f:var:`lambda` , :f:var:`item` % :f:var:`vg`
    !
    !
    integer              :: Nsym,ibath
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG allocate_dmft_bath"
#endif
    !
    !We don't check for Nbath>0 here: Fortran lets us allocate 0-dimensional arrays
    !but not access them. This avoids us filling the H part with ifs
    !
    if(dmft_bath%status)call deallocate_dmft_bath()
    !
    select case(bath_type)
    case default ;STOP "allocate_dmft_bath error: bath_type "//trim(bath_type)//" not valid"
    case('normal')
       !
       select case(ed_mode)
       case default                                 !normal [N,Sz]
          allocate(dmft_bath%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("superc")                              !superc [Sz]
          allocate(dmft_bath%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath%d(Nspin,Norb,Nbath))  !local SC order parameters the bath
          allocate(dmft_bath%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("nonsu2")                              !nonsu2 [N]
          allocate(dmft_bath%e(Nspin,Norb,Nbath))  !local energies of the bath
          allocate(dmft_bath%v(Nspin,Norb,Nbath))  !same-spin hybridization
          allocate(dmft_bath%u(Nspin,Norb,Nbath))  !spin-flip hybridization
       end select
       !
    case('hybrid')
       !
       select case(ed_mode)
       case default                                 !normal  [N,Sz]
          allocate(dmft_bath%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("superc")                              !superc  [Sz]
          allocate(dmft_bath%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath%d(Nspin,1,Nbath))     !local SC order parameters the bath
          allocate(dmft_bath%v(Nspin,Norb,Nbath))  !same-spin hybridization
       case ("nonsu2")                              !nonsu2 case [N] qn
          allocate(dmft_bath%e(Nspin,1,Nbath))     !local energies of the bath
          allocate(dmft_bath%v(Nspin,Norb,Nbath))  !same-spin hybridization
          allocate(dmft_bath%u(Nspin,Norb,Nbath))  !spin-flip hybridization
       end select
       !
    case('replica')
       !
       if(.not.Hb%status)stop "ERROR allocate_dmft_bath: Hb.basis not allocated"
       call deallocate_dmft_bath()     !
       Nsym=Hb%Nsym
       !
       allocate(dmft_bath%item(Nbath))
       dmft_bath%Nbasis=Nsym
       do ibath=1,Nbath
          allocate(dmft_bath%item(ibath)%lambda(Nsym))
       enddo
       !
    case('general')
       !
       if(.not.Hb%status)stop "ERROR allocate_dmft_bath: Hgeneral_basis not allocated"
       call deallocate_dmft_bath()     !
       Nsym=Hb%Nsym
       !
       allocate(dmft_bath%item(Nbath))
       dmft_bath%Nbasis=Nsym
       do ibath=1,Nbath
          allocate(dmft_bath%item(ibath)%lambda(Nsym))
          allocate(dmft_bath%item(ibath)%vg(Norb*Nspin))
       enddo
       !
    end select
    !
    dmft_bath%status=.true.
    !
  end subroutine allocate_dmft_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_dmft_bath()
    !
    ! Deallocate the  :f:var:`effective_bath` input  :f:var:`dmft_bath` 
    !
    integer              :: ibath,isym
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG deallocate_dmft_bath"
#endif
    !
    !We don't check for Nbath>0 here: Fortran lets us allocate 0-dimensional arrays
    !but not access them. This avoids us filling the H part with ifs
    !
    if(.not.dmft_bath%status)return
    if(allocated(dmft_bath%e))   deallocate(dmft_bath%e)
    if(allocated(dmft_bath%d))   deallocate(dmft_bath%d)
    if(allocated(dmft_bath%v))   deallocate(dmft_bath%v)
    if(allocated(dmft_bath%u))   deallocate(dmft_bath%u)
    if(bath_type=="replica")then
       dmft_bath%Nbasis= 0
       do ibath=1,Nbath
          deallocate(dmft_bath%item(ibath)%lambda)
       enddo
       deallocate(dmft_bath%item)
    endif
    if(bath_type=="general")then
       dmft_bath%Nbasis= 0
       do ibath=1,Nbath
          deallocate(dmft_bath%item(ibath)%lambda)
          deallocate(dmft_bath%item(ibath)%vg)
       enddo
       deallocate(dmft_bath%item)
    endif
    dmft_bath%status=.false.
  end subroutine deallocate_dmft_bath




  subroutine init_dmft_bath()
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    !
    ! Initialize the :f:var:`effective_bath` input :f:var:`dmft_bath`.
    ! The energy parameters :f:var:`e` are constructed using a centered :f:var:`nbath` discretization of the flat band of width 2 :f:var:`ed_hw_bath`.
    ! The hybridization parameters :f:var:`v` are set to :math:`\tfrac{1}{\sqrt{nbath}`.
    ! The superconducting parameters :f:var:`d`, if any, are set to :f:var:`deltasc`
    ! The spin-flip parameters :f:var:`u`, if any, are set equal to :f:var:`v`.
    ! The replica/general local bath Hamiltonians are built out of the initial values of :math:`\vec{\lambda}`, using the matrix decomposition.
    !
    !.. note::
    !
    !   If the file :f:var:`hfile` with the correct :f:var:`ed_file_suffix` and extension :f:var:`.restart` is found in the running directory, the bath is initilized by reading the from the file.
    !   If the file :f:var:`bfile` with the correct :f:var:`ed_file_suffix` and extension :f:var:`.restart` is found in the running directory, the replica/general bath matrix basi is re-initilized by reading the from the file.
    !
    !
    integer              :: Nbasis
    integer              :: i,ibath,isym,unit,flen,Nh,Nsym
    integer              :: io,jo,iorb,ispin,jorb,jspin
    logical              :: IOfile,diagonal_hsym,equal_lambdas
    real(8)              :: de
    real(8)              :: offset(Nbath)
    !
#ifdef _DEBUG
    write(Logfile,"(A)")"DEBUG init_dmft_bath"
#endif
    !
    if(Nbath.eq.0)return    
    !
    if(.not.dmft_bath%status)stop "ERROR init_dmft_bath error: bath not allocated"
    !
    select case(bath_type)
    case default
       !Get energies:
       dmft_bath%e(:,:,1)    =-ed_hw_bath
       dmft_bath%e(:,:,Nbath)= ed_hw_bath
       Nh=Nbath/2
       if(mod(Nbath,2)==0.and.Nbath>=4)then
          de=ed_hw_bath/max(Nh-1,1)
          dmft_bath%e(:,:,Nh)  = -1.d-1
          dmft_bath%e(:,:,Nh+1)=  1.d-1
          do i=2,Nh-1
             dmft_bath%e(:,:,i)   =-ed_hw_bath + (i-1)*de
             dmft_bath%e(:,:,Nbath-i+1)= ed_hw_bath - (i-1)*de
          enddo
       elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
          de=ed_hw_bath/Nh
          dmft_bath%e(:,:,Nh+1)= 0d0
          do i=2,Nh
             dmft_bath%e(:,:,i)        =-ed_hw_bath + (i-1)*de
             dmft_bath%e(:,:,Nbath-i+1)= ed_hw_bath - (i-1)*de
          enddo
       endif
       !Get spin-keep yhbridizations
       do i=1,Nbath
          dmft_bath%v(:,:,i)=max(0.1d0,1d0/sqrt(dble(Nbath)))
       enddo
       !Get SC amplitudes
       if(ed_mode=="superc")dmft_bath%d(:,:,:)=deltasc
       !Get spin-flip hybridizations
       if(ed_mode=="nonsu2")then
          do i=1,Nbath
             dmft_bath%u(:,:,i) = dmft_bath%v(:,:,i)  ! x ed_vsf_ratio
          enddo
       endif
       !
    case('replica','general')
       offset=0.d0
       if(Nbath>1) offset=linspace(-ed_offset_bath,ed_offset_bath,Nbath)
       !
       !BATH V INITIALIZATION
       do ibath=1,Nbath
          if(bath_type=='replica')then
             dmft_bath%item(ibath)%v    =max(0.1d0,1d0/sqrt(dble(Nbath)))
          elseif(bath_type=='general')then
             dmft_bath%item(ibath)%vg(:)=max(0.1d0,1d0/sqrt(dble(Nbath)))
          endif
       enddo
       !
       !BATH LAMBDAS INITIALIZATION
       !Do not need to check for Hreplica_basis: this is done at allocation time of the dmft_bath.
       Nsym = dmft_bath%Nbasis
       do isym=1,Nsym
          do ibath=1,Nbath
             dmft_bath%item(ibath)%lambda(isym) =  Hb%linit(ibath,isym)
          enddo
          !
          diagonal_hsym = is_diagonal(Hb%basis(isym)%O)
          equal_lambdas = all(Hb%linit(:,isym)==Hb%linit(Nbath,isym))
          !
          if(diagonal_hsym.AND.equal_lambdas.AND.Nbath>1)then
             offset=linspace(-ed_offset_bath,ed_offset_bath,Nbath)
             if(is_identity(Hb%basis(isym)%O).AND.mod(Nbath,2)==0)then
                offset(Nbath/2)     = max(-1.d-1,offset(Nbath/2))
                offset(Nbath/2 + 1) = min(1.d-1,offset(Nbath/2 + 1))
             endif
             do ibath=1,Nbath
                dmft_bath%item(ibath)%lambda(isym) =  Hb%linit(ibath,isym) + offset(ibath)
             enddo
             if(MpiMaster)then
               write(LOGfile,*) "                                                                    "
               write(LOGfile,*) "WARNING: some of your inital lambda have been internally changed    "
               write(LOGfile,*) "         while calling ed_init_solver. This happens whenever the    "
               write(LOGfile,*) "         corresponding Hsym is diagonal and all the replicas receive"
               write(LOGfile,*) "         the same initial lambda value, due to the deprecated legacy"
               write(LOGfile,*) "         practice of defining a unique lambda vector forall replicas"
               write(LOGfile,*) "         and let the solver decide how to handle these degeneracies."
               write(LOGfile,*) "         >>> If you really intend to have a degenerate diagonal term"
               write(LOGfile,*) "             in the bath you can define a suitable restart file.    "
               write(LOGfile,*) "         >>> If instead this is what you expected please consider to"
               write(LOGfile,*) "             move the desired rescaling in your driver, since this  "
               write(LOGfile,*) "             funcionality might be removed in a future update.      "
               write(LOGfile,*) "                                                                    "
             endif
          endif
       enddo
       !
    end select
    !
    !
    !
    !Read from file if exist:
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       call read_dmft_bath(used=.false.)
    endif
    
  end subroutine init_dmft_bath









  subroutine read_dmft_bath(used)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    !
    ! Read the :f:var:`effective_bath` from a file associated to the [optional] unit.
    !
    logical,optional              :: used
    logical                       :: used_    
    character(len=120)            :: file,bathfile
    logical                       :: IOfile
    integer                       :: i,io,jo,iorb,ispin,flen,unit
    character(len=20)             :: hsuffix
    character(len=64)             :: string_fmt,string_fmt_first
    character(len=1)              :: a
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG read_dmft_bath"
#endif
    !
    if(Nbath.eq.0)return      
    !
    used_   = .true.      ;if(present(used))used_=used
    hsuffix = ".restart"  ;if(used_)hsuffix=reg(".used") !default=used
    file    = str(Hfile)//str(ed_file_suffix)//str(hsuffix)
    bathfile   = str(Bfile)//str(ed_file_suffix)//str(hsuffix)
    !
    inquire(file=str(file),exist=IOfile)
    if(.not.IOfile)stop "read_dmft_bath ERROR: indicated file does not exist"
    if(ed_verbose>2)write(LOGfile,"(A)")'Reading bath from file '//str(file)
    !
    if(dmft_bath%status)call deallocate_dmft_bath()
    call allocate_dmft_bath()
    !
    flen = file_length(file,verbose=ed_verbose>2)
    if(flen==0)stop "Can not read file "//str(file)
    !
    open(free_unit(unit),file=str(file))
    !
    !Take care of the comment in the preamble:
    read(unit,*)
    !
    !BATH TYPE=normal,hybrid,replica,general
    select case(bath_type)
    case default
       !ED MODE=normal,superc,nonsu2
       select case(ed_mode)
       case default
          do i=1,min(flen,Nbath)
             read(unit,*)((&
                  dmft_bath%e(ispin,iorb,i),&
                  dmft_bath%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       case ("superc")
          do i=1,min(flen,Nbath)
             read(unit,*)((&
                  dmft_bath%e(ispin,iorb,i),&
                  dmft_bath%d(ispin,iorb,i),&
                  dmft_bath%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       case("nonsu2")
          do i=1,min(flen,Nbath)
             read(unit,*)((&
                  dmft_bath%e(ispin,iorb,i),&
                  dmft_bath%v(ispin,iorb,i),&
                  dmft_bath%u(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       end select
       !
    case('hybrid')
       !ED MODE=normal,superc,nonsu2
       select case(ed_mode)
       case default
          do i=1,min(flen,Nbath)
             read(unit,*)(&
                  dmft_bath%e(ispin,1,i),&
                  (&
                  dmft_bath%v(ispin,iorb,i),&
                  iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       case ("superc")
          do i=1,min(flen,Nbath)
             read(unit,*)(&
                  dmft_bath%e(ispin,1,i),&
                  dmft_bath%d(ispin,1,i),&
                  (&
                  dmft_bath%v(ispin,iorb,i),&
                  iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       case ("nonsu2")
          do i=1,min(flen,Nbath)
             read(unit,*)(&
                  dmft_bath%e(ispin,1,i),&
                  (&
                  dmft_bath%v(ispin,iorb,i),&
                  dmft_bath%u(ispin,iorb,i),&
                  iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       end select
       !
    case ('replica')
       read(unit,*)dmft_bath%Nbasis
       do i=1,Nbath
          read(unit,*)dmft_bath%item(i)%v,&
               (dmft_bath%item(i)%lambda(io),io=1,dmft_bath%Nbasis)
       enddo
       inquire(file=str(bathfile),exist=IOfile)
       if(IOfile)call read_Hreplica(str(bathfile))
       if(allocated(Hb%linit))then
          do i=1,Nbath
             Hb%linit(i,:) = dmft_bath%item(i)%lambda(:)
          enddo
       endif
       !
    case ('general')
       read(unit,*)dmft_bath%Nbasis
       do i=1,Nbath
          read(unit,*)dmft_bath%item(i)%vg(:),&
               (dmft_bath%item(i)%lambda(io),io=1,dmft_bath%Nbasis)
       enddo
       inquire(file=str(bathfile),exist=IOfile)
       if(IOfile)call read_Hgeneral(str(bathfile))
       if(allocated(Hb%linit))then
          do i=1,Nbath
             Hb%linit(i,:) = dmft_bath%item(i)%lambda(:)
          enddo
       endif
       !
    end select
    !
    close(unit)
    !
  end subroutine read_dmft_bath







  subroutine save_dmft_bath(used)
    !
    ! Save the :f:var:`effective_bath` to a file with .used or .restart extension according to input.
    !
    character(len=256)        :: file,bathfile
    logical,optional          :: used
    logical                   :: used_
    character(len=16)         :: extension
    integer                   :: unit_
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG save_dmft_bath"
#endif
    !
    if(Nbath.eq.0)return    
    !
    if(.not.dmft_bath%status)stop "save_dmft_bath error: bath is not allocated"
    !
    used_    =.false.    ;if(present(used))used_=used
    extension=".restart" ;if(used_)extension=".used"
    !
    file     =str(Hfile)//str(ed_file_suffix)//str(extension)
    bathfile =str(Bfile)//str(ed_file_suffix)//str(extension)
    !
    unit_=free_unit()
    if(MpiMaster)then
       open(unit_,file=str(file))
       call write_dmft_bath(unit_)
       close(unit_)
       select case (bath_type)
       case default;
       case ("replica","general")
          call save_Hreplica(str(bathfile))
       end select
    endif
  end subroutine save_dmft_bath


  subroutine write_dmft_bath(unit)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    !
    ! Write the :f:var:`effective_bath` on input to std.output or to a file associated to the [optional] unit.
    !
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i,Nsym
    integer              :: io,jo,iorb,ispin,isym
    real(8)              :: hybr_aux
    complex(8)           :: ho(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb)
    character(len=64)    :: string_fmt,string_fmt_first
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG write_dmft_bath"
#endif
    !
    if(Nbath.eq.0)return
    !
    unit_=LOGfile;if(present(unit))unit_=unit
    !
    if(.not.dmft_bath%status)stop "write_dmft_bath error: bath not allocated"
    !
    !BATH TYPE=normal,hybrid,replica,general
    select case(bath_type)
    case default
       !ED MODE=normal,superc,nonsu2
       select case(ed_mode)
       case default
          write(unit_,"(*(A21,1X))")&
               ((&
               "#Ek_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(*(ES21.12,1X))")((&
                  dmft_bath%e(ispin,iorb,i),&
                  dmft_bath%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       case ("superc")
          write(unit_,"(*(A21,1X))")&
               ((&
               "#Ek_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Dk_l"//reg(str(iorb))//"_s"//reg(str(ispin)) ,&
               "Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(*(ES21.12,1X))")((&
                  dmft_bath%e(ispin,iorb,i),&
                  dmft_bath%d(ispin,iorb,i),&
                  dmft_bath%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       case ("nonsu2")
          write(unit_,"(*(A21,1X))")&
               ((&
               "#Ek_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vak_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vbk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb), ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(*(ES21.12,1X))")((&
                  dmft_bath%e(ispin,iorb,i),&
                  dmft_bath%v(ispin,iorb,i),&
                  dmft_bath%u(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       end select
       !
    case('hybrid')
       !ED MODE=normal,superc,nonsu2
       select case(ed_mode)
       case default
          write(unit_,"(*(A21,1X))")(&
               "#Ek_s"//reg(str(ispin)),&
               ("Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),iorb=1,Norb),&
               ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(*(ES21.12,1X))")(&
                  dmft_bath%e(ispin,1,i),&
                  (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       case ("superc")
          write(unit_,"(*(A21,1X))")(&
               "#Ek_s"//reg(str(ispin)),&
               "Dk_s"//reg(str(ispin)) ,&
               ("Vk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),iorb=1,Norb),&
               ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(*(ES21.12,1X))")(&
                  dmft_bath%e(ispin,1,i),&
                  dmft_bath%d(ispin,1,i),&
                  (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       case ("nonsu2")
          write(unit_,"(*(A21,1X))")(&
               "#Ek_s"//reg(str(ispin)),&
               (&
               "Vak_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               "Vbk_l"//reg(str(iorb))//"_s"//reg(str(ispin)),&
               iorb=1,Norb),&
               ispin=1,Nspin)
          do i=1,Nbath
             write(unit_,"(*(ES21.12,1X))")(&
                  dmft_bath%e(ispin,1,i),    &
                  (dmft_bath%v(ispin,iorb,i),dmft_bath%u(ispin,iorb,i),iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
       end select
       !
    case ('replica')
       !
       write(unit_,"(*(A21,1X))")"#V_i",("Lambda_i"//reg(str(io)),io=1,dmft_bath%Nbasis)
       write(unit_,"(I3)")dmft_bath%Nbasis
       do i=1,Nbath
          write(unit_,"(*(ES21.12,1X))")dmft_bath%item(i)%v,&
               (dmft_bath%item(i)%lambda(io),io=1,dmft_bath%Nbasis)
       enddo
       !
       if(unit_/=LOGfile)then
          string_fmt = "("//str(Nnambu*Nspin*Norb)//"(A1,F5.2,A1,F5.2,A1,2x))"
          write(unit_,*)""
          do isym=1,Hb%Nsym
             Ho = nn2so_reshape( Hb%basis(isym)%O, Nnambu*Nspin,Norb)
             do io=1,Nnambu*Nspin*Norb
                write(unit_,string_fmt)&
                     ('(',dreal(Ho(io,jo)),',',dimag(Ho(io,jo)),')',jo =1,Nnambu*Nspin*Norb)
             enddo
             write(unit_,*)""
          enddo
       endif
    case ('general')
       !
       write(unit_,"(A1,*(A21,1X))")"#",("V_i"//reg(str(io)),io=1,Nspin*Norb),("Lambda_i"//reg(str(io)),io=1,dmft_bath%Nbasis)
       write(unit_,"(I3)")dmft_bath%Nbasis
       do i=1,Nbath
          write(unit_,"(*(ES21.12,1X))")(dmft_bath%item(i)%vg(io),io=1,Nspin*Norb),&
               (dmft_bath%item(i)%lambda(io),io=1,dmft_bath%Nbasis) 
       enddo
       !
       if(unit_/=LOGfile)then
          string_fmt = "("//str(Nnambu*Nspin*Norb)//"(A1,F5.2,A1,F5.2,A1,2x))"
          write(unit_,*)""
          do isym=1,Hb%Nsym
             Ho = nn2so_reshape( Hb%basis(isym)%O, Nnambu*Nspin,Norb)
             do io=1,Nnambu*Nspin*Norb
                write(unit_,string_fmt)&
                     ('(',dreal(Ho(io,jo)),',',dimag(Ho(io,jo)),')',jo =1,Nnambu*Nspin*Norb)
             enddo
             write(unit_,*)""
          enddo
       endif
    end select
  end subroutine write_dmft_bath




















  subroutine set_dmft_bath(bath_)
    !
    ! Set the :f:var:`effective_bath` components from the input user bath :f:var:`bath_` , i.e. it dumps the user bath to the internal data structure. 
    !
    real(8),dimension(:)   :: bath_
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath
    logical                :: check, diag_condition
    !
    diag_condition = check_diag_impHloc()
    if (bath_type=="normal" .and. .not. diag_condition) STOP "Normal bath cannot be used if impHloc(_anomalous) are off-diagonal"
    !
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG set_dmft_bath: dmft_bath <- user_bath"
#endif
    !
    if(Nbath.eq.0)return
    !
    if(.not.dmft_bath%status)stop "set_dmft_bath error: bath not allocated"
    !
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
          !
       case default
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%d(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%e(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   dmft_bath%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
       !
    case ('hybrid')
       !
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath%d(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = 2*Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                dmft_bath%e(ispin,1,i) = bath_(io)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath%v(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   dmft_bath%u(ispin,iorb,i) = bath_(io)
                enddo
             enddo
          enddo
       end select
       !
       !
    case ('replica')
       !
       stride = 1
       !Get Nbasis
       dmft_bath%Nbasis = NINT(bath_(stride))
       !get Lambdas
       do ibath=1,Nbath
          stride = stride + 1
          dmft_bath%item(ibath)%v = bath_(stride)
          dmft_bath%item(ibath)%lambda=bath_(stride+1 :stride+dmft_bath%Nbasis)
          stride=stride+dmft_bath%Nbasis
       enddo
    case ('general')
       !
       stride = 1
       !Get Nbasis
       dmft_bath%Nbasis = NINT(bath_(stride))
       !get Lambdas
       do ibath=1,Nbath
          dmft_bath%item(ibath)%vg(:) = bath_(stride+1:stride+Nspin*Norb)
          stride = stride + Nspin*Norb
          dmft_bath%item(ibath)%lambda=bath_(stride+1 :stride+dmft_bath%Nbasis)
          stride=stride+dmft_bath%Nbasis
       enddo
    end select
  end subroutine set_dmft_bath



  subroutine get_dmft_bath(bath_)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    !
    ! Set the user input bath :f:var:`bath_` from the components of the :f:var:`effective_bath` :f:var:`dmft_bath` , i.e. it dumps the internal data structure to the user bath. 
    !
    real(8),dimension(:)   :: bath_
    real(8)                :: hrep_aux(Nspin*Norb,Nspin*Norb)
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath,maxspin
    logical                :: check
#ifdef _DEBUG
    if(ed_verbose>1)write(Logfile,"(A)")"DEBUG get_dmft_bath: dmft_bath -> user_bath"
#endif
    !
    if(Nbath.eq.0)return  
    !
    if(.not.dmft_bath%status)stop "get_dmft_bath error: bath not allocated"
    !
    check=check_bath_dimension(bath_)
    if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
    !
    bath_ = 0d0
    !
    select case(bath_type)
    case default
       !
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath%d(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) =  dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath%e(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = 2*Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                   bath_(io) = dmft_bath%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
    case ('hybrid')
       !
       select case(ed_mode)
       case default
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) =  dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case ("superc")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) =  dmft_bath%d(ispin,1,i)
             enddo
          enddo
          stride = 2*Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) =  dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          !
       case("nonsu2")
          stride = 0
          do ispin=1,Nspin
             do i=1,Nbath
                io = stride + i + (ispin-1)*Nbath
                bath_(io) = dmft_bath%e(ispin,1,i)
             enddo
          enddo
          stride = Nspin*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath%v(ispin,iorb,i)
                enddo
             enddo
          enddo
          stride = Nspin*Nbath + Nspin*Norb*Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                do i=1,Nbath
                   io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                   bath_(io) = dmft_bath%u(ispin,iorb,i)
                enddo
             enddo
          enddo
       end select
       !
       !
    case ('replica')
       !
       stride = 1
       bath_(stride)=dmft_bath%Nbasis
       do ibath=1,Nbath
          stride = stride + 1
          bath_(stride)=dmft_bath%item(ibath)%v
          bath_(stride+1 : stride+dmft_bath%Nbasis)=dmft_bath%item(ibath)%lambda
          stride=stride+dmft_bath%Nbasis
       enddo
    case ('general')
       !
       stride = 1
       bath_(stride)=dmft_bath%Nbasis
       do ibath=1,Nbath
          bath_(stride+1:stride+Nspin*Norb)=dmft_bath%item(ibath)%vg(:)
          stride = stride + Nspin*Norb
          bath_(stride+1 : stride+dmft_bath%Nbasis)=dmft_bath%item(ibath)%lambda
          stride=stride+dmft_bath%Nbasis
       enddo
    end select
  end subroutine get_dmft_bath



END MODULE ED_BATH_DMFT









