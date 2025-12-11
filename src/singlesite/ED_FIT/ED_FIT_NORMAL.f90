MODULE ED_FIT_NORMAL
!:synopsis: Routines for bath fitting, :code:`NORMAL` case
  USE ED_FIT_COMMON
  USE SF_LINALG, only: diagonal

  implicit none
  private

  public :: chi2_fitgf_normal_normal
  public :: chi2_fitgf_normal_superc
  public :: chi2_fitgf_normal_nonsu2
contains



  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Irreducible bath normal phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_normal_normal(fg,bath_,ispin,iorb)
    complex(8),dimension(:,:,:)                 :: fg ![Norb][Norb][Lmats]
    real(8),dimension(:),intent(inout)          :: bath_
    integer                                     :: ispin
    integer,optional                            :: iorb
    real(8),dimension(:),allocatable            :: array_bath
    integer                                     :: iter,stride,jorb,myorb,i,io,j,Asize
    real(8)                                     :: chi
    logical                                     :: check
    character(len=256)                          :: suffix
    integer                                     :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_normal_normal: Fit"
#endif
    !
    if(size(fg,1)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,1]!=Norb"
    if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,2]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_normal_normal error: wrong bath dimensions"
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
    !
    allocate(Gdelta(1,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    !
    !Asize = get_chi2_bath_size()
    !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
    !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
    Asize = Nbath + Nbath
    allocate(array_bath(Asize))
    !
    do jorb=1,Norb
       if(present(iorb))then
          if(jorb/=iorb)cycle
       endif
       Orb_indx=jorb
       Spin_indx=ispin
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG chi2_fitgf_normal_normal: Fit orb"//str(Orb_indx)//", spin"//str(Spin_indx)
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG chi2_fitgf_normal_normal: cg_method:"//str(cg_method)//&
            ", cg_grad:"//str(cg_grad)//&
            ", cg_scheme:"//str(cg_scheme)
#endif
       !
       Gdelta(1,1:Ldelta) = fg(jorb,jorb,1:Ldelta)
       !
       !Nbath + Nbath
       stride = 0
       do i=1,Nbath
          io = stride + i
          array_bath(io) = dmft_bath%e(ispin,jorb,i)
       enddo
       stride = Nbath
       do i=1,Nbath
          io = stride + i
          array_bath(io) = dmft_bath%v(ispin,jorb,i)
       enddo
       !
       select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
       case default
          if(cg_grad==0)then
             select case (cg_scheme)
             case ("weiss")
                call fmin_cg(array_bath,chi2_weiss_normal_normal,grad_chi2_weiss_normal_normal,iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case ("delta")
                call fmin_cg(array_bath,chi2_delta_normal_normal,grad_chi2_delta_normal_normal,iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case default
                stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
             end select
          else
             select case (cg_scheme)
             case ("weiss")
                call fmin_cg(array_bath,chi2_weiss_normal_normal,iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case ("delta")
                call fmin_cg(array_bath,chi2_delta_normal_normal,iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case default
                stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
             end select
          endif
          !
          !
       case (1)
          if(cg_grad==0)then
             write(*,*) "                                                                                "
             write(*,*) "WARNING: analytic gradient not available with cg-method=1 (minimize f77 routine)"
             write(*,*) "         > we will force cg_grad=1 (so let the routine estimate the gradient)   "
             write(*,*) "                                                                                "
          endif
          select case (cg_scheme)
          case ("weiss")
             call fmin_cgminimize(array_bath,chi2_weiss_normal_normal,&
                  iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                  new_version=cg_minimize_ver,&
                  hh_par=cg_minimize_hh,&
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cgminimize(array_bath,chi2_delta_normal_normal,&
                  iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                  new_version=cg_minimize_ver,&
                  hh_par=cg_minimize_hh,&
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
          end select
          !
       end select
       !
       !
       write(LOGfile,"(A,ES18.9,A,I5,A)")&
            "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
            "  <--  Orb"//reg(txtfy(jorb))//" Spin"//reg(txtfy(ispin))
       !
       suffix="_orb"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       unit=free_unit()
       open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
       write(unit,"(ES18.9,1x,I5)") chi,iter
       close(unit)
       !
       !Nbath + Nbath
       stride = 0
       do i=1,Nbath
          io = stride + i
          dmft_bath%e(ispin,jorb,i) = array_bath(io) 
       enddo
       stride = Nbath
       do i=1,Nbath
          io = stride + i
          dmft_bath%v(ispin,jorb,i) = array_bath(io)
       enddo
       !
    enddo
    !
    call write_dmft_bath(LOGfile)
    !
    call save_dmft_bath()
    !
    allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand = g0and_bath_function(xi*Xdelta(:))
    else
       fgand = delta_bath_function(xi*Xdelta(:))
    endif
    call write_fit_result(ispin)
    deallocate(fgand)
    !
    call get_dmft_bath(bath_)
    call deallocate_dmft_bath()
    deallocate(Gdelta,Xdelta,Wdelta)
    !
  contains
    !
    subroutine write_fit_result(ispin)
      integer           :: i,jorb,ispin
      do jorb=1,Norb
         if(present(iorb))then
            if(jorb/=iorb)cycle
         endif
         suffix="_orb"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
         unit=free_unit()
         if(cg_scheme=='weiss')then
            open(unit,file="fit_weiss"//reg(suffix)//".ed")
         else
            open(unit,file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(jorb,jorb,i)),dimag(fgand(ispin,ispin,jorb,jorb,i)),&
                 dreal(fg(jorb,jorb,i)),dreal(fgand(ispin,ispin,jorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_normal_normal





  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_normal_superc(fg,bath_,ispin,iorb)
    complex(8),dimension(:,:,:,:)      :: fg ![2][Norb][Norb][Lmats]
    real(8),dimension(:),intent(inout) :: bath_
    integer                            :: ispin
    integer,optional                   :: iorb
    real(8),dimension(:),allocatable   :: array_bath
    integer                            :: iter,stride,i,io,j,jorb,Asize
    real(8)                            :: chi
    logical                            :: check
    character(len=256)                 :: suffix
    integer                            :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand,ffand ![Nspin][][Norb][][Ldelta]
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_normal_superc: Fit"
#endif
    if(size(fg,1)/=2)stop "chi2_fitgf_normal_superc error: size[fg,1]!=2"
    if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
    if(size(fg,3)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_normal_superc: wrong bath dimensions"
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,4))Ldelta=size(fg,4)
    !
    allocate(Gdelta(1,Ldelta))
    allocate(Fdelta(1,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*dble(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    !
    !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
    !D_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
    !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
    Asize = Nbath + Nbath + Nbath
    allocate(array_bath(Asize))
    !
    do jorb=1,Norb
       if(present(iorb))then
          if(jorb/=iorb)cycle
       endif
       Orb_indx=jorb
       Spin_indx=ispin
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG chi2_fitgf_normal_superc: Fit orb"//str(Orb_indx)//", spin"//str(Spin_indx)
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG chi2_fitgf_normal_superc: cg_method:"//str(cg_method)//&
            ", cg_grad:"//str(cg_grad)//&
            ", cg_scheme:"//str(cg_scheme)
#endif
       Gdelta(1,1:Ldelta) = fg(1,jorb,jorb,1:Ldelta)
       Fdelta(1,1:Ldelta) = fg(2,jorb,jorb,1:Ldelta)
       !
       !3*Nbath == Nbath + Nbath + Nbath
       stride = 0
       do i=1,Nbath
          io = stride + i
          array_bath(io) = dmft_bath%e(ispin,jorb,i)
       enddo
       stride = Nbath
       do i=1,Nbath
          io = stride + i
          array_bath(io) = dmft_bath%d(ispin,jorb,i)
       enddo
       stride = 2*Nbath
       do i=1,Nbath
          io = stride + i
          array_bath(io) = dmft_bath%v(ispin,jorb,i)
       enddo
       !
       !
       select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
       case default
          if(cg_grad==0)then
             select case (cg_scheme)
             case ("weiss")
                call fmin_cg(array_bath,chi2_weiss_normal_superc,grad_chi2_weiss_normal_superc,&
                                ! call fmin_cg(array_bath,chi2_weiss_normal_superc,&
                     iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case ("delta")
                call fmin_cg(array_bath,chi2_delta_normal_superc,grad_chi2_delta_normal_superc,&
                                ! call fmin_cg(array_bath,chi2_delta_normal_superc,&
                     iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case default
                stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
             end select
          else
             select case (cg_scheme)
             case ("weiss")
                call fmin_cg(array_bath,chi2_weiss_normal_superc,&
                     iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case ("delta")
                call fmin_cg(array_bath,chi2_delta_normal_superc,&
                     iter,chi,&
                     itmax=cg_niter,&
                     ftol=cg_Ftol,  &
                     istop=cg_stop, &
                     iverbose=(ed_verbose>3))
             case default
                stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
             end select
          endif
          !
       case (1)
          select case (cg_scheme)
          case ("weiss")
             call fmin_cgminimize(array_bath,chi2_weiss_normal_superc,&
                  iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                  new_version=cg_minimize_ver,&
                  hh_par=cg_minimize_hh,&
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cgminimize(array_bath,chi2_delta_normal_superc,&                
                  iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                  new_version=cg_minimize_ver,&
                  hh_par=cg_minimize_hh,&
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_normal_superc error: cg_scheme != [weiss,delta]"
          end select
          !
          !
       end select
       !
       write(LOGfile,"(A,ES18.9,A,I5,A)")&
            'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
            "  <--  Orb"//reg(txtfy(jorb))//" Spin"//reg(txtfy(ispin))
       !
       suffix="_orb"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       unit=free_unit()
       open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
       write(unit,"(ES18.9,1x,I5)") chi,iter
       close(unit)
       !
       stride = 0
       do i=1,Nbath
          io = stride + i
          dmft_bath%e(ispin,jorb,i) = array_bath(io)
       enddo
       stride = Nbath
       do i=1,Nbath
          io = stride + i
          dmft_bath%d(ispin,jorb,i) = array_bath(io) 
       enddo
       stride = 2*Nbath
       do i=1,Nbath
          io = stride + i
          dmft_bath%v(ispin,jorb,i) = array_bath(io)
       enddo
       !
    enddo
    call write_dmft_bath(LOGfile)
    !
    call save_dmft_bath()
    !
    allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
    allocate(ffand(Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand = g0and_bath_function(xi*Xdelta(:))
       ffand = f0and_bath_function(xi*Xdelta(:))
    else
       fgand = delta_bath_function(xi*Xdelta(:))
       ffand =fdelta_bath_function(xi*Xdelta(:))
    endif
    call write_fit_result(ispin)
    deallocate(fgand,ffand)
    call get_dmft_bath(bath_)
    call deallocate_dmft_bath()
    deallocate(Gdelta,Fdelta,Xdelta,Wdelta)
    !
  contains
    !
    subroutine write_fit_result(ispin)
      integer           :: jorb,ispin,gunit,funit
      do jorb=1,Norb
         if(present(iorb))then
            if(jorb/=iorb)cycle
         endif
         suffix="_l"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
         if(cg_scheme=='weiss')then
            open(free_unit(gunit),file="fit_weiss"//reg(suffix)//".ed")
            open(free_unit(funit),file="fit_fweiss"//reg(suffix)//".ed")
         else
            open(free_unit(gunit),file="fit_delta"//reg(suffix)//".ed")
            open(free_unit(funit),file="fit_fdelta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(gunit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(1,jorb,jorb,i)),dimag(fgand(ispin,ispin,jorb,jorb,i)),&
                 dreal(fg(1,jorb,jorb,i)),dreal(fgand(ispin,ispin,jorb,jorb,i))
            write(funit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(2,jorb,jorb,i)),dimag(ffand(ispin,ispin,jorb,jorb,i)),&
                 dreal(fg(2,jorb,jorb,i)),dreal(ffand(ispin,ispin,jorb,jorb,i))
         enddo
         close(gunit)
         close(funit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_normal_superc





  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Irreducible bath nonSU2 phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_normal_nonsu2(fg,bath_)
    complex(8),dimension(:,:,:,:,:)             :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
    real(8),dimension(:),intent(inout)          :: bath_
    real(8),dimension(:),allocatable            :: array_bath
    integer                                     :: iter,i,j,io,stride,iorb,ispin,jspin,cspin,Asize
    real(8)                                     :: chi
    logical                                     :: check
    character(len=20)                           :: suffix
    integer                                     :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_normal_nonsu2: Fit"
#endif
    if(size(fg,1)/=Nspin)stop "chi2_fitgf_normal_nonsu2 error: size[fg,1]!=Nspin"
    if(size(fg,2)/=Nspin)stop "chi2_fitgf_normal_nonsu2 error: size[fg,2]!=Nspin"
    if(size(fg,3)/=Norb)stop "chi2_fitgf_normal_nonsu2 error: size[fg,3]!=Norb"
    if(size(fg,4)/=Norb)stop "chi2_fitgf_normal_nonsu2 error: size[fg,4]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_normal_nonsu2 error: wrong bath dimensions"
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
    !
    totNspin = Nspin*(Nspin+1)/2
    allocate(getIspin(totNspin),getJspin(totNspin))
    cspin=0
    do ispin=1,Nspin
       do jspin=ispin,Nspin
          cspin=cspin+1
          getIspin(cspin)=ispin
          getJspin(cspin)=jspin
       enddo
    enddo
    if(cspin/=totNspin)stop "chi2_fitgf_normal_nonsu2: error counting the spins"
    !
    allocate(Gdelta(totNspin,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    !
    !Dimension for a given orbital
    Asize = Nbath + Nbath + Nbath
    Asize = Nspin*Asize
    ! if(.not.para_)Asize=Nspin*Asize !fit all spin components
    allocate(array_bath(Asize))
    do iorb=1,Norb
       Orb_indx=iorb
       !
#ifdef _DEBUG
       if(ed_verbose>3)write(Logfile,"(A)")"DEBUG chi2_fitgf_normal_nonsu2: Fit orb"//str(Orb_indx)
       if(ed_verbose>3)write(Logfile,"(A)")&
            "DEBUG chi2_fitgf_normal_nonsu2: cg_method:"//str(cg_method)//&
            ", cg_grad:"//str(cg_grad)//&
            ", cg_scheme:"//str(cg_scheme)
#endif
       !
       do i=1,totNspin
          Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),iorb,iorb,1:Ldelta)
       enddo
       !
       stride = 0
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             array_bath(io) = dmft_bath%e(ispin,iorb,i)
          enddo
       enddo
       stride = Nspin*Nbath
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             array_bath(io) = dmft_bath%v(ispin,iorb,i)
          enddo
       enddo
       stride = Nspin*Nbath + Nspin*Nbath
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             array_bath(io) = dmft_bath%u(ispin,iorb,i)
          enddo
       enddo
       !
       select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
       case default
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_normal_nonsu2,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_normal_nonsu2,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_normal_nonsu2 error: cg_scheme != [weiss,delta]"
          end select
          !
       case (1)
          select case (cg_scheme)
          case ("weiss")
             call fmin_cgminimize(array_bath,chi2_weiss_normal_nonsu2,&
                  iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                  new_version=cg_minimize_ver,&
                  hh_par=cg_minimize_hh,&
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cgminimize(array_bath,chi2_delta_normal_nonsu2,&
                  iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                  new_version=cg_minimize_ver,&
                  hh_par=cg_minimize_hh,&
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_normal_nonsu2 error: cg_scheme != [weiss,delta]"
          end select
       end select
       !
       !
       write(LOGfile,"(A,ES18.9,A,I5,A)")&
            "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
            "  <--  Orb"//reg(txtfy(iorb))//" All spins"
       !
       suffix="_orb"//reg(txtfy(iorb))//"_ALLspins_"//reg(ed_file_suffix)
       unit=free_unit()
       open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
       write(unit,"(ES18.9,1x,I5)") chi,iter
       close(unit)
       !
       stride = 0
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             dmft_bath%e(ispin,iorb,i) = array_bath(io)
          enddo
       enddo
       stride = Nspin*Nbath
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             dmft_bath%v(ispin,iorb,i)  = array_bath(io)
          enddo
       enddo
       stride = Nspin*Nbath + Nspin*Nbath
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             dmft_bath%u(ispin,iorb,i)  = array_bath(io)
          enddo
       enddo
       ! dmft_bath%e(Nspin,iorb,:) = dmft_bath%e(1,iorb,:)
       ! dmft_bath%v(Nspin,iorb,:) = dmft_bath%v(1,iorb,:)
       ! dmft_bath%u(Nspin,iorb,:) = dmft_bath%u(1,iorb,:)
       !
    enddo
    !
    call write_dmft_bath(LOGfile)
    !
    call save_dmft_bath()
    allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand = g0and_bath_function(xi*Xdelta(:))
    else
       fgand = delta_bath_function(xi*Xdelta(:))
    endif
    call write_fit_result()
    deallocate(fgand)
    !
    call get_dmft_bath(bath_)
    call deallocate_dmft_bath()
    deallocate(Gdelta,Xdelta,Wdelta)
    deallocate(getIspin,getJspin)
    !
  contains
    !
    subroutine write_fit_result()
      integer           :: i,j,s,iorb,ispin,jspin
      !
      do iorb=1,Norb
         do s=1,totNspin
            ispin = getIspin(s)
            jspin = getJspin(s)
            suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
            if(cg_scheme=='weiss')then
               open(free_unit(unit),file="fit_weiss"//reg(suffix)//".ed")
            else
               open(free_unit(unit),file="fit_delta"//reg(suffix)//".ed")
            endif
            !
            do i=1,Ldelta
               write(unit,"(5F24.15)")Xdelta(i),&
                    dimag(fg(ispin,jspin,iorb,iorb,i)),dimag(fgand(ispin,jspin,iorb,iorb,i)),&
                    dreal(fg(ispin,jspin,iorb,iorb,i)),dreal(fgand(ispin,jspin,iorb,iorb,i))
            enddo
            close(unit)
         enddo
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_normal_nonsu2














  !##################################################################
  ! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
  !##################################################################
  !> NORMAL
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function chi2_delta_normal_normal(a) result(chi2)
    real(8),dimension(:)         ::  a
    real(8)                      ::  chi2
    complex(8),dimension(Ldelta) ::  Delta
    real(8),dimension(Ldelta)    ::  Ctmp
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_normal. a:",a
#endif
    Delta = delta_normal_normal(a)
    !
    Ctmp = abs(Gdelta(1,:)-Delta(:))
    chi2=sum( Ctmp**cg_pow/Wdelta )
    chi2=chi2/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_normal_normal. Chi**2:",chi2
#endif
    !
  end function chi2_delta_normal_normal

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_delta_normal_normal(a) result(dchi2)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2
    real(8),dimension(size(a))           :: df
    complex(8),dimension(Ldelta)         :: Delta
    complex(8),dimension(Ldelta)         :: Ftmp
    real(8),dimension(Ldelta)            :: Ctmp
    complex(8),dimension(Ldelta,size(a)) :: dDelta
    integer                              :: j
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_normal. a:",a
#endif
    Delta   = delta_normal_normal(a)
    dDelta  = grad_delta_normal_normal(a)
    !
    Ftmp = Gdelta(1,:)-Delta(:)
    Ctmp = abs(Ftmp)**(cg_pow-2)
    do j=1,size(a)
       df(j)=sum( dreal(Ftmp)*dreal(dDelta(:,j))*Ctmp/Wdelta ) + &
            sum(  dimag(Ftmp)*dimag(dDelta(:,j))*Ctmp/Wdelta )
    enddo
    !
    dchi2 = -cg_pow*df/Ldelta
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_normal. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_delta_normal_normal

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  !+-------------------------------------------------------------+
  function chi2_weiss_normal_normal(a) result(chi2)
    real(8),dimension(:)         ::  a
    complex(8),dimension(Ldelta) ::  g0and
    real(8),dimension(Ldelta)    ::  Ctmp
    real(8)                      ::  chi2,w
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_normal_normal. a:",a
#endif
    !
    g0and  = g0and_normal_normal(a)
    !
    Ctmp = abs(Gdelta(1,:)-g0and(:))
    chi2 = sum( Ctmp**cg_pow/Wdelta )
    chi2 = chi2/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_normal_normal. Chi**2:",chi2
#endif
    !
  end function chi2_weiss_normal_normal

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_weiss_normal_normal(a) result(dchi2)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2
    real(8),dimension(size(a))           :: df
    complex(8),dimension(Ldelta)         :: g0and,Ftmp
    real(8),dimension(Ldelta)            :: Ctmp
    complex(8),dimension(Ldelta,size(a)) :: dg0and
    integer                              :: j
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_normal_normal. a:",a
#endif
    !
    g0and  = g0and_normal_normal(a)
    dg0and = grad_g0and_normal_normal(a)
    !
    Ftmp = Gdelta(1,:)-g0and(:)
    Ctmp = abs(Ftmp)**(cg_pow-2)
    do j=1,size(a)
       df(j)=sum( dreal(Ftmp)*dreal(dg0and(:,j))*Ctmp/Wdelta ) + &
            sum(  dimag(Ftmp)*dimag(dg0and(:,j))*Ctmp/Wdelta )
    enddo
    !
    dchi2 = -cg_pow*df/Ldelta
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_normal_normal. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_weiss_normal_normal





  !##################################################################
  ! THESE PROCEDURES EVALUATES THE 
  ! - \delta
  ! - \grad \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################

  function delta_normal_normal(a) result(Delta)
    real(8),dimension(:)         :: a
    complex(8),dimension(Ldelta) :: Delta
    integer                      :: i,io,stride
    real(8),dimension(Nbath)     :: eps,vps
    !
    !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       eps(i) = a(io) 
    enddo
    stride = Nbath
    do i=1,Nbath
       io = stride + i
       vps(i) = a(io)
    enddo
    !
    do i=1,Ldelta
       Delta(i) = sum( vps(:)*vps(:)/(xi*Xdelta(i) - eps(:)) )
    enddo
    !
  end function delta_normal_normal

  function grad_delta_normal_normal(a) result(dDelta)
    real(8),dimension(:)                 :: a
    complex(8),dimension(Ldelta,size(a)) :: dDelta
    integer                              :: i,k,ik,io,stride
    real(8),dimension(Nbath)             :: eps,vps
    complex(8)                           :: iw
    !
    !
    !\grad_{E_{a}(k)} \Delta_{bb}^{rr} = [ V_{a}(k)*V_{a}(k) / ( iw_n - E_{a}(k) )**2 ]
    !
    !\grad_{V_{a}(k)} \Delta_{bb}^{rr} = [ 2*V_{a}(k) / ( iw_n - E_{a}(k) ) ]
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       eps(i) = a(io) 
    enddo
    stride = Nbath
    do i=1,Nbath
       io = stride + i
       vps(i) = a(io)
    enddo
    !
    stride = 0
    do k=1,Nbath
       ik = stride + k
       dDelta(:,ik) = vps(k)*vps(k)/(xi*Xdelta(:) - eps(k))**2
    enddo
    stride = Nbath
    do k=1,Nbath
       ik = stride + k
       dDelta(:,ik) = 2d0*vps(k)/(xi*Xdelta(:) - eps(k))
    enddo
    !
  end function grad_delta_normal_normal


  function g0and_normal_normal(a) result(G0and)
    real(8),dimension(:)         :: a
    complex(8),dimension(Ldelta) :: G0and,Delta
    integer                      :: i,io,iorb,ispin
    !
    iorb   = Orb_indx
    ispin  = Spin_indx
    !
    Delta(:) = delta_normal_normal(a)
    G0and(:) = xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(:)
    G0and(:) = one/G0and(:)
    !
  end function g0and_normal_normal

  function grad_g0and_normal_normal(a) result(dG0and)
    real(8),dimension(:)                 :: a
    complex(8),dimension(Ldelta)         :: G0and,Delta
    complex(8),dimension(Ldelta,size(a)) :: dG0and,dDelta
    integer                              :: k,iorb,ispin
    !
    iorb   = Orb_indx
    ispin  = Spin_indx
    !
    Delta  = delta_normal_normal(a)
    dDelta = grad_delta_normal_normal(a)
    G0and  = xi*Xdelta + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta
    do k=1,size(a)
       dG0and(:,k) = one/G0and/G0and*dDelta(:,k)
    enddo
    !
  end function grad_g0and_normal_normal












  !>SUPERC
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
  !         in the SUPERCONDUCTING case.
  !+-------------------------------------------------------------+
  function chi2_delta_normal_superc(a) result(chi2)
    real(8),dimension(:)           ::  a
    real(8)                        ::  chi2
    complex(8),dimension(2,Ldelta) ::  Delta
    real(8),dimension(Ldelta)      ::  Ctmp,Ftmp
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_superc. a:",a
#endif
    !
    Delta(:,:) = delta_normal_superc(a)
    !
    Ctmp = abs(Gdelta(1,:)-Delta(1,:))
    Ftmp = abs(Fdelta(1,:)-Delta(2,:))
    chi2 = sum( Ctmp**cg_pow/Wdelta ) + sum( Ftmp**cg_pow/Wdelta )
    chi2 = chi2/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_normal_superc. Chi**2:",chi2
#endif
    !
  end function chi2_delta_normal_superc

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson 
  ! function in the SUPERCONDUCTING case.
  !+-------------------------------------------------------------+
  function grad_chi2_delta_normal_superc(a) result(dchi2)
    real(8),dimension(:)                   ::  a
    real(8),dimension(size(a))             ::  dchi2
    real(8),dimension(size(a))             ::  df
    complex(8),dimension(2,Ldelta)         ::  Delta
    complex(8),dimension(Ldelta)           ::  Gtmp,Ftmp
    real(8),dimension(Ldelta)              ::  Ctmp,Btmp
    complex(8),dimension(2,Ldelta,size(a)) ::  dDelta
    integer                                ::  j
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_superc. a:",a
#endif
    !
    Delta(:,:)    = delta_normal_superc(a)
    dDelta(:,:,:) = grad_delta_normal_superc(a)
    !
    Gtmp = Gdelta(1,:)-Delta(1,:)
    Ftmp = Fdelta(1,:)-Delta(2,:)
    !
    Ctmp = abs(Gtmp)**(cg_pow-2)
    Btmp = abs(Ftmp)**(cg_pow-2)
    !
    do j=1,size(a)
       df(j) = &
            sum( dreal(Gtmp)*dreal(dDelta(1,:,j))*Ctmp/Wdelta ) + &
            sum( dimag(Gtmp)*dimag(dDelta(1,:,j))*Ctmp/Wdelta ) + &
            sum( dreal(Ftmp)*dreal(dDelta(2,:,j))*Btmp/Wdelta ) + &
            sum( dimag(Ftmp)*dimag(dDelta(2,:,j))*Btmp/Wdelta )
    enddo
    !
    dchi2 = -cg_pow*df/Ldelta
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_superc. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_delta_normal_superc

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance for G_0 function 
  !         in the SUPERCONDUCTING case.
  !+-------------------------------------------------------------+
  function chi2_weiss_normal_superc(a) result(chi2)
    real(8),dimension(:)           ::  a
    complex(8),dimension(2,Ldelta) ::  g0and
    real(8),dimension(Ldelta)      ::  Gtmp,Ftmp
    real(8)                        ::  chi2
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_normal_superc. a:",a
#endif
    !
    g0and(:,:)  = g0and_normal_superc(a)
    !
    Gtmp = abs(Gdelta(1,:)-g0and(1,:))
    Ftmp = abs(Fdelta(1,:)-g0and(2,:))
    chi2 =        sum( Gtmp**cg_pow/Wdelta )
    chi2 = chi2 + sum( Ftmp**cg_pow/Wdelta )
    chi2 = chi2/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_normal_superc. chi**2:",chi2
#endif
    !
  end function chi2_weiss_normal_superc

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_weiss_normal_superc(a) result(dchi2)
    real(8),dimension(:)                   :: a
    real(8),dimension(size(a))             :: dchi2
    real(8),dimension(size(a))             :: df
    complex(8),dimension(2,Ldelta)         :: g0and
    complex(8),dimension(2,Ldelta,size(a)) :: dg0and
    complex(8),dimension(Ldelta)           :: Gtmp,Ftmp
    real(8),dimension(Ldelta)              :: Ctmp,Btmp
    integer                                :: j
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_normal_superc. a:",a
#endif
    !
    g0and  = g0and_normal_superc(a)
    dg0and = grad_g0and_normal_superc(a)
    !
    Gtmp = (Gdelta(1,:)-g0and(1,:))
    Ftmp = (Fdelta(1,:)-g0and(2,:))
    !
    Ctmp = abs(Gtmp)**(cg_pow-2)
    Btmp = abs(Ftmp)**(cg_pow-2)
    !
    do j=1,size(a)
       df(j) = &
            sum( dreal(Gtmp)*dreal(dg0and(1,:,j))*Ctmp/Wdelta ) + &
            sum( dimag(Gtmp)*dimag(dg0and(1,:,j))*Ctmp/Wdelta ) + &
            sum( dreal(Ftmp)*dreal(dg0and(2,:,j))*Btmp/Wdelta ) + &
            sum( dimag(Ftmp)*dimag(dg0and(2,:,j))*Btmp/Wdelta )
    enddo
    !
    dchi2 = -cg_pow*df/Ldelta
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_normal_superc. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_weiss_normal_superc







  !##################################################################
  ! THESE PROCEDURES EVALUATES THE 
  ! - \delta
  ! - \grad \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################
  function delta_normal_superc(a) result(Delta)
    real(8),dimension(:)            :: a
    complex(8),dimension(2,Ldelta)  :: Delta
    integer                         :: i,k,io,stride
    real(8),dimension(Nbath)        :: eps,vps,dps
    real(8),dimension(Nbath)        :: Den
    !
    !\Delta_{aa} = - \sum_k [ V_{a}(k) * V_{a}(k) * (iw_n + E_{a}(k)) / Den(k) ]
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       eps(i) = a(io)
    enddo
    stride = Nbath
    do i=1,Nbath
       io = stride + i
       dps(i) = a(io) 
    enddo
    stride = 2*Nbath
    do i=1,Nbath
       io = stride + i
       vps(i) = a(io)
    enddo
    !
    do i=1,Ldelta
       Delta(1,i) = -sum( vps(:)*vps(:)*( xi*Xdelta(i) + eps(:) )/(Xdelta(i)**2 + eps(:)**2 + dps(:)**2) )
       Delta(2,i) =  sum( dps(:)*vps(:)*vps(:)/(Xdelta(i)**2 + eps(:)**2 + dps(:)**2) )
    enddo
    !
  end function delta_normal_superc

  function grad_delta_normal_superc(a) result(dDelta)
    real(8),dimension(:)                   :: a
    complex(8),dimension(2,Ldelta,size(a)) :: dDelta
    integer                                :: i,k,ik,io,stride
    real(8),dimension(Nbath)               :: eps,vps,dps
    real(8),dimension(Ldelta,Nbath)        :: Den
    !
    !\grad_{E_{a}(k)} \Delta_{bb} = -V_{a}(k)*V_{a}(k)*[ 1/den(k) - 2*E_{a}(k)*(iw_n + E_{a}(k))/den(k)**2 ]
    !
    !\grad_{\D_{a}(k)} \Delta_{bb} = V_{a}(k)*V_{a}(k)*\D_{a}(k)*(iw_n + E_{a}(k)) /den(k)**2
    !
    !\grad_{ V_{a}(k)} \Delta_{bb} = -2*V_{a}(k)*(iw_n + E_{a}(k))/den(k)
    !
    !
    !
    !\grad_{E_{a}(k)} \FDelta_{aa} = -2 * V_{a}(k) * V_{a}(k) * E_{a}(k) * \Delta_{a}(k) / Den**2
    !
    !\grad_{\Delta_{a}(k)} \FDelta_{aa} = V_{a}(k) * V_{a}(k) * [ 1/den - 2* \Delta_{a}(k)*\Delta_{a}(k)/den**2 ]
    !
    !\grad_{ V_{a}(k)} \FDelta_{aa} =  2 * V_{a}(k) * \Delta_{a}(k) / den
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       eps(i) = a(io)
    enddo
    stride = Nbath
    do i=1,Nbath
       io = stride + i
       dps(i) = a(io) 
    enddo
    stride = 2*Nbath
    do i=1,Nbath
       io = stride + i
       vps(i) = a(io)
    enddo
    !
    !Den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
    forall(i=1:Ldelta,k=1:Nbath)Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
    !
    stride = 0
    do k=1,Nbath
       ik = stride + k
       dDelta(1,:,ik) = -vps(k)*vps(k)*(1d0/Den(:,k) - 2d0*eps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2)
    enddo
    stride = Nbath
    do k=1,Nbath
       ik = stride + k
       dDelta(1,:,ik) = 2d0*vps(k)*vps(k)*dps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2
    enddo
    stride = 2*Nbath
    do k=1,Nbath
       ik = stride + k
       dDelta(1,:,ik) = -2d0*vps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)
    enddo
    !
    !
    stride = 0
    do k=1,Nbath
       ik = stride + k
       dDelta(2,:,ik) = -2d0*vps(k)*vps(k)*eps(k)*dps(k)/Den(:,k)**2
    enddo
    stride = Nbath
    do k=1,Nbath
       ik = stride + k
       dDelta(2,:,ik) = vps(k)*vps(k)*(1d0/Den(:,k) - 2d0*dps(k)*dps(k)/Den(:,k)**2) 
    enddo
    stride = 2*Nbath
    do k=1,Nbath
       ik = stride + k
       dDelta(2,:,ik) = 2d0*vps(k)*dps(k)/Den(:,k)
    enddo
    !
  end function grad_delta_normal_superc

  function g0and_normal_superc(a) result(G0and)
    real(8),dimension(:)            :: a
    complex(8),dimension(2,Ldelta)  :: G0and,Delta
    real(8),dimension(Ldelta)       :: det
    complex(8),dimension(Ldelta)    :: fg,ff
    integer                         :: iorb,ispin
    !
    iorb   = Orb_indx
    ispin  = Spin_indx
    !
    Delta    = delta_normal_superc(a)
    !
    fg(:)    = xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) -  Delta(1,:)
    ff(:)    = -diagonal(impHloc_anomalous(1,1,:,:))               -  Delta(2,:)
    det(:)   = abs(fg(:))**2 + ff(:)**2
    G0and(1,:) = conjg(fg(:))/det(:)
    G0and(2,:) = ff(:)/det(:)
    !
  end function g0and_normal_superc

  function grad_g0and_normal_superc(a) result(dG0and)
    real(8),dimension(:)                   :: a
    complex(8),dimension(2,Ldelta)         :: G0and,Delta
    complex(8),dimension(2,Ldelta,size(a)) :: dG0and,dDelta
    integer                                :: i,k,ik,io,stride
    real(8),dimension(Nbath)               :: eps,vps,dps
    integer                                :: iorb,ispin
    real(8),dimension(Ldelta,Nbath)        :: Den
    complex(8),dimension(Ldelta)           :: g0,f0,dD,dC,dDet,zeta
    real(8),dimension(Ldelta)              :: det
    !
    iorb   = Orb_indx
    ispin  = Spin_indx
    !
    Delta  = delta_normal_superc(a)
    dDelta = grad_delta_normal_superc(a)
    !
    zeta= xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) 
    g0  = zeta - Delta(1,:)
    f0  = -pair_field(1:Norb) - Delta(2,:)
    Det = abs(g0)**2 + f0**2
    !
    do k=1,size(a)
       dD = conjg(dDelta(1,:,k))
       dC = dDelta(2,:,k)
       dDet = 2*dreal(g0*dD)+2*f0*dC
       dG0and(1,:,k) = -Det*dD + conjg(g0)*dDet
       dG0and(2,:,k) = -Det*dC + f0*dDet
       dG0and(1,:,k) = dG0and(1,:,k)/Det**2
       dG0and(2,:,k) = dG0and(2,:,k)/Det**2
    enddo
  end function grad_g0and_normal_superc










  !>NONSU2
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function chi2_delta_normal_nonsu2(a) result(chi2)
    real(8),dimension(:)                     ::  a
    real(8)                                  ::  chi2
    real(8),dimension(totNspin)              ::  chi2_spin
    complex(8),dimension(Nspin,Nspin,Ldelta) ::  Delta
    real(8),dimension(Ldelta)                ::  Ctmp
    integer                                  ::  s,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_normal_nonsu2. a:",a
#endif
    Delta = delta_normal_nonsu2(a)
    !
    do s=1,totNspin
       ispin=getIspin(s)
       jspin=getJspin(s)
       Ctmp = abs(Gdelta(s,:)-Delta(ispin,jspin,:))
       chi2_spin(s) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_spin)/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_normal_nonsu2. Chi**2:",chi2
#endif
    !
  end function chi2_delta_normal_nonsu2


  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  !+-------------------------------------------------------------+
  function chi2_weiss_normal_nonsu2(a) result(chi2)
    real(8),dimension(:)                     :: a
    real(8),dimension(totNspin)              :: chi2_spin
    complex(8),dimension(Nspin,Nspin,Ldelta) :: g0and
    real(8),dimension(Ldelta)                ::  Ctmp
    real(8)                                  :: chi2
    real(8)                                  :: w
    integer                                  :: i,s,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_normal_nonsu2. a:",a
#endif
    !
    g0and(:,:,:) = g0and_normal_nonsu2(a)
    !
    do s=1,totNspin
       ispin=getIspin(s)
       jspin=getJspin(s)
       Ctmp = abs(Gdelta(s,:)-g0and(ispin,jspin,:))
       chi2_spin(s) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_spin)/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_normal_nonsu2. chi**2:",chi2
#endif
    !
  end function chi2_weiss_normal_nonsu2






  !##################################################################
  ! THESE PROCEDURES EVALUATES THE 
  ! - \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################
  function delta_normal_nonsu2(a) result(Delta)
    real(8),dimension(:)                     :: a
    complex(8),dimension(Nspin,Nspin,Ldelta) :: Delta
    integer                                  :: ispin,jspin
    integer                                  :: i,io,stride,ih
    real(8),dimension(Nspin,Nbath)           :: hps
    real(8),dimension(Nspin,Nbath)           :: vps
    real(8),dimension(Nspin,Nbath)           :: ups
    real(8),dimension(Nhel,Nbath)            :: ehel
    real(8),dimension(Nhel,Nhel,Nbath)       :: whel

    stride = 0
    do ispin=1,Nspin
       do i=1,Nbath
          io = stride + i + (ispin-1)*Nbath
          hps(ispin,i) = a(io)
       enddo
    enddo
    stride = Nspin*Nbath
    do ispin=1,Nspin
       do i=1,Nbath
          io = stride + i + (ispin-1)*Nbath
          vps(ispin,i)  = a(io)
       enddo
    enddo
    stride = Nspin*Nbath + Nspin*Nbath
    do ispin=1,Nspin
       do i=1,Nbath
          io = stride + i + (ispin-1)*Nbath
          ups(ispin,i)  = a(io)
       enddo
    enddo
    ! hps(Nspin,:)=hps(ispin,:)
    ! vps(Nspin,:)=vps(ispin,:)
    ! ups(Nspin,:)=ups(ispin,:)
    !
    whel = get_Whyb_matrix(vps(1:Nspin,1:Nbath),ups(1:Nspin,1:Nbath))
    !
    Delta=zero
    do i=1,Ldelta
       do ispin=1,Nspin
          do jspin=1,Nspin
             do ih=1,Nspin
                Delta(ispin,jspin,i) = Delta(ispin,jspin,i) + sum( whel(ispin,ih,:)*whel(jspin,ih,:)/(xi*Xdelta(i) - hps(ih,:)) )
             enddo
          enddo
       enddo
    enddo
    !
  end function delta_normal_nonsu2

  function g0and_normal_nonsu2(a) result(G0and)
    integer                                  :: iorb,ispin,jspin
    real(8),dimension(:)                     :: a
    complex(8),dimension(Nspin,Nspin,Ldelta) :: G0and,Delta
    complex(8),dimension(Nspin,Nspin)        :: zeta,fgorb
    integer                                  :: i
    !
    iorb  = Orb_indx
    !
    Delta(:,:,:) = delta_normal_nonsu2(a)
    !
    do i=1,Ldelta
       zeta  = (xi*Xdelta(i)+xmu)*zeye(Nspin)
       fgorb = zeta - impHloc(:,:,iorb,iorb) - Delta(:,:,i)
       call inv(fgorb)
       G0and(:,:,i) = fgorb
    enddo
    !
  end function g0and_normal_nonsu2




END MODULE ED_FIT_NORMAL
