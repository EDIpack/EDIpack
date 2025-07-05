MODULE ED_FIT_HYBRID
!:synopsis: Routines for bath fitting, :code:`HYBRID` case
  USE ED_FIT_COMMON

  implicit none
  private


  public :: chi2_fitgf_hybrid_normal
  public :: chi2_fitgf_hybrid_superc
  public :: chi2_fitgf_hybrid_nonsu2

  
contains

  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Hybrid bath normal phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_hybrid_normal(fg,bath_,ispin)
    complex(8),dimension(:,:,:)                 :: fg ![Norb][Norb][Lmats]
    real(8),dimension(:),intent(inout)          :: bath_
    integer                                     :: ispin
    real(8),dimension(:),allocatable            :: array_bath
    integer                                     :: iter,stride,i,io,j,corb,l,Asize
    integer                                     :: iorb,jorb
    real(8)                                     :: chi
    logical                                     :: check
    character(len=256)                          :: suffix
    integer                                     :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_hybrid_normal: Fit"
#endif
    !
    if(size(fg,1)/=Norb)stop "chi2_fitgf_hybrid_normal error: size[fg,1]!=Norb"
    if(size(fg,2)/=Norb)stop "chi2_fitgf_hybrid_normal error: size[fg,2]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_hybrid_normal error: wrong bath dimensions"
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
    !
    allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
    corb=0
    do iorb=1,Norb
       do jorb=iorb,Norb
          corb=corb+1
          getIorb(corb)=iorb
          getJorb(corb)=jorb
       enddo
    enddo
    totNorb=corb
    if(totNorb/=(Norb*(Norb+1)/2))stop "chi2_fitgf_hybrid_normal error counting the orbitals"
    !
    allocate(Gdelta(totNorb,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    do i=1,totNorb
       Gdelta(i,1:Ldelta) = fg(getIorb(i),getJorb(i),1:Ldelta)
    enddo
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
    Spin_indx=ispin
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    !
    !E_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
    !V_{\s,:}(:)  [ 1 ][ Norb][Nbath]
    Asize = Nbath + Norb*Nbath
    allocate(array_bath(Asize))
    !
    !Nbath + Norb*Nbath
    stride = 0
    do i=1,Nbath
       io = stride + i
       array_bath(io) = dmft_bath%e(ispin,1,i)
    enddo
    stride = Nbath
    do iorb=1,Norb
       do i=1,Nbath
          io = stride + i + (iorb-1)*Nbath
          array_bath(io) = dmft_bath%v(ispin,iorb,i)
       enddo
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_normal_normal: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)
#endif
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
    case default
       if(cg_grad==0)then
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_hybrid_normal,grad_chi2_weiss_hybrid_normal,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_hybrid_normal,grad_chi2_delta_hybrid_normal,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
          end select
       else
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_hybrid_normal,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol, &
                  istop=cg_stop,&
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_hybrid_normal,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol, &
                  istop=cg_stop,&
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
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
          call fmin_cgminimize(array_bath,chi2_weiss_hybrid_normal,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_hybrid_normal,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5)")&
         'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
         "  <--  all Orbs, Spin"//reg(txtfy(ispin))
    !
    suffix="_ALLorb_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       dmft_bath%e(ispin,1,i) = array_bath(io)
    enddo
    stride = Nbath
    do iorb=1,Norb
       do i=1,Nbath
          io = stride + i + (iorb-1)*Nbath
          dmft_bath%v(ispin,iorb,i) = array_bath(io)
       enddo
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
    deallocate(getIorb,getJorb)
    !
  contains
    !
    subroutine write_fit_result(ispin)
      integer :: i,j,l,m,iorb,jorb,ispin,jspin
      !
      do l=1,totNorb
         iorb=getIorb(l)
         jorb=getJorb(l)
         suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
         unit=free_unit()
         if(cg_scheme=='weiss')then
            open(unit,file="fit_weiss"//reg(suffix)//".ed")
         else
            open(unit,file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(Gdelta(l,i)),dimag(fgand(ispin,ispin,iorb,jorb,i)),&
                 dreal(Gdelta(l,i)),dreal(fgand(ispin,ispin,iorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_hybrid_normal









  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_hybrid_superc(fg,bath_,ispin)
    complex(8),dimension(:,:,:,:)               :: fg ![2][Norb][Norb][Lmats]
    real(8),dimension(:),intent(inout)          :: bath_
    integer                                     :: ispin
    real(8),dimension(:),allocatable            :: array_bath
    integer                                     :: iter,stride,i,j,io,corb,l,Asize
    integer                                     :: iorb,jorb
    real(8)                                     :: chi
    logical                                     :: check
    character(len=256)                          :: suffix
    integer                                     :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand,ffand ![Nspin][][Norb][][Ldelta]
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_hybrid_superc: Fit"
#endif
    !
    if(size(fg,1)/=2)stop "chi2_fitgf_normal_superc error: size[fg,1]!=2"
    if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,2]!=Norb"
    if(size(fg,3)/=Norb)stop "chi2_fitgf_normal_superc error: size[fg,3]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_hybrid_superc: wrong bath dimensions"
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,4))Ldelta=size(fg,4)
    !
    allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
    corb=0
    do iorb=1,Norb
       do jorb=iorb,Norb
          corb=corb+1
          getIorb(corb)=iorb
          getJorb(corb)=jorb
       enddo
    enddo
    totNorb=corb
    if(totNorb/=(Norb*(Norb+1)/2))stop "chi2_fitgf_hybrid_superc error counting the orbitals"
    !
    allocate(Gdelta(totNorb,Ldelta))
    allocate(Fdelta(totNorb,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    do i=1,totNorb
       Gdelta(i,1:Ldelta) = fg(1,getIorb(i),getJorb(i),1:Ldelta)
       Fdelta(i,1:Ldelta) = fg(2,getIorb(i),getJorb(i),1:Ldelta)
    enddo
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
    Spin_indx=ispin 
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    !
    !E_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
    !D_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
    !V_{\s,:}(:)  [ 1 ][ Norb][Nbath]
    Asize = Nbath + Nbath + Norb*Nbath
    allocate(array_bath(Asize))
    !
    !Nbath + Nbath + Norb*Nbath
    stride = 0
    do i=1,Nbath 
       io = stride + i
       array_bath(io)   = dmft_bath%e(ispin,1,i)
    enddo
    stride = Nbath
    do i=1,Nbath 
       io = stride + i
       array_bath(io)    = dmft_bath%d(ispin,1,i)
    enddo
    stride = Nbath + Nbath
    do iorb=1,Norb
       do i=1,Nbath
          io = stride + i + (iorb-1)*Nbath
          array_bath(io) = dmft_bath%v(ispin,iorb,i)
       enddo
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_normal_superc: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)
#endif
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default

       if(cg_grad==0)then
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_hybrid_superc,&!grad_chi2_weiss_hybrid_superc,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_hybrid_superc,grad_chi2_delta_hybrid_superc,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
          end select
       else
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_hybrid_superc,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_hybrid_superc,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
          end select
       endif
       !
    case (1)
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_hybrid_superc,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_hybrid_superc,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_hybrid_superc error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5)")&
         'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
         "  <--  all Orbs, Spin"//reg(txtfy(ispin))
    !
    suffix="_ALLorb_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    stride = 0
    do i=1,Nbath 
       io = stride + i
       dmft_bath%e(ispin,1,i) = array_bath(io)   
    enddo
    stride = Nbath
    do i=1,Nbath 
       io = stride + i
       dmft_bath%d(ispin,1,i)   = array_bath(io)
    enddo
    stride = 2*Nbath
    do iorb=1,Norb
       do i=1,Nbath
          io = stride + i + (iorb-1)*Nbath
          dmft_bath%v(ispin,iorb,i) = array_bath(io)
       enddo
    enddo
    !
    call write_dmft_bath(LOGfile)
    !
    call save_dmft_bath()

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
    deallocate(getIorb,getJorb)
    !
  contains
    !
    subroutine write_fit_result(ispin)
      integer                                :: i,j,l,m,iorb,jorb,ispin,jspin
      integer                                :: gunit,funit
      real(8)                                :: w
      !
      do l=1,totNorb
         iorb=getIorb(l)
         jorb=getJorb(l)
         suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
         if(cg_scheme=='weiss')then
            open(free_unit(gunit),file="fit_weiss"//reg(suffix)//".ed")
            open(free_unit(funit),file="fit_fweiss"//reg(suffix)//".ed")
         else
            open(free_unit(gunit),file="fit_delta"//reg(suffix)//".ed")
            open(free_unit(funit),file="fit_fdelta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(gunit,"(10F24.15)")Xdelta(i),&
                 dimag(Gdelta(l,i)),dimag(fgand(ispin,ispin,iorb,jorb,i)),&
                 dreal(Gdelta(l,i)),dreal(fgand(ispin,ispin,iorb,jorb,i))
            write(funit,"(10F24.15)")Xdelta(i),&
                 dimag(Fdelta(l,i)),dimag(ffand(ispin,ispin,iorb,jorb,i)),&
                 dreal(Fdelta(l,i)),dreal(ffand(ispin,ispin,iorb,jorb,i))
         enddo
         close(gunit)
         close(funit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_hybrid_superc





  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Hybridized bath nonsu2 phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_hybrid_nonsu2(fg,bath_)
    complex(8),dimension(:,:,:,:,:)    :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
    real(8),dimension(:),intent(inout) :: bath_
    real(8),dimension(:),allocatable   :: array_bath
    integer                            :: i,j,iorb,jorb,ispin,jspin,io,jo
    integer                            :: iter,stride,count,Asize
    real(8)                            :: chi
    logical                            :: check
    character(len=256)                 :: suffix
    integer                            :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_hybrid_nonsu2: Fit"
#endif
    !
    if(size(fg,1)/=Nspin)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,1]!=Nspin"
    if(size(fg,2)/=Nspin)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,2]!=Nspin"
    if(size(fg,3)/=Norb)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,3]!=Norb"
    if(size(fg,4)/=Norb)stop "chi2_fitgf_hybrid_nonsu2 error: size[fg,4]!=Norb"
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_hybrid_nonsu2 error: wrong bath dimensions"
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
    !
    totNso = Nspin*(Nspin+1)/2 * Norb*(Norb+1)/2
    allocate(getIspin(totNso),getJspin(totNso))
    allocate(getIorb(totNso),getJorb(totNso))
    count=0
    do ispin=1,Nspin
       do jspin=ispin,Nspin
          do iorb=1,Norb
             do jorb=iorb,Norb
                count=count+1
                getIspin(count) = ispin
                getIorb(count)  = iorb
                getJspin(count) = jspin
                getJorb(count)  = jorb
             enddo
          enddo
       enddo
    enddo
    if(totNso/=count)stop "chi2_fitgf_hybrid_nonsu2: error counting the spin-orbitals"
    !
    allocate(Gdelta(totNso,Ldelta))
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
    do i=1,totNso
       Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
    enddo
    !
    Asize = Nbath + Norb*Nbath + Norb*Nbath
    Asize = Nspin*Asize
    ! if(.not.para_)Asize=Nspin*Asize !fit all spin components
    allocate(array_bath(Asize))
    !
    ! size = Nspin*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
    stride = 0
    do ispin=1,Nspin
       do i=1,Nbath
          io = stride + i + (ispin-1)*Nbath
          array_bath(io) = dmft_bath%e(ispin,1,i)
       enddo
    enddo
    stride = Nspin*Nbath
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nbath
             io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
             array_bath(io) = dmft_bath%v(ispin,iorb,i)
          enddo
       enddo
    enddo
    stride =  Nspin*Nbath + Nspin*Norb*Nbath
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nbath
             io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
             array_bath(io) = dmft_bath%u(ispin,iorb,i)
          enddo
       enddo
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_normal_nonsu2: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)
#endif
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
    case default
       select case (cg_scheme)
       case ("weiss")
          call fmin_cg(array_bath,chi2_weiss_hybrid_nonsu2,&
               iter,chi,&
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cg(array_bath,chi2_delta_hybrid_nonsu2,&
               iter,chi,&
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_hybrid_nonsu2 error: cg_scheme != [weiss,delta]"
       end select
       !
    case (1)
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_hybrid_nonsu2,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_hybrid_nonsu2,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,&
               hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_hybrid_nonsu2 error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")&
         "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
         "  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    ! size = Nspin*Nbath + Nspin*Norb*Nbath + Nspin*Norb*Nbath
    stride = 0
    do ispin=1,Nspin
       do i=1,Nbath
          io = stride + i + (ispin-1)*Nbath
          dmft_bath%e(ispin,1,i) = array_bath(io)
       enddo
    enddo
    stride = Nspin*Nbath
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nbath
             io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
             dmft_bath%v(ispin,iorb,i) = array_bath(io)
          enddo
       enddo
    enddo
    stride =  Nspin*Nbath + Nspin*Norb*Nbath
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nbath
             io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
             dmft_bath%u(ispin,iorb,i) = array_bath(io)
          enddo
       enddo
    enddo
    ! dmft_bath%e(Nspin,1,:) = dmft_bath%e(1,1,:)
    ! dmft_bath%v(Nspin,:,:) = dmft_bath%v(1,:,:)
    ! dmft_bath%u(Nspin,:,:) = dmft_bath%u(1,:,:)
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
    call write_fit_result()
    deallocate(fgand)
    !
    call get_dmft_bath(bath_)
    call deallocate_dmft_bath()
    deallocate(Gdelta,Xdelta,Wdelta)
    deallocate(getIspin,getJspin)
    deallocate(getIorb,getJorb)
    !
  contains
    !
    subroutine write_fit_result()
      integer           :: i,j,s,l,iorb,jorb,ispin,jspin
      !
      do l=1,totNso
         iorb = getIorb(l)
         jorb = getJorb(l)
         ispin = getIspin(l)
         jspin = getJspin(l)
         suffix="_l"//reg(txtfy(iorb))//&
              "_m"//reg(txtfy(jorb))//&
              "_s"//reg(txtfy(ispin))//&
              "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
         if(cg_scheme=='weiss')then
            open(free_unit(unit),file="fit_weiss"//reg(suffix)//".ed")
         else
            open(free_unit(unit),file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
                 dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_hybrid_nonsu2





  !##################################################################
  ! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
  !##################################################################
  !>NORMAL
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function chi2_delta_hybrid_normal(a) result(chi2)
    real(8),dimension(:)                   :: a
    real(8)                                :: chi2
    real(8),dimension(totNorb)             :: chi2_orb
    complex(8),dimension(Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)              :: Ctmp
    integer                                :: i,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_normal. a:",a
#endif
    !
    Delta = delta_hybrid_normal(a)
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       Ctmp=abs(Gdelta(l,:)-Delta(iorb,jorb,:))
       chi2_orb(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_orb)
    chi2=chi2/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_hybrid_normal. Chi**2:",chi2
#endif
    !
  end function chi2_delta_hybrid_normal

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_delta_hybrid_normal(a) result(dchi2)
    real(8),dimension(:)                           :: a
    real(8),dimension(size(a))                     :: dchi2
    real(8),dimension(totNorb,size(a))             :: df
    complex(8),dimension(Norb,Norb,Ldelta)         :: Delta
    complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Ldelta)                   :: Ftmp
    real(8),dimension(Ldelta)                      :: Ctmp
    integer                                        :: i,j,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_normal. a:",a
#endif
    !
    Delta  = delta_hybrid_normal(a)
    dDelta = grad_delta_hybrid_normal(a)
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       !
       Ftmp = Gdelta(l,:)-Delta(iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dDelta(iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dDelta(iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_normal. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_delta_hybrid_normal




  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  !+-------------------------------------------------------------+
  function chi2_weiss_hybrid_normal(a) result(chi2)
    real(8),dimension(:)                   :: a
    real(8),dimension(totNorb)             :: chi2_orb
    complex(8),dimension(Norb,Norb,Ldelta) :: g0and
    real(8),dimension(Ldelta)              :: Ctmp
    real(8)                                :: chi2
    integer                                :: i,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_hybrid_normal. a:",a
#endif
    !
    g0and(:,:,:) = g0and_hybrid_normal(a)
    !

    do l=1,totNorb
       iorb = getIorb(l)
       jorb = getJorb(l)
       Ctmp = abs(Gdelta(l,:)-g0and(iorb,jorb,:))
       chi2_orb(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_orb)/Norb
    chi2=chi2/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_hybrid_normal. Chi**2:",chi2
#endif
    !
    !
  end function chi2_weiss_hybrid_normal

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function grad_chi2_weiss_hybrid_normal(a) result(dchi2)
    real(8),dimension(:)                           :: a
    real(8),dimension(size(a))                     :: dchi2
    real(8),dimension(totNorb,size(a))             :: df
    complex(8),dimension(Norb,Norb,Ldelta)         :: g0and
    complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dg0and
    complex(8),dimension(Ldelta)                   :: Ftmp
    real(8),dimension(Ldelta)                      :: Ctmp
    integer                                        :: i,j,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_hybrid_normal. a:",a
#endif
    !
    g0and  = g0and_hybrid_normal(a)
    dg0and = grad_g0and_hybrid_normal(a)
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       !
       Ftmp = Gdelta(l,:)-g0and(iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dg0and(iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dg0and(iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_hybrid_normal. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_weiss_hybrid_normal





  !##################################################################
  ! THESE PROCEDURES EVALUATES THE 
  ! - \delta
  ! - \grad \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################
  function delta_hybrid_normal(a) result(Delta)
    real(8),dimension(:)                   :: a
    complex(8),dimension(Norb,Norb,Ldelta) :: Delta
    integer                                :: iorb,jorb
    integer                                :: i,io,l,stride
    real(8),dimension(Nbath)               :: eps
    real(8),dimension(Norb,Nbath)          :: vops
    !
    !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       eps(i)    = a(io)
    enddo
    stride = Nbath
    do l=1,Norb
       do i=1,Nbath
          io = stride + i + (l-1)*Nbath
          vops(l,i) = a(io)
       enddo
    enddo
    !
    do i=1,Ldelta
       do iorb=1,Norb
          Delta(iorb,iorb,i) = sum( vops(iorb,:)*vops(iorb,:)/(xi*Xdelta(i) - eps(:)) )
          do jorb=iorb+1,Norb
             Delta(iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(xi*Xdelta(i) - eps(:)) )
             Delta(jorb,iorb,i) = sum( vops(jorb,:)*vops(iorb,:)/(xi*Xdelta(i) - eps(:)) )
          enddo
       enddo
    enddo
    !
  end function delta_hybrid_normal

  function grad_delta_hybrid_normal(a) result(dDelta)
    real(8),dimension(:)                           :: a
    complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dDelta
    integer                                        :: iorb,jorb
    integer                                        :: i,k,ik,l,io,stride
    real(8),dimension(Nbath)                       :: eps
    real(8),dimension(Norb,Nbath)                  :: vops
    real(8),dimension(Norb,Norb)                   :: delta_orb
    !
    !\grad_{E_{1}(k)} \Delta_{ab} = [ V_{a}(k)*V_{b}(k) / ( iw_n - E_{1}(k) )**2 ]
    !
    !\grad_{V_{g}(k)} \Delta_{ab} = [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ] / ( iw_n - E_{1}(k) )
    !
    stride = 0
    do i=1,Nbath
       io = stride + i
       eps(i)    = a(io)
    enddo
    stride = Nbath
    do l=1,Norb
       do i=1,Nbath
          io = stride + i + (l-1)*Nbath
          vops(l,i) = a(io)
       enddo
    enddo
    !
    delta_orb = eye(Norb)
    !
    do iorb=1,Norb
       do jorb=1,Norb
          stride = 0
          do k=1,Nbath
             ik = stride + k
             dDelta(iorb,jorb,:,ik) = vops(iorb,k)*vops(jorb,k)/(xi*Xdelta(:) - eps(k))**2
          enddo
          stride = Nbath
          do l=1,Norb
             do k=1,Nbath
                ik = stride + k + (l-1)*Nbath
                dDelta(iorb,jorb,:,ik) = (delta_orb(l,iorb)*vops(jorb,k) + delta_orb(l,jorb)*vops(iorb,k))/(xi*Xdelta(:) - eps(k))
             enddo
          enddo
       enddo
    enddo
    !
  end function grad_delta_hybrid_normal




  function g0and_hybrid_normal(a) result(G0and)
    real(8),dimension(:)                   :: a
    complex(8),dimension(Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(Norb,Norb)        :: zeta,fgorb
    integer                                :: i,ispin
    !
    ispin  = Spin_indx
    !
    Delta = delta_hybrid_normal(a)
    !
    do i=1,Ldelta
       fgorb=zero
       zeta = (xi*Xdelta(i)+xmu)*eye(Norb)
       fgorb(:,:)   = zeta(:,:) - impHloc(ispin,ispin,:,:) - Delta(:,:,i)
       call inv(fgorb)
       G0and(:,:,i) = fgorb
    enddo
    !
  end function g0and_hybrid_normal

  function grad_g0and_hybrid_normal(a) result(dG0and)
    real(8),dimension(:)                           :: a
    complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dG0and
    complex(8),dimension(Norb,Norb,Ldelta)         :: G0and
    complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dDelta
    integer                                        :: iorb,jorb
    integer                                        :: ik
    !
    G0and  = g0and_hybrid_normal(a)
    dDelta = grad_delta_hybrid_normal(a)
    do iorb=1,Norb
       do jorb=1,Norb
          do ik=1,size(a)
             dG0and(iorb,jorb,:,ik) = G0and(iorb,jorb,:)*G0and(iorb,jorb,:)*dDelta(iorb,jorb,:,ik)
          enddo
       enddo
    enddo
    !
  end function grad_g0and_hybrid_normal









  !>SUPERC
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
  !         in the SUPERCONDUCTING case.
  !+-------------------------------------------------------------+
  function chi2_delta_hybrid_superc(a) result(chi2)
    real(8),dimension(:)                     :: a
    real(8)                                  :: chi2
    real(8),dimension(totNorb)               :: chi_orb
    complex(8),dimension(2,Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)                ::  Ctmp,Ftmp
    integer                                  ::  i,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_superc. a:",a
#endif
    !
    Delta = delta_hybrid_superc(a)
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       Ctmp = abs(Gdelta(l,:)-Delta(1,iorb,jorb,:))
       Ftmp = abs(Fdelta(l,:)-Delta(2,iorb,jorb,:))
       chi_orb(l) = sum( Ctmp**cg_pow/Wdelta ) + sum( Ftmp**cg_pow/Wdelta )
    enddo
    chi2=sum(chi_orb)/Ldelta/totNorb
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_hybrid_superc. Chi**2:",chi2
#endif
    !
  end function chi2_delta_hybrid_superc

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson 
  ! function in the SUPERCONDUCTING case.
  !+-------------------------------------------------------------+
  function grad_chi2_delta_hybrid_superc(a) result(dchi2)
    real(8),dimension(:)                             ::  a
    real(8),dimension(size(a))                       ::  dchi2
    real(8),dimension(totNorb,size(a))               ::  df
    complex(8),dimension(2,Norb,Norb,Ldelta)         ::  Delta
    complex(8),dimension(Ldelta)                     :: Gtmp,Ftmp
    real(8),dimension(Ldelta)                        :: Ctmp,Btmp
    complex(8),dimension(2,Norb,Norb,Ldelta,size(a)) ::  dDelta
    integer                                          ::  i,j,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_superc. a:",a
#endif
    !
    Delta  = delta_hybrid_superc(a)
    dDelta = grad_delta_hybrid_superc(a)
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       !
       Gtmp =Gdelta(l,:)-Delta(1,iorb,jorb,:)
       Ftmp =Fdelta(l,:)-Delta(2,iorb,jorb,:)
       !
       Ctmp = abs(Gtmp)**(cg_pow-2)
       Btmp = abs(Ftmp)**(cg_pow-2)
       !
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Gtmp)*dreal(dDelta(1,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Gtmp)*dimag(dDelta(1,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dreal(Ftmp)*dreal(dDelta(2,iorb,jorb,:,j))*Btmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dDelta(2,iorb,jorb,:,j))*Btmp/Wdelta )
       enddo
       !
    enddo
    !
    dchi2=-cg_pow*sum(df,dim=1)/Ldelta/totNorb
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_superc. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_delta_hybrid_superc

  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance for G_0 function 
  ! in the SUPERCONDUCTING case.
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  !+-------------------------------------------------------------+
  function chi2_weiss_hybrid_superc(a) result(chi2)
    real(8),dimension(:)                     ::  a
    real(8),dimension(totNorb)               ::  chi2_orb
    complex(8),dimension(2,Norb,Norb,Ldelta) ::  g0and
    real(8)                                  ::  chi2
    real(8),dimension(Ldelta)                ::  Gtmp,Ftmp
    integer                                  ::  i,l,iorb,jorb
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_hybrid_superc. a:",a
#endif
    !
    g0and = g0and_hybrid_superc(a)
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       Gtmp = abs(Gdelta(l,:)-g0and(1,iorb,jorb,:))
       Ftmp = abs(Fdelta(l,:)-g0and(2,iorb,jorb,:))
       chi2_orb(l) = sum( Gtmp**cg_pow/Wdelta )  + sum( Ftmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_orb)/Ldelta/totNorb
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_hybrid_superc. Chi**2:",chi2
#endif
    !
  end function chi2_weiss_hybrid_superc


  ! !+-------------------------------------------------------------+
  ! !PURPOSE: Evaluate the gradient \Grad\chi^2 of 
  ! ! \Delta_Anderson function.
  ! !+-------------------------------------------------------------+
  ! function grad_chi2_weiss_hybrid_superc(a) result(dchi2)
  !   real(8),dimension(:)                   :: a
  !   real(8),dimension(size(a))             :: dchi2
  !   real(8),dimension(size(a))             :: df
  !   complex(8),dimension(2,Ldelta)         :: g0and
  !   complex(8),dimension(2,Ldelta,size(a)) :: dg0and
  !   complex(8),dimension(Ldelta)           :: Gtmp,Ctmp
  !   real(8),dimension(Ldelta)              :: Ctmp,Btmp
  !   integer                                ::  i,j,l,iorb,jorb
  !   !
  ! #ifdef _DEBUG
  !   if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_hybrid_superc. a:",a
  ! #endif
  !   !
  !   g0and  = g0and_normal_normal(a)
  !   dg0and = grad_g0and_normal_normal(a)
  !   !
  !   do l=1,totNorb
  !      iorb=getIorb(l)
  !      jorb=getJorb(l)
  !      !
  !      Gtmp =Gdelta(l,:)-g0and(1,iorb,jorb,:)
  !      Ftmp =Fdelta(l,:)-g0and(2,iorb,jorb,:)
  !      !
  !      Ctmp = abs(Gtmp)**(cg_pow-2)
  !      Btmp = abs(Ftmp)**(cg_pow-2)
  !      !
  !      do j=1,size(a)
  !         df(l,j)=&
  !              sum( dreal(Gtmp)*dreal(dG0and(1,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
  !              sum( dimag(Gtmp)*dimag(dG0and(1,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
  !              sum( dreal(Ftmp)*dreal(dG0and(2,iorb,jorb,:,j))*Btmp/Wdelta ) + &
  !              sum( dimag(Ftmp)*dimag(dG0and(2,iorb,jorb,:,j))*Btmp/Wdelta )
  !      enddo
  !      !
  !   enddo
  !   !
  !   dchi2=-cg_pow*sum(df,dim=1)/Ldelta/totNorb
  ! #ifdef _DEBUG
  !   if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_hybrid_superc. dChi**2:",dchi2
  ! #endif
  !   !
  ! end function grad_chi2_weiss_hybrid_superc







  !##################################################################
  ! THESE PROCEDURES EVALUATES THE 
  ! - \delta
  ! - \grad \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################
  function delta_hybrid_superc(a) result(Delta)
    real(8),dimension(:)                     :: a
    complex(8),dimension(2,Norb,Norb,Ldelta) :: Delta
    integer                                  :: iorb,jorb
    integer                                  :: i,io,k,l,stride
    real(8),dimension(Nbath)                 :: eps,dps
    real(8),dimension(Norb,Nbath)            :: vops
    real(8),dimension(Ldelta,Nbath)          :: Den
    !
    !\Delta_{ab} = - \sum_k [ V_{a}(k) * V_{b}(k) * (iw_n + E(k)) / Den(k) ]
    !
    !\FDelta_{aa} = \sum_k [ \Delta_{a}(k) * V_{a}(k) * V_{a}(k)  / Den(k) ]
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
    do l=1,Norb
       do i=1,Nbath
          io = stride + i + (l-1)*Nbath
          vops(l,i) = a(io)
       enddo
    enddo
    !
    forall(i=1:Ldelta,k=1:Nbath)& !den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
         Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
    !
    do i=1,Ldelta
       do iorb=1,Norb
          Delta(1,iorb,iorb,i) = -sum( vops(iorb,:)*vops(iorb,:)*(xi*Xdelta(i) + eps(:))/Den(i,:) )
          !
          Delta(2,iorb,iorb,i) = sum( dps(:)*vops(iorb,:)*vops(iorb,:)/Den(i,:) )
          do jorb=iorb+1,Norb
             Delta(1,iorb,jorb,i) = -sum( vops(iorb,:)*vops(jorb,:)*(xi*Xdelta(i) + eps(:))/Den(i,:) )
             Delta(1,jorb,iorb,i) = -sum( vops(jorb,:)*vops(iorb,:)*(xi*Xdelta(i) + eps(:))/Den(i,:) )
             !
             Delta(2,iorb,jorb,i) = sum( dps(:)*vops(iorb,:)*vops(jorb,:)/Den(i,:) )
             Delta(2,jorb,iorb,i) = sum( dps(:)*vops(jorb,:)*vops(iorb,:)/Den(i,:) )
          enddo
       enddo
    enddo
    !
  end function delta_hybrid_superc

  function grad_delta_hybrid_superc(a) result(dDelta)
    real(8),dimension(:)                             :: a
    complex(8),dimension(2,Norb,Norb,Ldelta,size(a)) :: dDelta
    integer                                          :: iorb,jorb
    integer                                          :: i,k,io,ik,l,stride
    real(8),dimension(Nbath)                         :: eps,dps
    real(8),dimension(Norb,Nbath)                    :: vops
    real(8),dimension(Ldelta,Nbath)                  :: Den
    real(8),dimension(Norb,Norb)                     :: delta_orb
    !
    !\grad_{E_{1}(k)} \Delta_{ab} = -V_{a}(k)*V_{b}(k)*[ 1/den(k) - 2*E(k)*(iw_n + E(k))/den(k)**2 ]
    !
    !\grad_{\D_{1}(k)} \Delta_{ab} = V_{a}(k)*V_{b}(k)*\D(k)*(iw_n + E(k)) /den(k)**2
    !
    !\grad_{ V_{g}(k)} \Delta_{ab} =  [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ]*(iw_n + E_{a}(k))/den(k)
    !
    !
    !\grad_{E(k)} \FDelta_{ab} = -2 * V_{a}(k) * V_{b}(k) * E(k) * \D(k) / Den**2
    !
    !\grad_{\D(k)} \FDelta_{ab} = V_{a}(k) * V_{b}(k) * [ 1/den - 2* \D(k)*\D(k)/den**2 ]
    !
    !\grad_{ V_{g}(k)} \FDelta_{ab} =  [delta(g,a)*V_{b}(k) + delta(g,b)*V_{a}(k)] * \D(k) / den
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
    do l=1,Norb
       do i=1,Nbath
          io = stride + i + (l-1)*Nbath
          vops(l,i) = a(io)
       enddo
    enddo
    !
    forall(i=1:Ldelta,k=1:Nbath)& !den(k) = (w_n**2 + E_{a}(k)**2 + \D_{a}(k)**2
         Den(i,k) = Xdelta(i)**2 + eps(k)**2 + dps(k)**2 
    !
    delta_orb = eye(Norb)
    !
    do iorb=1,Norb
       do jorb=1,Norb
          !
          stride=0
          do k=1,Nbath
             ik = stride + k
             dDelta(1,iorb,jorb,:,ik) = -vops(iorb,k)*vops(jorb,k)*(1d0/Den(:,k) - 2d0*eps(k)*( xi*Xdelta(:) + eps(k))/Den(:,k)**2)
             dDelta(2,iorb,jorb,:,ik) = 2d0*vops(iorb,k)*vops(jorb,k)*eps(k)*dps(k)/Den(:,k)**2
          enddo
          stride=Nbath
          do k=1,Nbath
             ik = stride + k
             dDelta(1,iorb,jorb,:,ik) = 2d0*vops(iorb,k)*vops(jorb,k)*dps(k)*(xi*Xdelta(:) + eps(k))/Den(:,k)**2
             dDelta(2,iorb,jorb,:,ik) = -vops(iorb,k)*vops(jorb,k)*( 1d0/Den(:,k) - 2d0*dps(k)*dps(k)/Den(:,k)**2 )
          enddo
          stride=2*Nbath
          do l=1,Norb
             do k=1,Nbath
                ik = stride + k + (l-1)*Nbath
                dDelta(1,iorb,jorb,:,ik) = (delta_orb(l,iorb)*vops(jorb,k) + delta_orb(l,jorb)*vops(iorb,k))*(xi*Xdelta(:) + eps(k))/Den(:,k)
                dDelta(2,iorb,jorb,:,ik) = -(delta_orb(l,iorb)*vops(jorb,k) + delta_orb(l,jorb)*vops(iorb,k))*dps(k)/Den(:,k)
             enddo
          enddo
          !
       enddo
    enddo
    !
  end function grad_delta_hybrid_superc

  function g0and_hybrid_superc(a) result(G0and)
    real(8),dimension(:)                     :: a
    complex(8),dimension(2,Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(2,Norb,Norb)        :: zeta
    complex(8),dimension(2*Norb,2*Norb)      :: fgorb
    integer                                  :: i,k,ispin
    !
    ispin  = Spin_indx
    !
    Delta = delta_hybrid_superc(a)
    !
    do i=1,Ldelta
       fgorb=zero
       zeta(1,:,:) = (xi*Xdelta(i)+xmu)*eye(Norb) - impHloc(ispin,ispin,:,:)
       zeta(2,:,:) = (xi*Xdelta(i)-xmu)*eye(Norb) + impHloc(ispin,ispin,:,:)
       fgorb(1:Norb,1:Norb)                     = zeta(1,:,:) - Delta(1,:,:,i)
       fgorb(1:Norb,Norb+1:Norb+Norb)           =             - Delta(2,:,:,i)
       fgorb(Norb+1:Norb+Norb,1:Norb)           =             - Delta(2,:,:,i)
       fgorb(Norb+1:Norb+Norb,Norb+1:Norb+Norb) = zeta(2,:,:) + conjg( Delta(1,:,:,i) )
       call inv(fgorb)
       G0and(1,:,:,i) = fgorb(1:Norb,1:Norb)
       G0and(2,:,:,i) = fgorb(1:Norb,Norb+1:Norb+Norb)
    enddo
    !
  end function g0and_hybrid_superc

  ! function grad_g0and_hybrid_superc(a) result(dG0and)
  !   real(8),dimension(:)                   :: a
  !   complex(8),dimension(2,Ldelta)         :: G0and,Delta
  !   complex(8),dimension(2,Ldelta,size(a)) :: dG0and,dDelta
  !   integer                                :: i,k,ik,io,stride
  !   real(8),dimension(Nbath)               :: eps,vps,dps
  !   integer                                :: iorb,ispin
  !   real(8),dimension(Ldelta,Nbath)        :: Den
  !   complex(8),dimension(Ldelta)           :: g0,f0,dD,dC,dDet,zeta
  !   real(8),dimension(Ldelta)              :: det
  !   !
  !   iorb   = Orb_indx
  !   ispin  = Spin_indx
  !   !
  !   Delta  = delta_normal_superc(a)
  !   dDelta = grad_delta_normal_superc(a)
  !   !
  !   zeta= xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) 
  !   g0  = zeta - Delta(1,:)
  !   f0  =      - Delta(2,:)
  !   Det = abs(g0)**2 + f0**2
  !   !
  !   do k=1,size(a)
  !      dD = conjg(dDelta(1,:,k))
  !      dC = dDelta(2,:,k)
  !      dDet = 2*dreal(g0*dD)+2*f0*dC
  !      dG0and(1,:,k) = -Det*dD + conjg(g0)*dDet
  !      dG0and(2,:,k) = -Det*dC + f0*dDet
  !      dG0and(1,:,k) = dG0and(1,:,k)/Det**2
  !      dG0and(2,:,k) = dG0and(2,:,k)/Det**2
  !   enddo
  ! end function grad_g0and_hybrid_superc








  !>NONSU2
  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  !+-------------------------------------------------------------+
  function chi2_delta_hybrid_nonsu2(a) result(chi2)
    real(8),dimension(:)                               ::  a
    real(8)                                            ::  chi2
    real(8),dimension(totNso)                          ::  chi2_so
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) ::  Delta
    real(8),dimension(Ldelta)                          ::  Ctmp
    integer                                            ::  i,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_hybrid_nonsu2. a:",a
#endif
    !
    Delta = delta_hybrid_nonsu2(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       Ctmp=abs(Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_hybrid_nonsu2. Chi**2:",chi2
#endif
    !
  end function chi2_delta_hybrid_nonsu2


  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  !+-------------------------------------------------------------+
  function chi2_weiss_hybrid_nonsu2(a) result(chi2)
    real(8),dimension(:)                               :: a
    real(8)                                            :: chi2
    real(8),dimension(totNso)                          :: chi2_so
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    real(8),dimension(Ldelta)                          ::  Ctmp
    integer                                            :: i,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_hybrid_nonsu2. a:",a
#endif
    !
    g0and = g0and_hybrid_nonsu2(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       Ctmp=abs(Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)/Ldelta
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_hybrid_nonsu2. Chi**2:",chi2
#endif
    !
    !
  end function chi2_weiss_hybrid_nonsu2






  !##################################################################
  ! THESE PROCEDURES EVALUATES THE 
  ! - \delta
  ! - g0
  ! FUNCTIONS. 
  !##################################################################
  function delta_hybrid_nonsu2(a) result(Delta)
    real(8),dimension(:)                               :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    integer                                            :: ispin,jspin,iorb,jorb
    integer                                            :: i,io,stride,ih
    real(8),dimension(Nspin,Nbath)                     :: hps
    real(8),dimension(Nspin,Norb,Nbath)                :: vops
    real(8),dimension(Nspin,Norb,Nbath)                :: uops
    real(8),dimension(Nhel,Nbath)                      :: ehel
    real(8),dimension(Nhel,Nhel,Norb,Nbath)            :: wohel

    stride = 0
    do ispin=1,Nspin
       do i=1,Nbath
          io = stride + i + (ispin-1)*Nbath
          hps(ispin,i) = a(io)
       enddo
    enddo
    stride = Nspin*Nbath
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nbath
             io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
             vops(ispin,iorb,i) = a(io)
          enddo
       enddo
    enddo
    stride =  Nspin*Nbath + Nspin*Norb*Nbath
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nbath
             io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
             uops(ispin,iorb,i) = a(io)
          enddo
       enddo
    enddo
    ! hps(Nspin,:)    = hps(ispin,:)
    ! vops(Nspin,:,:) = vops(ispin,:,:)
    ! uops(Nspin,:,:) = uops(ispin,:,:)
    !
    wohel = get_Whyb_matrix(vops(1:Nspin,1:Norb,1:Nbath),uops(1:Nspin,1:Norb,1:Nbath))
    !
    Delta=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                do i=1,Ldelta
                   do ih=1,Nspin
                      Delta(ispin,jspin,iorb,jorb,i) = Delta(ispin,jspin,iorb,jorb,i) + &
                           sum( wohel(ispin,ih,iorb,:)*wohel(jspin,ih,jorb,:)/(xi*Xdelta(i) - hps(ih,:)) )
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !
  end function delta_hybrid_nonsu2

  function g0and_hybrid_nonsu2(a) result(G0and)
    real(8),dimension(:)                               :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: zeta,fgorb
    integer                                            :: i,Nso
    integer                                            :: iorb,jorb,ispin,jspin,io,jo
    !
    Nso = Norb*Nspin
    !
    Delta(:,:,:,:,:) = delta_hybrid_nonsu2(a)
    !
    do i=1,Ldelta
       fgorb=zero
       zeta = (xi*Xdelta(i)+xmu)*zeye(Nso)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   fgorb(io,jo)   = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
                enddo
             enddo
          enddo
       enddo
       call inv(fgorb)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
  end function g0and_hybrid_nonsu2





END MODULE ED_FIT_HYBRID
