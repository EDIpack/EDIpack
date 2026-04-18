
subroutine get_denmat_n4(fulldm,doprint)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
    use ED_VARS_GLOBAL, only: Ns, Nnambu
#endif
      complex(8),dimension(:,:,:,:),intent(inout) :: fulldm
      logical,intent(in),optional             :: doprint
      logical                                 :: doprint_
      !
      doprint_=.false.; if(present(doprint)) doprint_=doprint
      !
      if(.not.allocated(full_denmat))stop "ERROR: full_1body_density_matrix is not allocated"
      !
      call assert_shape(fulldm,shape(full_denmat),"get_full_denmat","fulldm")
      fulldm = full_denmat
      !
      if(doprint_)call print_rdm(fulldm,Ns)
      !
    contains
      !
      subroutine print_rdm(dm,Nrdm)
        integer                  ,intent(in)            :: Nrdm
        complex(8),dimension(:,:,:,:),intent(in)        :: dm
        integer                                         :: unit
        character(len=64)                               :: fname, ssname
        integer                                         :: io,jo, is,js
        !
        if(size(dm,1)/=Nrdm.OR.size(dm,2)/=Nrdm)then
           stop "ERROR: actual dm argument has incogruent size wrt explicitly passed Nrdm"
        endif
        !
        !
        do is=1,Nspin*Nnambu
            do js=1,Nspin*Nnambu
                write(ssname,'(I0,I0)') is, js
                fname = "full_1body_density_matrix"//str(ed_file_suffix)//"_s"//ssname//".ed"
                !
                unit = free_unit()
                open(unit,file=fname,action="write",position="rewind",status='unknown')
                do io=1,Nrdm
                    write(unit,"(*(F20.16,1X))") (dreal(dm(is,js,io,jo)),jo=1,Nrdm)
                enddo
                if(any(dimag(dm)/=0d0))then
                    write(unit,*)
                    do io=1,Nrdm
                        write(unit,"(*(F20.16,1X))") (dimag(dm(is,js,io,jo)),jo=1,Nrdm)
                    enddo
                endif
                close(unit)
            enddo
        enddo
        !
      end subroutine print_rdm
      !
end subroutine get_denmat_n4

subroutine get_denmat_n2(fulldm,doprint)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
    use ED_VARS_GLOBAL, only: Ns, Nnambu
#endif
          complex(8),dimension(:,:),intent(inout) :: fulldm
          logical,intent(in),optional             :: doprint
          logical                                 :: doprint_
          integer, dimension(2)                   :: target_shape
          integer                                         :: io,jo, is,js
          !
          doprint_=.false.; if(present(doprint)) doprint_=doprint
          !
          if(.not.allocated(full_denmat))stop "ERROR: full_1body_density_matrix is not allocated"
          !
          target_shape(1) = size(full_denmat,1)*size(full_denmat,3)
          target_shape(2) = size(full_denmat,2)*size(full_denmat,4)
          call assert_shape(fulldm,target_shape,"get_full_denmat","fulldm")
          !
          do is=1,Nspin*Nnambu
            do js=1,Nspin*Nnambu
              do io=1,Ns
                do jo=1,Ns
                  fulldm( io +(is-1)*Ns , jo +(js-1)*Ns ) = full_denmat(is,js,io,jo)
                enddo
              enddo
            enddo
          enddo
          !
          if(doprint_)call print_rdm(fulldm,Ns*Nspin*Nnambu)
          !
        contains
          !
          subroutine print_rdm(dm,Nrdm)
            integer                  ,intent(in)            :: Nrdm
            complex(8),dimension(:,:),intent(in)            :: dm
            integer                                         :: unit
            character(len=64)                               :: fname, ssname
            integer                                         :: io,jo
            !
            if(size(dm,1)/=Nrdm.OR.size(dm,2)/=Nrdm)then
               stop "ERROR: actual dm argument has incogruent size wrt explicitly passed Nrdm"
            endif
            !
            !
            write(ssname,'(I0,I0)') is, js
            fname = "full_1body_density_matrix"//str(ed_file_suffix)//".ed"
            !
            unit = free_unit()
            open(unit,file=fname,action="write",position="rewind",status='unknown')
            do io=1,Nrdm
                write(unit,"(*(F20.16,1X))") (dreal(dm(io,jo)),jo=1,Nrdm)
            enddo
            if(any(dimag(dm)/=0d0))then
                write(unit,*)
                do io=1,Nrdm
                    write(unit,"(*(F20.16,1X))") (dimag(dm(io,jo)),jo=1,Nrdm)
                enddo
            endif
            close(unit)
            !
          end subroutine print_rdm
          !
end subroutine get_denmat_n2
    
    