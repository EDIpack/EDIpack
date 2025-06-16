MODULE ED_PARSE_UMATRIX
  !:synopsis: Routines to parse files containing two-body operators
  USE SCIFOR, only: str,free_unit,file_length,to_lower,txtfy,eye
  USE ED_VARS_GLOBAL
  USE ED_INPUT_VARS
  implicit none
  private

  public :: set_umatrix
  public :: add_twobody_operator
  public :: reset_umatrix
  public :: read_umatrix_file
  public :: save_umatrix_file

contains

  subroutine add_twobody_operator(oi,si,oj,sj,ok,sk,ol,sl,Uijkl)
    !This subroutine lets the user add a two-body operator at runtime.
    !The convention is consistent with that of the umatrix file
    integer                       :: oi !Orbital index of :math:`c^{\dagger}_{i}`
    integer                       :: oj !Orbital index of :math:`c^{\dagger}_{j}`
    integer                       :: ok !Orbital index of :math:`c_{k}`
    integer                       :: ol !Orbital index of :math:`c_{l}`
    character(len=1)              :: si !Spin index of :math:`c^{\dagger}_{i}`
    character(len=1)              :: sj !Spin index of :math:`c^{\dagger}_{j}`
    character(len=1)              :: sk !Spin index of :math:`c_{k}`
    character(len=1)              :: sl !Spin index of :math:`c_{l}`  
    type(coulomb_matrix_element)  :: opline
    character(len=300)            :: dummy
    real(8)                       :: Uijkl !Interaction coefficient
    if(max(oi, oj, ok, ol)>Norb) stop "add_twobody_operator: too many orbitals" 
    if(min(oi, oj, ok, ol) <  1) stop "add_twobody_operator: orbital index < 1" 
    if (any(( [si, sj, sk, sl] /= "u" ) .AND. ( [si, sj, sk, sl] /= "d" ))) stop "add_twobody_operator: spin index malformed" 
    opline%cd_i = [oi, merge(1, 2, si == "u")]
    opline%cd_j = [oj, merge(1, 2, sj == "u")]
    opline%c_k  = [ok, merge(1, 2, sk == "u")]
    opline%c_l  = [ol, merge(1, 2, sl == "u")]
    opline%U    = Uijkl
    if(ed_verbose>3)then
       write(dummy, '(F10.6)') opline%U 
       write(LOGfile,"(A)")'Runtime two-body operator:     '//&
            trim(dummy)//&
            ' cd_['//str(oi)//str(si)//']'//&
            ' cd_['//str(oj)//str(sj)//']'//&
            ' c_['//str(ok)//str(sk)//']'//&
            ' c_['//str(ol)//str(sl)//']'
    endif    
    
    call grow_runtime_array(opline)  
  
  end subroutine add_twobody_operator

  subroutine set_umatrix()
    integer(8)         :: iline
    !This subroutine sets the internal interaction matrices and saves the list of 
    !operators to a file
    mfHloc        = zero
    Uloc_internal = zero
    Ust_internal  = zero
    Jh_internal   = zero
    Jx_internal   = zero
    Jp_internal   = zero
    if(allocated(coulomb_sundry))deallocate(coulomb_sundry)
  
    !If we need to read a umatrix file, we read a umatrix file. This is not at present 
    !compatible with ED_TOTAL_UD: in that case, just use the Kanamori coefficients.
    if(ed_read_umatrix)then
       if(.not. ED_TOTAL_UD) STOP "ED_TOTAL_UD = F and ED_READ_UMATRIX = T are incompatible"
       call read_umatrix_file(umatrix_file)
     endif
     
     !If the user has passed some extra operators, we read them.
     !This in principle doesn't care about ED_TOTAL_UD: we will check later.
     if(allocated(coulomb_runtime))then
       do iline=1,size(coulomb_runtime)
         call parse_umatrix_line(coulomb_runtime(iline))
       enddo
     endif

    !Here we need to operate on the various Uloc, Ust, Jh, Jx, Jp matrices
    !-the Hubbard terms are passed with an 1/2, but if the user made things
    !correctly they already have the two nup/ndw and ndw/nup terms. So nothing to do here.
    !
    !-Ust is symmetric with respect to orbital exchange: the elements we got should be ordered
    !in increasing orbital order. So the matrix should be upper triangular. In the routine that
    !creates Hloc, it does nup*ndw + ndw*nup, so here we divide by 2.0
    !
    Ust_internal = (Ust_internal + transpose(Ust_internal))/2.0
    !
    !-Jh needs to be rescaled. First, it also is symmetric w.r.t. orbital exchange.
    !Then, the coefficient of the terms in the H constructor is actually Ust-Jh. 
    !So, if the user passed this as Jh, we need to recast it as Ust - what the user passed
    !
    Jh_internal = (Jh_internal + transpose(Jh_internal))/2.0
    Jh_internal = Ust_internal - Jh_internal
    !
    !Jx and Jp have a summation that goes from 1 to Norb for both orbital indices, so no change there
 
    !If we use the default input variables, we set them  here. It is important to do that here, after
    !the reshufflings before, because these coefficients are not subject to those!
    if(ed_use_kanamori)then
      if(Norb > 5)STOP "ED_READ_UMATRIX = F: max 5 orbitals allowed"
      Uloc_internal = Uloc_internal + Uloc(1:Norb)
      Ust_internal  = Ust_internal  + Ust - Ust*eye(Norb)
      Jh_internal   = Jh_internal   + Jh  - Jh*eye(Norb)
      Jx_internal   = Jx_internal   + Jx  - Jx*eye(Norb)
      Jp_internal   = Jp_internal   + Jp  - Jp*eye(Norb)
    endif
 
 
    !Print interaction terms
    if(ed_verbose>2)then
      call print_umatrix()
    endif
    
    !Second check for ED_TOTAL_UD: it cannot be 'F' if 'Jx_internal, Jp_internal' or 'Coulomb_sundry'
    !are nonzero, with a caveat. Slightly different condition than the standard one: if only 
    !one orbital is present, Coulomb sundry terms are all fine. They are either the Hubbard U,
    !something impossible or something like n^2 = n. This latter case goes in sundry but does not
    !interfere with ed_total_ud. If more than one orbital is present, the standard check applies.
    if(.not.ed_total_ud)then
      if( (Norb>1 .AND. ( any((Jx_internal/=0d0)) .OR. any((Jp_internal/=0d0)) ) ) &
      .OR. (Norb>1 .AND. allocated(coulomb_sundry)) ) then
        STOP "ED ERROR: ed_total_ud=F cannot be used if non-density-density interaction terms are present"
      endif
    endif
    !
    !Save to file
    call save_umatrix_file()
  end subroutine set_umatrix
  
  subroutine reset_umatrix()
    !This subroutine resets ro zero the internal interaction matrices 
    
     !Uloc = zero
     !Ust = zero
     !Jh = zero
     !Jx = zero
     !Jp = zero
     if(ed_verbose>3)then
       write(LOGfile,"(A)")'Clearing all internal interaction matrix coefficients and user-prodived two-body orbitals.'
     endif
     if(allocated(mfHloc))mfHloc = zero
     if(allocated(Uloc_internal))Uloc_internal = zero
     if(allocated(Ust_internal))Ust_internal = zero
     if(allocated(Jh_internal))Jh_internal = zero
     if(allocated(Jx_internal))Jx_internal = zero
     if(allocated(Jp_internal))Jp_internal = zero
     if(allocated(coulomb_sundry))deallocate(coulomb_sundry)
     if(allocated(coulomb_runtime))deallocate(coulomb_runtime)
  end subroutine reset_umatrix
  

  subroutine print_umatrix()
     !This subroutine pretty-prints the interaction terms divided in Hubbard-Kanamori 
     !and sundry ones
     integer                       :: iorb,jorb,iline,ispin,jspin
     real(8)                       :: dummyU
     character(len=300)            :: dummy
     integer                       :: o1, o2, o3, o4
     character(len=1)              :: s1, s2, s3, s4
 
     write(LOGfile,"(A)")''
     write(LOGfile,"(A)")'Interaction coefficients:'
     write(LOGfile,"(A)")''
 
     write(LOGfile,"(A)")'ULOC:'
     write(LOGfile,"(90(F15.9,1X))") (Uloc_internal(iorb),iorb=1,Norb)
     write(LOGfile,"(A)")''
     !

     write(LOGfile,"(A)")'UST:'
     do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Ust_internal(iorb,jorb),jorb=1,Norb)
     enddo
     write(LOGfile,"(A)")''

     write(LOGfile,"(A)")'JH:'
     do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Jh_internal(iorb,jorb),jorb=1,Norb)
     enddo
     write(LOGfile,"(A)")''

     write(LOGfile,"(A)")'JX:'
     do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Jx_internal(iorb,jorb),jorb=1,Norb)
     enddo
     write(LOGfile,"(A)")'JP:'
     do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Jp_internal(iorb,jorb),jorb=1,Norb)
     enddo
     write(LOGfile,"(A)")''

     write(LOGfile,"(A)")'Mean-field terms from anticommutators:'
     do ispin=1,2
        do iorb=1,Norb
           write(LOGfile,"(100(A1,F8.4,A1,F8.4,A1,2x))")&
                (&
                (&
                '(',dreal(mfHloc(ispin,jspin,iorb,jorb)),',',dimag(mfHloc(ispin,jspin,iorb,jorb)),')',&
                jorb =1,Norb),&
                jspin=1,Nspin)
        enddo
     enddo
     write(LOGfile,"(A)")''

    if(allocated(coulomb_sundry))then
       write(LOGfile,"(A)")'There are '//str(size(coulomb_sundry))//' sundry terms.'
       do iline=1,size(coulomb_sundry)
          o1 = coulomb_sundry(iline)%cd_i(1)
          o2 = coulomb_sundry(iline)%cd_j(1)
          o3 = coulomb_sundry(iline)%c_k(1)
          o4 = coulomb_sundry(iline)%c_l(1)
          s1 = merge("u", "d", coulomb_sundry(iline)%cd_i(2)==1)
          s2 = merge("u", "d", coulomb_sundry(iline)%cd_j(2)==1)
          s3 = merge("u", "d", coulomb_sundry(iline)%c_k(2)==1)
          s4 = merge("u", "d", coulomb_sundry(iline)%c_l(2)==1)
          write(dummy, '(F10.6)') coulomb_sundry(iline)%U
          write(LOGfile,"(A)")'Sundry operator '//txtfy(iline,3)//':     '//&
               trim(dummy)//&
               ' cd_['//str(o1)//str(s1)//']'//&
               ' c_['//str(o2)//str(s2)//']'//&
               ' cd_['//str(o3)//str(s3)//']'//&
               ' c_['//str(o4)//str(s4)//']'
       enddo
       write(LOGfile,"(A)")''
    endif

  end subroutine print_umatrix



  subroutine save_umatrix_file(ufile)
    !This subroutine saves the currently used interaction Hamiltonian into a file. It takes
    !the elements from the internal :f:var:`ULOC`, :f:var:`UST`, :f:var:`JH`, :f:var:`JX`,
    !:f:var:`JP` and sundry terms
    character(len=*),optional     :: ufile  !File to which the interaction operators are saved. Default :f:var:`UMATRIX_FILE` :code:`.used`
    character(len=256)            :: outfile 
    integer                       :: iline,flen,unit_umatrix,rank, nops, ispin, jspin, iorb, jorb, ierr
    real(8)                       :: dummyU
    integer                       :: o1, o2, o3, o4
    character(len=1)              :: s1, s2, s3, s4

    outfile=trim(umatrix_file)//reg(ed_file_suffix)//".used"
    if(present(ufile))outfile=trim(ufile)

    unit_umatrix = free_unit()

    open(unit_umatrix,file=outfile)
    write(unit_umatrix, '(A)') "#Interaction two-body operators"
    write(unit_umatrix, '(I0,1X,A)') Norb, "BANDS"
    !First: write the Uloc terms
    do iorb = 1, Norb
       if(Uloc_internal(iorb)/=0)then
          write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'u', iorb, 'd', iorb, 'u', iorb, 'd', Uloc_internal(iorb)
          write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'd', iorb, 'u', iorb, 'd', iorb, 'u', Uloc_internal(iorb)
       endif
    enddo
    !Second: write the Ust terms
    do iorb = 1, Norb
       do jorb = 1,Norb
          if(Ust_internal(iorb,jorb)/=0.0)then
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'd', jorb, 'u', iorb, 'd', jorb, 'u', Ust_internal(iorb,jorb)
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'u', jorb, 'd', iorb, 'u', jorb, 'd', Ust_internal(iorb,jorb)
          endif
       enddo
    enddo
    !Third a: write the Ust-Jh terms
    do iorb = 1, Norb
       do jorb = 1,Norb
          if(Ust_internal(iorb,jorb)/=0.0)then
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'u', jorb, 'u', iorb, 'u', jorb, 'u', Ust_internal(iorb,jorb)
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'd', jorb, 'd', iorb, 'd', jorb, 'd', Ust_internal(iorb,jorb)
          endif
       enddo
    enddo
    !Third b: write the Ust-Jh terms
    do iorb = 1, Norb
       do jorb = 1,Norb
          if(Jh_internal(iorb,jorb)/=0.0)then
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'u', jorb, 'u', jorb, 'u', iorb, 'u', Jh_internal(iorb,jorb)
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'd', jorb, 'd', jorb, 'd', iorb, 'd', Jh_internal(iorb,jorb)
          endif
       enddo
    enddo
    !Fourth: write the Jx terms
    do iorb = 1, Norb
       do jorb = 1,Norb
          if(Jx_internal(iorb,jorb)/=0.0)then
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'd', jorb, 'u', jorb, 'd', iorb, 'u', Jx_internal(iorb,jorb)
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'u', jorb, 'd', jorb, 'u', iorb, 'd', Jx_internal(iorb,jorb)
          endif
       enddo
    enddo
    !Fifth: write the Jp terms
    do iorb = 1, Norb
       do jorb = 1,Norb
          if(Jp_internal(iorb,jorb)/=0.0)then
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'd', iorb, 'u', jorb, 'd', jorb, 'u', Jp_internal(iorb,jorb)
             write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') iorb, 'u', iorb, 'd', jorb, 'u', jorb, 'd', Jp_internal(iorb,jorb)
          endif
       enddo
    enddo
    !Sixth: write sundry terms
    if(allocated(coulomb_sundry))then
       do iline=1,size(coulomb_sundry)
          o1 = coulomb_sundry(iline)%cd_i(1)
          o2 = coulomb_sundry(iline)%cd_j(1)
          o3 = coulomb_sundry(iline)%c_k(1)
          o4 = coulomb_sundry(iline)%c_l(1)
          s1 = merge("u", "d", coulomb_sundry(iline)%cd_i(2)==1)
          s2 = merge("u", "d", coulomb_sundry(iline)%cd_j(2)==1)
          s3 = merge("u", "d", coulomb_sundry(iline)%c_k(2)==1)
          s4 = merge("u", "d", coulomb_sundry(iline)%c_l(2)==1)
          dummyU = coulomb_sundry(iline)%U
          !Don't forget the anticommutator!
          if(all((coulomb_sundry(iline)%cd_j==coulomb_sundry(iline)%c_k)))then
             dummyU = dummyU + dummyU
          endif
          write(unit_umatrix, '(4(I0,1X,A,1X),ES21.12)') o1, s1, o2, s2, o3, s3, o4, s4, dummyU 
       enddo
    endif
    close(unit_umatrix)

  end subroutine save_umatrix_file


  subroutine read_umatrix_file(ufile)
    !This subroutine reads the interaction Hamiltonian from a user-specified file.
    !Empty and commented lines are ignored. Each two-body operator is then parsed
    !to recognize whether it is of the form :f:var:`ULOC`, :f:var:`UST`, :f:var:`JH`,
    !:f:var:`JX`, :f:var:`JP`. If so, those coefficients are populated. If not, the
    !operator is saved in a dynamically allocated "sundry" array.
    !Since EDIpack applies operators from right to left as 
    !:math:`c\rightarrow c^{\dagger} \rightarrow c \rightarrow c^{\dagger}`
    !the read two-body operators need to be properly ordered. If mean-field terms arise
    !from the anticommutation, these are stored in an array :f:var:`mfHloc` of dimension
    ![:code:`2`, :code:`2`, :f:var:`Norb`, :f:var:`Norb`], which will be added to
    !:f:var:`impHloc` upon Fock space H construction.
    !If :f:var:`ED_VERBOSE` > :code:`2`, this routine will print extensive information
    !about the read file(s) and the type of the operators therein contained.
    !
    character(len=*)              :: ufile  !File containing a properly formatted interaction Hamiltonian. Default :f:var:`UMATRIX_FILE` :code:`.restart`
    character(len=300)            :: dummy
    type(coulomb_matrix_element)  :: opline
    logical                       :: master=.true.,ufile_exists, verbose, preamble
    integer                       :: iline,flen,unit_umatrix,rank, nops, ispin, jspin, iorb, jorb, ierr
    integer                       :: o1, o2, o3, o4
    character(len=1)              :: s1, s2, s3, s4
    !
#ifdef _MPI    
    if(check_MPI())then
       master=get_Master_MPI()
       rank  =get_Rank_MPI()
    endif
#endif
    inquire(file=str(ufile)//str(ed_file_suffix)//".restart",exist=ufile_exists)
    if(.not.ufile_exists)stop "read_umatrix_file ERROR: indicated file does not exist"
    if(ed_verbose>0)write(LOGfile,"(A)")'Reading interaction Hamiltonian from file '//trim(ufile)//reg(ed_file_suffix)//".restart"
    !
    !
    !
    open(free_unit(unit_umatrix),file=trim(ufile)//reg(ed_file_suffix)//".restart")
    !
    !Take care of the comment in the preamble:
    preamble = .true.
    do while(preamble)
       read(unit_umatrix,*) dummy
       if(dummy(1:1)/="#" .and. dummy(1:1)/="!" .and. dummy(1:1)/="%")then
          preamble = .false.
          read(dummy, *) o1
       endif
    enddo
    !
    if(o1/=Norb)STOP "Wrong number of orbitals in umatrix file header"
    !
    !Parse lines
    nops = 0
    do
       read(unit_umatrix, '(A)', iostat=ierr) dummy
       if (ierr /= 0) then
          exit  ! end-of-file error
       else
          read(dummy, *, iostat=ierr) o1, s1, o2, s2, o3, s3, o4, s4, opline%U
          if (ierr /= 0) then
             cycle  ! not a valid data line
          else
             nops = nops + 1
             if(max(o1, o2, o3, o4)>Norb) stop "read_umatrix_file: at line "//str(iline)//" too many orbitals" 
             if(min(o1, o2, o3, o4) <  1) stop "read_umatrix_file: at line "//str(iline)//" orbital index < 1" 
             if (any(( [s1, s2, s3, s4] /= "u" ) .AND. ( [s1, s2, s3, s4] /= "d" ))) stop "read_umatrix_file: at line "//str(iline)//" spin index malformed" 
             opline%cd_i = [o1, merge(1, 2, s1 == "u")]
             opline%cd_j = [o2, merge(1, 2, s2 == "u")]
             opline%c_k  = [o3, merge(1, 2, s3 == "u")]
             opline%c_l  = [o4, merge(1, 2, s4 == "u")]
             if(ed_verbose>2)then
                write(dummy, '(F10.6)') opline%U 
                write(LOGfile,"(A)")'Two-body operator '//txtfy(nops,3)//':     '//&
                     trim(dummy)//&
                     ' cd_['//str(o1)//str(s1)//']'//&
                     ' cd_['//str(o2)//str(s2)//']'//&
                     ' c_['//str(o3)//str(s3)//']'//&
                     ' c_['//str(o4)//str(s4)//']'
             endif
          endif
          call parse_umatrix_line(opline)
       endif
    enddo
    if(ed_verbose>1)write(LOGfile,"(A)")str(nops)//" two-body operators parsed."
    !
    close(unit_umatrix)
    !    
    !
  end subroutine read_umatrix_file


  subroutine parse_umatrix_line(line)
    !This subroutine switches the second and third elements of the coulomb matrix applying
    !commutation relations. If the resulting series of operators is of the type that would
    !couple to Uloc, Ust, Ust-Jh, Jx, Jp, store the resulting coefficient in the appropriate
    !matrix .
    !If some mean-field terms result from the anticommutation, store them in a mean-field 
    !matrix to be added along with impHloc at Fock space H creation time.
    !If a Coulomb matrix element comes along that would be inconsistent with the previously
    !determined value of Uloc, Ust, Ust-Jh, Jx, Jp then set the previously determined coeff
    !to zero and add this to the list of "everything_else".
    type(coulomb_matrix_element)        :: line
    character(len=10)                   :: operator_type
    integer,dimension(2)                :: dummy
    !
    !Zero: onsistently with w2dynamics, 1/2 prefactor is applied by the code.
    !It is also multiplied by -1, because the line in the w2dynamics file is the following:
    ![i j k l U_ijkl], but the operator is defined as cd_i cd_j U_ijkl c_l c_k.
    !Note that l and k are inverted.
    line%U = -0.5d0 * line%U
    !
    !First: order the two creation operators so that they
    !are set in an increasing order of (first) spin and (second) orbital from left to right
    !
    !orbital
    if (line%cd_i(1) > line%cd_j(1)) then
       dummy = line%cd_i
       line%cd_i = line%cd_j
       line%cd_j = dummy
       line%U    = -1.0*line%U
    endif
    !spin (overrides orbital)
    if (line%cd_i(2) > line%cd_j(2)) then
       dummy = line%cd_i
       line%cd_i = line%cd_j
       line%cd_j = dummy
       line%U    = -1.0*line%U
    endif
    !
    !Second: order the two annihilation operators so that they
    !are set in an increasing order of (first) spin and (second) orbital from left to right
    !
    !orbital
    if (line%c_k(1) > line%c_l(1)) then
       dummy = line%c_k
       line%c_k = line%c_l
       line%c_l = dummy
       line%U    = -1.0*line%U
    endif
    !spin (overrides orbital)
    if (line%c_k(2) > line%c_l(2)) then
       dummy = line%c_k
       line%c_k = line%c_l
       line%c_l = dummy
       line%U    = -1.0*line%U
    endif
    !
    !Third: are the indices of the swapped operators the same? If so, there is a mean-field
    !term coming out of the anticommutator
    !
    if(all((line%cd_j==line%c_k)))then
       mfHloc(line%cd_i(2),line%c_k(2),line%cd_i(1),line%c_k(1)) = mfHloc(line%cd_i(2),line%c_k(2),line%cd_i(1),line%c_k(1)) + line%U
    endif
    !
    !Fourth: multiply U by (-1) because operators will always be applied from right to left
    !as c->cd->c->cd, so second and third element are to be swapped
    !
    line%U = -1.0 * line%U
    !
    !Fifth: look at the four-operator term. Does it look like any of the ones we would put 
    !in an Hubbard-Kanamori density-density interaction?
    !
    if(line%cd_i(1)==line%c_k(1))then          !c_dag c: First two operators form a density
       if(line%cd_j(1)==line%c_l(1))then        !c_dag c: Second two operators form a density
          if(line%cd_i(2)/=line%cd_j(2))then          !c_dag spin are opposite: It is either Uloc or Ust
             if(line%cd_i(1)==line%cd_j(1))then        !It is Uloc
                Uloc_internal(line%cd_i(1)) = Uloc_internal(line%cd_i(1)) + line%U
                return
             else    !It is Ust
                Ust_internal(line%cd_i(1),line%cd_j(1)) = Ust_internal(line%cd_i(1),line%cd_j(1)) + line%U
                return
             endif
          else
             if(line%cd_i(1)/=line%cd_j(1))then        !It is Ust - Jh. Store this coefficient in Jh for now. Will need to do Jh = Ust - this in the end.
                Jh_internal(line%cd_i(1),line%cd_j(1)) = Jh_internal(line%cd_i(1),line%cd_j(1)) + line%U
                return
             endif
          endif
       endif
    endif
    !
    !Sixth: is it spin-exchange?
    !S-E: -J c^+_a_up c_a_dw c^+_b_dw c_b_up (i.ne.j)
    !Note that the sign change was already done at step 0
    !
    if(line%cd_i(1) /= line%cd_j(1) .and.& !iorb != jorb
         line%cd_i(1) /= line%c_k(1)  .and.& 
         line%cd_i(2) == line%c_k(2)  .and.&
         line%cd_j(1) /= line%c_l(1)  .and.&
         line%cd_j(2) == line%c_l(2))  then         
       Jx_internal(line%cd_i(1),line%c_k(1)) = Jx_internal(line%cd_i(1),line%c_k(1)) + line%U
       return
    endif
    !
    !Seventh: is it pair-hopping?
    !P-H: -J c^+_iorb_up c_jorb_up c^+_iorb_dw  c_jorb_dw (i.ne.j)
    !Note that the sign change was already done at step 0
    !
    if(line%cd_i(1) /= line%c_k(1)  .and.& !iorb != jorb
         line%cd_i(1) == line%cd_j(1) .and.& 
         line%cd_i(2) /= line%cd_j(2) .and.&
         line%c_k(1)  == line%c_l(1)  .and.&
         line%c_k(2)  /= line%c_l(2))  then          
       Jp_internal(line%cd_i(1),line%c_k(1)) = Jp_internal(line%cd_i(1),line%c_k(1)) + line%U
       return
    endif
    !
    !Eight: if it is none of the above, put this into coulomb_everything_else
    !
    call grow_sundry_array(line)


  end subroutine parse_umatrix_line


  subroutine grow_sundry_array(new_element)
    type(coulomb_matrix_element),intent(in)                 :: new_element
    type(coulomb_matrix_element),dimension(:),allocatable   :: temp
    integer                                                 :: dim_old
    !
    if(.not.allocated(coulomb_sundry))then
       allocate(coulomb_sundry(1))
       coulomb_sundry(1) = new_element
    else     
       dim_old = size(coulomb_sundry)
       allocate(temp(dim_old+1))
       temp(1:dim_old) = coulomb_sundry
       temp(dim_old+1) = new_element
       call move_alloc(temp,coulomb_sundry)
    endif
    !
  end subroutine grow_sundry_array
  
  
  subroutine grow_runtime_array(new_element)
    type(coulomb_matrix_element),intent(in)                 :: new_element
    type(coulomb_matrix_element),dimension(:),allocatable   :: temp
    integer                                                 :: dim_old
    !
    if(.not.allocated(coulomb_runtime))then
       allocate(coulomb_runtime(1))
       coulomb_runtime(1) = new_element
    else     
       dim_old = size(coulomb_runtime)
       allocate(temp(dim_old+1))
       temp(1:dim_old) = coulomb_runtime
       temp(dim_old+1) = new_element
       call move_alloc(temp,coulomb_runtime)
    endif
    !
  end subroutine grow_runtime_array


end module ED_PARSE_UMATRIX






