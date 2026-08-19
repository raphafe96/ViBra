module read_input_file
  implicit none

  contains

  subroutine read_inp(input_file, N_modes, N_expansion, constants_file, &
                      constants_mode, N_quanta, N_states, conv_scf,     &
                      N_threads, point_group, proj_cutoff, max_iter_sci, &
                      sci_mode, i_ref)
    implicit none
    character(len=*), intent(in)  :: input_file
    integer,          intent(out) :: N_modes, N_expansion, N_quanta
    integer,          intent(out) :: N_states, conv_scf, N_threads
    integer,          intent(out) :: max_iter_sci
    real*8,           intent(out) :: proj_cutoff
    character(len=*), intent(out) :: constants_file
    character(len=*), intent(out) :: constants_mode
    character(len=*), intent(out) :: point_group
    character(len=4), intent(out) :: sci_mode
    integer,          intent(out) :: i_ref
    character(len=1)  :: quanta_max_reference_sci
    integer           :: ios, unit
    integer           :: io_mode,    io_expan,   io_ctefile, io_ctmode
    integer           :: io_quanta,  io_nstates, io_conv,    io_threads
    integer           :: io_pg,      io_proj,    io_maxiter
    integer           :: io_scimode, io_quantaref
    character(len=256):: line, keyword
    logical           :: valid

    !--- Initialise outputs ---
    N_modes        = -99999
    N_expansion    = -99999
    N_quanta       = -99999
    N_states       = -99999
    conv_scf       = -99999
    N_threads      = -99999
    max_iter_sci   = -99999
    proj_cutoff    = -99999.d0
    constants_file = ""
    constants_mode = ""
    point_group    = ""
    sci_mode       = ""
    quanta_max_reference_sci = ""
    i_ref          = -99999

    !--- Sentinel values for io_* flags (mark as "not found") ---
    io_mode    = -1 ;  io_expan   = -1 ;  io_ctefile = -1
    io_ctmode  = -1 ;  io_quanta  = -1 ;  io_nstates = -1
    io_conv    = -1 ;  io_threads = -1 ;  io_pg      = -1
    io_proj    = -1 ;  io_maxiter = -1
    io_scimode = -1 ;  io_quantaref = -1

    valid = .true.

    !--- Banner ---
    write(*,'(A)') '========================================'
    write(*,'(A)') '     READING FORCE FIELD CONSTANTS     '
    write(*,'(A)') '========================================'

    open(newunit=unit, file=input_file, status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') '----------------------------------------'
      write(*,'(A)') '     ERROR: Cannot open input file      '
      write(*,'(A)') '----------------------------------------'
      write(*,'(A)') 'File: ' // trim(input_file)
      write(*,'(A)') '========================================'
      stop
    end if

    write(*,'(A)') ' Reading input file: ' // trim(input_file)
    write(*,'(A)') '----------------------------------------'

    !=========================================================================
    ! Parse loop
    !=========================================================================
    do
      read(unit,'(A)', iostat=ios) line
      if (ios /= 0) exit

      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

      read(line,*) keyword

      select case(trim(keyword))

        !--- existing keywords --------------------------------------------------
        case('NMODES','N_modes','number_of_modes')
          read(line,*,iostat=io_mode) keyword, N_modes
          write(*,'(A30,I9)') ' Number of modes: ', N_modes

        case('NEXPAN','N_expansion','number_of_expansion')
          read(line,*,iostat=io_expan) keyword, N_expansion
          write(*,'(A30,I9)') ' Number of Hermite functions: ', N_expansion

        case('FILECT','constants_file','file_name')
          read(line,*,iostat=io_ctefile) keyword, constants_file
          write(*,'(A30,A)') ' Constants file: ', trim(constants_file)

        case('CTEMOD')
          read(line,*,iostat=io_ctmode) keyword, constants_mode
          write(*,'(A30,A)') ' Constants mode: ', trim(constants_mode)

        case('NQUANT')
          read(line,*,iostat=io_quanta) keyword, N_quanta
          write(*,'(A30,I9)') ' Number of quanta: ', N_quanta

        case('CVGSCF')
          read(line,*,iostat=io_conv) keyword, conv_scf
          write(*,'(A30,I9)') ' Convergence VSCF: ', conv_scf

        case('NSTATE')
          read(line,*,iostat=io_nstates) keyword, N_states
          write(*,'(A30,I9)') ' Number of states: ', N_states

        case('THREAD')
          read(line,*,iostat=io_threads) keyword, N_threads
          write(*,'(A30,I9)') ' Number of threads: ', N_threads

        !--- NEW keywords -------------------------------------------------------

        ! Point group symbol, e.g.: PGROUP C2h
        case('PGROUP')
          read(line,*,iostat=io_pg) keyword, point_group
          write(*,'(A30,A)') ' Point group symbol: ', trim(point_group)

        ! Projection cutoff (Angstrom) for symmetry assignment,
        ! i.e. how much an atom is allowed to move under a symop
        ! and still be considered "unmoved".
        ! e.g.: PROJCT 0.05
        case('PROJCT')
          read(line,*,iostat=io_proj) keyword, proj_cutoff
          write(*,'(A30,F12.6)') ' Projection cutoff (Ang): ', proj_cutoff

        ! Maximum number of iterations for Selected CI, optionally followed
        ! by the SCI mode (auto/list) and, if mode is 'auto', a 1-char
        ! lowercase reference excitation level (s, d, t, or q).
        ! e.g.: MAXSCI 200
        !       MAXSCI 200 auto d
        !       MAXSCI 200 list
        case('MAXSCI')
          read(line,*,iostat=io_maxiter) keyword, max_iter_sci
          write(*,'(A30,I9)') ' Max iter Selected CI: ', max_iter_sci

          read(line,*,iostat=io_scimode) keyword, max_iter_sci, sci_mode
          if (io_scimode == 0) then
            write(*,'(A30,A)') ' SCI mode: ', trim(sci_mode)
          end if

          read(line,*,iostat=io_quantaref) keyword, max_iter_sci, sci_mode, &
                                            quanta_max_reference_sci
          if (io_quantaref == 0) then
            write(*,'(A30,A)') ' Quanta max reference SCI: ', &
                               quanta_max_reference_sci
          end if

        !--- unknown ------------------------------------------------------------
        case default
          write(*,*)
          write(*,'(A30,A)') ' >> WARNING: Unknown keyword: ', trim(keyword)
          write(*,'(A)')     '    Available keywords:'
          write(*,'(A)')     '      NMODES  NEXPAN  FILECT  CTEMOD'
          write(*,'(A)')     '      NQUANT  CVGSCF  NSTATE  THREAD'
          write(*,'(A)')     '      PGROUP  PROJCT  MAXSCI'

      end select
    end do

    close(unit)

    !=========================================================================
    ! Validation
    !=========================================================================
    write(*,'(A)') '----------------------------------------'
    write(*,'(A)') '      Validating input parameters       '
    write(*,'(A)') '----------------------------------------'

    !--- NMODES ---
    if (N_modes <= 0 .or. io_mode /= 0) then
      write(*,'(A)') ' ERROR: NMODES must be a positive integer'
      valid = .false.
    else
      write(*,'(A,I6)') ' N_modes    validation : PASSED  => ', N_modes
    end if

    !--- NEXPAN ---
    if (N_expansion <= 0 .or. io_expan /= 0) then
      write(*,'(A)') ' ERROR: NEXPAN must be a positive integer'
      valid = .false.
    else
      write(*,'(A,I6)') ' NEXPAN     validation : PASSED  => ', N_expansion
    end if

    !--- FILECT ---
    if (len_trim(constants_file) == 0 .or. io_ctefile /= 0) then
      write(*,'(A)') ' ERROR: FILECT (constants file) not specified'
      valid = .false.
    else
      open(101, file=constants_file, status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
        write(*,'(A)')   ' ERROR: Constants file not found or unreadable'
        write(*,'(A,A)') '        File: ', trim(constants_file)
        valid = .false.
      else
        write(*,'(A)')   ' FILECT     validation : PASSED'
        close(101)
      end if
    end if

    !--- CTEMOD ---
    if (len_trim(constants_mode) == 0 .or. io_ctmode /= 0) then
      write(*,'(A)') ' ERROR: CTEMOD not specified'
      valid = .false.
    else
      select case(trim(constants_mode))
        case('orca_vpt2','tba')
          write(*,'(A,A)') ' CTEMOD     validation : PASSED  => ', &
                           trim(constants_mode)
        case default
          write(*,'(A)')   ' ERROR: CTEMOD must be orca_vpt2 (or tba)'
          write(*,'(A,A)') '        Current value: ', trim(constants_mode)
          valid = .false.
      end select
    end if

    !--- NQUANT ---
    if (io_quanta /= 0 .or. N_quanta <= -10000) then
      write(*,'(A)') ' ERROR: NQUANT must be an integer'
      valid = .false.
    else
      write(*,'(A,I6)') ' NQUANT     validation : PASSED  => ', N_quanta
    end if

    !--- NSTATE ---
    if (io_nstates /= 0 .or. N_states <= -10000) then
      write(*,'(A)') ' ERROR: NSTATE must be an integer'
      valid = .false.
    else
      write(*,'(A,I6)') ' NSTATE     validation : PASSED  => ', N_states
    end if

    !--- CVGSCF ---
    if (conv_scf <= 0 .or. io_conv /= 0) then
      write(*,'(A)') ' ERROR: CVGSCF must be a positive integer'
      valid = .false.
    else
      write(*,'(A,I6)') ' CVGSCF     validation : PASSED  => ', conv_scf
    end if

    !--- THREAD ---
    if (N_threads <= 0 .or. io_threads /= 0) then
      write(*,'(A)') ' ERROR: THREAD must be a positive integer'
      valid = .false.
    else
      write(*,'(A,I6)') ' THREAD     validation : PASSED  => ', N_threads
    end if

    !--- PGROUP ---
    ! Optional: if absent default to C1 (no symmetry).
    if (io_pg /= 0 .or. len_trim(point_group) == 0) then
      write(*,'(A)') ' PGROUP     not found  : defaulting to C1'
      point_group = 'C1'
      io_pg       = 0          ! clear error flag so validation passes
    else
      write(*,'(A,A)') ' PGROUP     validation : PASSED  => ', &
                       trim(point_group)
    end if

    !--- PROJCT ---
    ! Optional: default 0.05 Angstrom.
    if (io_proj /= 0 .or. proj_cutoff < 0.d0) then
      write(*,'(A)') ' PROJCT     not found  : defaulting to 0.05 Ang'
      proj_cutoff = 0.05d0
      io_proj     = 0
    else
      write(*,'(A,F10.5)') ' PROJCT     validation : PASSED  => ', proj_cutoff
    end if

    !--- MAXSCI ---
    ! Optional: default 100.
    if (io_maxiter /= 0 .or. max_iter_sci < 0) then
      write(*,'(A)') ' MAXSCI     not found  : defaulting to 100'
      max_iter_sci = 100
      io_maxiter   = 0
    else
      write(*,'(A,I6)') ' MAXSCI     validation : PASSED  => ', max_iter_sci
    end if

    !--- SCI mode (auto/list) and quanta reference ---
    ! Always try to read/validate sci_mode if present, but only enforce
    ! strict rules when SCI is actually active (max_iter_sci > 0).
    if (io_scimode == 0) then
      ! Mode was provided on the MAXSCI line.
      select case(trim(sci_mode))
        case('auto','list')
          ! Valid.
          write(*,'(A,A)') ' SCIMODE    validation : PASSED  => ', &
                           trim(sci_mode)
        case default
          if (max_iter_sci > 0) then
            write(*,'(A)') ' ERROR: SCIMODE must be auto or list (lowercase).'
            write(*,'(A,A)') '        Got: ', trim(sci_mode)
            valid = .false.
          else
            write(*,'(A)') ' WARNING: Invalid SCIMODE - defaulting to auto.'
            sci_mode = 'auto'
          end if
      end select
    else
      ! Mode was not given.
      if (max_iter_sci > 0) then
        write(*,'(A)') ' ERROR: SCIMODE not specified. Must be auto or list.'
        valid = .false.
      else
        sci_mode = 'auto'
        write(*,'(A)') ' SCIMODE    not found  : defaulting to auto'
      end if
    end if

    if (max_iter_sci > 0) then
      ! SCI is active – process quanta reference if mode is 'auto'.
      if (trim(sci_mode) == 'auto') then
        if (io_quantaref /= 0) then
          write(*,'(A)') ' SCIQREF    not found  : defaulting to d'
          quanta_max_reference_sci = 'd'
        else
          write(*,'(A,A)') ' SCIQREF    validation : PASSED  => ', &
                           quanta_max_reference_sci
        end if

        ! Convert reference letter to excitation level.
        select case(trim(quanta_max_reference_sci))
          case('s'); i_ref = 1
          case('d'); i_ref = 2
          case('t'); i_ref = 3
          case('q'); i_ref = 4
          case default
            write(*,'(A)') '----------------------------------------'
            write(*,'(A)') ' ERROR: SCIQREF must be one of: s, d, t, q'
            write(*,'(A)') '        (lowercase only)'
            write(*,'(A,A)') '        Got: ', quanta_max_reference_sci
            write(*,'(A)') '----------------------------------------'
            stop
        end select

        if (N_quanta < i_ref + 1) then
          write(*,'(A)') '----------------------------------------'
          write(*,'(A)') ' ERROR: NQUANT is too small for the requested'
          write(*,'(A)') '        Selected CI reference excitation level'
          write(*,'(A,A,A,I3)') '        SCIQREF = ', &
                                trim(quanta_max_reference_sci), &
                                '  requires NQUANT >= ', i_ref + 1
          write(*,'(A,I3)') '        Current NQUANT = ', N_quanta
          write(*,'(A)') '        Please increase NQUANT in the input file.'
          write(*,'(A)') '----------------------------------------'
          stop
        end if
      else   ! sci_mode == 'list'
        i_ref = 2
      end if
    else
      ! SCI not used; i_ref is irrelevant.
      i_ref = 2
    end if

    !--- Final result ---
    write(*,'(A)') '----------------------------------------'
    if (valid) then
      write(*,'(A)') ' INPUT VALIDATION      : SUCCESS'
      write(*,'(A)') '----------------------------------------'
    else
      write(*,'(A)') ' INPUT VALIDATION      : FAILED'
      write(*,'(A)') '----------------------------------------'
      stop
    end if

  end subroutine read_inp

  subroutine read_intg(keyword, val, default_val)
  implicit none
  character(len=6), intent(in)  :: keyword
  integer,          intent(out) :: val
  integer,          intent(in)  :: default_val

  integer            :: iounit, ios, ios2
  logical            :: exists
  character(len=100) :: line

  ! Check if file exists
  inquire(file='extra_input.txt', exist=exists)
  if (.not. exists) then
    val = default_val
    return
  end if
  
  ! Open file and search for keyword
  open(newunit=iounit, file='extra_input.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    val = default_val
    return
  end if

  do
    read(iounit, '(A)', iostat=ios) line
    if (ios /= 0) exit   ! end of file or error

    ! Skip lines shorter than 6 characters
    if (len_trim(line) < 6) cycle

    ! Compare first 6 characters with the given keyword
    if (line(1:6) == keyword) then
      ! Attempt to read integer from the remainder of the line
      read(line(7:), *, iostat=ios2) val
      if (ios2 == 0) then
        close(iounit)
        return
      else
        ! Read error – fall back to default
        exit
      end if
    end if
  end do

  close(iounit)
  val = default_val
end subroutine read_intg


subroutine read_real(keyword, val, default_val)
  implicit none
  character(len=6), intent(in)  :: keyword
  real*8,             intent(out) :: val
  real*8,             intent(in)  :: default_val

  integer            :: iounit, ios, ios2
  logical            :: exists
  character(len=100) :: line

  ! Check if file exists
  inquire(file='extra_input.txt', exist=exists)
  if (.not. exists) then
    val = default_val
    return
  end if

  ! Open file and search for keyword
  open(newunit=iounit, file='extra_input.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    val = default_val
    return
  end if

  do
    read(iounit, '(A)', iostat=ios) line
    if (ios /= 0) exit

    if (len_trim(line) < 6) cycle

    if (line(1:6) == keyword) then
      ! Attempt to read real number from the remainder of the line
      read(line(7:), *, iostat=ios2) val
      if (ios2 == 0) then
        close(iounit)
        return
      else
        exit
      end if
    end if
  end do

  close(iounit)
  val = default_val
end subroutine read_real



end module read_input_file