module symmetry_module
  implicit none

  real*8, parameter :: pi_sym = 3.14159265358979323846d0

  !===========================================================================
  ! Point group and irrep data
  !===========================================================================
  character(len=10) :: point_group_name = 'C1'
  integer           :: n_irreps  = 1
  integer           :: n_symops  = 1

  character(len=8), allocatable :: irrep_names(:)
  integer,          allocatable :: mode_irrep(:)
  integer,          allocatable :: irrep_product(:,:)

  real*8 :: proj_cutoff_sym

  !===========================================================================
  ! Geometry / normal mode data
  !===========================================================================
  integer                        :: n_atoms      = 0
  integer                        :: n_modes_file = 0
  real*8,  allocatable           :: atom_xyz(:,:)
  real*8,  allocatable           :: atom_mass(:)
  character(len=4), allocatable  :: atom_label(:)
  integer,          allocatable  :: atom_Z(:)
  real*8,  allocatable           :: mode_freq(:)
  real*8,  allocatable           :: mode_disp(:,:)

  !===========================================================================
  ! Symmetry operations and character table
  !===========================================================================
  real*8,  allocatable           :: symop_mat(:,:,:)      ! PA-frame ops
  character(len=8), allocatable  :: symop_name(:)
  real*8,  allocatable           :: char_table(:,:)

  !===========================================================================
  ! Principal axis frame
  !===========================================================================
  real*8              :: pa_axes(3,3)
  real*8, allocatable :: xyz_pa_global(:,:)

  !===========================================================================
  ! Atom permutation table
  !===========================================================================
  integer, allocatable :: atom_image(:,:)

  real*8, allocatable :: mode_char(:,:)

contains

  !==========================================================================
  ! PUBLIC: init_symmetry
  !==========================================================================
  subroutine init_symmetry(N_modes_vci, irrep_of_mode, proj_cutoff_in, warning)
    integer, intent(in)  :: N_modes_vci
    integer, intent(out) :: irrep_of_mode(N_modes_vci)
    real*8,  intent(in)  :: proj_cutoff_in
    integer :: i, n_vib, warning
    integer, parameter :: n_skip = 6
    character(len=10)  :: pg_input

    proj_cutoff_sym = proj_cutoff_in

    call read_normal_modes()
    call read_point_group_from_user(pg_input)
    call setup_group_from_name(pg_input)

    call build_atom_image_table()

    call debug_mode_displacements(1, n_modes_file)

    call assign_mode_irreps_all()

    warning = 0
    irrep_of_mode = 1
    n_vib = 0
    do i = n_skip + 1, n_modes_file
      n_vib = n_vib + 1
      if (n_vib <= N_modes_vci) irrep_of_mode(n_vib) = mode_irrep(i)
    end do

    if (n_vib < N_modes_vci) &
      write(*,'(A,I4,A,I4)') &
        ' WARNING sym: only ', n_vib, &
        ' modes after skipping 6, expected ', N_modes_vci

    call print_symmetry_summary(101, N_modes_vci, irrep_of_mode, warning)
    call print_symmetry_summary(6,   N_modes_vci, irrep_of_mode, warning)

    call rewrite_normal_mode_pa_frame()
  end subroutine init_symmetry


  !==========================================================================
  ! Debug: show how normal mode displacements behave under each symop
  !==========================================================================
  subroutine debug_mode_displacements(mode_first, mode_last)
    integer, intent(in) :: mode_first, mode_last

    real*8  :: L_pa(3,n_atoms), disp_lab(3)
    real*8  :: RL(3)
    real*8  :: numerator, denominator, char_val
    real*8  :: contrib
    integer :: imode, iat, iop, jat


    write(102,'(A)') ' '
    write(102,'(A)') ' ============================================================'
    write(102,'(A)') '  NORMAL MODE DISPLACEMENT DIAGNOSTIC'
    write(102,'(A)') ' ============================================================'

    do imode = mode_first, mode_last

      write(102,'(A)') ' '
      write(102,'(A,I4,A,F12.3,A)') &
        ' --- Mode ', imode, '   freq = ', mode_freq(imode), ' cm-1'

      !--- Lab-frame displacements
      write(102,'(A)') '   Lab-frame displacements (dx, dy, dz):'
      write(102,'(A6,A4,3A14)') '  idx','sym','    dx','    dy','    dz'
      do iat = 1, n_atoms
        disp_lab(1) = mode_disp(3*(iat-1)+1, imode)
        disp_lab(2) = mode_disp(3*(iat-1)+2, imode)
        disp_lab(3) = mode_disp(3*(iat-1)+3, imode)
        write(102,'(2X,I4,1X,A4,3F14.7)') &
          iat, trim(atom_label(iat)), disp_lab(1), disp_lab(2), disp_lab(3)
      end do

      !--- Rotate to PA frame
      do iat = 1, n_atoms
        disp_lab(1) = mode_disp(3*(iat-1)+1, imode)
        disp_lab(2) = mode_disp(3*(iat-1)+2, imode)
        disp_lab(3) = mode_disp(3*(iat-1)+3, imode)
        L_pa(1,iat) = dot_product(pa_axes(:,1), disp_lab)
        L_pa(2,iat) = dot_product(pa_axes(:,2), disp_lab)
        L_pa(3,iat) = dot_product(pa_axes(:,3), disp_lab)
      end do

      write(102,'(A)') ' '
      write(102,'(A)') '   PA-frame displacements (da, db, dc):'
      write(102,'(A6,A4,3A14,A14)') &
        '  idx','sym','    da','    db','    dc','  |d|^2'
      do iat = 1, n_atoms
        write(102,'(2X,I4,1X,A4,3F14.7,F14.7)') &
          iat, trim(atom_label(iat)), &
          L_pa(1,iat), L_pa(2,iat), L_pa(3,iat), &
          L_pa(1,iat)**2 + L_pa(2,iat)**2 + L_pa(3,iat)**2
      end do

      denominator = sum(L_pa**2)
      write(102,'(A,F14.7)') '   Total |L|^2 = ', denominator

      !--- Character under each symop (already in PA frame)
      do iop = 1, n_symops

        numerator = 0.d0

        write(102,'(A)') ' '
        write(102,'(A,I2,A,A)') &
          '   >>> Symop ', iop, ': ', trim(symop_name(iop))

        write(102,'(A)') '       All atom contributions:'
        write(102,'(A6,A4,A6,3A12,3A12,A12,A8)') &
          '  idx','sym','->img', &
          '  da','  db','  dc', &
          '  Rda','  Rdb','  Rdc', &
          '  d_j.Rd_i','  type'

        do iat = 1, n_atoms
          jat = atom_image(iat, iop)
          if (jat == 0) then
            write(102,'(2X,I4,1X,A4,A6,A)') &
              iat, trim(atom_label(iat)), ' NONE', '  NO IMAGE - skipped'
            cycle
          end if

          ! Apply PA-frame symop to PA-frame displacement
          RL(1) = symop_mat(1,1,iop)*L_pa(1,iat) &
                + symop_mat(1,2,iop)*L_pa(2,iat) &
                + symop_mat(1,3,iop)*L_pa(3,iat)
          RL(2) = symop_mat(2,1,iop)*L_pa(1,iat) &
                + symop_mat(2,2,iop)*L_pa(2,iat) &
                + symop_mat(2,3,iop)*L_pa(3,iat)
          RL(3) = symop_mat(3,1,iop)*L_pa(1,iat) &
                + symop_mat(3,2,iop)*L_pa(2,iat) &
                + symop_mat(3,3,iop)*L_pa(3,iat)

          contrib = L_pa(1,jat)*RL(1) &
                  + L_pa(2,jat)*RL(2) &
                  + L_pa(3,jat)*RL(3)

          numerator = numerator + contrib

          if (jat == iat) then
            write(102,'(2X,I4,1X,A4,I6,3F12.6,3F12.6,F12.6,A8)') &
              iat, trim(atom_label(iat)), jat, &
              L_pa(1,iat), L_pa(2,iat), L_pa(3,iat), &
              RL(1), RL(2), RL(3), contrib, '  SELF'
          else
            write(102,'(2X,I4,1X,A4,I6,3F12.6,3F12.6,F12.6,A8)') &
              iat, trim(atom_label(iat)), jat, &
              L_pa(1,iat), L_pa(2,iat), L_pa(3,iat), &
              RL(1), RL(2), RL(3), contrib, '  PERM'
          end if

        end do

        if (denominator > 1.d-20) then
          char_val = numerator / denominator
        else
          char_val = 1.d0
        end if

        write(102,'(A,F12.6,A,F12.6,A,F10.5)') &
          '       numerator=', numerator, &
          '  denominator=', denominator, &
          '  chi=', char_val

      end do ! iop

      write(102,'(A)') ' '
      write(102,'(A)') ' ============================================================'

    end do ! imode

    

  end subroutine debug_mode_displacements


  !==========================================================================
  ! Print symmetry summary
  !==========================================================================
  subroutine print_symmetry_summary(iunit, N_modes_vci, irrep_of_mode, warning)
    integer, intent(in)    :: iunit, N_modes_vci
    integer, intent(in)    :: irrep_of_mode(N_modes_vci)
    integer, intent(inout) :: warning
    integer :: i, j, n_vib
    integer, parameter :: n_skip = 6

    write(iunit,'(A)') '========================================'
    write(iunit,'(A,A)')      '  Point group : ', trim(point_group_name)
    write(iunit,'(A,I4)')     '  1D irreps   : ', n_irreps
    write(iunit,'(A,F8.4,A)') '  Proj. cutoff: ', proj_cutoff_sym, ' Ang'
    write(iunit,'(A)') '========================================'
    write(iunit,'(A)') '  Mode  Freq(cm-1)   Irrep'
    write(iunit,'(A)') '  ----  ----------   -----'
    n_vib = 0
    do i = n_skip + 1, n_modes_file
      n_vib = n_vib + 1
      if (n_vib <= N_modes_vci) &
        write(iunit,'(I6,F12.3,3X,A)') n_vib, mode_freq(i), &
             trim(irrep_names(irrep_of_mode(n_vib)))
    end do
    write(iunit,'(A)') '========================================'
    write(iunit,'(A)') '  Characters per mode:'
    write(iunit,'(A)', advance='no') '    Mode   '
    do j = 1, n_symops
      write(iunit,'(A8)', advance='no') symop_name(j)
    end do
    write(iunit,*)
    n_vib = 0
    do i = n_skip + 1, n_modes_file
      n_vib = n_vib + 1
      if (n_vib <= N_modes_vci) then
        write(iunit,'(I8)', advance='no') n_vib
        do j = 1, n_symops
          write(iunit,'(F8.3)', advance='no') mode_char(i, j)
          if (iunit == 101) then
            if (abs(mode_char(i,j)) < 0.95d0 .and. &
                abs(mode_char(i,j)) > 0.05d0) warning = 1
          end if
        end do
        write(iunit,*)
      end if
    end do
    write(iunit,'(A)') '========================================'
    write(iunit,'(A)') ' '

    if (warning == 1) then
      write(iunit,'(A)') '<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>'
      write(iunit,'(A)') '                WARNING                '
      write(iunit,'(A)') '        Characters badly assigned      '
      write(iunit,'(A)') 'Consider using different point group/C1'
      write(iunit,'(A)') '<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>'
      write(iunit,'(A)') ' '
    end if

  end subroutine print_symmetry_summary


  !==========================================================================
  ! Read point group
  !==========================================================================
  subroutine read_point_group_from_user(pg)
    character(len=10), intent(out) :: pg
    integer :: ios
    character(len=80) :: line

    open(unit=56, file='point_group.txt', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') ' ERROR: point_group.txt not found.'; stop
    end if
    read(56,'(A)', iostat=ios) line; close(56)
    if (ios /= 0) then
      write(*,'(A)') ' ERROR: could not read point_group.txt'; stop
    end if
    line = adjustl(trim(line))
    call canonicalise_pg_name(line, pg)
    write(*,'(A,A)') ' Point group: ', trim(pg)

  end subroutine read_point_group_from_user


  !==========================================================================
  ! Canonicalise
  !==========================================================================
  subroutine canonicalise_pg_name(raw, canonical)
    character(len=*),  intent(in)  :: raw
    character(len=10), intent(out) :: canonical
    character(len=80) :: upper
    integer :: i, n

    n = len_trim(raw); upper = ''
    do i = 1, n
      upper(i:i) = raw(i:i)
      if (raw(i:i) >= 'a' .and. raw(i:i) <= 'z') &
        upper(i:i) = char(ichar(raw(i:i)) - 32)
    end do
    upper = trim(adjustl(upper))

    select case (trim(upper))
      case ('C1');  canonical = 'C1'
      case ('CS');  canonical = 'Cs'
      case ('CI');  canonical = 'Ci'
      case ('C2');  canonical = 'C2'
      case ('C2H'); canonical = 'C2h'
      case ('C2V'); canonical = 'C2v'
      case ('D2');  canonical = 'D2'
      case ('D2H'); canonical = 'D2h'
      case default
        write(*,'(A,A)') ' ERROR: unsupported point group: ', trim(raw); stop
    end select

  end subroutine canonicalise_pg_name


  !==========================================================================
  ! Dispatch
  !==========================================================================
  subroutine setup_group_from_name(pg)
    character(len=10), intent(in) :: pg

    call compute_principal_frame()

    select case (trim(pg))
      case ('C1');  call setup_C1()
      case ('Cs');  call setup_Cs()
      case ('Ci');  call setup_Ci()
      case ('C2');  call setup_C2()
      case ('C2h'); call setup_C2h()
      case ('C2v'); call setup_C2v()
      case ('D2');  call setup_D2()
      case ('D2h'); call setup_D2h()
      case default
        write(*,'(A,A)') ' ERROR: unknown group: ', trim(pg); stop
    end select

    write(*,'(A,A,A,I2,A,I2,A)') ' Group ', trim(point_group_name), &
      ': ', n_irreps, ' irreps, ', n_symops, ' symops.'

  end subroutine setup_group_from_name


  !==========================================================================
  ! PUBLIC: state_irrep_func
  !==========================================================================
  function state_irrep_func(config, N_modes) result(irr)
    integer, intent(in) :: N_modes, config(N_modes)
    integer :: irr, i
    irr = 1
    do i = 1, N_modes
      if (mod(config(i), 2) == 1) irr = irrep_product(irr, mode_irrep(i))
    end do
  end function state_irrep_func


  !==========================================================================
  ! PUBLIC: irrep_label
  !==========================================================================
  function irrep_label(irr) result(nm)
    integer, intent(in) :: irr
    character(len=8)    :: nm
    if (irr >= 1 .and. irr <= n_irreps) then
      nm = irrep_names(irr)
    else
      nm = '?'
    end if
  end function irrep_label


  !**************************************************************************
  !  PRIVATE ROUTINES
  !**************************************************************************

  !==========================================================================
  ! Read normal_mode.txt
  !==========================================================================
  subroutine read_normal_modes()
    integer           :: ios, imode, iat, i
    character(len=10) :: word1
    integer           :: dummy_int

    open(unit=55, file='normal_mode.txt', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: cannot open normal_mode.txt'; stop
    end if

    read(55,*) n_atoms
    n_modes_file = 3 * n_atoms

    if (allocated(atom_label)) deallocate(atom_label)
    if (allocated(atom_xyz))   deallocate(atom_xyz)
    if (allocated(atom_mass))  deallocate(atom_mass)
    if (allocated(atom_Z))     deallocate(atom_Z)
    allocate(atom_label(n_atoms), atom_xyz(3,n_atoms), &
             atom_mass(n_atoms),  atom_Z(n_atoms))

    do iat = 1, n_atoms
      read(55,*) atom_label(iat), &
                 atom_xyz(1,iat), atom_xyz(2,iat), atom_xyz(3,iat)
      atom_mass(iat) = element_mass(atom_label(iat))
      atom_Z(iat)    = element_Z(atom_label(iat))
    end do

    if (allocated(mode_freq)) deallocate(mode_freq)
    if (allocated(mode_disp)) deallocate(mode_disp)
    allocate(mode_freq(n_modes_file), mode_disp(n_modes_file, n_modes_file))
    mode_freq = 0.d0; mode_disp = 0.d0

    do imode = 1, n_modes_file
      read(55,*) word1, dummy_int
      read(55,*) mode_freq(imode)
      do i = 1, n_modes_file
        read(55,*) mode_disp(i, imode)
      end do
    end do

    close(55)

    write(*,'(A,I4,A,I4,A)') ' Read ', n_atoms, &
      ' atoms, ', n_modes_file, ' modes from normal_mode.txt'
    write(*,'(A)') ' Atoms:'
    do iat = 1, n_atoms
      write(*,'(2X,I3,1X,A4,A,I3,A,F8.4,A,3F10.4)') &
        iat, trim(atom_label(iat)), &
        '  Z=', atom_Z(iat), '  m=', atom_mass(iat), &
        '  xyz=', atom_xyz(:,iat)
    end do

  end subroutine read_normal_modes


  !==========================================================================
  ! Atomic number
  !==========================================================================
  function element_Z(sym) result(Z)
    character(len=*), intent(in) :: sym
    integer :: Z
    select case (trim(adjustl(sym)))
      case('H','D'); Z=1;  case('He'); Z=2;  case('Li'); Z=3
      case('Be');    Z=4;  case('B');  Z=5;  case('C');  Z=6
      case('N');     Z=7;  case('O');  Z=8;  case('F');  Z=9
      case('Ne');    Z=10; case('Na'); Z=11; case('Mg'); Z=12
      case('Al');    Z=13; case('Si'); Z=14; case('P');  Z=15
      case('S');     Z=16; case('Cl'); Z=17; case('Ar'); Z=18
      case('K');     Z=19; case('Ca'); Z=20; case('Fe'); Z=26
      case('Zn');    Z=30; case('Br'); Z=35; case('I');  Z=53
      case default
        write(*,*) 'WARNING: unknown element "'//trim(sym)//'", Z=0'; Z=0
    end select
  end function element_Z


  !==========================================================================
  ! Atomic mass
  !==========================================================================
  function element_mass(sym) result(mass)
    character(len=*), intent(in) :: sym
    real*8 :: mass
    select case (trim(adjustl(sym)))
      case('H');  mass=1.00794d0;  case('D');  mass=2.01410d0
      case('He'); mass=4.00260d0;  case('Li'); mass=6.941d0
      case('Be'); mass=9.01218d0;  case('B');  mass=10.811d0
      case('C');  mass=12.0107d0;  case('N');  mass=14.0067d0
      case('O');  mass=15.9994d0;  case('F');  mass=18.9984d0
      case('Ne'); mass=20.1797d0;  case('Na'); mass=22.9898d0
      case('Mg'); mass=24.305d0;   case('Al'); mass=26.9815d0
      case('Si'); mass=28.0855d0;  case('P');  mass=30.9738d0
      case('S');  mass=32.065d0;   case('Cl'); mass=35.453d0
      case('Ar'); mass=39.948d0;   case('K');  mass=39.0983d0
      case('Ca'); mass=40.078d0;   case('Fe'); mass=55.845d0
      case('Zn'); mass=65.38d0;    case('Br'); mass=79.904d0
      case('I');  mass=126.904d0
      case default
        write(*,*) 'WARNING: unknown element "'//trim(sym)//'", mass=1'
        mass=1.0d0
    end select
  end function element_mass


  !==========================================================================
  ! Center of mass
  !==========================================================================
  subroutine center_of_mass(xyz_in, xyz_out)
    real*8, intent(in)  :: xyz_in(3, n_atoms)
    real*8, intent(out) :: xyz_out(3, n_atoms)
    real*8  :: xcm(3), total_mass
    integer :: iat
    total_mass = sum(atom_mass(1:n_atoms))
    xcm = 0.d0
    do iat = 1, n_atoms
      xcm = xcm + atom_mass(iat) * xyz_in(:,iat)
    end do
    xcm = xcm / total_mass
    do iat = 1, n_atoms
      xyz_out(:,iat) = xyz_in(:,iat) - xcm
    end do
  end subroutine center_of_mass


  !==========================================================================
  ! Principal axes via LAPACK dsyev
  !==========================================================================
  subroutine principal_axes(xyz_com, Ia, Ib, Ic, axes)
    real*8, intent(in)  :: xyz_com(3, n_atoms)
    real*8, intent(out) :: Ia, Ib, Ic, axes(3,3)
    real*8  :: IT(3,3), m, x, y, z, work(50), evals(3), A(3,3)
    real*8  :: mass_mod(n_atoms)
    integer :: iat, info

    mass_mod = 1.0d0

    IT = 0.d0
    do iat = 1, n_atoms
      m=atom_mass(iat)*mass_mod(iat)
      x=xyz_com(1,iat); y=xyz_com(2,iat); z=xyz_com(3,iat)
      IT(1,1)=IT(1,1)+m*(y*y+z*z); IT(2,2)=IT(2,2)+m*(x*x+z*z)
      IT(3,3)=IT(3,3)+m*(x*x+y*y); IT(1,2)=IT(1,2)-m*x*y
      IT(1,3)=IT(1,3)-m*x*z;       IT(2,3)=IT(2,3)-m*y*z
    end do
    IT(2,1)=IT(1,2); IT(3,1)=IT(1,3); IT(3,2)=IT(2,3)
    A = IT
    call dsyev('V','U',3,A,3,evals,work,50,info)
    if (info /= 0) write(*,*) 'WARNING principal_axes: dsyev info=',info
    Ia=evals(1); Ib=evals(2); Ic=evals(3); axes=A

  end subroutine principal_axes


  !==========================================================================
  ! Compute principal axis frame
  !==========================================================================
  subroutine compute_principal_frame()
    real*8  :: xyz_com(3,n_atoms), Ia, Ib, Ic
    integer :: iat

    if (allocated(xyz_pa_global)) deallocate(xyz_pa_global)
    allocate(xyz_pa_global(3,n_atoms))

    call center_of_mass(atom_xyz, xyz_com)
    call principal_axes(xyz_com, Ia, Ib, Ic, pa_axes)

    write(*,'(A,3F12.4)') ' Principal moments (amu*Ang^2): ', Ia, Ib, Ic
    write(*,'(A)') ' PA axes (cols=a,b,c; rows=lab x,y,z):'
    write(*,'(A,3F10.5)') '   x: ', pa_axes(1,:)
    write(*,'(A,3F10.5)') '   y: ', pa_axes(2,:)
    write(*,'(A,3F10.5)') '   z: ', pa_axes(3,:)

    do iat = 1, n_atoms
      xyz_pa_global(1,iat) = dot_product(pa_axes(:,1), xyz_com(:,iat))
      xyz_pa_global(2,iat) = dot_product(pa_axes(:,2), xyz_com(:,iat))
      xyz_pa_global(3,iat) = dot_product(pa_axes(:,3), xyz_com(:,iat))
    end do

    write(*,'(A)') ' Atoms in PA frame:'
    write(*,'(4A12)') 'symbol','    a','    b','    c'
    do iat = 1, n_atoms
      write(*,'(1A12, 3F12.6)') trim(atom_label(iat)), &
        xyz_pa_global(1,iat), xyz_pa_global(2,iat), xyz_pa_global(3,iat)
    end do

  end subroutine compute_principal_frame


  !==========================================================================
  ! Build atom image table  (symop_mat is already in PA frame)
  !==========================================================================
  subroutine build_atom_image_table()
    integer :: iat, iop, kat, best_kat, n_missing
    real*8  :: new_pos(3), d, best_d

    if (allocated(atom_image)) deallocate(atom_image)
    allocate(atom_image(n_atoms, n_symops))
    atom_image = 0

    write(*,'(A)') ' '
    write(*,'(A)') ' === Atom permutation table ==='
    write(*,'(A,F8.4,A)') ' (proj_cutoff = ', proj_cutoff_sym, ' Ang)'

    do iop = 1, n_symops
      n_missing = 0

      write(*,'(A)') ' '
      write(*,'(A,I2,A,A)') ' Symop ', iop, ': ', trim(symop_name(iop))
      write(*,'(A)') '   PA-frame matrix:'
      write(*,'(A,3F9.4)') '     ', symop_mat(1,:,iop)
      write(*,'(A,3F9.4)') '     ', symop_mat(2,:,iop)
      write(*,'(A,3F9.4)') '     ', symop_mat(3,:,iop)
      write(*,'(A)') ' '
      write(*,'(A)') &
        '   atom  sym    new_a        new_b        new_c     ' &
        //'  ->img   dist       status'

      do iat = 1, n_atoms
        ! Apply PA-frame symop to PA-frame coordinates
        new_pos(1) = symop_mat(1,1,iop)*xyz_pa_global(1,iat) &
                   + symop_mat(1,2,iop)*xyz_pa_global(2,iat) &
                   + symop_mat(1,3,iop)*xyz_pa_global(3,iat)
        new_pos(2) = symop_mat(2,1,iop)*xyz_pa_global(1,iat) &
                   + symop_mat(2,2,iop)*xyz_pa_global(2,iat) &
                   + symop_mat(2,3,iop)*xyz_pa_global(3,iat)
        new_pos(3) = symop_mat(3,1,iop)*xyz_pa_global(1,iat) &
                   + symop_mat(3,2,iop)*xyz_pa_global(2,iat) &
                   + symop_mat(3,3,iop)*xyz_pa_global(3,iat)

        best_d = 1.d10; best_kat = 0
        do kat = 1, n_atoms
          if (atom_Z(kat) /= atom_Z(iat)) cycle
          d = sqrt( (new_pos(1)-xyz_pa_global(1,kat))**2 &
                  + (new_pos(2)-xyz_pa_global(2,kat))**2 &
                  + (new_pos(3)-xyz_pa_global(3,kat))**2 )
          if (d < best_d) then; best_d=d; best_kat=kat; end if
        end do

        if (best_d < proj_cutoff_sym) then
          atom_image(iat,iop) = best_kat
        else
          atom_image(iat,iop) = 0
          n_missing = n_missing + 1
        end if

        if (atom_image(iat,iop) == iat) then
          write(*,'(A,I4,1X,A4,3F13.6,I6,F11.6,A)') &
            '   ', iat, trim(atom_label(iat)), &
            new_pos, best_kat, best_d, '  SELF'
        else if (atom_image(iat,iop) > 0) then
          write(*,'(A,I4,1X,A4,3F13.6,I6,F11.6,A)') &
            '   ', iat, trim(atom_label(iat)), &
            new_pos, best_kat, best_d, '  ->perm'
        else
          write(*,'(A,I4,1X,A4,3F13.6,I6,F11.6,A)') &
            '   ', iat, trim(atom_label(iat)), &
            new_pos, best_kat, best_d, '  *** NO IMAGE ***'
        end if

      end do

      if (n_missing > 0) then
        write(*,'(A,I3,A,A,A)') &
          '   *** WARNING: ', n_missing, &
          ' atoms have no image under ', trim(symop_name(iop)), &
          ' ***'
      end if

    end do

    write(*,'(A)') ' '

  end subroutine build_atom_image_table


  !==========================================================================
  ! Count mapped atoms  (used by orientation-detection routines)
  ! The trial matrix R is already in PA coordinates.
  !==========================================================================
  function count_mapped(R) result(n)
    real*8, intent(in) :: R(3,3)
    integer :: n, iat, kat
    real*8  :: new_pos(3), d

    n = 0
    do iat = 1, n_atoms
      new_pos(1)=R(1,1)*xyz_pa_global(1,iat)+R(1,2)*xyz_pa_global(2,iat) &
                +R(1,3)*xyz_pa_global(3,iat)
      new_pos(2)=R(2,1)*xyz_pa_global(1,iat)+R(2,2)*xyz_pa_global(2,iat) &
                +R(2,3)*xyz_pa_global(3,iat)
      new_pos(3)=R(3,1)*xyz_pa_global(1,iat)+R(3,2)*xyz_pa_global(2,iat) &
                +R(3,3)*xyz_pa_global(3,iat)
      do kat = 1, n_atoms
        if (atom_Z(kat) /= atom_Z(iat)) cycle
        d = sqrt( (new_pos(1)-xyz_pa_global(1,kat))**2 &
                + (new_pos(2)-xyz_pa_global(2,kat))**2 &
                + (new_pos(3)-xyz_pa_global(3,kat))**2 )
        if (d < proj_cutoff_sym) then; n=n+1; exit; end if
      end do
    end do

  end function count_mapped


  subroutine make_C2(k,R)
    integer,intent(in)::k; real*8,intent(out)::R(3,3); real*8::u(3)
    u=0.d0;u(k)=1.d0;call rot_mat(u,pi_sym,R)
  end subroutine make_C2

  subroutine make_refl(k,R)
    integer,intent(in)::k; real*8,intent(out)::R(3,3); real*8::nv(3)
    nv=0.d0;nv(k)=1.d0;call refl_mat(nv,R)
  end subroutine make_refl

  subroutine rot_mat(u,theta,R)
    real*8,intent(in)::u(3),theta; real*8,intent(out)::R(3,3)
    real*8::c,s,t,ux,uy,uz
    c=cos(theta);s=sin(theta);t=1.d0-c
    ux=u(1);uy=u(2);uz=u(3)
    R(1,1)=t*ux*ux+c;    R(1,2)=t*ux*uy-s*uz; R(1,3)=t*ux*uz+s*uy
    R(2,1)=t*ux*uy+s*uz; R(2,2)=t*uy*uy+c;    R(2,3)=t*uy*uz-s*ux
    R(3,1)=t*ux*uz-s*uy; R(3,2)=t*uy*uz+s*ux; R(3,3)=t*uz*uz+c
  end subroutine rot_mat

  subroutine refl_mat(nv,R)
    real*8,intent(in)::nv(3); real*8,intent(out)::R(3,3); integer::i,j
    do i=1,3;do j=1,3
      R(i,j)=-2.d0*nv(i)*nv(j);if(i==j)R(i,j)=R(i,j)+1.d0
    end do;end do
  end subroutine refl_mat

  subroutine alloc_group(nirr,nsym)
    integer,intent(in)::nirr,nsym
    if(allocated(irrep_names))   deallocate(irrep_names)
    if(allocated(char_table))    deallocate(char_table)
    if(allocated(symop_mat))     deallocate(symop_mat)
    if(allocated(symop_name))    deallocate(symop_name)
    if(allocated(irrep_product)) deallocate(irrep_product)
    allocate(irrep_names(nirr), char_table(nirr,nsym))
    allocate(symop_mat(3,3,nsym), symop_name(nsym))
    char_table=0.d0; symop_mat=0.d0
    symop_name=' '; irrep_names=' '
  end subroutine alloc_group

  subroutine set_E(k)
    integer,intent(in)::k
    symop_mat(:,:,k)=0.d0
    symop_mat(1,1,k)=1.d0; symop_mat(2,2,k)=1.d0; symop_mat(3,3,k)=1.d0
  end subroutine set_E

  function ax_label(k) result(c)
    integer,intent(in)::k; character(len=1)::c
    select case(k)
      case(1); c='a'; case(2); c='b'; case(3); c='c'
      case default; c='?'
    end select
  end function ax_label

  subroutine other_two_axes(main_ax,ax1,ax2)
    integer,intent(in)::main_ax; integer,intent(out)::ax1,ax2
    select case(main_ax)
      case(1); ax1=2; ax2=3
      case(2); ax1=1; ax2=3
      case(3); ax1=1; ax2=2
    end select
  end subroutine other_two_axes


  !==========================================================================
  ! Orientation detection
  !==========================================================================
  subroutine detect_c2v_orientation(c2ax,sv1ax,sv2ax)
    integer,intent(out)::c2ax,sv1ax,sv2ax
    integer::try_c2,ax1,ax2,sc_c2,sc_sv1,sc_sv2,total,best
    real*8::R_c2(3,3),R_sv1(3,3),R_sv2(3,3)
    best=-1; c2ax=3; sv1ax=1; sv2ax=2
    do try_c2=1,3
      call other_two_axes(try_c2,ax1,ax2)
      call make_C2(try_c2,R_c2); call make_refl(ax1,R_sv1); call make_refl(ax2,R_sv2)
      sc_c2=count_mapped(R_c2); sc_sv1=count_mapped(R_sv1); sc_sv2=count_mapped(R_sv2)
      total=sc_c2+sc_sv1+sc_sv2
      write(*,'(A,I2,A,I2,A,I2,A,3I4,A,I5)') &
        '  C2v trial C2=ax',try_c2,' sv_n=ax',ax1,' sv_n=ax',ax2, &
        '  mapped(C2,sv1,sv2)=',sc_c2,sc_sv1,sc_sv2,'  total=',total
      if(total>best)then; best=total; c2ax=try_c2; sv1ax=ax1; sv2ax=ax2; end if
    end do
    write(*,'(A,I2,A,I2,A,I2)') &
      ' C2v chosen: C2=ax',c2ax,'  sv normals: ax',sv1ax,' ax',sv2ax
  end subroutine detect_c2v_orientation

  subroutine detect_cs_orientation(norm_ax)
    integer,intent(out)::norm_ax; integer::try_ax,sc,best; real*8::R(3,3)
    best=-1; norm_ax=3
    do try_ax=1,3
      call make_refl(try_ax,R); sc=count_mapped(R)
      write(*,'(A,I2,A,I4)') '  Cs norm_ax=',try_ax,' mapped=',sc
      if(sc>best)then; best=sc; norm_ax=try_ax; end if
    end do
    write(*,'(A,I2)') ' Cs: mirror normal=ax',norm_ax
  end subroutine detect_cs_orientation

  subroutine detect_c2_orientation(c2ax)
    integer,intent(out)::c2ax; integer::try_ax,sc,best; real*8::R(3,3)
    best=-1; c2ax=3
    do try_ax=1,3
      call make_C2(try_ax,R); sc=count_mapped(R)
      write(*,'(A,I2,A,I4)') '  C2 ax=',try_ax,' mapped=',sc
      if(sc>best)then; best=sc; c2ax=try_ax; end if
    end do
    write(*,'(A,I2)') ' C2: axis=ax',c2ax
  end subroutine detect_c2_orientation

  subroutine detect_c2h_orientation(c2ax)
    integer,intent(out)::c2ax
    integer::try_ax,sc_c2,sc_sh,total,best; real*8::R_c2(3,3),R_sh(3,3)
    best=-1; c2ax=3
    do try_ax=1,3
      call make_C2(try_ax,R_c2); call make_refl(try_ax,R_sh)
      sc_c2=count_mapped(R_c2); sc_sh=count_mapped(R_sh); total=sc_c2+sc_sh
      write(*,'(A,I2,A,2I4,A,I5)') &
        '  C2h ax=',try_ax,' mapped(C2,sh)=',sc_c2,sc_sh,' total=',total
      if(total>best)then; best=total; c2ax=try_ax; end if
    end do
    write(*,'(A,I2)') ' C2h: axis=ax',c2ax
  end subroutine detect_c2h_orientation


  !==========================================================================
  ! GROUP SETUP  (all matrices are directly in PA frame)
  !==========================================================================
  subroutine setup_C1()
    point_group_name='C1'; n_irreps=1; n_symops=1; call alloc_group(1,1)
    irrep_names(1)='A'; char_table(1,1)=1.d0
    call set_E(1); symop_name(1)='E'; call build_product_table()
  end subroutine setup_C1

  subroutine setup_Cs()
    integer::norm_ax; real*8::R(3,3)
    call detect_cs_orientation(norm_ax)
    point_group_name='Cs'; n_irreps=2; n_symops=2; call alloc_group(2,2)
    irrep_names(1)="A'"; irrep_names(2)='A"'
    char_table(1,:)=[1.d0,1.d0]; char_table(2,:)=[1.d0,-1.d0]
    call set_E(1); symop_name(1)='E'
    call make_refl(norm_ax,R); symop_mat(:,:,2)=R
    symop_name(2)='sh_'//ax_label(norm_ax); call build_product_table()
  end subroutine setup_Cs

  subroutine setup_Ci()
    real*8::Ri(3,3)
    point_group_name='Ci'; n_irreps=2; n_symops=2; call alloc_group(2,2)
    irrep_names(1)='Ag'; irrep_names(2)='Au'
    char_table(1,:)=[1.d0,1.d0]; char_table(2,:)=[1.d0,-1.d0]
    call set_E(1); symop_name(1)='E'
    Ri=0.d0; Ri(1,1)=-1.d0; Ri(2,2)=-1.d0; Ri(3,3)=-1.d0
    symop_mat(:,:,2)=Ri; symop_name(2)='i'; call build_product_table()
  end subroutine setup_Ci

  subroutine setup_C2()
    integer::c2ax; real*8::R(3,3)
    call detect_c2_orientation(c2ax)
    point_group_name='C2'; n_irreps=2; n_symops=2; call alloc_group(2,2)
    irrep_names(1)='A'; irrep_names(2)='B'
    char_table(1,:)=[1.d0,1.d0]; char_table(2,:)=[1.d0,-1.d0]
    call set_E(1); symop_name(1)='E'
    call make_C2(c2ax,R); symop_mat(:,:,2)=R
    symop_name(2)='C2_'//ax_label(c2ax); call build_product_table()
  end subroutine setup_C2

  subroutine setup_C2h()
    integer::c2ax; real*8::R(3,3),Ri(3,3),Rsh(3,3)
    call detect_c2h_orientation(c2ax)
    point_group_name='C2h'; n_irreps=4; n_symops=4; call alloc_group(4,4)
    irrep_names(1)='Ag'; irrep_names(2)='Bg'
    irrep_names(3)='Au'; irrep_names(4)='Bu'
    char_table(1,:)=[1.d0, 1.d0, 1.d0, 1.d0]
    char_table(2,:)=[1.d0,-1.d0, 1.d0,-1.d0]
    char_table(3,:)=[1.d0, 1.d0,-1.d0,-1.d0]
    char_table(4,:)=[1.d0,-1.d0,-1.d0, 1.d0]
    call set_E(1); symop_name(1)='E'
    call make_C2(c2ax,R); symop_mat(:,:,2)=R; symop_name(2)='C2_'//ax_label(c2ax)
    Ri=0.d0; Ri(1,1)=-1.d0; Ri(2,2)=-1.d0; Ri(3,3)=-1.d0
    symop_mat(:,:,3)=Ri; symop_name(3)='i'
    call make_refl(c2ax,Rsh); symop_mat(:,:,4)=Rsh
    symop_name(4)='sh_'//ax_label(c2ax); call build_product_table()
  end subroutine setup_C2h

  subroutine setup_C2v()
    integer::c2ax,sv1ax,sv2ax; real*8::R(3,3)
    call detect_c2v_orientation(c2ax,sv1ax,sv2ax)
    point_group_name='C2v'; n_irreps=4; n_symops=4; call alloc_group(4,4)
    irrep_names(1)='A1'; irrep_names(2)='A2'
    irrep_names(3)='B1'; irrep_names(4)='B2'
    char_table(1,:)=[1.d0, 1.d0, 1.d0, 1.d0]
    char_table(2,:)=[1.d0, 1.d0,-1.d0,-1.d0]
    char_table(3,:)=[1.d0,-1.d0, 1.d0,-1.d0]
    char_table(4,:)=[1.d0,-1.d0,-1.d0, 1.d0]
    call set_E(1); symop_name(1)='E'
    call make_C2(c2ax,R); symop_mat(:,:,2)=R; symop_name(2)='C2_'//ax_label(c2ax)
    call make_refl(sv1ax,R); symop_mat(:,:,3)=R; symop_name(3)='sv_'//ax_label(sv1ax)
    call make_refl(sv2ax,R); symop_mat(:,:,4)=R; symop_name(4)='sv_'//ax_label(sv2ax)
    call build_product_table()
  end subroutine setup_C2v

  subroutine setup_D2()
    real*8::R(3,3),u(3)
    point_group_name='D2'; n_irreps=4; n_symops=4; call alloc_group(4,4)
    irrep_names(1)='A '; irrep_names(2)='B1'
    irrep_names(3)='B2'; irrep_names(4)='B3'
    char_table(1,:)=[1.d0, 1.d0, 1.d0, 1.d0]
    char_table(2,:)=[1.d0, 1.d0,-1.d0,-1.d0]
    char_table(3,:)=[1.d0,-1.d0, 1.d0,-1.d0]
    char_table(4,:)=[1.d0,-1.d0,-1.d0, 1.d0]
    call set_E(1); symop_name(1)='E'
    u=[0.d0,0.d0,1.d0]; call rot_mat(u,pi_sym,R)
    symop_mat(:,:,2)=R; symop_name(2)='C2c'
    u=[0.d0,1.d0,0.d0]; call rot_mat(u,pi_sym,R)
    symop_mat(:,:,3)=R; symop_name(3)='C2b'
    u=[1.d0,0.d0,0.d0]; call rot_mat(u,pi_sym,R)
    symop_mat(:,:,4)=R; symop_name(4)='C2a'
    call build_product_table()
  end subroutine setup_D2

  subroutine setup_D2h()
    real*8::R(3,3),Ri(3,3),nv(3),u(3)
    point_group_name='D2h'; n_irreps=8; n_symops=8; call alloc_group(8,8)
    irrep_names(1)='Ag '; irrep_names(2)='B1g'
    irrep_names(3)='B2g'; irrep_names(4)='B3g'
    irrep_names(5)='Au '; irrep_names(6)='B1u'
    irrep_names(7)='B2u'; irrep_names(8)='B3u'
    char_table(1,:)=[1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0]
    char_table(2,:)=[1.d0, 1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0]
    char_table(3,:)=[1.d0,-1.d0, 1.d0,-1.d0, 1.d0,-1.d0, 1.d0,-1.d0]
    char_table(4,:)=[1.d0,-1.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0, 1.d0]
    char_table(5,:)=[1.d0, 1.d0, 1.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0]
    char_table(6,:)=[1.d0, 1.d0,-1.d0,-1.d0,-1.d0,-1.d0, 1.d0, 1.d0]
    char_table(7,:)=[1.d0,-1.d0, 1.d0,-1.d0,-1.d0, 1.d0,-1.d0, 1.d0]
    char_table(8,:)=[1.d0,-1.d0,-1.d0, 1.d0,-1.d0, 1.d0, 1.d0,-1.d0]
    call set_E(1); symop_name(1)='E'
    u=[0.d0,0.d0,1.d0]; call rot_mat(u,pi_sym,R)
    symop_mat(:,:,2)=R; symop_name(2)='C2c'
    u=[0.d0,1.d0,0.d0]; call rot_mat(u,pi_sym,R)
    symop_mat(:,:,3)=R; symop_name(3)='C2b'
    u=[1.d0,0.d0,0.d0]; call rot_mat(u,pi_sym,R)
    symop_mat(:,:,4)=R; symop_name(4)='C2a'
    Ri=0.d0; Ri(1,1)=-1.d0; Ri(2,2)=-1.d0; Ri(3,3)=-1.d0
    symop_mat(:,:,5)=Ri; symop_name(5)='i'
    nv=[0.d0,0.d0,1.d0]; call refl_mat(nv,R)
    symop_mat(:,:,6)=R; symop_name(6)='shc'
    nv=[0.d0,1.d0,0.d0]; call refl_mat(nv,R)
    symop_mat(:,:,7)=R; symop_name(7)='shb'
    nv=[1.d0,0.d0,0.d0]; call refl_mat(nv,R)
    symop_mat(:,:,8)=R; symop_name(8)='sha'
    call build_product_table()
  end subroutine setup_D2h


  !==========================================================================
  ! Product table
  !==========================================================================
  subroutine build_product_table()
    integer::i,j,k,iop; real*8::prod_char(n_symops),diff; logical::found
    if(allocated(irrep_product)) deallocate(irrep_product)
    allocate(irrep_product(n_irreps,n_irreps)); irrep_product=1
    do i=1,n_irreps; do j=1,n_irreps
      do iop=1,n_symops
        prod_char(iop)=char_table(i,iop)*char_table(j,iop)
      end do
      found=.false.
      do k=1,n_irreps
        diff=sum(abs(prod_char(1:n_symops)-char_table(k,1:n_symops)))
        if(diff<1.d-6)then; irrep_product(i,j)=k; found=.true.; exit; end if
      end do
      if(.not.found) irrep_product(i,j)=1
    end do; end do
  end subroutine build_product_table


  !==========================================================================
  ! Assign irreps  — symop_mat is already in PA frame
  !==========================================================================
  subroutine assign_mode_irreps_all()
    real*8  :: L_pa(3,n_atoms), disp_lab(3), RL(3)
    real*8  :: numerator, denominator, diff_char, diff_pos
    real*8  :: char_vec(n_symops)
    integer :: imode, iat, iop, irr, best_irr, jat

    if (allocated(mode_irrep)) deallocate(mode_irrep)
    allocate(mode_irrep(n_modes_file))
    if (allocated(mode_char)) deallocate(mode_char)
    allocate(mode_char(n_modes_file, n_symops))
    mode_char = 0.d0; mode_irrep = 1
    if (n_irreps == 1) return

    write(*,'(A,F8.4,A)') ' assign_mode_irreps: proj_cutoff=', &
      proj_cutoff_sym, ' Ang'

    do imode = 1, n_modes_file

      ! Rotate displacements to PA frame
      do iat = 1, n_atoms
        disp_lab(1) = mode_disp(3*(iat-1)+1, imode)
        disp_lab(2) = mode_disp(3*(iat-1)+2, imode)
        disp_lab(3) = mode_disp(3*(iat-1)+3, imode)
        L_pa(1,iat) = dot_product(pa_axes(:,1), disp_lab)
        L_pa(2,iat) = dot_product(pa_axes(:,2), disp_lab)
        L_pa(3,iat) = dot_product(pa_axes(:,3), disp_lab)
      end do

      denominator = sum(L_pa**2)

      if (denominator < 1.d-20) then
        mode_irrep(imode) = 1
        mode_char(imode, :) = 1.d0
        cycle
      end if

      do iop = 1, n_symops
        numerator = 0.d0

        do iat = 1, n_atoms
          jat = atom_image(iat, iop)
          if (jat == 0) cycle

          ! PA-frame symop matrix applied to PA-frame displacement
          RL(1) = symop_mat(1,1,iop)*L_pa(1,iat) &
                + symop_mat(1,2,iop)*L_pa(2,iat) &
                + symop_mat(1,3,iop)*L_pa(3,iat)
          RL(2) = symop_mat(2,1,iop)*L_pa(1,iat) &
                + symop_mat(2,2,iop)*L_pa(2,iat) &
                + symop_mat(2,3,iop)*L_pa(3,iat)
          RL(3) = symop_mat(3,1,iop)*L_pa(1,iat) &
                + symop_mat(3,2,iop)*L_pa(2,iat) &
                + symop_mat(3,3,iop)*L_pa(3,iat)

          numerator = numerator &
            + L_pa(1,jat)*RL(1) &
            + L_pa(2,jat)*RL(2) &
            + L_pa(3,jat)*RL(3)
        end do

        char_vec(iop) = numerator / denominator

      end do

      mode_char(imode, 1:n_symops) = char_vec(1:n_symops)

      best_irr = 1; diff_char = 1.d10
      do irr = 1, n_irreps
        diff_pos = sum(abs(char_vec(1:n_symops) - char_table(irr, 1:n_symops)))
        if (diff_pos < diff_char) then
          diff_char = diff_pos; best_irr = irr
        end if
      end do
      mode_irrep(imode) = best_irr

    end do

  end subroutine assign_mode_irreps_all

    !==========================================================================
  ! Replace geometry in normal_mode.txt with PA-frame coordinates
  !==========================================================================
  subroutine rewrite_normal_mode_pa_frame()
    integer :: ios, imode, iat, i
    real*8  :: disp_lab(3), disp_pa(3)

    ! First, rotate ALL displacement vectors to PA frame
    do imode = 1, n_modes_file
      do iat = 1, n_atoms
        disp_lab(1) = mode_disp(3*(iat-1)+1, imode)
        disp_lab(2) = mode_disp(3*(iat-1)+2, imode)
        disp_lab(3) = mode_disp(3*(iat-1)+3, imode)
        
        ! Apply the SAME rotation used for coordinates
        disp_pa(1) = dot_product(pa_axes(:,1), disp_lab)
        disp_pa(2) = dot_product(pa_axes(:,2), disp_lab)
        disp_pa(3) = dot_product(pa_axes(:,3), disp_lab)
        
        mode_disp(3*(iat-1)+1, imode) = disp_pa(1)
        mode_disp(3*(iat-1)+2, imode) = disp_pa(2)
        mode_disp(3*(iat-1)+3, imode) = disp_pa(3)
      end do
    end do

    ! Now write the file with BOTH rotated
    open(unit=57, file='normal_mode.txt', status='replace', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: cannot open normal_mode.txt for writing'; return
    end if

    write(57,'(I12)') n_atoms

    do iat = 1, n_atoms
      write(57,'(10X,A4,3F13.6)') &
        atom_label(iat), &
        xyz_pa_global(1,iat), xyz_pa_global(2,iat), xyz_pa_global(3,iat)
    end do

    do imode = 1, n_modes_file
      write(57,'(A,I13)') ' Mode', imode
      write(57,'(F14.6)') mode_freq(imode)
      do i = 1, n_modes_file
        write(57,'(F14.8)') mode_disp(i, imode)
      end do
    end do

    close(57)
    write(*,'(A)') ' normal_mode.txt rewritten with PA-frame geometry AND displacements.'

end subroutine rewrite_normal_mode_pa_frame

end module symmetry_module