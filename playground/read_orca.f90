
module read_orca_file

use jacobi_diagonalization

contains

subroutine read_orca(file, linear, vibrations, cubic, quartic, N_modes, save_dip, save_second_dipole)
  
  implicit none
  
  ! Input
  character(len=200), intent(in) :: file
  
  integer:: n_atoms, n_coords
  real*8, allocatable :: atomic_masses(:)
  character(len=2), allocatable :: atomic_symbols(:)
  real*8, allocatable :: atomic_coords(:,:), dipole_derivatives(:,:), dipole_derivatives_rotate(:,:), dipole_derivatives_transpose(:,:), transp(:,:)
  real*8, allocatable :: second_dipole(:, :, :), second_dipole_rotate(:, :, :)
  real*8, allocatable :: hessian(:,:)
  real*8, allocatable :: cubic(:,:,:)
  real*8, allocatable :: quartic(:,:,:,:)
  real*8, allocatable :: eigenvectors_mass(:,:)
  real*8, allocatable :: intensities(:)

  
  integer :: elements_quadratic, elements_cubic, i, j, k, l, io, idx, linear
  integer :: index1, index2, index3, index4, N_modes, atom_k, mode_idx
  real*8 :: temp_val, intensity_value
  character(len=200) :: line
  logical :: found_hessian
  real*8, allocatable :: frequencies(:), vibrations(:)
  real*8, allocatable :: eigenvectors(:,:), sqrt_mass(:,:)
  logical :: is_linear
  real*8 :: vc, eps_0, na, cte, h_bar, cm_to_si, dip_au_to_si
  real*8 :: save_dip(N_modes,3), save_second_dipole(N_modes,N_modes,3)

vc = 299792458.d0
eps_0 = 8.8541878128/(10**(12.d0))
na = 6.02214076 * 10**23.d0
h_bar = 1.054571817/(10**(34.d0))
cm_to_si = 1.986445857/(10**23.d0)
dip_au_to_si =  3.335640951/(10**(30.d0))


!vc = 299792458.0 m/s (speed of light)
!eps_0 = 8.8541878128 × 10⁻¹² F/m (vacuum permittivity)
!na = 6.02214076 × 10²³ mol⁻¹ (Avogadro's number)
!h_bar = 1.054571817 × 10⁻³⁴ J·s (reduced Planck constant)
!cm_to_si = 1.986445857 × 10⁻²³ J (converts cm⁻¹ to Joules)
!dip_au_to_si = 3.335640951 × 10⁻³⁰ C·m (converts atomic units of dipole to SI)

cte = na*cm_to_si*dip_au_to_si**2/(1000.*6.*(vc**2)*eps_0*h_bar**2)
cte = 1.d0
  is_linear = .false.

  if (linear == 1) then
    is_linear = .true.
  end if

  ! Initialize
  found_hessian = .false.
  n_atoms = 0
  n_coords = 0

  write(*, *) 
  write(*, '(A)') '----------------------------------------'
  write(*, '(A)') '        Reading ORCA .vpt2 file         '
  write(*, '(A)') '----------------------------------------'

  open(101, file=trim(file), status='old', iostat=io)
  if (io /= 0) then
    write(*,*) 'ERROR: Cannot open file ', trim(file)
    stop
  end if

  do 
    read(101, '(A)', iostat=io) line
    if (io /= 0) exit

    ! Read atomic coordinates and masses
    if (index(line, '# Atomic coordinates in Angstroem') > 0) then
      read(101, *) n_atoms
      allocate(atomic_coords(3, n_atoms))
      allocate(atomic_masses(n_atoms))
      allocate(atomic_symbols(n_atoms))
      allocate(dipole_derivatives(n_atoms*3,3))
      allocate(dipole_derivatives_rotate(n_atoms*3,3))
      allocate(intensities(n_atoms*3))
      allocate(transp(3,n_atoms*3))
      allocate(sqrt_mass(3*n_atoms, 3*n_atoms))
      allocate(eigenvectors_mass(3*n_atoms, 3*n_atoms))
      allocate(second_dipole(N_modes, 3*n_atoms, 3), second_dipole_rotate(N_modes, 3*n_atoms, 3))
      dipole_derivatives = 0.d0

      do i = 1, n_atoms
        read(101, *) atomic_symbols(i), idx, atomic_masses(i), &
                    atomic_coords(1,i), atomic_coords(2,i), atomic_coords(3,i)
        if (linear == 1) then
          if (N_modes /= 3*n_atoms - 5) then
            write(*,*)
            write(*,*) 'The number of modes in the input file is inconsistent. Expected: ', 3*n_atoms - 5
            write(*,*)
            stop
          end if
      else
          if (N_modes /= 3*n_atoms - 6) then
            write(*,*)
            write(*,*) 'The number of modes in the input file is inconsistent. Expected: ', 3*n_atoms - 6
            write(*,*)
            stop
          end if
      end if
      end do
    !  write(*,'(1A30, 1I6, 1A20)') 'Found ', n_atoms, ' atoms'
    end if

    ! Read Hessian
    if (trim(line) .eq. '# Hessian[i][j] in Eh/(bohr**2)') then
      read(101,*) elements_quadratic
      n_coords = elements_quadratic
      
      if (.not. allocated(hessian)) then
        allocate(hessian(n_coords, n_coords))
        hessian = 0.0d0
      end if
      
      do i = 1, elements_quadratic * elements_quadratic
        read(101, *) index1, index2, hessian(index1+1, index2+1)
      end do
      found_hessian = .true.
      write(*,'(1A30)') 'Reading Hessian matrix'
    end if
    
    ! Read Cubic force constants
    if (index(line, '# Cubic[i][j][k] force field in 1/cm') > 0) then
      read(101,*) elements_cubic
      if (.not. allocated(cubic)) then
        allocate(cubic(elements_cubic, elements_cubic, elements_cubic))
        cubic = 0.0d0
      end if
      
      do i = 1, elements_cubic**3
        read(101, *) index1, index2, index3, temp_val
        cubic(index1+1, index2+1, index3+1) = temp_val
      end do
      write(*,'(1A30)') 'Reading cubic force constants'
    end if
    
    ! Read Semi-quartic force constants - with 4 indices
    if (index(line, '# Semi-quartic[i][j][k][k] force field in 1/cm') > 0 .or. &
        index(line, '# Quartic[i][j][k][l] force field in 1/cm') > 0) then
      read(101,*) elements_cubic !just to recycle variable 
      if (.not. allocated(quartic)) then
        allocate(quartic(elements_cubic, elements_cubic, elements_cubic, elements_cubic))
        quartic = 0.0d0
      end if
      
      ! Read until we've read all elements
      do i = 1, elements_cubic**4
        read(101, *, iostat=io) index1, index2, index3, temp_val
        if (io /= 0) exit
        quartic(index1+1, index2+1, index3+1, index3+1) = temp_val
      end do
      write(*,'(1A30)') 'Reading quartic force constants'
    end if

    ! Read dipole derivatives
    if (index(line, '# Dipole Derivatives[i][j] in (Eh*bohr)^1/2') > 0) then
      read(101,*) 
      
      do i = 1, n_atoms*3*3
        read(101, *) index1, index2, temp_val
        dipole_derivatives(index1+1, index2+1) = temp_val
      end do
      write(*,'(1A30)') 'Reading dipole derivatives'
    end if

    ! Read dipole second derivatives
    if (index(line, '# 2nd Dipole Derivatives[NVib][threeN][xyz] in (Eh*bohr)^1/2') > 0) then
      read(101,*) 
      
      do i = 1, n_atoms*3*3*N_modes
        read(101, *) index1, index2, index3, temp_val
        second_dipole(index1+1, index2+1, index3+1) = temp_val
      end do
      write(*,'(1A30)') 'Reading second dipole derivat.'
    end if

  end do

  close(101)

  
  if (.not. found_hessian) then
    write(*,*) '>> ERROR: Hessian not found in file!'
    stop
  end if
  
  write(*,*)
  write(*,*) '  >> File reading complete <<'
 ! write(*,*)

  call compute_frequencies(hessian, atomic_masses, n_atoms, frequencies, eigenvectors)

  call print_frequencies(frequencies, n_atoms, is_linear)

  if (linear == 1) then
    do i = 1, 3*n_atoms-5
      if (N_modes /= 3*n_atoms - 5) then
        write(*,*)
        write(*,*) 'The number of modes in the input file is inconsistent. Expected: ', 3*n_atoms - 5
        write(*,*)
        stop
      end if
      vibrations(i) = frequencies(i+5)
    end do
  else
    do i = 1, 3*n_atoms-6
      if (N_modes /= 3*n_atoms - 6) then
        write(*,*)
        write(*,*) 'The number of modes in the input file is inconsistent. Expected: ', 3*n_atoms - 6
        write(*,*)
        stop
      end if
      vibrations(i) = frequencies(i+6)
    end do
  end if


  ! Build mass-weighting matrix
  sqrt_mass = 0.d0
  do k = 1, n_atoms*3
    atom_k = (k-1)/3 + 1
    sqrt_mass(k,k) = 1.d0/sqrt(atomic_masses(atom_k))
  end do

  ! Build mass-weighted eigenvectors
  eigenvectors_mass = 0.d0
  do i = 1, n_atoms*3
    do k = 1, n_atoms*3
      eigenvectors_mass(k,i) = eigenvectors(k,i) * sqrt_mass(k,k)
    end do
  end do

  ! Transform dipole derivatives to normal coordinates
  dipole_derivatives_rotate = 0.d0
  second_dipole_rotate = 0.d0
  do i = 1, n_atoms*3  ! normal mode index
    do j = 1, 3        ! Cartesian direction (x,y,z)
      do k = 1, n_atoms*3  ! Cartesian coordinate index
        dipole_derivatives_rotate(i,j) = dipole_derivatives_rotate(i,j) + &
             dipole_derivatives(k,j) * eigenvectors_mass(k,i)
             do l = 1, N_modes
                second_dipole_rotate(l,i,j) = second_dipole_rotate(l,i,j) + &
                  second_dipole(l,k,j)*eigenvectors_mass(k,i)
             end do
      end do
    end do
  end do

  !saving normal modes

  open(1010, file = 'normal_mode.txt') 

  write(1010,*) n_atoms
  do i = 1, n_atoms

    write(1010,'(1A12, 3F12.6)') atomic_symbols(i), atomic_coords(1:3, i)


  end do


  do i = 1, n_atoms*3
    write(1010,*) 'Mode ', i
      write(1010,'(F12.6)') frequencies(i)
      do j = 1, n_atoms*3
        write(1010,'(F14.8)') eigenvectors(j,i)
      end do
  end do

  close(1010)

  save_dip = 0.d0
  write(*,*)
  write(*,'(A)') '--------------------------------------------'
  write(*,'(A)') '      IR Intensities (double harmonic)      '
  write(*,'(A)') '--------------------------------------------'
  write(*,'(A)') ' Mode   Frequency(cm^-1)  Intensity(arb. u.)'
  write(*,'(A)') '--------------------------------------------'
  
  ! Only calculate intensities for vibrational modes (skip trans/rot)
  if (linear == 1) then
    do i = 1, 3*n_atoms-5
      ! mode_idx = i + 5 (vibrational mode)
      mode_idx = i + 5
      ! Square of dipole derivative (sum over x,y,z components)
      intensity_value = dipole_derivatives_rotate(mode_idx,1)**2 + &
                        dipole_derivatives_rotate(mode_idx,2)**2 + &
                        dipole_derivatives_rotate(mode_idx,3)**2
      ! Multiply by conversion factor
      intensities(mode_idx) = intensity_value * cte
      save_dip(i, 1:3) = dipole_derivatives_rotate(mode_idx,1:3)
      write(*,'(I5, F15.3, 1F15.3)') i, vibrations(i), intensities(mode_idx)*974.88 !this is from orca... I got by scalling, better than guess the correct numbers.
    end do
  else
    do i = 1, 3*n_atoms-6
      ! mode_idx = i + 6 (vibrational mode)
      mode_idx = i + 6
      ! Square of dipole derivative (sum over x,y,z components)
      intensity_value = dipole_derivatives_rotate(mode_idx,1)**2 + &
                        dipole_derivatives_rotate(mode_idx,2)**2 + &
                        dipole_derivatives_rotate(mode_idx,3)**2
      ! Multiply by conversion factor
      intensities(mode_idx) = intensity_value * cte
      write(*,'(I5, F15.3, 1F15.3)') i, vibrations(i), intensities(mode_idx)*974.88
      save_dip(i, 1:3) = dipole_derivatives_rotate(mode_idx,1:3)
    end do
  end if


   ! Only calculate intensities for vibrational modes (skip trans/rot)
  if (linear == 1) then
    do i = 1, N_modes
      do j = 1, N_modes
        mode_idx = j + 5
        save_second_dipole(i, j, 1:3) = second_dipole_rotate(i,mode_idx,1:3)
      end do
    end do
  else
    do i = 1, N_modes
      do j = 1, N_modes
        mode_idx = j + 6
        save_second_dipole(i, j, 1:3) = second_dipole_rotate(i,mode_idx,1:3)
      end do
    end do
  end if

  write(200,*) 'HO'
  write(200,*) N_modes
  do i = 1, N_modes
    if (linear == 1) then
      write(200,'(2F18.6)') vibrations(i), 100*intensities(i+5)/sqrt(sum(intensities(:)**2))
    else 
      write(200,'(2F18.6)') vibrations(i), 100*intensities(i+6)/sqrt(sum(intensities(:)**2))
    end if
  end do
  
  write(*,'(A)') '--------------------------------------------'
  write(*,*)

end subroutine read_orca



subroutine compute_frequencies(hessian, masses, n_atoms, frequencies, eigenvectors)
  implicit none
  
  ! Input
  real*8, intent(in) :: hessian(:,:)
  real*8, intent(in) :: masses(:)
  integer, intent(in) :: n_atoms
  
  ! Output
  real*8, allocatable, intent(out) :: frequencies(:)
  real*8, allocatable, intent(out), optional :: eigenvectors(:,:)
  
  ! Local variables
  integer :: n_coords, i, j, k, l, iter
  real*8, allocatable :: mass_weighted_hessian(:,:)
  real*8, allocatable :: eigvals(:), eigvecs(:,:)
  real*8, parameter :: au_to_cm1 = 219474.63d0
  real*8, parameter :: amu_to_au = 1822.8884862d0
  real*8, parameter :: bohr_to_angstrom = 0.529177210903d0
  real*8 :: conversion_factor
  
  n_coords = size(hessian, 1)
  
  ! Allocate arrays
  allocate(mass_weighted_hessian(n_coords, n_coords))
  allocate(eigvals(n_coords), eigvecs(n_coords, n_coords))
  allocate(frequencies(n_coords))
  
  ! Conversion factor: from Eh/bohr² to cm⁻¹ through mass-weighting
  ! In atomic units: 
  !   1 Eh = 4.3597447222071e-18 J
  !   1 bohr = 5.29177210903e-11 m
  !   1 amu = 1.66053906660e-27 kg
  !   1 cm⁻¹ = 1.98630e-23 J
  ! But we can use the direct conversion factor:
  conversion_factor = 5140.4866d0  ! sqrt(1/(1amu in au)) * au_to_cm1 / bohr_to_angstrom
  
  
  ! Construct mass-weighted Hessian
  do i = 1, n_coords
    do j = 1, n_coords
      k = (i-1)/3 + 1  ! atom index for coordinate i
      l = (j-1)/3 + 1  ! atom index for coordinate j
      
      ! Correct mass-weighting: 
      ! H_mw(i,j) = H(i,j) / sqrt(m_i * m_j) * conversion_factor²
      mass_weighted_hessian(i,j) = hessian(i,j) / sqrt(masses(k) * masses(l)) * conversion_factor**2
    end do
  end do
  
  ! Diagonalize
  call dsyevd_A(mass_weighted_hessian, eigvals, eigvecs)
  
  ! The eigenvalues are now in (cm⁻¹)², so take square root
  do i = 1, n_coords
    if (eigvals(i) > 0.0d0) then
      frequencies(i) = sqrt(eigvals(i))
    else if (eigvals(i) < 0.0d0) then
      frequencies(i) = -sqrt(abs(eigvals(i)))
    else
      frequencies(i) = 0.0d0
    end if
  end do
  
  ! Return eigenvectors if requested
  if (present(eigenvectors)) then
    allocate(eigenvectors(n_coords, n_coords))
    eigenvectors = eigvecs
  end if
  
  deallocate(mass_weighted_hessian, eigvals, eigvecs)
  
end subroutine compute_frequencies

subroutine print_frequencies(frequencies, n_atoms, is_linear, unit)
  implicit none
  
  ! Input
  real*8, intent(in) :: frequencies(:)
  integer, intent(in) :: n_atoms
  logical, intent(in), optional :: is_linear
  integer, intent(in), optional :: unit
  
  ! Local variables
  integer :: i, n_coords, out_unit, start_idx
  logical :: linear
  
  n_coords = size(frequencies)
  out_unit = 6
  if (present(unit)) out_unit = unit
  
  linear = .false.
  if (present(is_linear)) linear = is_linear
  
  write(out_unit, *)
  write(out_unit, '(A)') '----------------------------------------'
  write(out_unit, '(A)') '     Normal Mode Frequencies (cm^-1)    '
  write(out_unit, '(A)') '----------------------------------------'
  write(out_unit, *)
  
  if (linear) then
    start_idx = 6
    write(out_unit, '(A)') 'Vibrational frequencies (excluding 5 trans/rot modes):'
  else
    start_idx = 7
    write(out_unit, '(A)') 'Vibrational frequencies (excluding 6 trans/rot modes):'
  end if
  write(out_unit, *)
  
  do i = start_idx, n_coords
    if (linear) then
      write(out_unit, '(I4, F12.3)') i-5, frequencies(i)
    else
      write(out_unit, '(I4, F12.3)') i-6, frequencies(i)
    end if
  end do
  
  write(out_unit, *)
  write(out_unit, '(A)') 'All eigenvalues (including trans/rot):'
  write(out_unit, *)
  do i = 1, n_coords
    if ((.not. linear .and. i <= 6) .or. (linear .and. i <= 5)) then
      write(out_unit, '(I4, F12.3, A)') i, frequencies(i), ' (trans/rot)'
    else
      write(out_unit, '(I4, F12.3)') i, frequencies(i)
    end if
  end do
  
  write(*, *) '------'
  write(*, '(A, F12.3)') ' HO ZPE: ', sum(frequencies(start_idx:n_coords))/2.d0

end subroutine print_frequencies


end module