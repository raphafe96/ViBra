program main_vscf
  use jacobi_diagonalization
  use compute_integrals
  use one_mode
  use read_input_file
  use read_orca_file
  use combination
  use vib_ci
  use symmetry_module
  use mkl_service
  use omp_lib
  implicit none

  !===========================================================================
  ! Scalar integers
  !===========================================================================
  integer :: N_expansion, N_modes, i, j, k, l, mu, nu, working_mode, p
  integer :: total_combinations, N_quanta, N_states, m, n, contar, test, total_combinations2
  integer :: state_number
  integer :: conv_scf, N_threads, calculate_sci, max_iter_sci        
  integer :: use_symmetry, warning
  integer :: quanta_max_reference_sci

  !===========================================================================
  ! Mode irrep array
  !===========================================================================
  integer, allocatable :: mode_irrep_arr(:)

  !===========================================================================
  ! Combination / state arrays
  !===========================================================================
  integer, allocatable :: mode_excite(:)
  integer, allocatable :: combination_vec(:,:),combination_vec2(:,:)

  !===========================================================================
  ! VSCF arrays
  !===========================================================================
  real*8, allocatable :: HO_freq(:)
  real*8, allocatable :: Integral_sobrepos(:,:)
  real*8, allocatable :: Coeff(:,:), new_coeff(:,:)
  real*8, allocatable :: Potential_3(:,:,:)
  real*8, allocatable :: Potential_4(:,:,:,:)
  real*8, allocatable :: store_integrals(:,:,:,:)
  real*8, allocatable :: full_coef(:,:,:)
  real*8, allocatable :: dipole_derivatives(:,:)
  real*8, allocatable :: intensities_vscf(:)
  real*8, allocatable :: transition_energy_vscf(:)
  real*8, allocatable :: second_dipole_derivatives(:,:,:)
  integer, allocatable :: converged_state(:)

  !===========================================================================
  ! Inverted index arrays
  !===========================================================================
  integer :: total_3, total_4, n_cubic_max, n_quartic_max
  real*8,  allocatable :: Potential_3_vec(:), Potential_4_vec(:)
  integer, allocatable :: final_index_3(:,:), count_index_3(:,:)
  integer, allocatable :: n_unique_3(:), unique_modes_3(:,:)
  integer, allocatable :: check3(:)
  integer, allocatable :: final_index_4(:,:), count_index_4(:,:)
  integer, allocatable :: n_unique_4(:), unique_modes_4(:,:)
  integer, allocatable :: check4(:)
  integer, allocatable :: cubic_for_mode(:,:),   n_cubic_for_mode(:)
  integer, allocatable :: quartic_for_mode(:,:),  n_quartic_for_mode(:)

  !===========================================================================
  ! Scalars
  !===========================================================================
  real*8  :: cm_to_hartree, energy_ground, end_time, start_time
  real*8  :: convergence, trash, energy, new_energy, delta_e
  real*8  :: eps_sci
  real*8  :: proj_cutoff                                              

  character(len=200) :: input_file, constants_file, constants_mode
  character(len=10)  :: point_group_input 
  character(len=4) :: sci_mode                           
  integer(8) :: t_start, t_end, count_rate, count_max
  real(8)    :: elapsed_time


  !===========================================================================
  ! list of states
  !===========================================================================
  integer :: max_states


  !===========================================================================
  ! Timing
  !===========================================================================
  call system_clock(count_rate=count_rate, count_max=count_max)
  call system_clock(t_start)
  call cpu_time(start_time)



  !===========================================================================
  ! Banner
  !===========================================================================
  write(*,'(A)') '========================================'
  write(*,'(A)') ' Centro Brasileiro de Pesquisas Fisicas '
  write(*,'(A)') '       CBPF - Rio de Janeiro, Brasil'
  write(*,'(A)') '========================================'
  write(*,'(A)') '                VCI@VSCF'
  write(*,'(A)') '----------------------------------------'
  write(*,'(A)') '  Vibrational Configuration Interaction'
  write(*,'(A)') '  at'
  write(*,'(A)') '  Vibrational Self-Consistent Field'
  write(*,'(A)') '----------------------------------------'
  write(*,'(A)') '========================================'
  write(*,*)

  !===========================================================================
  ! Read input                                                        CHANGED
  !===========================================================================
  input_file = 'input_vscf.txt'
  test = 0 !this will compare the results with implementation from CRYSTAL for H2O found on their webpage (tutorials) https://tutorials.crystalsolutions.eu/tutorial.html?td=anharmonicity&tf=anharm!
  sci_mode = 'auto'

  if (test == 0) then
    call read_inp(input_file, N_modes, N_expansion, constants_file,   &
                  constants_mode, N_quanta, N_states, conv_scf,       &
                  N_threads, point_group_input, proj_cutoff,           &
                  max_iter_sci, sci_mode, quanta_max_reference_sci)                                        
  else
    N_modes         = 3
    N_expansion     = 10
    N_quanta        = 5
    N_states        = -1
    conv_scf        = 3
    N_threads       = 1
    point_group_input = 'C1'                                          
    proj_cutoff     = 0.05d0                                          
    max_iter_sci    = 0                                             
  end if

  conv_scf    = 10**conv_scf
  convergence = 1.d0 / dble(conv_scf)
  eps_sci     = 0.00005d0

  call mkl_set_num_threads(N_threads)
  call omp_set_num_threads(N_threads)

  !===========================================================================
  ! Allocate main arrays
  !===========================================================================
  allocate(HO_freq(N_modes))
  allocate(Integral_sobrepos(N_modes, N_expansion))
  allocate(Coeff(N_modes, N_expansion), new_coeff(N_modes, N_expansion))
  allocate(Potential_3(N_modes, N_modes, N_modes))
  allocate(Potential_4(N_modes, N_modes, N_modes, N_modes))
  allocate(full_coef(N_modes, N_expansion, N_expansion))
  allocate(dipole_derivatives(N_modes, 3))
  allocate(intensities_vscf(2*N_modes + N_modes*(N_modes-1)/2))
  allocate(transition_energy_vscf(2*N_modes + N_modes*(N_modes-1)/2))
  allocate(second_dipole_derivatives(N_modes, N_modes, 3))
  allocate(converged_state(2*N_modes + N_modes*(N_modes-1)/2))
  allocate(mode_irrep_arr(N_modes))

  converged_state = 0
  Potential_3     = 0.d0
  Potential_4     = 0.d0
  HO_freq         = 0.d0
  mode_irrep_arr  = 1     ! default: all totally symmetric

  open(200, file='intensities.txt')
  dipole_derivatives = 0.d0

  !===========================================================================
  ! Read force constants
  !===========================================================================
  if (test == 0) then
    call read_orca(constants_file, 0, HO_freq, Potential_3, Potential_4, &
                   N_modes, dipole_derivatives, second_dipole_derivatives)

    write(*,*) 'DIPOLES'
    write(*,*) '-------'
    do i = 1, N_modes
      write(*,'(3F12.6)') dipole_derivatives(i, 1:3)
    end do
  else
    HO_freq(1) = 1639.07d0
    HO_freq(2) = 3798.22d0
    HO_freq(3) = 3898.95d0
    Potential_3(1,1,1) = -275.7935d0
    Potential_3(2,2,2) = -1823.2787d0
    Potential_3(3,3,3) = -0.0026d0
    Potential_3(1,2,1) =  316.5947d0
    Potential_3(1,2,2) =  74.2228d0
    Potential_3(1,1,3) = -0.0001d0
    Potential_3(1,3,3) =  276.1912d0
    Potential_3(2,2,3) = -0.0001d0
    Potential_3(3,3,2) = -1827.7982d0
    Potential_3(1,2,3) = -11.7001d0
    Potential_4(1,1,1,1) = -45.1547d0
    Potential_4(2,2,2,2) =  765.7648d0
    Potential_4(3,3,3,3) =  775.3768d0
    Potential_4(2,2,1,1) = -303.1598d0
    Potential_4(1,1,3,3) = -364.7396d0
    Potential_4(2,2,3,3) =  793.6590d0
    Potential_4(1,1,1,2) =  156.5137d0
    Potential_4(1,1,1,3) =  0.0002d0
    Potential_4(1,2,2,2) = -67.5310d0
    Potential_4(2,2,2,3) =  0.0011d0
    Potential_4(1,3,3,3) =  0.0297d0
    Potential_4(2,3,3,3) = -0.0545d0
    Potential_4(1,1,2,3) =  5.3400d0
    Potential_4(1,2,2,3) = -5.3404d0
    Potential_4(1,2,3,3) = -123.7959d0
  end if

  call apply_degeneracy_factors(N_modes, Potential_3, Potential_4)

  !===========================================================================
  ! Unit conversion
  !===========================================================================
  cm_to_hartree = 0.0000045563350d0
  HO_freq     = HO_freq     * cm_to_hartree
  Potential_3 = Potential_3 * cm_to_hartree / 6.d0
  Potential_4 = Potential_4 * cm_to_hartree / 24.d0

  full_coef  = 0.d0
  Coeff(:,:) = 0.001d0

  !===========================================================================
  ! Build inverted index
  !===========================================================================
  total_3 = (N_modes+2)*(N_modes+1)*N_modes/6
  total_4 = (N_modes+3)*(N_modes+2)*(N_modes+1)*N_modes/24

  allocate(Potential_3_vec(total_3), Potential_4_vec(total_4))
  allocate(final_index_3(total_3,3),  count_index_3(total_3,3))
  allocate(n_unique_3(total_3),        unique_modes_3(total_3,3))
  allocate(check3(total_3))
  allocate(final_index_4(total_4,4),  count_index_4(total_4,4))
  allocate(n_unique_4(total_4),        unique_modes_4(total_4,4))
  allocate(check4(total_4))

  call build_inverted_index(N_modes, Potential_3, Potential_4,         &
      total_3, total_4,                                                 &
      Potential_3_vec, Potential_4_vec,                                 &
      final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
      final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
      cubic_for_mode, n_cubic_for_mode, n_cubic_max,                    &
      quartic_for_mode, n_quartic_for_mode, n_quartic_max)

  write(*,'(A)') '========================================'
  write(*,'(A,I8,A,I8,A)') ' Inverted index built: ', count(check3==1), &
      ' cubic, ', count(check4==1), ' quartic nonzero terms'
  write(*,'(A,I6,A,I6)') ' max cubic/mode=', n_cubic_max, &
      ', max quartic/mode=', n_quartic_max
  write(*,'(A)') '========================================'

  !===========================================================================
  ! Precompute 1-D integrals
  !===========================================================================
  allocate(store_integrals(N_modes, N_expansion, N_expansion, 0:5))
  allocate(mode_excite(N_modes))

  write(*,'(A)') '========================================'
  write(*,'(A)') '             STARTING VSCF              '
  write(*,'(A)') '========================================'

  do working_mode = 1, N_modes
    do nu = 1, N_expansion
      do mu = 1, N_expansion
        do p = 0, 5
          call sobrepos(nu-1, mu-1, 1.0d0, 1.0d0, &
              store_integrals(working_mode,nu,mu,p), HO_freq(working_mode), p)
        end do
      end do
    end do
  end do

  open(101, file='vscf.out')

  !===========================================================================
  ! VSCF: ground state
  !===========================================================================
  mode_excite(:) = 0
 ! write(*,'(A)') '----------------------------------------'
 ! write(*,'(1A,100I6)') 'WORKING ON STATE', mode_excite(:)
 ! write(*,'(A)') '----------------------------------------'

  write(*,'(A)') '----------------------------------------'
  write(*,'(1A)') 'WORKING ON GROUND STATE: '
  write(*,'(A)') '----------------------------------------'

  write(101,'(A)') '----------------------------------------'
  write(101,'(1A,100I6)') 'WORKING ON STATE', mode_excite(:)
  write(101,'(A)') '----------------------------------------'

  new_energy = 0.d0
  energy     = 0.d0
  Coeff(:,:) = 0.001d0

  do
    energy = new_energy
    call constant_one_mode(Coeff, Potential_3, Potential_4,            &
        N_modes, N_expansion, HO_freq, new_coeff, store_integrals,     &
        new_energy, 0, mode_excite, full_coef,                         &
        0, dipole_derivatives, trash, N_threads,                       &
        total_3, total_4,                                               &
        Potential_3_vec, Potential_4_vec,                               &
        final_index_3, count_index_3, check3,                          &
        final_index_4, count_index_4, check4)

    Coeff(:,:) = new_coeff(:,:)
    delta_e    = new_energy - energy
    write(*,'(2F20.6)') new_energy, delta_e
    write(101,'(2F20.6)') new_energy, delta_e

    if (abs(delta_e) < convergence) then
      call constant_one_mode(Coeff, Potential_3, Potential_4,          &
          N_modes, N_expansion, HO_freq, new_coeff, store_integrals,   &
          new_energy, 1, mode_excite, full_coef,                       &
          0, dipole_derivatives, trash, N_threads,                     &
          total_3, total_4,                                             &
          Potential_3_vec, Potential_4_vec,                             &
          final_index_3, count_index_3, check3,                        &
          final_index_4, count_index_4, check4)
      energy_ground = new_energy
      write(*,'(2F20.6)') new_energy, new_energy - energy
      write(*,'(1A)') '>>>>>>>> CONVERGENCE REACHED'
      write(101,'(2F20.6)') new_energy, new_energy - energy
      write(101,'(1A)') '>>>>>>>> CONVERGENCE REACHED'
      write(101,'(1A28,1F12.6)') 'STATE ENERGY (cm-1): ',      new_energy
      write(101,'(1A28,1F12.6)') 'TRANSITION ENERGY (cm-1): ', &
                                   energy - energy_ground
      exit
    end if
  end do

  write(*,'(A)') '>> Calculating VSCF fundamentals '
  !===========================================================================
  ! VSCF: single excitations
  !===========================================================================
  contar = 0
  intensities_vscf = 0.d0

  do j = 1, 1
    do i = 1, N_modes
      contar = contar + 1
      Coeff(:,:) = 0.001d0
      mode_excite    = 0
      mode_excite(i) = j

     ! write(*,*)
     ! write(*,'(A)') '----------------------------------------'
     ! write(*,'(1A,100I6)') 'WORKING ON STATE', mode_excite(:)
     ! write(*,'(A)') '----------------------------------------'
      write(101,*)
      write(101,'(A)') '----------------------------------------'
      write(101,'(1A,100I6)') 'WORKING ON STATE', mode_excite(:)
      write(101,'(A)') '----------------------------------------'

      new_energy = 0.d0
      energy     = 0.d0
      write(101,'(1A20,1A20)') 'ENERGY(cm-1)', 'DELTA E'

      do
        energy = new_energy
        call constant_one_mode(Coeff, Potential_3, Potential_4,        &
            N_modes, N_expansion, HO_freq, new_coeff, store_integrals, &
            new_energy, 0, mode_excite, full_coef,                     &
            0, dipole_derivatives, trash, N_threads,                   &
            total_3, total_4,                                           &
            Potential_3_vec, Potential_4_vec,                           &
            final_index_3, count_index_3, check3,                      &
            final_index_4, count_index_4, check4)

        Coeff(:,:) = new_coeff(:,:)
        delta_e    = new_energy - energy
        write(101,'(2F20.6)') new_energy, delta_e

        if (abs(delta_e) < convergence) then
          call constant_one_mode(Coeff, Potential_3, Potential_4,      &
              N_modes, N_expansion, HO_freq, new_coeff, store_integrals,&
              new_energy, 0, mode_excite, full_coef,                   &
              1, dipole_derivatives, intensities_vscf(contar),         &
              N_threads,                                                &
              total_3, total_4,                                         &
              Potential_3_vec, Potential_4_vec,                         &
              final_index_3, count_index_3, check3,                    &
              final_index_4, count_index_4, check4)
          converged_state(contar) = 1
          transition_energy_vscf(contar) = new_energy - energy_ground
          write(101,'(2F20.6)') new_energy, delta_e
          write(101,'(1A)') '>>>>>>>> CONVERGENCE REACHED'
         ! write(101,'(1A28,1F12.6)')  'TRANSITION ENERGY: ', &
         !                            new_energy - energy_ground
          write(101,'(1A28,1F12.6)') 'STATE ENERGY (cm-1): ',   new_energy
          write(101,'(1A28,1F12.6)') 'TRANSITION ENERGY (cm-1): ', &
                                      new_energy - energy_ground
          write(101,'(1A)') '---'
          write(101,'(1A23,1I3,1A2,1F12.6)') 'HARMONIC ENERGY', i, ': ', &
              HO_freq(i)/cm_to_hartree
          exit
        end if
      end do
    end do
  end do

  !--- Write VSCF intensities ---
  n = 0
  do i = 1, contar
    if (converged_state(i) == 1) n = n + 1
  end do

  write(200,*) 'VSCF'
  write(200,*) n
  do i = 1, contar
    if (converged_state(i) == 1) then
      write(200,'(2F18.6)') transition_energy_vscf(i), &
          100*intensities_vscf(i)/sqrt(sum(intensities_vscf(:)**2))
    end if
  end do


  write(*,'(A)') '========================================'
  write(*,'(A)') '              FINISHED VSCF             '
  write(*,'(A)') '========================================'

  write(*,'(A)') '========================================'
  write(*,'(A)') '        STARTING SYMMETRY MODULE        '
  write(*,'(A)') '========================================'
  !===========================================================================
  ! Symmetry initialisation                                           
  ! PGROUP is now read from input_vscf.txt via the PGROUP keyword.
  ! point_group.txt is written here so that the symmetry_module (which
  ! still reads that file internally) gets the user-supplied value.
  ! The projection cutoff (proj_cutoff) is passed to init_symmetry.
  !===========================================================================
  use_symmetry = 0

  block
    logical :: nm_exists

    !--- Write the point group the user specified in the input file ---
    open(12, file='point_group.txt', status='replace')
      write(12,'(A)') trim(point_group_input)
    close(12)

    inquire(file='normal_mode.txt', exist=nm_exists)
    !nm_exists = .true.
    if (nm_exists .and. point_group_input /= 'C1') then
      write(*,'(A)') '========================================'
      write(*,'(A)') ' Found normal_mode.txt'
      write(*,'(A,A)') ' Point group from input: ', trim(point_group_input)
     ! write(*,'(A)') ' Initialising symmetry module...        '
      write(*,'(A)') '========================================'
      write(101,'(A)') '========================================'
      write(101,'(A)') ' Symmetry module initialised            '
      write(101,'(A)') '========================================'
      warning = 0
      open(102, file='debug_sym.txt')
      call init_symmetry(N_modes, mode_irrep_arr, proj_cutoff, warning)         
      close(102)
      use_symmetry = 1
    else
      write(*,'(A)') ' normal_mode.txt not found: symmetry disabled.'
      write(101,'(A)') ' Symmetry disabled (normal_mode.txt missing).'
    end if
  end block

  !===========================================================================
  ! VCI
  !===========================================================================
  ! max_iter_sci and eps_sci are now fully controlled from the input file.  

  if (N_quanta > 0) then

    call mkl_set_num_threads(N_threads)
    call omp_set_num_threads(N_threads)


    write(*,'(A)') '========================================'
    write(*,'(A)') '            FINISHED SYMMETRY           '
    write(*,'(A)') '========================================'

    write(*,'(A)') '========================================'
    write(*,'(A)') '              STARTING VCI              '
    write(*,'(A)') '========================================'

    write(101,*)
    write(101,'(A)') '========================================'
    write(101,*)
    write(101,'(A)') 'VCI performed using the VSCF ground state'
    write(101,*)
    write(101,'(A)') '========================================'
    write(101,'(A)') '            STARTING  VCI               '
    write(101,'(A)') '========================================'
    
    !--- Generate configurations ---
    j = 0
    total_combinations = 0
    total_combinations2 = 0

    !if (max_iter_sci == 0 .and. sci_mode == 'auto') sci_mode = 'auto'

    if (sci_mode .eq. 'auto') then
      if (N_quanta >= N_expansion) N_quanta = N_expansion - 1

      do i = 0, N_quanta
        call count_combinations(N_modes, i, j)
        total_combinations = total_combinations + j
      end do

      if (N_states > total_combinations .or. N_states < 1) &
          N_states = total_combinations

      allocate(combination_vec(total_combinations, N_modes))
      combination_vec = 0
      state_number    = 0
      j               = 0

      do i = 0, N_quanta
        state_number = state_number + j
        call generate_combinations(N_modes, i, combination_vec, &
                                  total_combinations, j, state_number)
      end do

      write(101,*)
      write(101,'(A,I8)') ' States generated: ', total_combinations
      write(*,*)
      write(*,'(A,I8)')   ' States generated: ', total_combinations

    end if

    if (sci_mode .eq. 'list') then
      if (N_quanta >= N_expansion) N_quanta = N_expansion - 1

      do i = 0, N_quanta
        call count_combinations(N_modes, i, j)
        total_combinations2 = total_combinations2 + j
      end do

      if (N_states > total_combinations2 .or. N_states < 1) &
          N_states = total_combinations2

      allocate(combination_vec2(total_combinations2, N_modes))
      combination_vec2 = 0
      state_number    = 0
      j               = 0

      do i = 0, N_quanta
        state_number = state_number + j
        call generate_combinations(N_modes, i, combination_vec2, &
                                  total_combinations2, j, state_number)
      end do

      write(101,*)
      write(101,'(A,I8)') ' States generated: ', total_combinations2
      write(*,*)
      write(*,'(A,I8)')   ' States generated: ', total_combinations2

      open(199, file = 'list_states.txt') 

      read(199, *) max_states
      total_combinations = max_states
      N_states = max_states
      
      allocate(combination_vec(max_states, N_modes))
      write(101,'(A,I8)') ' States read: ', max_states
      write(*,'(A,I8)') ' States read: ', max_states
      
      do i = 1, max_states
        read(199,*) combination_vec(i, 1:N_modes)
        write(*, '(100000I8)') i, combination_vec(i, 1:N_modes)
        write(101, '(100000I8)') i, combination_vec(i, 1:N_modes)
      end do

      close(199)

      write(*,*) 
      write(*,*) '>>> MAX. QUANTA PER STATE: ', maxval(sum(combination_vec(:, 1:N_modes), dim=2))
      write(*,*) 

      if(maxval(sum(combination_vec(:, 1:N_modes), dim=2)) .gt. N_quanta) then

        write(*,*) 'ERROR, min. value for  NQUANT: ', maxval(sum(combination_vec(:, 1:N_modes), dim=2)), ' FOUND: ', N_quanta
        stop

      end if

    end if

    second_dipole_derivatives = second_dipole_derivatives / 2.d0

    !--- Determine calculate_sci from max_iter_sci:                    
    !    If the user set MAXSCI in the input, run Selected CI;
    !    otherwise run Standard VCI (or Symmetry-adapted if applicable).
    !    A value > 0 that was explicitly set triggers SCI.
    !    (proj_cutoff is used inside init_symmetry, not here directly.)
    calculate_sci = 0

    if (max_iter_sci > 0 .and. point_group_name == 'C1') calculate_sci = 1
    if (max_iter_sci == 0 .and. sci_mode == 'list') calculate_sci = 1

    !=========================================================================
    ! Dispatch
    !=========================================================================
    if (use_symmetry == 1 .and. point_group_input /= 'C1') then

      write(*,'(A)')   ' >>> Using SYMMETRY-ADAPTED VCI (VCI-SYM)'
      write(101,'(A)') ' >>> Using SYMMETRY-ADAPTED VCI (VCI-SYM)'

      call vibrational_ci_sym(                                         &
          Potential_3, Potential_4, N_modes, N_expansion,              &
          HO_freq, store_integrals, full_coef, combination_vec,        &
          total_combinations, N_states, N_threads,                     &
          dipole_derivatives, second_dipole_derivatives, N_quanta,     &
          total_3, total_4,                                             &
          Potential_3_vec, Potential_4_vec,                             &
          final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
          final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
          cubic_for_mode, n_cubic_for_mode, n_cubic_max,               &
          quartic_for_mode, n_quartic_for_mode, n_quartic_max,         &
          mode_irrep_arr)

    else if (calculate_sci == 1 .and. point_group_input == 'C1') then

      write(*,'(A)')   ' >>> Using SELECTED VCI (SCI@VSCF)'
      write(*,'(A,I6)')   '     max_iter = ', max_iter_sci             
      write(101,'(A)') ' >>> Using SELECTED VCI (SCI@VSCF)'
      write(101,'(A,I6)') '     max_iter = ', max_iter_sci             

      call selected_vibrational_ci(                                    &
          Potential_3, Potential_4, N_modes, N_expansion,              &
          HO_freq, store_integrals, full_coef, combination_vec,        &
          total_combinations, N_states, N_threads,                     &
          dipole_derivatives, second_dipole_derivatives,               &
          N_quanta, max_iter_sci,                                       &
          total_3, total_4,                                             &
          Potential_3_vec, Potential_4_vec,                             &
          final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
          final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
          cubic_for_mode, n_cubic_for_mode, n_cubic_max,               &
          quartic_for_mode, n_quartic_for_mode, n_quartic_max,         &
          combination_vec2, total_combinations2, sci_mode, quanta_max_reference_sci)

    else

      write(*,'(A)')   ' >>> Using STANDARD VCI'
      write(101,'(A)') ' >>> Using STANDARD VCI'

      call vibrational_ci(                                             &
          Potential_3, Potential_4, N_modes, N_expansion,              &
          HO_freq, store_integrals, full_coef, combination_vec,        &
          total_combinations, N_states, N_threads,                     &
          dipole_derivatives, second_dipole_derivatives, N_quanta,     &
          total_3, total_4,                                             &
          Potential_3_vec, Potential_4_vec,                             &
          final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
          final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
          cubic_for_mode, n_cubic_for_mode, n_cubic_max,               &
          quartic_for_mode, n_quartic_for_mode, n_quartic_max)

    end if

    call cpu_time(end_time)
    
    write(*,'(A)') '========================================'
    write(*,'(A)') '             FINISHED VCI               '
    write(*,'(A)') '========================================'
    
 
    write(*,'(1A,1F11.2,1A)') " Total CPU time:     ", &
                               end_time - start_time, " seconds"

  end if  ! N_quanta > 0

  !===========================================================================
  ! Final timing
  !===========================================================================
  call system_clock(t_end)
  elapsed_time = real(t_end - t_start) / real(count_rate)
  write(*,'(1A,1F11.2,1A)') " Total elapsed time: ", elapsed_time, " seconds"
  write(*,'(A)') '========================================'

  !===========================================================================
  ! Cleanup
  !===========================================================================
  deallocate(Potential_3_vec, Potential_4_vec)
  deallocate(final_index_3, count_index_3, n_unique_3, unique_modes_3, check3)
  deallocate(final_index_4, count_index_4, n_unique_4, unique_modes_4, check4)
  deallocate(cubic_for_mode,   n_cubic_for_mode)
  deallocate(quartic_for_mode, n_quartic_for_mode)
  deallocate(mode_irrep_arr)
  

    if(warning == 1) then
        
        write(*,'(1A)') '<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>'
        write(*,'(1A)')
        write(*,'(1A)') '                WARNING                '
        write(*,'(1A)') '        Characters badly assigned      '
        write(*,'(1A)') 'Consider using different point group/C1'
        write(*,'(1A)')
        write(*,'(1A)') '<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>'
         write(*,'(1A)') ''

        write(101,'(1A)') '<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>'
        write(101,'(1A)')
        write(101,'(1A)') '                WARNING                '
        write(101,'(1A)') '        Characters badly assigned      '
        write(101,'(1A)') 'Consider using different point group/C1'
        write(101,'(1A)')
        write(101,'(1A)') '<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>'
        write(101,'(1A)') ''
       
    end if
    warning = 0


  close(101)
  close(200)
  write(*,'(A)')
  write(*,'(A)') " <:> Normal termination."

end program main_vscf