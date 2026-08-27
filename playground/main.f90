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
  use v_pt2
  implicit none

  !Please report any problem to raphafe96@gmail.com. Collaboration is welcome!

  ! An important note
  !
  ! Everything that is performed, assumed dimensionless coordinates. Why? Because the force constants are given in such a way. If the force constants are not given in such a way
  ! you either convert the constants or change: the integral module needs to take the frequencies: instead of setting the omega=1 parameter- and also when p = 5, the kinect term must be adapated.
  ! To keep the code unifrom when p = 2 I removed the HO_freq from the integral block (omega=1) and moved this contribution in the harmonic terms outside. So check how the harmonic contributions are computed: ...+ modal_int(ii, vm(ii), vn(ii), 2) * HO_freq(ii) * 0.5d0 
  ! I need p=2 for other loops rather than only evaluating harmonic potentials - for examples in cubic and quartic terms with 2 modes equivalent. So this was the best option
  !
  ! This also has implications for the dipoles: when I took the constant that conversts from the atomic units to km/mol, I just used the ratio from the one I obtained and the one from ORCA.
  ! It works perfectly, and no matter the molecule and the mode, there is an unified convertion cosntant.
  !
  ! However, in VCI: (VSCF will be not accounted because it normalizes anayway and only uses the linear dipole):
  ! There is a factor of two in the intensiteis due to this scalling. 
  ! The problem is that when the operator <Psi_HO_0|Qi|Psi_HO_1> is applied, it returns, in dimensionless coordinates: 1/sqrt(2). So you just take this number and multiply by the value of the dipole. Easy.
  ! And for the HO approaximation, this is done! And I went to ORCA, got the value they sinalize as km/mol (it involves other constants to match an absorption) when reescaled the INTENSITY values: 974.88.
  ! Note that I never touched the integral, I just took as it is, and got the scaling. The probkem is that the dipoles are not the transition dipoles, they are just the ones that I got from ORCA. To get the transition values, I would have to dived by sqrt(2)
  ! When the transition dipoles are calculated via CI, they incorporate this constant. Morever, it incorportes terms with p=2. The consequence is that when the dipole is calculated and scaled by 974.88
  ! So when two CI states are used to get the dipole it calculates the integrals, get the values and multiply by the dipoles, either with p=1 or p=2 for quadratic dipoles.
  ! The consequence is that if we take the conversiton for the intensities from ORCA as it is, it will be missing a term of two. Because the intensity would be the square of the dipole and voila a factor of 2 in the intensites.
  ! It still holds true for terms going in the p = 2 kernel, because when this factor was built, only term with p = 1 was present (double harmonic approx.) and it is a mere converting factor. The calculation of anharmonic in the VCI accounts for the terms that are diveded by 2 instead of sqrt2 (p=2)
  ! A way of seeing that is I could get this scalling using the transition dipoles. intead of the pure dipoles.
  ! So the consequence is: for the HO this is the scaling factor. for VCI when the intensity is required as km/mol, it is multiplied by 2. simple. When the intensity is normalized, like the VSCF, this factor is irrelevant.

  !
  !
  !To be aware:
  !
  !I really did not find much about these dipoles from ORCA. I did not transform anything using frequencies (like I assumed they are in diomensionless coordinates already, like the constants) - It seems to work.
  !I tested for ethylene and water and it is pretty much the same value as reported in: J Chem Phys. 2007 May 28;126(20):204101. doi: 10.1063/1.2734970.
  !For ethylene theres one combination band (011000000000) I got 8.9km/mol in VCI and ORCA VPT2 reports 10.33! 


  !TO DO
  ! ------ work on EXCLUD and SA-VCI to work together. RN it is disabled.
  ! when using SA-VCI, the list of normal modes and the external file with the displacement, contain all normal modes (3N including translation and rotation)
  ! the program skip the first 6, and read the rest, but in a full list in order. The list has to adapt
  ! init_symmetry precisa receber a lista de modos excluídos (ou a lista list_new_modes) para filtrar mode_irrep na mesma ordem usada pelo resto do programa:

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
  integer :: use_vci_at_vscf, write_vscf_ref_energy

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
  ! DAV
  !===========================================================================
  integer :: use_davidson_int
  logical :: use_davidson

  !===========================================================================
  ! To exclude
  !===========================================================================
  integer, allocatable :: modes_to_exclude(:), list_new_modes(:), N_modes_new, number_to_exclude
  real*8, allocatable :: Potential_3_to_exclude(:,:,:)
  real*8, allocatable :: Potential_4_to_exclude(:,:,:,:)
  real*8, allocatable :: dipole_derivatives_to_exclude(:,:)
  real*8, allocatable :: second_dipole_derivatives_to_exclude(:,:,:)
  real*8, allocatable :: HO_freq_to_exclude(:)
  logical :: exclude_mode
  real*8 :: energy_excluded

  !
  !
  integer :: run_vpt2


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
  write(*,'(A)') '                ViBra'
  write(*,'(A)') '========================================'
  write(*,'(A)') ' Centro Brasileiro de Pesquisas Fisicas '
  write(*,'(A)') '       CBPF - Rio de Janeiro, Brasil'
  write(*,'(A)') '----------------------------------------'
  write(*,*)

  !===========================================================================
  ! Read input                                                        
  !===========================================================================
 
  input_file = 'input_vscf.txt'
  test = 0 !this will compare the results with implementation from CRYSTAL for H2O found on their webpage (tutorials) https://tutorials.crystalsolutions.eu/tutorial.html?td=anharmonicity&tf=anharm!
  sci_mode = 'auto'
  use_vci_at_vscf = 1
  warning = 0

  call read_intg('RUNSCF', use_vci_at_vscf, 1)
  call read_intg('RUNH2O', test, 0)
  call read_intg('RUNPT2', run_vpt2, 0)
  if (run_vpt2 /= 0) use_vci_at_vscf = 0 !we make the vpt2 theory based on harmonic oscilator...

  N_states = -1 !If less roots is needed, the code must run with davidson.

 ! if (test == 1) use_vci_at_vscf = 1
 


  if (test == 0) then
    call read_inp(input_file, N_modes, N_expansion, constants_file,   &
                  constants_mode, N_quanta, N_states, conv_scf,       &
                  N_threads, point_group_input, proj_cutoff,           &
                  max_iter_sci, sci_mode, quanta_max_reference_sci) 

   use_davidson = .false.
   call read_intg('RUNDAV', use_davidson_int, 0)
   if (use_davidson_int == 1) use_davidson = .true.
   if (use_davidson) N_states = -1        
                                  
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

  call mkl_set_num_threads(INT(N_threads, KIND=4))
  call omp_set_num_threads(INT(N_threads, KIND=4))

  !===========================================================================
  ! Allocate main arrays
  !===========================================================================
  allocate(HO_freq(N_modes))
  allocate(Potential_3(N_modes, N_modes, N_modes))
  allocate(Potential_4(N_modes, N_modes, N_modes, N_modes))
  allocate(dipole_derivatives(N_modes, 3))
  allocate(second_dipole_derivatives(N_modes, N_modes, 3))

  
  Potential_3     = 0.d0
  Potential_4     = 0.d0
  HO_freq         = 0.d0
  

  open(200, file='intensities.txt')
  dipole_derivatives = 0.d0
  second_dipole_derivatives = 0.d0

  !===========================================================================
  ! Read force constants
  !===========================================================================
  if (test == 0) then
    call read_orca(constants_file, 0, HO_freq, Potential_3, Potential_4, &
                   N_modes, dipole_derivatives, second_dipole_derivatives)

    write(*,*) 'DIPOLES (not transition)'
    write(*,*) '------------------------'
    do i = 1, N_modes
      write(*,'(3F12.6)') dipole_derivatives(i, 1:3)
    end do
    write(*,*) '------------------------'

  !!!!!!!!!!!!!!!!!!!
  ! EXCLUDING MODES !
  !!!!!!!!!!!!!!!!!!!

  

  allocate(modes_to_exclude(N_modes))

  exclude_mode = .false.
  modes_to_exclude = -1
  number_to_exclude = 0
  
  call read_exclude(N_modes, HO_freq, exclude_mode, modes_to_exclude, number_to_exclude)
  open(101, file='vscf.out')
  write(101,'(A)') '========================================'
  write(101,'(A)') '                ViBra'
  write(101,'(A)') '========================================'
  write(101,'(A)') ' Centro Brasileiro de Pesquisas Fisicas '
  write(101,'(A)') '       CBPF - Rio de Janeiro, Brasil'
  write(101,'(A)') '----------------------------------------'
  write(101,*)

  if(exclude_mode) then

    write(*,*) '------------------------'
    write(*,*) '    EXCLUDING MODES'
    write(*,*) '------------------------'
    write(101,*) '------------------------'
    write(101,*) '    EXCLUDING MODES'
    write(101,*) '------------------------'

    N_modes_new = N_modes - number_to_exclude
    allocate(list_new_modes(N_modes_new))

    allocate( Potential_3_to_exclude(N_modes_new, N_modes_new, N_modes_new), &
              Potential_4_to_exclude(N_modes_new,N_modes_new,N_modes_new,N_modes_new), &
              dipole_derivatives_to_exclude(N_modes_new, 3), &
              second_dipole_derivatives_to_exclude(N_modes_new,N_modes_new,3), &
              HO_freq_to_exclude(N_modes_new))

    list_new_modes = 0
    energy_excluded = 0.d0
    k = 0
    do i = 1, N_modes
        if (.not. any(modes_to_exclude(1:number_to_exclude) == i)) then
          k = k + 1
          list_new_modes(k) = i
        else
          write(*,'(1A, 1I4, 1A, 1F12.3)') 'Mode excluded (index before exclusion): ', i, '    HO freq. (cm-1): ', HO_freq(i)  
          write(101,'(1A, 1I4, 1A, 1F12.3)') 'Mode excluded (index before exclusion): ', i, '    HO freq. (cm-1): ', HO_freq(i)  
          energy_excluded = energy_excluded + HO_freq(i)
        end if
    end do

    Potential_3_to_exclude = 0.d0
    Potential_4_to_exclude = 0.d0
    dipole_derivatives_to_exclude = 0.d0
    second_dipole_derivatives_to_exclude = 0.d0

    do i = 1, N_modes_new
    dipole_derivatives_to_exclude(i,:) = dipole_derivatives(list_new_modes(i),:)
    HO_freq_to_exclude(i) = HO_freq(list_new_modes(i))
      do j = 1, N_modes_new
      second_dipole_derivatives_to_exclude(i,j,:) = second_dipole_derivatives(list_new_modes(i), list_new_modes(j), :)
        do k = 1, N_modes_new
        Potential_3_to_exclude(i,j,k) = Potential_3(list_new_modes(i),list_new_modes(j),list_new_modes(k))
          do l = 1, N_modes_new
             Potential_4_to_exclude(i,j,k,l) = Potential_4(list_new_modes(i),list_new_modes(j),list_new_modes(k),list_new_modes(l))
          end do
        end do
      end do
    end do

    write(*,*) ' >>> New modes list '
    write(101,*) ' >>> New modes list '

    do i = 1, N_modes_new
      write(*,'(1A, 1I4, 1A, 1F12.3)') 'Mode (index after exclusion): ', i, '    HO freq. (cm-1): ', HO_freq_to_exclude(i)  
      write(101,'(1A, 1I4, 1A, 1F12.3)') 'Mode (index after exclusion): ', i, '    HO freq. (cm-1): ', HO_freq_to_exclude(i)  
    end do


    if (allocated(HO_freq))               deallocate(HO_freq) 
    if (allocated(Potential_3))           deallocate(Potential_3)
    if (allocated(Potential_4))           deallocate(Potential_4)
    if (allocated(dipole_derivatives))    deallocate(dipole_derivatives)
    if (allocated(second_dipole_derivatives)) deallocate(second_dipole_derivatives)

    N_modes = N_modes_new

    allocate(HO_freq(N_modes))
    allocate(Potential_3(N_modes, N_modes, N_modes))
    allocate(Potential_4(N_modes, N_modes, N_modes, N_modes))
    allocate(dipole_derivatives(N_modes, 3))
    allocate(second_dipole_derivatives(N_modes, N_modes, 3))

    
    Potential_3     = 0.d0
    Potential_4     = 0.d0
    HO_freq         = 0.d0
    dipole_derivatives = 0.d0
    second_dipole_derivatives = 0.d0

    Potential_3 = Potential_3_to_exclude
    Potential_4 = Potential_4_to_exclude
    dipole_derivatives = dipole_derivatives_to_exclude
    second_dipole_derivatives = second_dipole_derivatives_to_exclude
    HO_freq = HO_freq_to_exclude

    deallocate( Potential_3_to_exclude, &
              Potential_4_to_exclude, &
              dipole_derivatives_to_exclude, &
              second_dipole_derivatives_to_exclude, &
              HO_freq_to_exclude, &
              modes_to_exclude, &
              list_new_modes)

  write(*,'(1A)') 'WARNING: absolute energies do not include contributions from excluded modes. These must be added manually where required (e.g., for zero-point energy).'
  write(*,'(1A, 1F18.4)') 'Harmonic contribution from excluded modes (1/2 times sum of excluded frequencies, cm-1): ', energy_excluded/2

  write(101,'(1A)') 'WARNING: absolute energies do not include contributions from excluded modes. These must be added manually where required (e.g., for zero-point energy).'
  write(101,'(1A, 1F18.4)') 'Harmonic contribution from excluded modes (1/2 times sum of excluded frequencies, cm-1): ', energy_excluded/2
  end if !this is the exclude mode if

  else !this is the test water if
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


  allocate(full_coef(N_modes, N_expansion, N_expansion))
  allocate(Integral_sobrepos(N_modes, N_expansion))
  allocate(Coeff(N_modes, N_expansion), new_coeff(N_modes, N_expansion))
  allocate(intensities_vscf(2*N_modes + N_modes*(N_modes-1)/2))
  allocate(transition_energy_vscf(2*N_modes + N_modes*(N_modes-1)/2))
  allocate(converged_state(2*N_modes + N_modes*(N_modes-1)/2))
  allocate(mode_irrep_arr(N_modes))

  mode_irrep_arr  = 1     
  converged_state = 0

  call apply_degeneracy_factors(N_modes, Potential_3, Potential_4)

  !===========================================================================
  ! Unit conversion
  !===========================================================================
  cm_to_hartree = 0.0000045563350d0
  HO_freq     = HO_freq     * cm_to_hartree
  Potential_3 = Potential_3 * cm_to_hartree / 6.d0
  Potential_4 = Potential_4 * cm_to_hartree / 24.d0

  !  if (use_vci_at_vscf == 0) then !this builds the perfect vectors from an HO basis.
    full_coef = 0.d0
    Coeff     = 0.d0
    do i = 1, N_modes
    Coeff(i,1) = 1.d0
      do j = 1, N_expansion
        full_coef(i,j,j) = 1.0d0   ! only diagonal elements
      end do
    end do
  !end if

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

if(use_vci_at_vscf == 1) then

  write(*,'(A)') '========================================'
  write(*,'(A)') '             STARTING VSCF              '
  write(*,'(A)') '========================================'

  !===========================================================================
  ! VSCF: ground state
  !===========================================================================
  mode_excite(:) = 0
  write_vscf_ref_energy = 0
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

  do
    energy = new_energy
    write_vscf_ref_energy = write_vscf_ref_energy + 1
    call constant_one_mode(Coeff, Potential_3, Potential_4,            &
        N_modes, N_expansion, HO_freq, new_coeff, store_integrals,     &
        new_energy, 0, mode_excite, full_coef,                         &
        0, dipole_derivatives, trash, N_threads,                       &
        total_3, total_4,                                               &
        Potential_3_vec, Potential_4_vec,                               &
        final_index_3, count_index_3, check3,                          &
        final_index_4, count_index_4, check4, write_vscf_ref_energy)

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
          final_index_4, count_index_4, check4, write_vscf_ref_energy)
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
    write_vscf_ref_energy = 0
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
            final_index_4, count_index_4, check4, write_vscf_ref_energy)

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
              final_index_4, count_index_4, check4, write_vscf_ref_energy)
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
end if !end checking if runs HO or VSCF



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

    call mkl_set_num_threads(INT(N_threads, KIND=4))
    call omp_set_num_threads(INT(N_threads, KIND=4))


    write(*,'(A)') '========================================'
    write(*,'(A)') '            FINISHED SYMMETRY           '
    write(*,'(A)') '========================================'
    if(run_vpt2 == 0) then
      write(*,'(A)') '========================================'
      write(*,'(A)') '              STARTING VCI              '
      write(*,'(A)') '========================================'

      write(101,*)
      write(101,'(A)') '========================================'
      write(101,*)
      if (use_vci_at_vscf == 1) then 
        write(101,'(A)') 'VCI performed using the VSCF ground state'
      else
        write(101,'(A)') 'VCI performed using the HO ground state'
      end if 
      write(101,*)
      write(101,'(A)') '========================================'
      write(101,'(A)') '            STARTING  VPT2              '
      write(101,'(A)') '========================================'
    else
      write(*,'(A)') '========================================'
      write(*,'(A)') '              STARTING VPT2             '
      write(*,'(A)') '========================================'

      write(101,*)
      write(101,'(A)') '========================================'
      write(101,*)
      if (use_vci_at_vscf == 1) then 
        write(101,'(A)') 'VPT2 performed using the VSCF ground state'
      else
        write(101,'(A)') 'VPT2 performed using the HO ground state'
      end if 
      write(101,*)
      write(101,'(A)') '========================================'
      write(101,'(A)') '            STARTING  VPT2              '
      write(101,'(A)') '========================================'
    end if
  !just for debug

if (test == 1) then
  do i = 1, N_modes
    ! Print the matrix for mode i
    write(*,*) 'MODE: ', i
    do j = 1, N_expansion
      write(*,'(*(F9.5))') (full_coef(i,k,j), k=1, N_expansion)
    end do
    write(*,*)
  end do
end if

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

      !if (N_states > total_combinations .or. N_states < 1) &
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

      !if (N_states > total_combinations2 .or. N_states < 1) &
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

    if (exclude_mode == .true. .and. point_group_input /= 'C1') then
      write(*,*)
      write(*,'(A)') 'ERROR: NOT YET IMPLEMENTED TO USE SA-VCI and EXCLUD, change to normal VCI or Selected VCI'
      stop
    end if
    !=========================================================================
    ! Dispatch
    !=========================================================================
    if(run_vpt2 == 1) then
      call vpt2(                                             &
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
      write(*,'(A)') '========================================'
      write(*,'(A)') '             FINISHED VPT2              '
      write(*,'(A)') '========================================'
      call cpu_time(end_time)
    else
    
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
    
  end if
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


  if(exclude_mode) write(*,'(1A)') 'WARNING: absolute energies do not include contributions from excluded modes. These must be added manually where required (e.g., for zero-point energy).'
  if(exclude_mode) write(*,'(1A, 1F18.4)') 'Harmonic contribution from excluded modes (1/2 times sum of excluded frequencies, cm-1): ', energy_excluded/2

  if(exclude_mode) write(101,'(1A)') 'WARNING: absolute energies do not include contributions from excluded modes. These must be added manually where required (e.g., for zero-point energy).'
  if(exclude_mode) write(101,'(1A, 1F18.4)') 'Harmonic contribution from excluded modes (1/2 times sum of excluded frequencies, cm-1): ', energy_excluded/2
  write(*,'(A)')
  write(*,'(A)') " <:> Normal termination."

  
  close(101)
  close(200)

end program main_vscf
