module vib_ci
use combination  ! Provides count_unique_elements subroutine for processing unique mode indices
use jacobi_diagonalization  ! Provides diagonalization routines (dsyevr_A, dsyevd_A)
use omp_lib  ! Provides OpenMP parallelization functions
use symmetry_module
use read_input_file
contains




!TO DO - remove sparse_diff_modes from everything- it is an old trash variable, that may increase memory but does nothing ;D I did it from the fullCI module, since it is the one with the largest hamiltonians.
    !Update: removed for selected vci as well
    !Missing: need to remove from SA VCI

!TO DO - add the restriction to ndiff .lt. 3 to SA-VCI and S-VCI
    !Update: done for S-VCI
    !MISSING: need to do for SA

!TO DO - add absolute dipoles to SA and S-VCI
    !Update: done for S-VCI
    !MISSING: need to do for SA

!TO DO - add davidson to SA and S-VCI. (specially S-VCI)
    !Update: done for S-VCI
    !MISSING: need to do for SA

!TO DO - Work on parallelization of computing H elements (got like 4x speedup when moving from 1 to 8 threads) for the vib_ci routine (full VCI) when calling davidson. Part of this less speedup I believe it is because my pc has other things running as well, so it takes cpu time
!TO DO - There is a critical statement in the sparse pair list, this dumps a lot the perfomance, but in the entire process of build H. the computing elements part is the one that takes longer. This part improved only 2x when 8 threads were used.
!The overall speedup for etileno with (from 1 to) 8 threads was around 3x only. Being the overall - I used 7 quantas to test. (20 and 65 secs - 1 and 8 cores)

!for the davison solver, same etileno, overall speeup was 5x (27s to 135s)


!TO DO - add extra timer counter for specif parts, klike I did for davidson solver... 

!good news :D
!Attention: orca uses a semi quartic force field, meaning all modes can differ up to 3, not 4! it is because it must have at least two indices the same in the quartic terms, which simplifies a LOT!


!=============================================================================!
! Subroutine: compute_H_element                                               !
!                                                                             !
! Computes a single Hamiltonian matrix element H(m,n) between two             !
! configurations identified by their indices into vec_combinations.           !
! This kernel is shared by both vibrational_ci and selected_vibrational_ci    !
! to avoid code duplication.                                                  !
!                                                                             !
! The element is computed as the sum of:                                      !
!   - One-body terms (harmonic + diagonal cubic/quartic corrections)          !
!   - Cubic coupling terms V3(i,j,k)                                          !
!   - Quartic coupling terms V4(i,j,k,l)                                      !
!                                                                             !
! All expensive mu/nu basis function loops are avoided by using the           !
! precomputed modal_int lookup table.                                         !
!=============================================================================!
subroutine compute_H_element(m_cfg, n_cfg, vm, vn, N_modes, max_quanta, &
    modal_int, &
    Potential_3, Potential_4, HO_freq, &
    Potential_3_vec, Potential_4_vec, &
    check3, check4, total_3, total_4, &
    final_index_3, count_index_3, n_unique_3, unique_modes_3, &
    final_index_4, count_index_4, n_unique_4, unique_modes_4, &
    cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    H_val)
implicit none

! Configuration quantum number vectors for the two configurations
! vm(i) = quantum number of mode i in configuration m
! vn(i) = quantum number of mode i in configuration n
integer, intent(in) :: m_cfg, n_cfg  
integer, intent(in) :: N_modes, max_quanta, total_3, total_4
integer, intent(in) :: n_cubic_max, n_quartic_max
integer, intent(in) :: vm(N_modes), vn(N_modes)
integer, intent(in) :: check3(total_3), check4(total_4)
integer, intent(in) :: final_index_3(total_3, 3), count_index_3(total_3, 3)
integer, intent(in) :: n_unique_3(total_3), unique_modes_3(total_3, 3)
integer, intent(in) :: final_index_4(total_4, 4), count_index_4(total_4, 4)
integer, intent(in) :: n_unique_4(total_4), unique_modes_4(total_4, 4)
integer, intent(in) :: cubic_for_mode(n_cubic_max, N_modes)
integer, intent(in) :: n_cubic_for_mode(N_modes)
integer, intent(in) :: quartic_for_mode(n_quartic_max, N_modes)
integer, intent(in) :: n_quartic_for_mode(N_modes)
real*8,  intent(in) :: modal_int(N_modes, 0:max_quanta, 0:max_quanta, 0:5)
real*8,  intent(in) :: Potential_3(N_modes, N_modes, N_modes)
real*8,  intent(in) :: Potential_4(N_modes, N_modes, N_modes, N_modes)
real*8,  intent(in) :: HO_freq(N_modes)
real*8,  intent(in) :: Potential_3_vec(total_3)
real*8,  intent(in) :: Potential_4_vec(total_4)
real*8,  intent(out) :: H_val

! Local variables
integer :: ii, jj, pp, mode_idx, Vc_check, n_diff, cc
integer :: diff_modes(4), to_sparse_cut
real*8  :: ovlp(N_modes), prefix_prod(0:N_modes), suffix_prod(N_modes+1)
real*8  :: prod_except_i(N_modes), prod_all, step_prod
logical :: term_valid
logical :: is_semiquartic_ff

H_val = 0.d0

!call classify_quartic_terms(total_4, final_index_4, n_unique_4, &
!      check4, Potential_4_vec, is_semiquartic_ff)

!to_sparse_cut = merge(3, 4, is_semiquartic_ff)
to_sparse_cut = 4


! --- Compute overlaps and prefix/suffix products ---
! ovlp(ii) = <vm(ii)|vn(ii)> (one-mode overlap)
do ii = 1, N_modes
    ovlp(ii) = modal_int(ii, vm(ii), vn(ii), 0)
end do
ovlp = ovlp + 1.d-16 !If we call the vscf, fine, this will not be 0. However, if this is exactly 0 like for some couplings when calling the VCI@HO, there is a problem when dividing by 0.
prefix_prod(0) = 1.0d0
do ii = 1, N_modes
    prefix_prod(ii) = prefix_prod(ii-1) * ovlp(ii)
end do
prod_all = prefix_prod(N_modes)

suffix_prod(N_modes+1) = 1.0d0
do ii = N_modes, 1, -1
    suffix_prod(ii) = suffix_prod(ii+1) * ovlp(ii)
end do

do ii = 1, N_modes
    prod_except_i(ii) = prefix_prod(ii-1) * suffix_prod(ii+1)
end do

! --- Count differing modes ---
n_diff = 0
diff_modes = 0
do ii = 1, N_modes
    if (vm(ii) /= vn(ii)) then
        n_diff = n_diff + 1
        if (n_diff <= to_sparse_cut) diff_modes(n_diff) = ii
        if (n_diff > to_sparse_cut) return  ! No coupling beyond 4-mode terms
    end if
end do

! --- One-body terms (n_diff <= 1) ---
! Includes: harmonic h_i, diagonal cubic V3(i,i,i), diagonal quartic V4(i,i,i,i)
if (n_diff <= 1) then !these checks were needes before
    do ii = 1, N_modes
       ! if (abs(prod_except_i(ii)) < 1.d-30) cycle
        H_val = H_val + ( &
            modal_int(ii, vm(ii), vn(ii), 5) &
            + modal_int(ii, vm(ii), vn(ii), 4) * Potential_4(ii,ii,ii,ii) &
            + modal_int(ii, vm(ii), vn(ii), 3) * Potential_3(ii,ii,ii) &
            + modal_int(ii, vm(ii), vn(ii), 2) * HO_freq(ii) * 0.5d0 &
            ) * prod_except_i(ii)
    end do
end if

! --- Cubic terms (n_diff <= 3) ---
if (n_diff <= 3) then
    if (n_diff == 0) then
        ! Diagonal: scan all nonzero cubic terms
        do jj = 1, total_3
            if (check3(jj) == 0) cycle
            step_prod = prod_all
            Vc_check = 1
            do pp = 1, n_unique_3(jj)
                mode_idx = unique_modes_3(jj, pp)
                !if (abs(ovlp(mode_idx)) < 1.d-300) then
                !    step_prod = 0.d0; exit
                !end if
                step_prod = step_prod / ovlp(mode_idx)
            end do
          !  if (abs(step_prod) < 1.d-30) cycle
            do pp = 1, 3
                if (final_index_3(jj, pp) > 0 .and. count_index_3(jj, pp) /= 3) then
                    Vc_check = 0
                    mode_idx = final_index_3(jj, pp)
                    step_prod = step_prod * modal_int(mode_idx, vm(mode_idx), vn(mode_idx), count_index_3(jj, pp))
                end if
            end do
            if (Vc_check == 0) H_val = H_val + Potential_3_vec(jj) * step_prod
        end do
    else
        ! Off-diagonal: use inverted index from first differing mode
        do cc = 1, n_cubic_for_mode(diff_modes(1))
            jj = cubic_for_mode(cc, diff_modes(1))
            ! All differing modes must appear in this cubic term
            term_valid = .true.
            do pp = 2, n_diff
                term_valid = .false.
                do ii = 1, n_unique_3(jj)
                    if (unique_modes_3(jj, ii) == diff_modes(pp)) then
                        term_valid = .true.; exit
                    end if
                end do
                if (.not. term_valid) exit
            end do
            if (.not. term_valid) cycle
            step_prod = prod_all
            Vc_check = 1
            do pp = 1, n_unique_3(jj)
                mode_idx = unique_modes_3(jj, pp)
              !  if (abs(ovlp(mode_idx)) < 1.d-300) then
              !      step_prod = 0.d0; exit
              !  end if
                step_prod = step_prod / ovlp(mode_idx)
            end do
          !  if (abs(step_prod) < 1.d-30) cycle
            do pp = 1, 3
                if (final_index_3(jj, pp) > 0 .and. count_index_3(jj, pp) /= 3) then
                    Vc_check = 0
                    mode_idx = final_index_3(jj, pp)
                    step_prod = step_prod * modal_int(mode_idx, vm(mode_idx), vn(mode_idx), count_index_3(jj, pp))
                end if
            end do
            if (Vc_check == 0) H_val = H_val + Potential_3_vec(jj) * step_prod
        end do
    end if
end if

! --- Quartic terms (n_diff <= 4) --- or 3, if ff is semiquartic like orca
if (n_diff <= to_sparse_cut) then
    if (n_diff == 0) then
        ! Diagonal: scan all nonzero quartic terms
        do jj = 1, total_4
            if (check4(jj) == 0) cycle
            step_prod = prod_all
            Vc_check = 1
            do pp = 1, n_unique_4(jj)
                mode_idx = unique_modes_4(jj, pp)
                !if (abs(ovlp(mode_idx)) < 1.d-300) then
                !    step_prod = 0.d0; exit
                !end if
                step_prod = step_prod / ovlp(mode_idx)
            end do
           ! if (abs(step_prod) < 1.d-30) cycle
            do pp = 1, 4
                if (final_index_4(jj, pp) > 0 .and. count_index_4(jj, pp) /= 4) then
                    Vc_check = 0
                    mode_idx = final_index_4(jj, pp)
                    step_prod = step_prod * modal_int(mode_idx, vm(mode_idx), vn(mode_idx), count_index_4(jj, pp))
                end if
            end do
            if (Vc_check == 0) H_val = H_val + Potential_4_vec(jj) * step_prod
        end do
    else
        ! Off-diagonal: use inverted index from first differing mode
        do cc = 1, n_quartic_for_mode(diff_modes(1))
            jj = quartic_for_mode(cc, diff_modes(1))
            term_valid = .true.
            do pp = 2, n_diff
                term_valid = .false.
                do ii = 1, n_unique_4(jj)
                    if (unique_modes_4(jj, ii) == diff_modes(pp)) then
                        term_valid = .true.; exit
                    end if
                end do
                if (.not. term_valid) exit
            end do
            if (.not. term_valid) cycle
            step_prod = prod_all
            Vc_check = 1
            do pp = 1, n_unique_4(jj)
                mode_idx = unique_modes_4(jj, pp)
                !if (abs(ovlp(mode_idx)) < 1.d-300) then
                !    step_prod = 0.d0; exit
                !end if
                step_prod = step_prod / ovlp(mode_idx)
            end do
           ! if (abs(step_prod) < 1.d-30) cycle
            do pp = 1, 4
                if (final_index_4(jj, pp) > 0 .and. count_index_4(jj, pp) /= 4) then
                    Vc_check = 0
                    mode_idx = final_index_4(jj, pp)
                    step_prod = step_prod * modal_int(mode_idx, vm(mode_idx), vn(mode_idx), count_index_4(jj, pp))
                end if
            end do
            if (Vc_check == 0) H_val = H_val + Potential_4_vec(jj) * step_prod
        end do
    end if
end if

end subroutine compute_H_element

!=============================================================================!
! Main VCI routine – This is from my old version,                             !
!I did not properly check every variable after changes... the old ones are    !
!correct, but new ones may be misisng, sorry.                                 !
!I dropped more cooments, so it is fine                                       !
!=============================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N_modes              : number of vibrational normal modes
! N_expansion          : number of modal basis functions per mode  
! N_quanta             : maximum number of quanta (determines size of modal_int table)
! N_states             : number of eigenstates requested from diagonalization
! N_threads            : number of OpenMP threads for parallel execution
! total_combinations   : total number of VCI configurations (size of basis set)
! Potential_3(i,j,k)   : cubic force constants for three-mode coupling
! Potential_4(i,j,k,l) : quartic force constants for four-mode coupling
! HO_freq(i)           : harmonic frequencies for each mode (hartree units)
! store_integrals(mode, mu, nu, power) : precomputed 1D integrals of harmonic oscillator basis
!     power=0: overlap integral <phi_mu|phi_nu>
!     power=1: <phi_mu|q|phi_nu> (position operator)
!     power=2: <phi_mu|q^2|phi_nu> (squared position)
!     power=3: <phi_mu|q^3|phi_nu> (cubic position)
!     power=4: <phi_mu|q^4|phi_nu> (quartic position)
!     power=5: <phi_mu|h_i|phi_nu> (one-mode Hamiltonian)
! full_coef(mode, quanta+1, mu) : VSCF modal coefficients
!     full_coef(i, v+1, mu) = coefficient of basis function mu for quantum number v of mode i
! vec_combinations(config, mode) : quantum number assignments for each configuration
!     vec_combinations(m, i) = quantum number of mode i in configuration m
! dipole_derivatives(mode, xyz) : first derivatives of the dipole moment with respect to mode displacement
! second_dipole_derivatives(i, j, xyz) : second derivatives of the dipole moment
! total_3, total_4             : number of unique cubic/quartic triplets/quartets
! Potential_3_vec, Potential_4_vec : vectorized potentials (precomputed)
! final_index_3/4, count_index_3/4 : unique mode indices and multiplicities (precomputed)
! n_unique_3/4, unique_modes_3/4   : number and list of unique modes per term (precomputed)
! check3/4                          : nonzero flags (precomputed)
! cubic_for_mode, n_cubic_for_mode  : inverted index for cubic terms (precomputed)
! quartic_for_mode, n_quartic_for_mode : inverted index for quartic terms (precomputed)
! n_cubic_max, n_quartic_max        : max entries per mode in inverted index (precomputed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vibrational_ci(Potential_3, Potential_4, N_modes, N_expansion, &
    HO_freq, store_integrals, full_coef, vec_combinations, &
    total_combinations, N_states, N_threads, dipole_derivatives, &
    second_dipole_derivatives, N_quanta, &
    total_3, total_4, &
    Potential_3_vec, Potential_4_vec, &
    final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
    final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
    cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    quartic_for_mode, n_quartic_for_mode, n_quartic_max)
  use omp_lib
  implicit none

 


  logical :: use_davidson 
  integer :: Nfirst_davidson, use_davidson_int, estimate_states
  real*8  :: davidson_tol, max_freq
  integer :: davidson_max_iter 
  real(8) :: H_threshold 
  
  integer :: pos, old_n_sparse
  integer, allocatable :: tmp_sparse_int(:)
  real(8), allocatable :: tmp_sparse_real(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT VARIABLES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, intent(in) :: N_modes, N_expansion, N_quanta, N_states, N_threads
  integer, intent(in) :: total_combinations
  integer, intent(in) :: total_3, total_4, n_cubic_max, n_quartic_max
  integer, intent(in) :: vec_combinations(total_combinations, N_modes)
  integer, intent(in) :: final_index_3(total_3, 3), count_index_3(total_3, 3)
  integer, intent(in) :: n_unique_3(total_3), unique_modes_3(total_3, 3)
  integer, intent(in) :: check3(total_3)
  integer, intent(in) :: final_index_4(total_4, 4), count_index_4(total_4, 4)
  integer, intent(in) :: n_unique_4(total_4), unique_modes_4(total_4, 4)
  integer, intent(in) :: check4(total_4)
  integer, intent(in) :: cubic_for_mode(n_cubic_max, N_modes)
  integer, intent(in) :: n_cubic_for_mode(N_modes)
  integer, intent(in) :: quartic_for_mode(n_quartic_max, N_modes)
  integer, intent(in) :: n_quartic_for_mode(N_modes)
  real*8,  intent(in) :: Potential_3(N_modes, N_modes, N_modes)
  real*8,  intent(in) :: Potential_4(N_modes, N_modes, N_modes, N_modes)
  real*8,  intent(in) :: HO_freq(N_modes)
  real*8,  intent(in) :: full_coef(N_modes, N_expansion, N_expansion)
  real*8,  intent(in) :: store_integrals(N_modes, N_expansion, N_expansion, 0:5)
  real*8,  intent(in) :: dipole_derivatives(N_modes, 3)
  real*8,  intent(in) :: second_dipole_derivatives(N_modes, N_modes, 3)
  real*8,  intent(in) :: Potential_3_vec(total_3)
  real*8,  intent(in) :: Potential_4_vec(total_4)

  ! Local variables
  integer :: l, i, j, k, p, pp, mu, nu, ii, jj, kk, indices(3)
  integer :: N_states_loc, Nfirst_eff, m, n, n_diff
  real*8  :: cm_to_hartree, values(3), sum_energy_HO
  integer :: max_quanta_actual, to_sparse_cut
  logical :: is_semiquartic_ff

  ! Modal integrals
  real*8, allocatable :: modal_int(:,:,:,:)

  ! Hamiltonian (only for LAPACK) and sparse arrays
  real*8, allocatable :: H(:,:), H_sparse(:)
  integer :: n_sparse, idx
  integer, allocatable :: diff_modes(:)

  integer(kind=4), allocatable :: sparse_m(:), sparse_n(:), sparse_ndiff(:)
  integer(kind=4), allocatable :: sparse_m_dav(:), sparse_n_dav(:)
  integer :: n_sparse_dav

  !--------------------------- NEW: Dipole-specific sparse list ---------------------------
  integer :: n_sparse_dip
  integer(kind=4), allocatable :: sparse_m_dip(:), sparse_n_dip(:), sparse_ndiff_dip(:)
  !---------------------------------------------------------------------------------------

  ! Eigenvalues / eigenvectors  (only needed states)
  real*8, allocatable :: eigenvalues(:)
  real*8, allocatable :: eigenvectors(:,:)   ! (total_combinations, N_states_loc)

  ! Dipole intermediates – no full matrix
  real*8, allocatable :: dipole_final(:,:), intensity(:)
  real*8, allocatable :: d_x(:), d_y(:), d_z(:)   ! D*v0 vectors
  integer :: number_to_print_int

  ! Overlap / product arrays for dipole (on‑the‑fly)
  real*8 :: ovlp(N_modes), prefix_prod(0:N_modes), suffix_prod(N_modes+1)
  real*8 :: prod_except_i(N_modes), step_prod, prod_all
  real*8 :: dip_local(3)
  integer :: final_index_dipole((N_modes*(N_modes+1))/2, 2)
  real*8  :: dipole_vec((N_modes*(N_modes+1))/2, 3)
  logical :: term_valid, single_mode_diff
  integer :: diff_count, diff_mode

  ! Work arrays
  real*8 :: H_val
  integer :: vm(N_modes), vn(N_modes)
  real*8, external :: ddot

  ! Timing
  integer(8) :: t_start, t_end, count_rate, count_max
  real(8) :: elapsed_time, start_time, end_time

  integer(8) :: t_start_dav, t_end_dav, count_rate_dav, count_max_dav
  real(8) :: elapsed_time_dav, start_time_dav, end_time_dav

  cm_to_hartree = 0.0000045563350d0
  call system_clock(count_rate=count_rate, count_max=count_max)
  call system_clock(t_start)
  call cpu_time(start_time)

  max_quanta_actual = N_quanta

  call classify_quartic_terms(total_4, final_index_4, n_unique_4, &
      check4, Potential_4_vec, is_semiquartic_ff)
  to_sparse_cut = merge(3, 4, is_semiquartic_ff)
  write(*,'(A,I2)') ' Sparse pair-list cutoff set to n_diff <= ', to_sparse_cut
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use_davidson = .false.
  call read_real('DAVCVG', davidson_tol, 0.0001d0)
  call read_real('DAVCUT', H_threshold, 0.0000005d0)
  call read_real('MAXFRQ', max_freq, 4500.d0)

  call read_intg('DAVMAX', davidson_max_iter, 25)
  call read_intg('DAVSTA', Nfirst_davidson, 10)
  call read_intg('RUNDAV', use_davidson_int, 0)
  call read_intg('RUNENR', estimate_states, 0)
  

  if (use_davidson_int == 1) use_davidson = .true.
  if (Nfirst_davidson .gt. total_combinations) Nfirst_davidson = total_combinations

  
  i = 0
  if(estimate_states .gt. 0) then
  write(*,'(A,F12.2)') ' >>> Estimating (based on HO) the number of states with energy bellow the cutoff: ',  max_freq
    do j = 1, total_combinations
        sum_energy_HO = 0
        do k = 1, N_modes
            sum_energy_HO = sum_energy_HO + HO_freq(k)*real(vec_combinations(j, k)) 
        end do
        if(sum_energy_HO / cm_to_hartree .lt. max_freq) i = i + 1
        if (estimate_states == 2 .and. sum_energy_HO / cm_to_hartree .lt. max_freq) then
            do l = 1, N_modes
                write(*, '(I4)', advance='no') vec_combinations(j, l)
            end do
            write(*, '(F18.2)') sum_energy_HO / cm_to_hartree
        end if
    end do
  write(*,'(A,I8)') ' States found: ', i
  end if


  allocate(diff_modes(to_sparse_cut))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 1: Precompute modal integrals
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(modal_int(N_modes, 0:max_quanta_actual, 0:max_quanta_actual, 0:5))
  modal_int = 0.d0
  write(*,'(A,I4,A,I4)') ' Precomputing modal integrals: ', N_modes, &
       ' modes, max_quanta=', max_quanta_actual
  call omp_set_num_threads(INT(N_threads, KIND=4))

!$OMP PARALLEL DO PRIVATE(ii, i, j, p, mu, nu) COLLAPSE(2) SCHEDULE(static) !This parallel section is not very efficient... loops are not so large.
  do ii = 1, N_modes
     do i = 0, max_quanta_actual
        do j = i, max_quanta_actual
           do p = 0, 5
            !  modal_int(ii, i, j, p) = store_integrals(ii, i+1, j+1, p)

              do mu = 1, N_expansion
                 do nu = 1, N_expansion
                    modal_int(ii, i, j, p) = modal_int(ii, i, j, p) &
                         + full_coef(ii, i+1, mu) * full_coef(ii, j+1, nu) &
                         * store_integrals(ii, mu, nu, p)
                 end do
              end do
           end do
           do p = 0, 5
              modal_int(ii, j, i, p) = modal_int(ii, i, j, p)
           end do
        end do
     end do
  end do
 !$OMP END PARALLEL DO
  write(*,'(A)') ' Modal integrals precomputed.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 2: Build dipole pair vector (for on‑the‑fly use)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dipole_vec = 0.d0
  kk = 0
  do i = 1, N_modes
     do j = i, N_modes
        kk = kk + 1
        dipole_vec(kk, 1:3) = second_dipole_derivatives(i, j, 1:3)
        final_index_dipole(kk, 1) = i
        final_index_dipole(kk, 2) = j
     end do
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 3: Build sparse pair list (Hamiltonian)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,'(A)') ' Building sparse pair list ...'
  n_sparse = 0
  !$omp parallel do default(none) &
  !$omp shared(total_combinations, N_modes, vec_combinations,to_sparse_cut) &
  !$omp private(n, m, n_diff, ii) &
  !$omp reduction(+:n_sparse)
  do n = 1, total_combinations
     do m = n, total_combinations
        n_diff = 0
        do ii = 1, N_modes
           if (vec_combinations(m, ii) /= vec_combinations(n, ii)) then
              n_diff = n_diff + 1
              if (n_diff > to_sparse_cut) exit
           end if
        end do
        if (n_diff <= to_sparse_cut) n_sparse = n_sparse + 1
     end do
  end do
  !$omp end parallel do

  write(*,'(A,I12,A,I12,A,F6.2,A)') ' Non zero pairs: ', n_sparse, ' out of ', &
       total_combinations*(total_combinations+1)/2, ' (', &
       100.d0*dble(n_sparse)/dble(total_combinations*(total_combinations+1)/2), '%)'

  allocate(sparse_m(n_sparse), sparse_n(n_sparse))
  allocate(sparse_ndiff(n_sparse))
  idx = 0
   !$omp parallel do private(m, n_diff, ii, diff_modes) &
   !$omp shared(idx, sparse_m, sparse_n, sparse_ndiff, vec_combinations, N_modes, total_combinations,to_sparse_cut) &
   !$omp default(none)
  do n = 1, total_combinations
     do m = n, total_combinations
        n_diff = 0
        diff_modes = 0
        do ii = 1, N_modes
           if (vec_combinations(m, ii) /= vec_combinations(n, ii)) then
              n_diff = n_diff + 1
              if (n_diff > to_sparse_cut) exit
              diff_modes(n_diff) = ii
           end if
        end do
        if (n_diff <= to_sparse_cut) then
           !$omp critical
           idx = idx + 1
           sparse_m(idx) = m
           sparse_n(idx) = n
           sparse_ndiff(idx) = n_diff
           !$omp end critical
        end if
     end do
  end do
   !$omp end parallel do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Build dipole-specific pair list (n_diff <= 2 only)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n_sparse_dip = 0

  do idx = 1, n_sparse
     if (sparse_ndiff(idx) <= 2) n_sparse_dip = n_sparse_dip + 1
  end do


  allocate(sparse_m_dip(n_sparse_dip), sparse_n_dip(n_sparse_dip), sparse_ndiff_dip(n_sparse_dip))
  n_sparse_dip = 0
  do idx = 1, n_sparse
     if (sparse_ndiff(idx) <= 2) then
        n_sparse_dip = n_sparse_dip + 1
        sparse_m_dip(n_sparse_dip) = sparse_m(idx)
        sparse_n_dip(n_sparse_dip) = sparse_n(idx)
        sparse_ndiff_dip(n_sparse_dip) = sparse_ndiff(idx)
     end if
  end do
  write(*,'(A,I0,A)') ' Dipole-specific pair list contains ', n_sparse_dip, ' elements.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 4: Build Hamiltonian / Precompute H_sparse (Davidson)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (.not. use_davidson) then
     allocate(H(total_combinations, total_combinations))
     H = 0.d0
     !$OMP PARALLEL DO PRIVATE(idx, m, n, H_val, vm, vn) SCHEDULE(dynamic,256) DEFAULT(NONE) &
     !$OMP& SHARED(n_sparse, sparse_m, sparse_n, N_modes, max_quanta_actual, vec_combinations, &
     !$OMP&        modal_int, HO_freq, Potential_3, Potential_4, Potential_3_vec, Potential_4_vec, &
     !$OMP&        check3, check4, total_3, total_4, final_index_3, count_index_3, n_unique_3, &
     !$OMP&        unique_modes_3, final_index_4, count_index_4, n_unique_4, unique_modes_4, &
     !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
     !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, H)
     do idx = 1, n_sparse
        m = sparse_m(idx)
        n = sparse_n(idx)
        vm(:) = vec_combinations(m, 1:N_modes)
        vn(:) = vec_combinations(n, 1:N_modes)
        call compute_H_element(m, n, vm, vn, N_modes, max_quanta_actual, &
             modal_int, Potential_3, Potential_4, HO_freq, &
             Potential_3_vec, Potential_4_vec, &
             check3, check4, total_3, total_4, &
             final_index_3, count_index_3, n_unique_3, unique_modes_3, &
             final_index_4, count_index_4, n_unique_4, unique_modes_4, &
             cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
             quartic_for_mode, n_quartic_for_mode, n_quartic_max, H_val)
        H(m,n) = H_val
     end do
     !$OMP END PARALLEL DO
     do m = 1, total_combinations
        do n = m+1, total_combinations
           H(m,n) = H(n,m)
        end do
     end do
  else
        allocate(H_sparse(n_sparse))
        H_sparse = 0.d0
        write(*,'(A)') ' Precomputing sparse Hamiltonian matrix (Davidson)...'
        !$OMP PARALLEL DO PRIVATE(idx, m, n, H_val, vm, vn) SCHEDULE(dynamic,256) DEFAULT(NONE) &
        !$OMP& SHARED(n_sparse, sparse_m, sparse_n, N_modes, max_quanta_actual, vec_combinations, &
        !$OMP&        modal_int, HO_freq, Potential_3, Potential_4, Potential_3_vec, Potential_4_vec, &
        !$OMP&        check3, check4, total_3, total_4, final_index_3, count_index_3, n_unique_3, &
        !$OMP&        unique_modes_3, final_index_4, count_index_4, n_unique_4, unique_modes_4, &
        !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
        !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, H_sparse)
        do idx = 1, n_sparse
           m = sparse_m(idx)
           n = sparse_n(idx)
           vm(:) = vec_combinations(m, 1:N_modes)
           vn(:) = vec_combinations(n, 1:N_modes)
           call compute_H_element(m, n, vm, vn, N_modes, max_quanta_actual, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, H_val)
           H_sparse(idx) = H_val
        end do
        !$OMP END PARALLEL DO
        write(*,'(A)') ' Sparse Hamiltonian precomputed.'

        ! Build a SEPARATE, compacted copy for the Davidson matvec only
        write(*,'(A,ES10.3)') ' Filtering elements below threshold (Davidson H only): ', H_threshold
        allocate(sparse_m_dav(n_sparse), sparse_n_dav(n_sparse))
        n_sparse_dav = 0
        do idx = 1, n_sparse
           if (abs(H_sparse(idx)) >= H_threshold) then
              n_sparse_dav = n_sparse_dav + 1
              sparse_m_dav(n_sparse_dav) = sparse_m(idx)
              sparse_n_dav(n_sparse_dav) = sparse_n(idx)
              H_sparse(n_sparse_dav)     = H_sparse(idx)
           end if
        end do
        write(*,'(A,I0,A,I0,A,F6.2,A)') &
             ' Kept ', n_sparse_dav, ' / ', n_sparse, &
             ' elements (', 100.0d0 * real(n_sparse_dav,8) / real(n_sparse,8), '% retained)'

        ! Reallocate to exact size
        allocate(tmp_sparse_int(n_sparse_dav))
        tmp_sparse_int(:) = sparse_m_dav(1:n_sparse_dav)
        deallocate(sparse_m_dav)
        allocate(sparse_m_dav(n_sparse_dav))
        sparse_m_dav(:) = tmp_sparse_int(:)
        deallocate(tmp_sparse_int)

        allocate(tmp_sparse_int(n_sparse_dav))
        tmp_sparse_int(:) = sparse_n_dav(1:n_sparse_dav)
        deallocate(sparse_n_dav)
        allocate(sparse_n_dav(n_sparse_dav))
        sparse_n_dav(:) = tmp_sparse_int(:)
        deallocate(tmp_sparse_int)

        allocate(tmp_sparse_real(n_sparse_dav))
        tmp_sparse_real(:) = H_sparse(1:n_sparse_dav)
        deallocate(H_sparse)
        allocate(H_sparse(n_sparse_dav))
        H_sparse(:) = tmp_sparse_real(:)
        deallocate(tmp_sparse_real)

        write(*,'(A,I0)') ' Final number of Davidson sparse elements: ', n_sparse_dav
  end if

  call cpu_time(end_time)
  call system_clock(t_end)
  elapsed_time = real(t_end - t_start) / real(count_rate)
  write(*,'(1A, 1F12.2, 1A)') " H build CPU time:     ", end_time - start_time, " seconds"
  write(*,'(1A, 1F12.2, 1A)') " H build elapsed time: ", elapsed_time, " seconds"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 5: Diagonalisation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (use_davidson) then
     N_states_loc = Nfirst_davidson
  else
     N_states_loc = N_states
  end if

  allocate(eigenvalues(N_states_loc))
  allocate(eigenvectors(total_combinations, N_states_loc))
  eigenvalues = 0.d0
  eigenvectors = 0.d0
  
  ! Print eigenvalues and CI decomposition
  if (use_davidson .and. N_threads == 1) write(*,'(1A,1F18.4)') ' H(1,1) = ', H_sparse(1)/cm_to_hartree !this necessarely not produce the H(1,1) at H_sparse(1) due threading. Only in single thread.
  
!  if (use_davidson .and. N_threads == 1) then
!    do i = 1, n_sparse_dav
!        if(sparse_m_dav(i) >= sparse_n_dav(i) .and. sparse_m_dav(i) .lt. 20) write(*,'(2I4, 1F18.2)') sparse_n_dav(i), sparse_m_dav(i), H_sparse(i)/cm_to_hartree
!    end do
!  end if

  if (use_davidson) then
     call system_clock(count_rate=count_rate_dav, count_max=count_max_dav)
     call system_clock(t_start_dav)
     call cpu_time(start_time_dav)

     write(*,'(1A,I4)') ' >>>>> Using Davidson eigensolver. Roots requested: ', Nfirst_davidson
     call jacobi_davidson_eigensolver(total_combinations, Nfirst_davidson, &
          n_sparse_dav, sparse_m_dav, sparse_n_dav, H_sparse, &
          davidson_tol, davidson_max_iter, eigenvalues, eigenvectors)
     write(*,'(1A)') ' Davidson eigensolver finished.'
     deallocate(H_sparse, sparse_m_dav, sparse_n_dav)

    call cpu_time(end_time_dav)
    call system_clock(t_end_dav)
    elapsed_time_dav = real(t_end_dav - t_start_dav) / real(count_rate_dav)
    write(*,'(1A, 1F12.2, 1A)') " DAV CPU time:     ", end_time_dav - start_time_dav, " seconds"
    write(*,'(1A, 1F12.2, 1A)') " DAV elapsed time: ", elapsed_time_dav, " seconds"

  else
     write(*,'(1A)') ' >>>>> H matrix assembled. Entering diagonalisation...'
     call system_clock(t_start)
     call cpu_time(start_time)
     if (N_states_loc /= total_combinations) then
        call dsyevr_A(H, N_states_loc, eigenvalues, eigenvectors, 'N')
     else
        call dsyevd_A(H, eigenvalues, eigenvectors)
     end if
     call cpu_time(end_time)
     call system_clock(t_end)
     elapsed_time = real(t_end - t_start) / real(count_rate)
     write(*,'(1A, 1F12.2, 1A)') " H diago CPU time:     ", end_time - start_time, " seconds"
     write(*,'(1A, 1F12.2, 1A)') " H diago elapsed time: ", elapsed_time, " seconds"
     deallocate(H)
  end if
  
  write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
  write(*,'(1A19, 1A22, 1A14, 6A12)') 'E (cm-1)', 'E - ZPE', &
       '              ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C'
  write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
  write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
  write(101,'(1A19, 1A22, 1A14, 6A12)') 'E (cm-1)', 'E - ZPE', &
       '              ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C'
  write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'

  do i = 1, N_states_loc
     values = 0.
     indices = 0
     call top_three_components(eigenvectors(:,i), total_combinations, values, indices)
     if (i <= 10) then
        write(*,'(1F19.4, 1F22.4, 1A14, 3F12.6, 3I12)') &
             eigenvalues(i)/cm_to_hartree, &
             eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree, &
             '      ||      ', values, indices
     end if
     write(101,'(1F19.4, 1F22.4, 1A14, 3F12.6, 3I12)') &
          eigenvalues(i)/cm_to_hartree, &
          eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree, &
          '      ||      ', values, indices
  end do
  write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
  write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 6: Transition dipoles – memory‑lean + energy filter
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*)
  write(*,*) '-Calculating transition dipoles (onthefly, no full D matrix)'

  ! Determine how many transitions lie below 4500 cm⁻¹
  number_to_print_int = 0
  do i = 1, N_states_loc - 1
     if ( (eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree < max_freq ) then
        number_to_print_int = number_to_print_int + 1
     end if
  end do
  if (number_to_print_int == 0) number_to_print_int = N_states_loc - 1   ! fallback

  if (number_to_print_int < 1) then
     deallocate(modal_int, sparse_m, sparse_n, sparse_ndiff)
     if (allocated(sparse_m_dip)) deallocate(sparse_m_dip, sparse_n_dip, sparse_ndiff_dip)
     return
  end if

  allocate(dipole_final(number_to_print_int, 3))
  allocate(intensity(number_to_print_int))
  allocate(d_x(total_combinations), d_y(total_combinations), d_z(total_combinations))
  d_x = 0.d0; d_y = 0.d0; d_z = 0.d0

  do idx = 1, n_sparse_dip
     n_diff = sparse_ndiff_dip(idx)          ! always ≤ 2
     m = sparse_m_dip(idx)
     n = sparse_n_dip(idx)

     ! overlaps
     do ii = 1, N_modes
        ovlp(ii) = modal_int(ii, vec_combinations(m,ii), vec_combinations(n,ii), 0)
     end do
     prefix_prod(0) = 1.0d0
     do ii = 1, N_modes
        prefix_prod(ii) = prefix_prod(ii-1) * ovlp(ii)
     end do
     prod_all = prefix_prod(N_modes)
     suffix_prod(N_modes+1) = 1.0d0
     do ii = N_modes, 1, -1
        suffix_prod(ii) = suffix_prod(ii+1) * ovlp(ii)
     end do
     do ii = 1, N_modes
        prod_except_i(ii) = prefix_prod(ii-1) * suffix_prod(ii+1)
     end do

     dip_local = 0.d0

     ! Determine if the change is a single-mode Δv = 1 or Δv = 2 (for diagonal second derivative)
     diff_count = 0
     diff_mode = 0
     do ii = 1, N_modes
        if (vec_combinations(m,ii) /= vec_combinations(n,ii)) then
           diff_count = diff_count + 1
           diff_mode = ii
        end if
     end do
     single_mode_diff = (diff_count == 1)

     ! one-body dipole (including Δv=2 diagonal)
     if (n_diff <= 1 .or. (n_diff == 2 .and. single_mode_diff)) then
        do ii = 1, N_modes
           if (abs(prod_except_i(ii)) < 1.d-30) cycle
           dip_local(1:3) = dip_local(1:3) + ( &
                modal_int(ii, vec_combinations(m,ii), vec_combinations(n,ii), 1)*dipole_derivatives(ii, 1:3)&
              + modal_int(ii, vec_combinations(m,ii), vec_combinations(n,ii), 2)*second_dipole_derivatives(ii, ii, 1:3) &
                ) * prod_except_i(ii)
        end do
     end if

     ! off-diagonal dipole (i ≠ j)
     do jj = 1, size(dipole_vec,1)
        if (final_index_dipole(jj,1) == final_index_dipole(jj,2)) cycle
        if (all(abs(dipole_vec(jj,1:3)) < 1.d-30)) cycle
        ii = final_index_dipole(jj, 1)
        kk = final_index_dipole(jj, 2)
        if (abs(ovlp(ii)) < 1.d-300 .or. abs(ovlp(kk)) < 1.d-300) then
           step_prod = 1.0d0
           do pp = 1, N_modes
              if (pp /= ii .and. pp /= kk) step_prod = step_prod * ovlp(pp)
           end do
        else
           step_prod = prod_all / (ovlp(ii) * ovlp(kk))
        end if
        if (abs(step_prod) < 1.d-30) cycle
        step_prod = step_prod &
             * modal_int(ii, vec_combinations(m,ii), vec_combinations(n,ii), 1)*modal_int(kk, vec_combinations(m,kk), vec_combinations(n,kk), 1)
        dip_local(1:3) = dip_local(1:3) + dipole_vec(jj, 1:3) * step_prod * 2.0d0
     end do

     ! D · v0 accumulation
     d_x(m) = d_x(m) + dip_local(1) * eigenvectors(n,1)
     d_y(m) = d_y(m) + dip_local(2) * eigenvectors(n,1)
     d_z(m) = d_z(m) + dip_local(3) * eigenvectors(n,1)
     if (m /= n) then
        d_x(n) = d_x(n) + dip_local(1) * eigenvectors(m,1)
        d_y(n) = d_y(n) + dip_local(2) * eigenvectors(m,1)
        d_z(n) = d_z(n) + dip_local(3) * eigenvectors(m,1)
     end if
  end do

  write(*,*) 'Dipole matrix multiplication finished.'

  ! Now build the transition dipoles for the filtered states
  kk = 0
  do i = 1, N_states_loc - 1
     if ( (eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree < max_freq ) then
        kk = kk + 1
        dipole_final(kk,1) = ddot(total_combinations, eigenvectors(1,i+1), 1, d_x, 1)
        dipole_final(kk,2) = ddot(total_combinations, eigenvectors(1,i+1), 1, d_y, 1)
        dipole_final(kk,3) = ddot(total_combinations, eigenvectors(1,i+1), 1, d_z, 1)
        intensity(kk) = sum(dipole_final(kk,:)**2)
     end if
  end do

    
    !---------------------------------------------------------------------
  ! OLD output: unit 200 (original code – kept unchanged)
  !---------------------------------------------------------------------
  write(200,*) 'VCI'
  write(200,*) number_to_print_int
  do i = 1, number_to_print_int
     write(200,'(2F18.6)') (eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree, &
          100.d0 * intensity(i) / sqrt(sum(intensity(1:number_to_print_int)**2))
  end do



    open(unit=300, file='dipoles_intensity_vci.txt', status='replace')
    write(300,'(A)') 'Transition dipoles and intensities (km/mol)'
    write(300,'(A)') 'State  Energy (cm-1)   Dip_x (a.u.)   Dip_y (a.u.)   Dip_z (a.u.)   Intensity (km/mol)'
    do i = 1, number_to_print_int
    write(300,'(I6,F18.6,3F14.8,F14.4)') i+1,                               &
            (eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree,                                                        &
            dipole_final(i,1), dipole_final(i,2), dipole_final(i,3), intensity(i)*974.88d0*2.d0
    end do

  close(300)
  !---------------------------------------------------------------------

  ! Cleanup
  deallocate(d_x, d_y, d_z, dipole_final, intensity)
  deallocate(modal_int, sparse_m, sparse_n, sparse_ndiff)
  if (allocated(sparse_m_dip)) deallocate(sparse_m_dip, sparse_n_dip, sparse_ndiff_dip)
  deallocate(eigenvalues, eigenvectors)

  write(*,*) '-DONE'
  write(*,*)
end subroutine vibrational_ci



!=============================================================================!
! Subroutine: selected_vibrational_ci                                         !
!                                                                             !
! Performs Selected Vibrational CI (SCI) starting from a CISD reference space.!
!                                                                             !
! Algorithm:                                                                  !
!   1. Identify the CISD reference: all configurations in vec_combinations    !
!      that differ from the ground state (configuration 1, all zeros) by at   !
!      most 2 modes (singles and doubles).                                    !
!   2. Build and diagonalize H in the reference space.                        !
!   3. For every configuration outside the reference (the "external" space),  !
!     compute the Epstein-Nesbet PT2 correction for each reference eigenstate:!
!         e_alpha = |<alpha|H|Psi_I>|^2 / (E_I - H_aa)                        !
!      where H_aa = <alpha|H|alpha> is the diagonal element of the external   !
!      configuration.                                                         !
!   4. For each CISD eigenstate I, the N_sel_per_state configurations with    !
!      the largest |e_alpha| are selected.                                    !
!   5. The union of all selected configurations across all states is added    !
!      to the CISD reference, the full H is built and diagonalized once.      !
!   6. Final energies and eigenvectors are printed.                           !
!                                                                             !
! Inputs mirror those of vibrational_ci plus N_sel_per_state (number of       !
! external configurations to select per CISD state).                          !
!=============================================================================!
subroutine selected_vibrational_ci(Potential_3, Potential_4, N_modes, N_expansion, &
    HO_freq, store_integrals, full_coef, vec_combinations, &
    total_combinations, N_states, N_threads, dipole_derivatives, &
    second_dipole_derivatives, N_quanta, N_sel_per_state, &
    total_3, total_4, &
    Potential_3_vec, Potential_4_vec, &
    final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
    final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
    cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    quartic_for_mode, n_quartic_for_mode, n_quartic_max, vec_combinations2, &
    total_combinations2, mode_sci, reference_max)
implicit none

!---------------------------------------------------------------------------
! INPUT VARIABLES
!---------------------------------------------------------------------------
integer,  intent(in) :: N_modes, N_expansion, N_quanta, N_threads, reference_max
integer,  intent(in) :: total_combinations, total_combinations2
integer,  intent(in) :: vec_combinations(total_combinations, N_modes), vec_combinations2(total_combinations2, N_modes)
integer,  intent(in) :: total_3, total_4, n_cubic_max, n_quartic_max
integer,  intent(in) :: final_index_3(total_3, 3), count_index_3(total_3, 3)
integer,  intent(in) :: n_unique_3(total_3), unique_modes_3(total_3, 3)
integer,  intent(in) :: check3(total_3)
integer,  intent(in) :: final_index_4(total_4, 4), count_index_4(total_4, 4)
integer,  intent(in) :: n_unique_4(total_4), unique_modes_4(total_4, 4)
integer,  intent(in) :: check4(total_4)
integer,  intent(in) :: cubic_for_mode(n_cubic_max, N_modes)
integer,  intent(in) :: n_cubic_for_mode(N_modes)
integer,  intent(in) :: quartic_for_mode(n_quartic_max, N_modes)
integer,  intent(in) :: n_quartic_for_mode(N_modes)
!integer,  intent(in) :: N_sel_per_state
real*8,   intent(in) :: Potential_3(N_modes, N_modes, N_modes)
real*8,   intent(in) :: Potential_4(N_modes, N_modes, N_modes, N_modes)
real*8,   intent(in) :: HO_freq(N_modes)
real*8,   intent(in) :: full_coef(N_modes, N_expansion, N_expansion)
real*8,   intent(in) :: store_integrals(N_modes, N_expansion, N_expansion, 0:5)
real*8,   intent(in) :: dipole_derivatives(N_modes, 3)
real*8,   intent(in) :: second_dipole_derivatives(N_modes, N_modes, 3)
real*8,   intent(in) :: Potential_3_vec(total_3)
real*8,   intent(in) :: Potential_4_vec(total_4)
character (len=4), intent(in) :: mode_sci 

!---------------------------------------------------------------------------
! SCALAR LOCALS
!---------------------------------------------------------------------------
integer :: i, j, k, ii, jj, kk, pp, m, n, n_states_list, N_sel_per_state, check_list, N_states
integer :: max_quanta
integer :: mu, nu, p  
integer :: mode_idx, Vc_check, n_diff
integer :: n_ref, n_sel, n_ext, n_new
integer :: n_sparse, idx, cc
integer :: n_cisd_states, i_state
integer :: number_to_print_int
integer, allocatable :: diff_modes(:)
real*8  :: cm_to_hartree
real*8  :: H_val, H_aa, e_alpha, denom
real*8  :: step_prod, prod_all
real*8  :: values(3)
integer :: indices(3)

!---------------------------------------------------------------------------
! MODAL INTEGRAL TABLE
!---------------------------------------------------------------------------
real*8, allocatable :: modal_int(:,:,:,:)

!---------------------------------------------------------------------------
! CONFIGURATION SELECTION ARRAYS
!---------------------------------------------------------------------------
logical, allocatable :: is_selected(:), is_ref(:)
integer, allocatable :: ref_list(:), sel_list(:), ext_list(:)

!---------------------------------------------------------------------------
! HAMILTONIAN AND EIGENSYSTEM FOR ACTIVE SPACE
!---------------------------------------------------------------------------
real*8, allocatable :: H_sel(:,:)
real*8, allocatable :: eigenvalues(:), eigenvectors(:,:)

!---------------------------------------------------------------------------
! SPARSE PAIR LIST FOR ACTIVE SPACE
!---------------------------------------------------------------------------
integer(kind=4), allocatable :: sparse_m(:), sparse_n(:)
integer, allocatable :: sparse_ndiff(:)

!---------------------------------------------------------------------------
! DIPOLE ARRAYS
!---------------------------------------------------------------------------
real*8, allocatable :: dipole_sel(:,:,:)
real*8, allocatable :: dipole_final(:,:)
real*8, allocatable :: dipole_vec(:,:)
integer, allocatable :: final_index_dipole(:,:)
real*8, allocatable :: intensity(:)
real*8, allocatable :: tmp_vec(:)
real*8, external    :: ddot

!---------------------------------------------------------------------------
! EN-PT2 CONTRIBUTION ARRAYS FOR TOP-N SELECTION
!---------------------------------------------------------------------------
real*8, allocatable :: pt2_contrib(:,:)
integer, allocatable :: top_n_indices(:,:)
real*8, allocatable :: top_n_values(:,:)

!---------------------------------------------------------------------------
! COUPLING VECTOR FOR PT2 SCREENING
!---------------------------------------------------------------------------
real*8, allocatable :: coupling_vec(:)
real*8, allocatable :: hrow(:)

!---------------------------------------------------------------------------
! SCRATCH LOCALS
!---------------------------------------------------------------------------
integer :: vm(N_modes), vn(N_modes)
real*8  :: dip_local(3), ovlp_loc(N_modes)
real*8  :: prefix_prod(0:N_modes), suffix_prod_loc(N_modes+1)
real*8  :: prod_except_i(N_modes), step_prod_loc

!---------------------------------------------------------------------------
! TIMING
!---------------------------------------------------------------------------
integer(8) :: t_start, t_end, count_rate_val, count_max_val
real(8)    :: elapsed_time, start_time, end_time

!---------------------------------------------------------------------------
! FOR DAV
!---------------------------------------------------------------------------
real*8 :: davidson_tol, H_threshold, max_freq
integer :: davidson_max_iter, Nfirst_davidson, use_davidson_int, estimate_states
logical :: use_davidson, is_semiquartic_ff
real*8, allocatable :: H_sparse(:)
integer(kind=4), allocatable :: sparse_m_dav(:), sparse_n_dav(:)
integer :: n_sparse_dav
integer :: n_sparse_dip
integer :: to_sparse_cut
integer(kind=4), allocatable :: sparse_m_dip(:), sparse_n_dip(:), sparse_ndiff_dip(:)
integer(8) :: t_start_dav, t_end_dav, count_rate_dav, count_max_dav
real(8) :: elapsed_time_dav, start_time_dav, end_time_dav
integer, allocatable :: tmp_sparse_int(:)
real(8), allocatable :: tmp_sparse_real(:)
integer :: max_quanta_actual
integer :: to_write_output
real*8, allocatable :: d_x(:), d_y(:), d_z(:)
integer :: n_loc

!==========================================================================
! INITIALIZATION
!==========================================================================
cm_to_hartree = 0.0000045563350d0
max_quanta    = N_quanta
max_quanta_actual = N_quanta

check_list = 0
if (mode_sci == 'list') then
    check_list = 1
    N_states = total_combinations2
    !N_sel_per_state = 8
end if

  !NEW INPUT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use_davidson = .false.
  call read_real('DAVCVG', davidson_tol, 0.0001d0)
  call read_real('DAVCUT', H_threshold, 0.0000005d0)
  call read_real('MAXFRQ', max_freq, 4500.d0)

  call read_intg('DAVMAX', davidson_max_iter, 25)
  call read_intg('DAVSTA', Nfirst_davidson, 10)
  call read_intg('RUNDAV', use_davidson_int, 0)
  call read_intg('RUNENR', estimate_states, 0)
  

  if (use_davidson_int == 1) use_davidson = .true.


call system_clock(count_rate=count_rate_val, count_max=count_max_val)
call system_clock(t_start)
call cpu_time(start_time)

write(*,'(A)') ' ============================================='
write(*,'(A)') '        SELECTED VCI (SCI@VSCF)              '
write(*,'(A)') '        Reference + EN-PT2 top-N select  '
write(*,'(A,I8,A)') '        N_sel_per_state = ', N_sel_per_state, ' configs/state'
write(*,'(A)') ' ============================================='

is_semiquartic_ff = .false.
call classify_quartic_terms(total_4, final_index_4, n_unique_4, &
      check4, Potential_4_vec, is_semiquartic_ff)
  to_sparse_cut = merge(3, 4, is_semiquartic_ff)
  write(*,'(A,I2)') ' Sparse pair-list cutoff set to n_diff <= ', to_sparse_cut

allocate(diff_modes(to_sparse_cut))

!==========================================================================
! STEP 1: Precompute modal integrals
!==========================================================================
allocate(modal_int(N_modes, 0:max_quanta, 0:max_quanta, 0:5))
modal_int = 0.d0

call omp_set_num_threads(INT(N_threads, KIND=4))

!$OMP PARALLEL DO PRIVATE(ii,i,j,p,mu,nu) COLLAPSE(2) SCHEDULE(static)
do ii = 1, N_modes
    do i = 0, max_quanta
        do j = i, max_quanta
            do p = 0, 5
                do mu = 1, N_expansion
                    do nu = 1, N_expansion
                        modal_int(ii,i,j,p) = modal_int(ii,i,j,p) &
                            + full_coef(ii,i+1,mu) * full_coef(ii,j+1,nu) &
                            * store_integrals(ii,mu,nu,p)
                    end do
                end do
            end do
            do p = 0, 5
                modal_int(ii,j,i,p) = modal_int(ii,i,j,p)
            end do
        end do
    end do
end do
!$OMP END PARALLEL DO

write(*,'(A)') ' Modal integrals precomputed.'

!==========================================================================
! STEP 2: Build dipole pair vector (second derivatives)
!==========================================================================
allocate(dipole_vec(N_modes*(N_modes+1)/2, 3))
allocate(final_index_dipole(N_modes*(N_modes+1)/2, 2))
dipole_vec = 0.d0

kk = 0
do i = 1, N_modes
    do j = i, N_modes
        kk = kk + 1
        dipole_vec(kk, 1:3)      = second_dipole_derivatives(i,j,1:3)
        final_index_dipole(kk,1) = i
        final_index_dipole(kk,2) = j
    end do
end do

write(*,'(A,I8,A,I8,A)') ' Potential terms: ', count(check3==1), &
    ' cubic, ', count(check4==1), ' quartic (nonzero)'
write(*,'(A,I6,A,I6)') ' Inverted index: max cubic/mode=', n_cubic_max, &
    ', max quartic/mode=', n_quartic_max

!==========================================================================
! STEP 3: Build CI reference space
!==========================================================================
if (check_list == 1) then

    allocate(is_ref(total_combinations2))
    allocate(is_selected(total_combinations2))
    is_ref      = .false.
    is_selected = .false.
    n_ref       = 0
    open (11, file='list_states.txt')

        read(11, *) n_states_list

    close(11)

    do m = 1, total_combinations ! Total_combinations in this case is just the number of combinations from the list, for the reference CI state, combination2 is for the all combinations
        do n = 1, n_states_list
            if (all(vec_combinations(m, :) == vec_combinations(n, :))) then
                is_ref(m)      = .true.
                is_selected(m) = .true.
                n_ref          = n_ref + 1
            end if
        end do
    end do

    allocate(ref_list(n_ref))
    j = 0
    do m = 1, total_combinations2
        if (is_ref(m)) then
            j = j + 1
            ref_list(j) = m
        end if
    end do


else 
    allocate(is_ref(total_combinations))
    allocate(is_selected(total_combinations))
    is_ref      = .false.
    is_selected = .false.
    n_ref       = 0

    do m = 1, total_combinations
        if (sum(vec_combinations(m, :)) <= reference_max) then
            is_ref(m)      = .true.
            is_selected(m) = .true.
            n_ref          = n_ref + 1
        end if
    end do

    allocate(ref_list(n_ref))
    j = 0
    do m = 1, total_combinations
        if (is_ref(m)) then
            j = j + 1
            ref_list(j) = m
        end if
    end do


end if

if (check_list == 1) then
    write(*,'(A,I8,A,I8,A)') ' LIST reference: ', n_ref, &
        ' configs out of ', total_combinations2, ' total'
    write(101,'(A,I8,A,I8,A)') ' LIST reference: ', n_ref, &
        ' configs out of ', total_combinations2, ' total'
else
    write(*,'(A,I8,A,I8,A)') ' CI   reference: ', n_ref, &
        ' configs out of ', total_combinations, ' total'
    write(101,'(A,I8,A,I8,A)') ' CI   reference: ', n_ref, &
        ' configs out of ', total_combinations, ' total'
end if 

!==========================================================================
! STEP 4: Build and diagonalize CISD Hamiltonian
!==========================================================================

!----------------------------------------------------------------------
! 4a: Build sel_list (= ref_list at this stage) and ext_list
!----------------------------------------------------------------------

if (check_list == 1) then
    n_sel = n_ref
    n_ext = total_combinations2 - n_sel
else 
    n_sel = n_ref
    n_ext = total_combinations - n_sel
end if

allocate(sel_list(n_sel))
allocate(ext_list(n_ext))

if (check_list == 1) then
    j = 0
    do m = 1, total_combinations2
        if (is_selected(m)) then
            j = j + 1
            sel_list(j) = m
        end if
    end do

    j = 0
    do m = 1, total_combinations2
        if (.not. is_selected(m)) then
            j = j + 1
            ext_list(j) = m
        end if
    end do

else 
    j = 0
    do m = 1, total_combinations
        if (is_selected(m)) then
            j = j + 1
            sel_list(j) = m
        end if
    end do

    j = 0
    do m = 1, total_combinations
        if (.not. is_selected(m)) then
            j = j + 1
            ext_list(j) = m
        end if
    end do

end if

write(*,'(A,I8,A,I8)') ' Reference =', n_sel, ', external=', n_ext
write(101,'(A,I8,A,I8)') ' Reference =', n_sel, ', external=', n_ext

!----------------------------------------------------------------------
! 4b: Build sparse pair list for the CISD space
!----------------------------------------------------------------------

if(check_list == 1) then 
    ! --- Count pass ---
    n_sparse = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff = 0
            do ii = 1, N_modes
                if (vec_combinations2(m,ii) /= vec_combinations2(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) n_sparse = n_sparse + 1
        end do
    end do

    allocate(sparse_m(n_sparse), sparse_n(n_sparse))
    allocate(sparse_ndiff(n_sparse))

    ! --- Fill pass ---
    idx = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff    = 0
            diff_modes = 0
            do ii = 1, N_modes
                if (vec_combinations2(m,ii) /= vec_combinations2(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff <= to_sparse_cut) diff_modes(n_diff) = ii
                    if (n_diff >  to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) then
                idx = idx + 1
                sparse_m(idx)          = j
                sparse_n(idx)          = i
                sparse_ndiff(idx)      = n_diff
            end if
        end do
    end do

    !----------------------------------------------------------------------
    ! 4c: Fill CISD H_sel
    !----------------------------------------------------------------------
    allocate(H_sel(n_sel, n_sel))
    H_sel = 0.d0

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
    !$OMP&        n_sel, sel_list, vec_combinations2, N_modes, max_quanta, &
    !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
    !$OMP&        Potential_3_vec, Potential_4_vec, &
    !$OMP&        check3, check4, total_3, total_4, &
    !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
    !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
    !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    !$OMP&        H_sel) &
    !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val) &
    !$OMP& SCHEDULE(dynamic, 64)
    do k = 1, n_sparse
        i = sparse_m(k)
        j = sparse_n(k)
        m = sel_list(i)
        n = sel_list(j)
        vm(1:N_modes) = vec_combinations2(m, 1:N_modes)
        vn(1:N_modes) = vec_combinations2(n, 1:N_modes)

        call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
            modal_int, Potential_3, Potential_4, HO_freq, &
            Potential_3_vec, Potential_4_vec, &
            check3, check4, total_3, total_4, &
            final_index_3, count_index_3, n_unique_3, unique_modes_3, &
            final_index_4, count_index_4, n_unique_4, unique_modes_4, &
            cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
            quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
            H_val)

        H_sel(i,j) = H_val
        H_sel(j,i) = H_val
    end do
    !$OMP END PARALLEL DO

    deallocate(sparse_m, sparse_n, sparse_ndiff)
end if

if (check_list == 0) then 
    !----------------------------------------------------------------------
    ! 4b: Build sparse pair list for the CISD space
    !----------------------------------------------------------------------

    ! --- Count pass ---
    n_sparse = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff = 0
            do ii = 1, N_modes
                if (vec_combinations(m,ii) /= vec_combinations(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) n_sparse = n_sparse + 1
        end do
    end do

    allocate(sparse_m(n_sparse), sparse_n(n_sparse))
    allocate(sparse_ndiff(n_sparse))

    ! --- Fill pass ---
    idx = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff    = 0
            diff_modes = 0
            do ii = 1, N_modes
                if (vec_combinations(m,ii) /= vec_combinations(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff <= to_sparse_cut) diff_modes(n_diff) = ii
                    if (n_diff >  to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) then
                idx = idx + 1
                sparse_m(idx)          = j
                sparse_n(idx)          = i
                sparse_ndiff(idx)      = n_diff
            end if
        end do
    end do

    !----------------------------------------------------------------------
    ! 4c: Fill CISD H_sel
    !----------------------------------------------------------------------
    allocate(H_sel(n_sel, n_sel))
    H_sel = 0.d0

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
    !$OMP&        n_sel, sel_list, vec_combinations, N_modes, max_quanta, &
    !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
    !$OMP&        Potential_3_vec, Potential_4_vec, &
    !$OMP&        check3, check4, total_3, total_4, &
    !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
    !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
    !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    !$OMP&        H_sel) &
    !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val) &
    !$OMP& SCHEDULE(dynamic, 64)
    do k = 1, n_sparse
        i = sparse_m(k)
        j = sparse_n(k)
        m = sel_list(i)
        n = sel_list(j)
        vm(1:N_modes) = vec_combinations(m, 1:N_modes)
        vn(1:N_modes) = vec_combinations(n, 1:N_modes)

        call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
            modal_int, Potential_3, Potential_4, HO_freq, &
            Potential_3_vec, Potential_4_vec, &
            check3, check4, total_3, total_4, &
            final_index_3, count_index_3, n_unique_3, unique_modes_3, &
            final_index_4, count_index_4, n_unique_4, unique_modes_4, &
            cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
            quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
            H_val)

        H_sel(i,j) = H_val
        H_sel(j,i) = H_val
    end do
    !$OMP END PARALLEL DO

    deallocate(sparse_m, sparse_n, sparse_ndiff)

end if

!----------------------------------------------------------------------
! 4d: Diagonalize CISD H_sel
!----------------------------------------------------------------------
allocate(eigenvalues(n_sel))
allocate(eigenvectors(n_sel, n_sel))
eigenvalues  = 0.d0
eigenvectors = 0.d0

if (n_sel <= N_states) then
    call dsyevd_A(H_sel, eigenvalues, eigenvectors)
else
    call dsyevr_A(H_sel, N_states, eigenvalues, eigenvectors, 'N')
end if

! Number of CISD eigenstates available for screening
n_cisd_states = min(N_states, n_sel)

! Print CISD energies (first 10 states)
write(*,'(A)') ' --- REF. energies (cm-1, relative to ZPE):'
write(101,'(A)') ' --- REF. energies (cm-1, relative to ZPE):'
do i = 1, min(10, n_cisd_states)
    write(*,'(I6,2F16.4)') i, &
        eigenvalues(i)/cm_to_hartree, &
        (eigenvalues(i)-eigenvalues(1))/cm_to_hartree
    write(101,'(I6,2F16.4)') i, &
        eigenvalues(i)/cm_to_hartree, &
        (eigenvalues(i)-eigenvalues(1))/cm_to_hartree
end do

!==========================================================================
! STEP 5: EN-PT2 screening — select top N_sel_per_state per CISD state
!==========================================================================

if (check_list == 0) then
    write(*,'(A,I8,A,I4,A)') ' Computing EN-PT2 contributions for ', n_ext, &
        ' external configs across ', n_cisd_states, ' REF. states'
    write(101,'(A,I8,A,I4,A)') ' Computing EN-PT2 contributions for ', n_ext, &
        ' external configs across ', n_cisd_states, ' REF. states'

    ! pt2_contrib(k, i_state) = |<alpha_k|H|Psi_I>|^2 / |E_I - H_aa|
    ! for external config k and CISD state i_state
    allocate(pt2_contrib(n_ext, n_cisd_states))
    pt2_contrib = 0.d0

    allocate(coupling_vec(n_cisd_states))

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& SHARED(n_ext, ext_list, n_sel, sel_list, vec_combinations, to_sparse_cut, &
    !$OMP&        N_modes, max_quanta, modal_int, &
    !$OMP&        Potential_3, Potential_4, HO_freq, &
    !$OMP&        Potential_3_vec, Potential_4_vec, &
    !$OMP&        check3, check4, total_3, total_4, &
    !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
    !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
    !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    !$OMP&        eigenvalues, eigenvectors, n_cisd_states, &
    !$OMP&        pt2_contrib) &
    !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val, H_aa, &
    !$OMP&         hrow, coupling_vec, e_alpha, denom, n_diff) &
    !$OMP& SCHEDULE(dynamic, 16)
    do k = 1, n_ext
        m  = ext_list(k)
        vm(1:N_modes) = vec_combinations(m, 1:N_modes)

        ! Compute diagonal element H_aa = <alpha|H|alpha>
        call compute_H_element(m, m, vm, vm, N_modes, max_quanta, &
            modal_int, Potential_3, Potential_4, HO_freq, &
            Potential_3_vec, Potential_4_vec, &
            check3, check4, total_3, total_4, &
            final_index_3, count_index_3, n_unique_3, unique_modes_3, &
            final_index_4, count_index_4, n_unique_4, unique_modes_4, &
            cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
            quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
            H_aa)

        ! Compute coupling row: hrow(j) = <alpha|H|phi_j> for all CISD configs j
        allocate(hrow(n_sel))
        hrow = 0.d0
        do j = 1, n_sel
            n  = sel_list(j)
            vn(1:N_modes) = vec_combinations(n, 1:N_modes)

            n_diff = 0
            do ii = 1, N_modes
                if (vm(ii) /= vn(ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > to_sparse_cut) exit
                end if
            end do
            if (n_diff > to_sparse_cut) cycle

            call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
                H_val)
            hrow(j) = H_val
        end do

        ! Compute coupling to each CISD eigenstate: <alpha|H|Psi_I> = sum_j c_jI * hrow(j)
        coupling_vec(1:n_cisd_states) = 0.d0
        do j = 1, n_sel
            if (abs(hrow(j)) < 1.d-30) cycle
            do i = 1, n_cisd_states
                coupling_vec(i) = coupling_vec(i) + eigenvectors(j,i) * hrow(j) !h_row = Hmn = state inside the CISD and outisde, eigenvectoirs = Cn the CISD coefficient
            end do
        end do

        ! Compute EN-PT2 contribution for each CISD state
        do i = 1, n_cisd_states
            denom = eigenvalues(i) - H_aa
            if (abs(denom) < 1.d-12) then
                pt2_contrib(k, i) = 0.d0
            else
                pt2_contrib(k, i) = abs(coupling_vec(i) * coupling_vec(i) / denom)
            end if
        end do

        deallocate(hrow)

    end do
    !$OMP END PARALLEL DO

    deallocate(coupling_vec)

end if

if (check_list == 1) then
    write(*,'(A,I8,A,I4,A)') ' Computing EN-PT2 contributions for ', n_ext, &
        ' external configs across ', n_cisd_states, ' LIST states'
    write(101,'(A,I8,A,I4,A)') ' Computing EN-PT2 contributions for ', n_ext, &
        ' external configs across ', n_cisd_states, ' LIST states'

    ! pt2_contrib(k, i_state) = |<alpha_k|H|Psi_I>|^2 / |E_I - H_aa|
    ! for external config k and CISD state i_state
    allocate(pt2_contrib(n_ext, n_cisd_states))
    pt2_contrib = 0.d0

    allocate(coupling_vec(n_cisd_states))

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& SHARED(n_ext, ext_list, n_sel, sel_list, vec_combinations2, to_sparse_cut, &
    !$OMP&        N_modes, max_quanta, modal_int, &
    !$OMP&        Potential_3, Potential_4, HO_freq, &
    !$OMP&        Potential_3_vec, Potential_4_vec, &
    !$OMP&        check3, check4, total_3, total_4, &
    !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
    !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
    !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    !$OMP&        eigenvalues, eigenvectors, n_cisd_states, &
    !$OMP&        pt2_contrib) &
    !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val, H_aa, &
    !$OMP&         hrow, coupling_vec, e_alpha, denom, n_diff) &
    !$OMP& SCHEDULE(dynamic, 16)
    do k = 1, n_ext
        m  = ext_list(k)
        vm(1:N_modes) = vec_combinations2(m, 1:N_modes)

        ! Compute diagonal element H_aa = <alpha|H|alpha>
        call compute_H_element(m, m, vm, vm, N_modes, max_quanta, &
            modal_int, Potential_3, Potential_4, HO_freq, &
            Potential_3_vec, Potential_4_vec, &
            check3, check4, total_3, total_4, &
            final_index_3, count_index_3, n_unique_3, unique_modes_3, &
            final_index_4, count_index_4, n_unique_4, unique_modes_4, &
            cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
            quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
            H_aa)

        ! Compute coupling row: hrow(j) = <alpha|H|phi_j> for all CISD configs j
        allocate(hrow(n_sel))
        hrow = 0.d0
        do j = 1, n_sel
            n  = sel_list(j)
            vn(1:N_modes) = vec_combinations2(n, 1:N_modes)

            n_diff = 0
            do ii = 1, N_modes
                if (vm(ii) /= vn(ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > to_sparse_cut) exit
                end if
            end do
            if (n_diff > to_sparse_cut) cycle

            call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
                H_val)
            hrow(j) = H_val
        end do

        ! Compute coupling to each CISD eigenstate: <alpha|H|Psi_I> = sum_j c_jI * hrow(j)
        coupling_vec(1:n_cisd_states) = 0.d0
        do j = 1, n_sel
            if (abs(hrow(j)) < 1.d-30) cycle
            do i = 1, n_cisd_states
                coupling_vec(i) = coupling_vec(i) + eigenvectors(j,i) * hrow(j) !h_row = Hmn = state inside the CISD and outisde, eigenvectoirs = Cn the CISD coefficient
            end do
        end do

        ! Compute EN-PT2 contribution for each CISD state
        do i = 1, n_cisd_states
            denom = eigenvalues(i) - H_aa
            if (abs(denom) < 1.d-12) then
                pt2_contrib(k, i) = 0.d0
            else
                pt2_contrib(k, i) = abs(coupling_vec(i) * coupling_vec(i) / denom)
            end if
        end do

        deallocate(hrow)

    end do
    !$OMP END PARALLEL DO

    deallocate(coupling_vec)

end if


!----------------------------------------------------------------------
! 5b: For each CISD state, find the top N_sel_per_state external configs
!     and mark them as selected (union across all states)
!----------------------------------------------------------------------
! top_n_indices(rank, state) = index into ext_list of the rank-th most important config
! top_n_values(rank, state)  = corresponding pt2 contribution value


allocate(top_n_indices(N_sel_per_state, n_cisd_states))
allocate(top_n_values(N_sel_per_state, n_cisd_states))
top_n_indices = 0
top_n_values  = 0.d0

! For each CISD state, find the N_sel_per_state largest pt2_contrib entries
do i_state = 1, n_cisd_states
    ! Initialize with first N_sel_per_state entries (or fewer if n_ext is small)
    do k = 1, min(N_sel_per_state, n_ext)
        top_n_values(k, i_state)  = pt2_contrib(k, i_state)
        top_n_indices(k, i_state) = k
    end do

    ! Simple insertion sort to maintain top-N: sort initial entries descending
    do i = 1, min(N_sel_per_state, n_ext)
        do j = i + 1, min(N_sel_per_state, n_ext)
            if (top_n_values(j, i_state) > top_n_values(i, i_state)) then
                e_alpha = top_n_values(i, i_state)
                top_n_values(i, i_state) = top_n_values(j, i_state)
                top_n_values(j, i_state) = e_alpha
                idx = top_n_indices(i, i_state)
                top_n_indices(i, i_state) = top_n_indices(j, i_state)
                top_n_indices(j, i_state) = idx
            end if
        end do
    end do

    ! Scan remaining external configs and insert if larger than current minimum
    do k = N_sel_per_state + 1, n_ext
        ! The minimum of the top-N is at position N_sel_per_state (sorted descending)
        if (pt2_contrib(k, i_state) > top_n_values(N_sel_per_state, i_state)) then
            ! Replace the smallest entry
            top_n_values(N_sel_per_state, i_state)  = pt2_contrib(k, i_state)
            top_n_indices(N_sel_per_state, i_state) = k

            ! Bubble up to maintain descending order
            do i = N_sel_per_state - 1, 1, -1
                if (top_n_values(i+1, i_state) > top_n_values(i, i_state)) then
                    e_alpha = top_n_values(i, i_state)
                    top_n_values(i, i_state) = top_n_values(i+1, i_state)
                    top_n_values(i+1, i_state) = e_alpha
                    idx = top_n_indices(i, i_state)
                    top_n_indices(i, i_state) = top_n_indices(i+1, i_state)
                    top_n_indices(i+1, i_state) = idx
                else
                    exit
                end if
            end do
        end if
    end do
end do

! Mark the union of all top-N configurations as selected
n_new = 0
do i_state = 1, n_cisd_states
    do k = 1, min(N_sel_per_state, n_ext)
        if (top_n_indices(k, i_state) > 0) then
            m = ext_list(top_n_indices(k, i_state))
            if (.not. is_selected(m)) then
                is_selected(m) = .true.
                n_new = n_new + 1
            end if
        end if
    end do
end do

deallocate(pt2_contrib, top_n_indices, top_n_values)

write(*,'(A,I8,A,I4,A)') '  --> ', n_new, &
    ' new configurations selected (union across ', n_cisd_states, ' states)'
write(101,'(A,I8,A,I4,A)') '  --> ', n_new, &
    ' new configurations selected (union across ', n_cisd_states, ' states)'

! Free CISD arrays no longer needed
deallocate(sel_list, ext_list, H_sel, eigenvalues, eigenvectors)

!==========================================================================
! STEP 6: Build and diagonalize the final expanded Hamiltonian
!==========================================================================

if (check_list == 0) then
    ! Build the final sel_list from the union of CISD + selected externals
    n_sel = count(is_selected)

    allocate(sel_list(n_sel))
    j = 0
    do m = 1, total_combinations
        if (is_selected(m)) then
            j = j + 1
            sel_list(j) = m
        end if
    end do

    write(*,'(A,I8,A)') ' Final active space: ', n_sel, ' configurations'
    write(*,'(A,I8,A,I8,A)') '   (REF.: ', n_ref, ' + selected: ', n_new, ')'
    write(101,'(A,I8,A)') ' Final active space: ', n_sel, ' configurations'
    write(101,'(A,I8,A,I8,A)') '   (REF.: ', n_ref, ' + selected: ', n_new, ')'

    ! Build sparse pair list for the final active space
    n_sparse = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff = 0
            do ii = 1, N_modes
                if (vec_combinations(m,ii) /= vec_combinations(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) n_sparse = n_sparse + 1
        end do
    end do

    allocate(sparse_m(n_sparse), sparse_n(n_sparse))
    allocate(sparse_ndiff(n_sparse))

    idx = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff    = 0
            diff_modes = 0
            do ii = 1, N_modes
                if (vec_combinations(m,ii) /= vec_combinations(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff <= to_sparse_cut) diff_modes(n_diff) = ii
                    if (n_diff >  to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) then
                idx = idx + 1
                sparse_m(idx)              = j
                sparse_n(idx)              = i
                sparse_ndiff(idx)          = n_diff
            end if
        end do
    end do
    

    if (.not. use_davidson) then
        ! Fill final H_sel
        allocate(H_sel(n_sel, n_sel))
        H_sel = 0.d0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
        !$OMP&        n_sel, sel_list, vec_combinations, N_modes, max_quanta, &
        !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
        !$OMP&        Potential_3_vec, Potential_4_vec, &
        !$OMP&        check3, check4, total_3, total_4, &
        !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
        !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
        !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
        !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
        !$OMP&        H_sel) &
        !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val) &
        !$OMP& SCHEDULE(dynamic, 64)
        do k = 1, n_sparse
            i = sparse_m(k)
            j = sparse_n(k)
            m = sel_list(i)
            n = sel_list(j)
            vm(1:N_modes) = vec_combinations(m, 1:N_modes)
            vn(1:N_modes) = vec_combinations(n, 1:N_modes)

            call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
                H_val)
            H_sel(i,j) = H_val
            H_sel(j,i) = H_val
        end do
        !$OMP END PARALLEL DO

    else !conditional to use DAV
        allocate(H_sparse(n_sparse))
        H_sparse = 0.d0
        write(*,'(A)') ' Precomputing sparse Hamiltonian matrix (Davidson)...'
                !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
        !$OMP&        n_sel, sel_list, vec_combinations, N_modes, max_quanta, &
        !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
        !$OMP&        Potential_3_vec, Potential_4_vec, &
        !$OMP&        check3, check4, total_3, total_4, &
        !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
        !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
        !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
        !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
        !$OMP&        H_sparse) &
        !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val) &
        !$OMP& SCHEDULE(dynamic, 64)
        do k = 1, n_sparse
            i = sparse_m(k)
            j = sparse_n(k)
            m = sel_list(i)
            n = sel_list(j)
            vm(1:N_modes) = vec_combinations(m, 1:N_modes)
            vn(1:N_modes) = vec_combinations(n, 1:N_modes)

            call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
                H_val)
            H_sparse(k) = H_val
        end do
        !$OMP END PARALLEL DO
        write(*,'(A)') ' Sparse Hamiltonian precomputed.'

        ! Build a SEPARATE, compacted copy for the Davidson matvec only
        write(*,'(A,ES10.3)') ' Filtering elements below threshold (Davidson H only): ', H_threshold
        allocate(sparse_m_dav(n_sparse), sparse_n_dav(n_sparse))
        n_sparse_dav = 0
        do idx = 1, n_sparse
           if (abs(H_sparse(idx)) >= H_threshold) then
              n_sparse_dav = n_sparse_dav + 1
              sparse_m_dav(n_sparse_dav) = sparse_m(idx)
              sparse_n_dav(n_sparse_dav) = sparse_n(idx)
              H_sparse(n_sparse_dav)     = H_sparse(idx)
           end if
        end do
        write(*,'(A,I0,A,I0,A,F6.2,A)') &
             ' Kept ', n_sparse_dav, ' / ', n_sparse, &
             ' elements (', 100.0d0 * real(n_sparse_dav,8) / real(n_sparse,8), '% retained)'

        ! Reallocate to exact size
        allocate(tmp_sparse_int(n_sparse_dav))
        tmp_sparse_int(:) = sparse_m_dav(1:n_sparse_dav)
        deallocate(sparse_m_dav)
        allocate(sparse_m_dav(n_sparse_dav))
        sparse_m_dav(:) = tmp_sparse_int(:)
        deallocate(tmp_sparse_int)

        allocate(tmp_sparse_int(n_sparse_dav))
        tmp_sparse_int(:) = sparse_n_dav(1:n_sparse_dav)
        deallocate(sparse_n_dav)
        allocate(sparse_n_dav(n_sparse_dav))
        sparse_n_dav(:) = tmp_sparse_int(:)
        deallocate(tmp_sparse_int)

        allocate(tmp_sparse_real(n_sparse_dav))
        tmp_sparse_real(:) = H_sparse(1:n_sparse_dav)
        deallocate(H_sparse)
        allocate(H_sparse(n_sparse_dav))
        H_sparse(:) = tmp_sparse_real(:)
        deallocate(tmp_sparse_real)

        write(*,'(A,I0)') ' Final number of Davidson sparse elements: ', n_sparse_dav

    end if

    if (.not. use_davidson) then
        ! Diagonalize the final Hamiltonian
        allocate(eigenvalues(n_sel), eigenvectors(n_sel, n_sel))
        eigenvalues  = 0.d0
        eigenvectors = 0.d0

        if (n_sel <= N_states) then
            call dsyevd_A(H_sel, eigenvalues, eigenvectors)
        else
            call dsyevr_A(H_sel, N_states, eigenvalues, eigenvectors, 'N')
        end if
        deallocate(H_sel)
    else 
        if (Nfirst_davidson .gt. n_sel) Nfirst_davidson = n_sel
        call system_clock(count_rate=count_rate_dav, count_max=count_max_dav)
        call system_clock(t_start_dav)
        call cpu_time(start_time_dav)
        allocate(eigenvalues(Nfirst_davidson), eigenvectors(n_sel, Nfirst_davidson))
        write(*,'(1A,I4)') ' >>>>> Using Davidson eigensolver. Roots requested: ', Nfirst_davidson
        call jacobi_davidson_eigensolver(n_sel, Nfirst_davidson, &
            n_sparse_dav, sparse_m_dav, sparse_n_dav, H_sparse, &
            davidson_tol, davidson_max_iter, eigenvalues, eigenvectors)
        write(*,'(1A)') ' Davidson eigensolver finished.'
        deallocate(H_sparse, sparse_m_dav, sparse_n_dav)

        call cpu_time(end_time_dav)
        call system_clock(t_end_dav)
        elapsed_time_dav = real(t_end_dav - t_start_dav) / real(count_rate_dav)
        write(*,'(1A, 1F12.2, 1A)') " DAV CPU time:     ", end_time_dav - start_time_dav, " seconds"
        write(*,'(1A, 1F12.2, 1A)') " DAV elapsed time: ", elapsed_time_dav, " seconds"
    end if

end if !conditional to check list

if (check_list == 1) then
    ! Build the final sel_list from the union of CISD + selected externals
    n_sel = count(is_selected)

    allocate(sel_list(n_sel))
    j = 0
    do m = 1, total_combinations2
        if (is_selected(m)) then
            j = j + 1
            sel_list(j) = m
        end if
    end do

    write(*,'(A,I8,A)') ' Final active space: ', n_sel, ' configurations'
    write(*,'(A,I8,A,I8,A)') '   (REF.: ', n_ref, ' + selected: ', n_new, ')'
    write(101,'(A,I8,A)') ' Final active space: ', n_sel, ' configurations'
    write(101,'(A,I8,A,I8,A)') '   (REF.: ', n_ref, ' + selected: ', n_new, ')'

    ! Build sparse pair list for the final active space
    n_sparse = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff = 0
            do ii = 1, N_modes
                if (vec_combinations2(m,ii) /= vec_combinations2(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) n_sparse = n_sparse + 1
        end do
    end do

    allocate(sparse_m(n_sparse), sparse_n(n_sparse))
    allocate(sparse_ndiff(n_sparse))

    idx = 0
    do j = 1, n_sel
        do i = j, n_sel
            m = sel_list(j) ; n = sel_list(i)
            n_diff    = 0
            diff_modes = 0
            do ii = 1, N_modes
                if (vec_combinations2(m,ii) /= vec_combinations2(n,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff <= to_sparse_cut) diff_modes(n_diff) = ii
                    if (n_diff >  to_sparse_cut) exit
                end if
            end do
            if (n_diff <= to_sparse_cut) then
                idx = idx + 1
                sparse_m(idx)              = j
                sparse_n(idx)              = i
                sparse_ndiff(idx)          = n_diff
            end if
        end do
    end do
    

    if (.not. use_davidson) then
        ! Fill final H_sel
        allocate(H_sel(n_sel, n_sel))
        H_sel = 0.d0

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
        !$OMP&        n_sel, sel_list, vec_combinations2, N_modes, max_quanta, &
        !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
        !$OMP&        Potential_3_vec, Potential_4_vec, &
        !$OMP&        check3, check4, total_3, total_4, &
        !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
        !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
        !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
        !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
        !$OMP&        H_sel) &
        !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val) &
        !$OMP& SCHEDULE(dynamic, 64)
        do k = 1, n_sparse
            i = sparse_m(k)
            j = sparse_n(k)
            m = sel_list(i)
            n = sel_list(j)
            vm(1:N_modes) = vec_combinations2(m, 1:N_modes)
            vn(1:N_modes) = vec_combinations2(n, 1:N_modes)

            call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
                H_val)
            H_sel(i,j) = H_val
            H_sel(j,i) = H_val
        end do
        !$OMP END PARALLEL DO

    else !conditional to use DAV
        allocate(H_sparse(n_sparse))
        H_sparse = 0.d0
        write(*,'(A)') ' Precomputing sparse Hamiltonian matrix (Davidson)...'
                !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
        !$OMP&        n_sel, sel_list, vec_combinations2, N_modes, max_quanta, &
        !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
        !$OMP&        Potential_3_vec, Potential_4_vec, &
        !$OMP&        check3, check4, total_3, total_4, &
        !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
        !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
        !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
        !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
        !$OMP&        H_sparse) &
        !$OMP& PRIVATE(k, i, j, m, n, ii, vm, vn, H_val) &
        !$OMP& SCHEDULE(dynamic, 64)
        do k = 1, n_sparse
            i = sparse_m(k)
            j = sparse_n(k)
            m = sel_list(i)
            n = sel_list(j)
            vm(1:N_modes) = vec_combinations2(m, 1:N_modes)
            vn(1:N_modes) = vec_combinations2(n, 1:N_modes)

            call compute_H_element(m, n, vm, vn, N_modes, max_quanta, &
                modal_int, Potential_3, Potential_4, HO_freq, &
                Potential_3_vec, Potential_4_vec, &
                check3, check4, total_3, total_4, &
                final_index_3, count_index_3, n_unique_3, unique_modes_3, &
                final_index_4, count_index_4, n_unique_4, unique_modes_4, &
                cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
                quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
                H_val)
            H_sparse(k) = H_val
        end do
        !$OMP END PARALLEL DO
        write(*,'(A)') ' Sparse Hamiltonian precomputed.'

        ! Build a SEPARATE, compacted copy for the Davidson matvec only
        write(*,'(A,ES10.3)') ' Filtering elements below threshold (Davidson H only): ', H_threshold
        allocate(sparse_m_dav(n_sparse), sparse_n_dav(n_sparse))
        n_sparse_dav = 0
        do idx = 1, n_sparse
           if (abs(H_sparse(idx)) >= H_threshold) then
              n_sparse_dav = n_sparse_dav + 1
              sparse_m_dav(n_sparse_dav) = sparse_m(idx)
              sparse_n_dav(n_sparse_dav) = sparse_n(idx)
              H_sparse(n_sparse_dav)     = H_sparse(idx)
           end if
        end do
        write(*,'(A,I0,A,I0,A,F6.2,A)') &
             ' Kept ', n_sparse_dav, ' / ', n_sparse, &
             ' elements (', 100.0d0 * real(n_sparse_dav,8) / real(n_sparse,8), '% retained)'

        ! Reallocate to exact size
        allocate(tmp_sparse_int(n_sparse_dav))
        tmp_sparse_int(:) = sparse_m_dav(1:n_sparse_dav)
        deallocate(sparse_m_dav)
        allocate(sparse_m_dav(n_sparse_dav))
        sparse_m_dav(:) = tmp_sparse_int(:)
        deallocate(tmp_sparse_int)

        allocate(tmp_sparse_int(n_sparse_dav))
        tmp_sparse_int(:) = sparse_n_dav(1:n_sparse_dav)
        deallocate(sparse_n_dav)
        allocate(sparse_n_dav(n_sparse_dav))
        sparse_n_dav(:) = tmp_sparse_int(:)
        deallocate(tmp_sparse_int)

        allocate(tmp_sparse_real(n_sparse_dav))
        tmp_sparse_real(:) = H_sparse(1:n_sparse_dav)
        deallocate(H_sparse)
        allocate(H_sparse(n_sparse_dav))
        H_sparse(:) = tmp_sparse_real(:)
        deallocate(tmp_sparse_real)

        write(*,'(A,I0)') ' Final number of Davidson sparse elements: ', n_sparse_dav

    end if

    if (.not. use_davidson) then
        ! Diagonalize the final Hamiltonian
        allocate(eigenvalues(n_sel), eigenvectors(n_sel, n_sel))
        eigenvalues  = 0.d0
        eigenvectors = 0.d0

        if (n_sel <= N_states) then
            call dsyevd_A(H_sel, eigenvalues, eigenvectors)
        else
            call dsyevr_A(H_sel, N_states, eigenvalues, eigenvectors, 'N')
        end if
        deallocate(H_sel)
    else 
        if (Nfirst_davidson .gt. n_sel) Nfirst_davidson = n_sel
        call system_clock(count_rate=count_rate_dav, count_max=count_max_dav)
        call system_clock(t_start_dav)
        call cpu_time(start_time_dav)
        allocate(eigenvalues(Nfirst_davidson), eigenvectors(n_sel, Nfirst_davidson))
        write(*,'(1A,I4)') ' >>>>> Using Davidson eigensolver. Roots requested: ', Nfirst_davidson
        call jacobi_davidson_eigensolver(n_sel, Nfirst_davidson, &
            n_sparse_dav, sparse_m_dav, sparse_n_dav, H_sparse, &
            davidson_tol, davidson_max_iter, eigenvalues, eigenvectors)
        write(*,'(1A)') ' Davidson eigensolver finished.'
        deallocate(H_sparse, sparse_m_dav, sparse_n_dav)

        call cpu_time(end_time_dav)
        call system_clock(t_end_dav)
        elapsed_time_dav = real(t_end_dav - t_start_dav) / real(count_rate_dav)
        write(*,'(1A, 1F12.2, 1A)') " DAV CPU time:     ", end_time_dav - start_time_dav, " seconds"
        write(*,'(1A, 1F12.2, 1A)') " DAV elapsed time: ", elapsed_time_dav, " seconds"
    end if

end if !conditional to check list

!==========================================================================
! STEP 7: Final output — mirror VCI format exactly
!==========================================================================
to_write_output = n_sel
if (use_davidson) to_write_output = Nfirst_davidson

if(check_list == 0) then

    write(*,'(A)') ' ============================================='
    write(*,'(A)') '   FINAL VCI ENERGIES'
    write(*,'(A)') ' ============================================='
    write(*,'(A,I8)') ' Final active space size: ', n_sel
    write(101,'(A)') ' ============================================='
    write(101,'(A)') '   FINAL VCI ENERGIES'
    write(101,'(A)') ' ============================================='
    write(101,'(A,I8)') ' Final active space size: ', n_sel

   ! write(*,*)
   ! write(*,*) '>> Final List' 
   ! do i = 1, n_sel
   !     write(*,'(99999I5)') sel_list(i), vec_combinations(sel_list(i),:)
   ! end do
   ! write(*,*)

    write(101,*)
    write(101,*) '>> Final List' 
    do i = 1, n_sel
        write(101,'(99999I5)') sel_list(i), vec_combinations(sel_list(i),:)
    end do
    write(101,*)
    if (to_write_output <= N_states) then
        write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(*,'(1A19, 1A22, 1A14, 6A12)') 'E (cm-1)', 'E - ZPE', &
            '              ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C'
        write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(101,'(1A19, 1A22, 1A14, 6A12)') 'E (cm-1)', 'E - ZPE', &
            '              ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C'
        write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'

        do i = 1, min(N_states, to_write_output)
            values  = 0.d0
            indices = 0
            call top_three_components(eigenvectors(:,i), n_sel, values, indices)
            if (i < 10) then
                write(*,'(1F19.4, 1F22.4, 1A14, 3F12.6, 3I12)') &
                    eigenvalues(i)/cm_to_hartree, &
                    eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree, &
                    '      ||      ', values, sel_list(indices(1)), sel_list(indices(2)), sel_list(indices(3))
            end if
            write(101,'(1F19.4, 1F22.4, 1A14, 3F12.6, 3I12)') &
                eigenvalues(i)/cm_to_hartree, &
                eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree, &
                '      ||      ', values, sel_list(indices(1)), sel_list(indices(2)), sel_list(indices(3))
        end do
        write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'

    else
        write(*,'(1A)') '-----------------------------------------'
        write(*,'(1A19, 1A22)') 'E (cm-1)', 'E - ZPE'
        write(*,'(1A)') '-----------------------------------------'
        write(101,'(1A)') '-----------------------------------------'
        write(101,'(1A19, 1A22)') 'E (cm-1)', 'E - ZPE'
        write(101,'(1A)') '-----------------------------------------'

        do i = 1, min(N_states, to_write_output)
            values  = 0.d0
            indices = 0
            call top_three_components(eigenvectors(:,i), min(N_states,to_write_output), values, indices)
            if (i < 10) then
                write(*,'(1F19.4, 1F22.4)') &
                    eigenvalues(i)/cm_to_hartree, &
                    eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree
            end if
            write(101,'(1F19.4, 1F22.4)') &
                eigenvalues(i)/cm_to_hartree, &
                eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree
        end do
        write(*,'(1A)') '-----------------------------------------'
        write(101,'(1A)') '-----------------------------------------'

    end if

end if

if (check_list == 1) then 

    write(*,'(A)') ' ============================================='
    write(*,'(A)') '   FINAL VCI ENERGIES'
    write(*,'(A)') ' ============================================='
    write(*,'(A,I8)') ' Final active space size: ', n_sel
    write(101,'(A)') ' ============================================='
    write(101,'(A)') '   FINAL VCI ENERGIES'
    write(101,'(A)') ' ============================================='
    write(101,'(A,I8)') ' Final active space size: ', n_sel

    !write(*,*)
    !write(*,*) '>> Final List' 
    !do i = 1, n_sel
    !    write(*,'(99999I5)') sel_list(i), vec_combinations2(sel_list(i),:)
    !end do
    !write(*,*)

    write(101,*)
    write(101,*) '>> Final List' 
    do i = 1, n_sel
        write(101,'(99999I5)') sel_list(i), vec_combinations2(sel_list(i),:)
    end do
    write(101,*)

    if (to_write_output <= N_states) then
        write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(*,'(1A19, 1A22, 1A14, 6A12)') 'E (cm-1)', 'E - ZPE', &
            '              ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C'
        write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(101,'(1A19, 1A22, 1A14, 6A12)') 'E (cm-1)', 'E - ZPE', &
            '              ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C'
        write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'

        do i = 1, min(N_states, to_write_output)
            values  = 0.d0
            indices = 0
            call top_three_components(eigenvectors(:,i), n_sel, values, indices)
            if (i < 10) then
                write(*,'(1F19.4, 1F22.4, 1A14, 3F12.6, 3I12)') &
                    eigenvalues(i)/cm_to_hartree, &
                    eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree, &
                    '      ||      ', values, sel_list(indices(1)), sel_list(indices(2)), sel_list(indices(3))
            end if
            write(101,'(1F19.4, 1F22.4, 1A14, 3F12.6, 3I12)') &
                eigenvalues(i)/cm_to_hartree, &
                eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree, &
                '      ||      ', values, sel_list(indices(1)), sel_list(indices(2)), sel_list(indices(3))
        end do
        write(*,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'
        write(101,'(1A)') '-------------------------------------------------------------------------------------------------------------------------------'

    else
        write(*,'(1A)') '-----------------------------------------'
        write(*,'(1A19, 1A22)') 'E (cm-1)', 'E - ZPE'
        write(*,'(1A)') '-----------------------------------------'
        write(101,'(1A)') '-----------------------------------------'
        write(101,'(1A19, 1A22)') 'E (cm-1)', 'E - ZPE'
        write(101,'(1A)') '-----------------------------------------'

        do i = 1, min(N_states, to_write_output)
            values  = 0.d0
            indices = 0
            call top_three_components(eigenvectors(:,i), min(N_states,to_write_output), values, indices)
            if (i < 10) then
                write(*,'(1F19.4, 1F22.4)') &
                    eigenvalues(i)/cm_to_hartree, &
                    eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree
            end if
            write(101,'(1F19.4, 1F22.4)') &
                eigenvalues(i)/cm_to_hartree, &
                eigenvalues(i)/cm_to_hartree - eigenvalues(1)/cm_to_hartree
        end do
        write(*,'(1A)') '-----------------------------------------'
        write(101,'(1A)') '-----------------------------------------'

    end if

end if

!==========================================================================
! STEP 8: Transition dipoles for the final selected space
!==========================================================================

if (check_list == 0) then
    write(*,*) ' -Calculating SCI transition dipoles (onthefly, no full D matrix)'

    !----------------------------------------------------------------------
    ! 8a: Build dipole-specific sparse list (n_diff <= 2) from the pair
    !     list already built for H in STEP 6
    !----------------------------------------------------------------------
    n_sparse_dip = 0
    do idx = 1, n_sparse
        if (sparse_ndiff(idx) <= 2) n_sparse_dip = n_sparse_dip + 1
    end do

    allocate(sparse_m_dip(n_sparse_dip), sparse_n_dip(n_sparse_dip), sparse_ndiff_dip(n_sparse_dip))
    n_sparse_dip = 0
    do idx = 1, n_sparse
        if (sparse_ndiff(idx) <= 2) then
            n_sparse_dip = n_sparse_dip + 1
            sparse_m_dip(n_sparse_dip) = sparse_m(idx)   ! LOCAL indices (1..n_sel)
            sparse_n_dip(n_sparse_dip) = sparse_n(idx)   ! LOCAL indices (1..n_sel)
            sparse_ndiff_dip(n_sparse_dip) = sparse_ndiff(idx)
        end if
    end do
    write(*,'(A,I0,A)') ' Dipole-specific pair list (SCI/CI space) contains ', n_sparse_dip, ' elements.'

    !----------------------------------------------------------------------
    ! 8b: On-the-fly D*v0 accumulation
    !----------------------------------------------------------------------
    allocate(d_x(n_sel), d_y(n_sel), d_z(n_sel))
    d_x = 0.d0; d_y = 0.d0; d_z = 0.d0

    do idx = 1, n_sparse_dip
        i     = sparse_m_dip(idx)          
        j     = sparse_n_dip(idx)          
        n_diff = sparse_ndiff_dip(idx)
        m     = sel_list(i)               
        n_loc = sel_list(j)                
        vm(1:N_modes) = vec_combinations(m,     1:N_modes)
        vn(1:N_modes) = vec_combinations(n_loc, 1:N_modes)

        do ii = 1, N_modes
            ovlp_loc(ii) = modal_int(ii, vm(ii), vn(ii), 0)
        end do
        prefix_prod(0) = 1.0d0
        do ii = 1, N_modes
            prefix_prod(ii) = prefix_prod(ii-1) * ovlp_loc(ii)
        end do
        prod_all = prefix_prod(N_modes)
        suffix_prod_loc(N_modes+1) = 1.0d0
        do ii = N_modes, 1, -1
            suffix_prod_loc(ii) = suffix_prod_loc(ii+1) * ovlp_loc(ii)
        end do
        do ii = 1, N_modes
            prod_except_i(ii) = prefix_prod(ii-1) * suffix_prod_loc(ii+1)
        end do

        dip_local = 0.d0

        if (n_diff <= 1) then
            do ii = 1, N_modes
                if (abs(prod_except_i(ii)) < 1.d-30) cycle
                dip_local(1:3) = dip_local(1:3) + ( &
                    modal_int(ii, vm(ii), vn(ii), 1) * dipole_derivatives(ii, 1:3) &
                  + modal_int(ii, vm(ii), vn(ii), 2) * second_dipole_derivatives(ii,ii,1:3) &
                    ) * prod_except_i(ii)
            end do
        end if

        do kk = 1, N_modes*(N_modes+1)/2
            if (final_index_dipole(kk,1) == final_index_dipole(kk,2)) cycle
            if (abs(dipole_vec(kk,1)) < 1.d-30 .and. &
                abs(dipole_vec(kk,2)) < 1.d-30 .and. &
                abs(dipole_vec(kk,3)) < 1.d-30) cycle
            ii = final_index_dipole(kk,1)
            jj = final_index_dipole(kk,2)
            if (abs(ovlp_loc(ii)) < 1.d-300 .or. abs(ovlp_loc(jj)) < 1.d-300) then
                step_prod_loc = 1.0d0
                do pp = 1, N_modes
                    if (pp /= ii .and. pp /= jj) step_prod_loc = step_prod_loc * ovlp_loc(pp)
                end do
            else
                step_prod_loc = prod_all / (ovlp_loc(ii) * ovlp_loc(jj))
            end if
            if (abs(step_prod_loc) < 1.d-30) cycle
            step_prod_loc = step_prod_loc &
                * modal_int(ii, vm(ii), vn(ii), 1) * modal_int(jj, vm(jj), vn(jj), 1)
            dip_local(1:3) = dip_local(1:3) + dipole_vec(kk,1:3) * step_prod_loc * 2.0d0
        end do

        d_x(i) = d_x(i) + dip_local(1) * eigenvectors(j,1)
        d_y(i) = d_y(i) + dip_local(2) * eigenvectors(j,1)
        d_z(i) = d_z(i) + dip_local(3) * eigenvectors(j,1)
        if (i /= j) then
            d_x(j) = d_x(j) + dip_local(1) * eigenvectors(i,1)
            d_y(j) = d_y(j) + dip_local(2) * eigenvectors(i,1)
            d_z(j) = d_z(j) + dip_local(3) * eigenvectors(i,1)
        end if
    end do

    write(*,*) 'SCI dipole accumulation finished.'

    number_to_print_int = 0
    do i = 1, min(N_states, to_write_output) - 1
        if ((eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree < max_freq) &
            number_to_print_int = number_to_print_int + 1
    end do
    if (number_to_print_int < 1) number_to_print_int = 1
    if (number_to_print_int >= to_write_output) number_to_print_int = to_write_output - 1

    allocate(dipole_final(number_to_print_int, 3))
    allocate(intensity(number_to_print_int))

    do i = 1, number_to_print_int
        dipole_final(i,1) = ddot(n_sel, eigenvectors(1,i+1), 1, d_x, 1)
        dipole_final(i,2) = ddot(n_sel, eigenvectors(1,i+1), 1, d_y, 1)
        dipole_final(i,3) = ddot(n_sel, eigenvectors(1,i+1), 1, d_z, 1)
        intensity(i) = sum(dipole_final(i,:)**2)
    end do

    deallocate(d_x, d_y, d_z)
    deallocate(sparse_m_dip, sparse_n_dip, sparse_ndiff_dip)
end if


if (check_list == 1) then
    write(*,*) ' -Calculating SCI transition dipoles (onthefly, no full D matrix)'

    !----------------------------------------------------------------------
    ! 8a: Build dipole-specific sparse list (n_diff <= 2) from the pair
    !     list already built for H in STEP 6
    !----------------------------------------------------------------------
    n_sparse_dip = 0
    do idx = 1, n_sparse
        if (sparse_ndiff(idx) <= 2) n_sparse_dip = n_sparse_dip + 1
    end do

    allocate(sparse_m_dip(n_sparse_dip), sparse_n_dip(n_sparse_dip), sparse_ndiff_dip(n_sparse_dip))
    n_sparse_dip = 0
    do idx = 1, n_sparse
        if (sparse_ndiff(idx) <= 2) then
            n_sparse_dip = n_sparse_dip + 1
            sparse_m_dip(n_sparse_dip) = sparse_m(idx)   ! LOCAL indices (1..n_sel)
            sparse_n_dip(n_sparse_dip) = sparse_n(idx)   ! LOCAL indices (1..n_sel)
            sparse_ndiff_dip(n_sparse_dip) = sparse_ndiff(idx)
        end if
    end do
    write(*,'(A,I0,A)') ' Dipole-specific pair list (SCI/LIST space) contains ', n_sparse_dip, ' elements.'

    !----------------------------------------------------------------------
    ! 8b: On-the-fly D*v0 accumulation
    !----------------------------------------------------------------------
    allocate(d_x(n_sel), d_y(n_sel), d_z(n_sel))
    d_x = 0.d0; d_y = 0.d0; d_z = 0.d0

    do idx = 1, n_sparse_dip
        i     = sparse_m_dip(idx)          
        j     = sparse_n_dip(idx)          
        n_diff = sparse_ndiff_dip(idx)
        m     = sel_list(i)                
        n_loc = sel_list(j)                
        vm(1:N_modes) = vec_combinations2(m,     1:N_modes)
        vn(1:N_modes) = vec_combinations2(n_loc, 1:N_modes)

        do ii = 1, N_modes
            ovlp_loc(ii) = modal_int(ii, vm(ii), vn(ii), 0)
        end do
        prefix_prod(0) = 1.0d0
        do ii = 1, N_modes
            prefix_prod(ii) = prefix_prod(ii-1) * ovlp_loc(ii)
        end do
        prod_all = prefix_prod(N_modes)
        suffix_prod_loc(N_modes+1) = 1.0d0
        do ii = N_modes, 1, -1
            suffix_prod_loc(ii) = suffix_prod_loc(ii+1) * ovlp_loc(ii)
        end do
        do ii = 1, N_modes
            prod_except_i(ii) = prefix_prod(ii-1) * suffix_prod_loc(ii+1)
        end do

        dip_local = 0.d0

        if (n_diff <= 1) then
            do ii = 1, N_modes
                if (abs(prod_except_i(ii)) < 1.d-30) cycle
                dip_local(1:3) = dip_local(1:3) + ( &
                    modal_int(ii, vm(ii), vn(ii), 1) * dipole_derivatives(ii, 1:3) &
                  + modal_int(ii, vm(ii), vn(ii), 2) * second_dipole_derivatives(ii,ii,1:3) &
                    ) * prod_except_i(ii)
            end do
        end if

        do kk = 1, N_modes*(N_modes+1)/2
            if (final_index_dipole(kk,1) == final_index_dipole(kk,2)) cycle
            if (abs(dipole_vec(kk,1)) < 1.d-30 .and. &
                abs(dipole_vec(kk,2)) < 1.d-30 .and. &
                abs(dipole_vec(kk,3)) < 1.d-30) cycle
            ii = final_index_dipole(kk,1)
            jj = final_index_dipole(kk,2)
            if (abs(ovlp_loc(ii)) < 1.d-300 .or. abs(ovlp_loc(jj)) < 1.d-300) then
                step_prod_loc = 1.0d0
                do pp = 1, N_modes
                    if (pp /= ii .and. pp /= jj) step_prod_loc = step_prod_loc * ovlp_loc(pp)
                end do
            else
                step_prod_loc = prod_all / (ovlp_loc(ii) * ovlp_loc(jj))
            end if
            if (abs(step_prod_loc) < 1.d-30) cycle
            step_prod_loc = step_prod_loc &
                * modal_int(ii, vm(ii), vn(ii), 1) * modal_int(jj, vm(jj), vn(jj), 1)
            dip_local(1:3) = dip_local(1:3) + dipole_vec(kk,1:3) * step_prod_loc * 2.0d0
        end do

        d_x(i) = d_x(i) + dip_local(1) * eigenvectors(j,1)
        d_y(i) = d_y(i) + dip_local(2) * eigenvectors(j,1)
        d_z(i) = d_z(i) + dip_local(3) * eigenvectors(j,1)
        if (i /= j) then
            d_x(j) = d_x(j) + dip_local(1) * eigenvectors(i,1)
            d_y(j) = d_y(j) + dip_local(2) * eigenvectors(i,1)
            d_z(j) = d_z(j) + dip_local(3) * eigenvectors(i,1)
        end if
    end do

    write(*,*) 'SCI dipole accumulation finished.'

    number_to_print_int = 0
    do i = 1, min(N_states, to_write_output) - 1
        if ((eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree < max_freq) &
            number_to_print_int = number_to_print_int + 1
    end do
    if (number_to_print_int < 1) number_to_print_int = 1
    if (number_to_print_int >= to_write_output) number_to_print_int = to_write_output - 1

    allocate(dipole_final(number_to_print_int, 3))
    allocate(intensity(number_to_print_int))

    do i = 1, number_to_print_int
        dipole_final(i,1) = ddot(n_sel, eigenvectors(1,i+1), 1, d_x, 1)
        dipole_final(i,2) = ddot(n_sel, eigenvectors(1,i+1), 1, d_y, 1)
        dipole_final(i,3) = ddot(n_sel, eigenvectors(1,i+1), 1, d_z, 1)
        intensity(i) = sum(dipole_final(i,:)**2)
    end do

    deallocate(d_x, d_y, d_z)
    deallocate(sparse_m_dip, sparse_n_dip, sparse_ndiff_dip)
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FINAL

write(200,*) 'VCI'
write(200,*) number_to_print_int
do i = 1, number_to_print_int
    write(200,'(2F18.6)') &
        (eigenvalues(i+1)-eigenvalues(1))/cm_to_hartree, &
        100.d0*intensity(i)/sqrt(sum(intensity(1:number_to_print_int)**2))
end do

!---------------------------------------------------------------------
! Transition dipoles + intensities (km/mol) — SCI analogue of the
! unit-300 output in vibrational_ci
!---------------------------------------------------------------------
open(unit=300, file='dipoles_intensity_sci.txt', status='replace')
write(300,'(A)') 'Transition dipoles and intensities (km/mol) - Selected VCI'
write(300,'(A)') 'State  Energy (cm-1)   Dip_x (a.u.)   Dip_y (a.u.)   Dip_z (a.u.)   Intensity (km/mol)'
do i = 1, number_to_print_int
    write(300,'(I6,F18.6,3F14.8,F14.4)') i+1, &
        (eigenvalues(i+1) - eigenvalues(1))/cm_to_hartree, &
        dipole_final(i,1), dipole_final(i,2), dipole_final(i,3), &
        intensity(i) * 974.88d0 * 2.d0
end do
close(300)
!---------------------------------------------------------------------

call cpu_time(end_time)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(count_rate_val)
write(*,'(A,F12.2,A)') ' SCI CPU time:     ', end_time - start_time, ' seconds'
write(*,'(A,F12.2,A)') ' SCI elapsed time: ', elapsed_time,          ' seconds'

!==========================================================================
! CLEANUP
!==========================================================================
! deallocate(modal_int)
! deallocate(dipole_vec, final_index_dipole)
! deallocate(is_ref, is_selected)
! deallocate(ref_list, sel_list)
! deallocate(eigenvalues, eigenvectors)
! deallocate(dipole_sel, dipole_final, intensity, tmp_vec)

end subroutine selected_vibrational_ci

!=============================================================================!
! vibrational_ci_sym                                                           !
! Symmetry-adapted VCI: blocks H by irrep, diagonalizes each block,           !
! merges and sorts all eigenvalues, prints with symmetry label column.         !
! Add this subroutine inside the vib_ci module, before end module vib_ci      !
!=============================================================================!
subroutine vibrational_ci_sym(Potential_3, Potential_4, N_modes, N_expansion, &
    HO_freq, store_integrals, full_coef, vec_combinations, &
    total_combinations, N_states, N_threads, dipole_derivatives, &
    second_dipole_derivatives, N_quanta, &
    total_3, total_4, &
    Potential_3_vec, Potential_4_vec, &
    final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
    final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
    cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    mode_irrep_in)
implicit none

!--- Input ---
integer,  intent(in) :: N_modes, N_expansion, N_quanta, N_states, N_threads
integer,  intent(in) :: total_combinations
integer,  intent(in) :: total_3, total_4, n_cubic_max, n_quartic_max
integer,  intent(in) :: vec_combinations(total_combinations, N_modes)
integer,  intent(in) :: final_index_3(total_3, 3), count_index_3(total_3, 3)
integer,  intent(in) :: n_unique_3(total_3), unique_modes_3(total_3, 3)
integer,  intent(in) :: check3(total_3)
integer,  intent(in) :: final_index_4(total_4, 4), count_index_4(total_4, 4)
integer,  intent(in) :: n_unique_4(total_4), unique_modes_4(total_4, 4)
integer,  intent(in) :: check4(total_4)
integer,  intent(in) :: cubic_for_mode(n_cubic_max, N_modes)
integer,  intent(in) :: n_cubic_for_mode(N_modes)
integer,  intent(in) :: quartic_for_mode(n_quartic_max, N_modes)
integer,  intent(in) :: n_quartic_for_mode(N_modes)
integer,  intent(in) :: mode_irrep_in(N_modes)
real*8,   intent(in) :: Potential_3(N_modes, N_modes, N_modes)
real*8,   intent(in) :: Potential_4(N_modes, N_modes, N_modes, N_modes)
real*8,   intent(in) :: HO_freq(N_modes)
real*8,   intent(in) :: full_coef(N_modes, N_expansion, N_expansion)
real*8,   intent(in) :: store_integrals(N_modes, N_expansion, N_expansion, 0:5)
real*8,   intent(in) :: dipole_derivatives(N_modes, 3)
real*8,   intent(in) :: second_dipole_derivatives(N_modes, N_modes, 3)
real*8,   intent(in) :: Potential_3_vec(total_3)
real*8,   intent(in) :: Potential_4_vec(total_4)

!--- Local scalars ---
integer :: i, j, k, ii, jj, kk, pp, m, n_loc
integer :: max_quanta, iirr, N_states_loc
integer :: mu, nu, p
integer :: n_diff, diff_modes(4)
integer :: n_sparse, idx, cc
integer :: n_total_states, global_idx
integer :: number_to_print_int
integer :: local_gs, local_i, n_blk
real*8  :: cm_to_hartree, zpe, H_val
logical :: term_valid
integer :: Vc_check, mode_idx
real*8  :: ovlp(N_modes), prefix_prod(0:N_modes), suffix_prod_loc(N_modes+1)
real*8  :: prod_except_i(N_modes), prod_all, step_prod
real*8  :: dip_local(3)
integer :: vm(N_modes), vn(N_modes)
real*8  :: values(3)
integer :: indices_3(3)

!--- Modal integrals ---
real*8, allocatable :: modal_int(:,:,:,:)

!--- Dipole pair data ---
real*8,  allocatable :: dipole_vec_2d(:,:)
integer, allocatable :: final_index_dipole(:,:)

!--- State symmetry ---
integer, allocatable :: state_irrep_vec(:)
integer, allocatable :: blk_size(:)
integer, allocatable :: blk_list(:,:)
integer, allocatable :: blk_pos(:)

!--- Sparse pair list ---
integer, allocatable :: sparse_m(:), sparse_n(:)
integer, allocatable :: sparse_ndiff(:)
integer, allocatable :: sparse_diff_modes(:,:)

!--- Per-block Hamiltonian/eigensystem ---
real*8, allocatable :: H_block(:,:)
real*8, allocatable :: evals_block(:), evecs_block(:,:)

!--- Storage across all blocks via derived type ---
type :: blk_evec_t
    real*8, allocatable :: ev(:,:)
    real*8, allocatable :: en(:)
end type blk_evec_t
type(blk_evec_t), allocatable :: blk_data(:)

!--- Merged sorted results ---
real*8,  allocatable :: all_en(:)
integer, allocatable :: all_irr(:)
integer, allocatable :: all_st(:,:)
real*8,  allocatable :: all_co(:,:)
integer, allocatable :: sort_idx(:)
integer, allocatable :: blk_local_idx(:)

!--- Dipole / intensity arrays ---
real*8, allocatable :: dip_mat(:,:,:)
real*8, allocatable :: dip_final(:,:)
real*8, allocatable :: intensity(:)
real*8, allocatable :: psi0(:), psiI(:), mu_psi0(:)
real*8, external    :: ddot

integer(8) :: t_start, t_end, count_rate_v, count_max_v
real(8)    :: elapsed_time, start_time, end_time

!==========================================================================
! INIT
!==========================================================================
cm_to_hartree = 0.0000045563350d0
max_quanta    = N_quanta
N_states_loc  = N_states

! Push mode_irrep_in into symmetry_module so state_irrep_func works
if (.not. allocated(mode_irrep)) allocate(mode_irrep(N_modes))
mode_irrep(1:N_modes) = mode_irrep_in(1:N_modes)

call system_clock(count_rate=count_rate_v, count_max=count_max_v)
call system_clock(t_start)
call cpu_time(start_time)

write(*,'(A)') '========================================'
write(*,'(A)') '   SYMMETRY-ADAPTED VCI (VCI-SYM)      '
write(*,'(A,A)') '   Point group: ', trim(point_group_name)
write(*,'(A,I4,A)') '   ', n_irreps, ' symmetry blocks'
write(*,'(A)') '========================================'
write(101,'(A)') '========================================'
write(101,'(A)') '   SYMMETRY-ADAPTED VCI (VCI-SYM)      '
write(101,'(A,A)') '   Point group: ', trim(point_group_name)
write(101,'(A,I4,A)') '   ', n_irreps, ' symmetry blocks'
write(101,'(A)') '========================================'

!==========================================================================
! STEP 1: Modal integrals
!==========================================================================
allocate(modal_int(N_modes, 0:max_quanta, 0:max_quanta, 0:5))
modal_int = 0.d0
call omp_set_num_threads(INT(N_threads, KIND=4))

!$OMP PARALLEL DO PRIVATE(ii,i,j,p,mu,nu) COLLAPSE(2) SCHEDULE(static)
do ii = 1, N_modes
    do i = 0, max_quanta
        do j = i, max_quanta
            do p = 0, 5
                do mu = 1, N_expansion
                    do nu = 1, N_expansion
                        modal_int(ii,i,j,p) = modal_int(ii,i,j,p) &
                            + full_coef(ii,i+1,mu)*full_coef(ii,j+1,nu) &
                            * store_integrals(ii,mu,nu,p)
                    end do
                end do
            end do
            do p = 0, 5
                modal_int(ii,j,i,p) = modal_int(ii,i,j,p)
            end do
        end do
    end do
end do
!$OMP END PARALLEL DO

write(*,'(A)') ' Modal integrals precomputed.'
write(101,'(A)') ' Modal integrals precomputed.'

!==========================================================================
! STEP 2: Dipole pair vector from second derivatives
!==========================================================================
allocate(dipole_vec_2d(N_modes*(N_modes+1)/2, 3))
allocate(final_index_dipole(N_modes*(N_modes+1)/2, 2))
dipole_vec_2d = 0.d0
kk = 0
do i = 1, N_modes
    do j = i, N_modes
        kk = kk + 1
        dipole_vec_2d(kk,1:3)    = second_dipole_derivatives(i,j,1:3)
        final_index_dipole(kk,1) = i
        final_index_dipole(kk,2) = j
    end do
end do

!==========================================================================
! STEP 3: Assign irrep to each VCI configuration
!==========================================================================
allocate(state_irrep_vec(total_combinations))
do m = 1, total_combinations
    state_irrep_vec(m) = state_irrep_func(vec_combinations(m,:), N_modes)
end do

! Count block sizes and report
allocate(blk_size(n_irreps))
blk_size = 0
do m = 1, total_combinations
    blk_size(state_irrep_vec(m)) = blk_size(state_irrep_vec(m)) + 1
end do

write(*,'(A)') ' Block sizes by irrep:'
write(101,'(A)') ' Block sizes by irrep:'
do iirr = 1, n_irreps
    write(*,'(4X,A8,A,I8)')  trim(irrep_names(iirr)), ': ', blk_size(iirr)
    write(101,'(4X,A8,A,I8)') trim(irrep_names(iirr)), ': ', blk_size(iirr)
end do

! Build blk_list(irrep, position_in_block) = global_config_index
! Build blk_pos(global_config_index)       = position in its block
allocate(blk_list(n_irreps, maxval(blk_size)+1))
allocate(blk_pos(total_combinations))
blk_list = 0
blk_pos  = 0
blk_size = 0   ! reuse as fill counter
do m = 1, total_combinations
    iirr = state_irrep_vec(m)
    blk_size(iirr) = blk_size(iirr) + 1
    blk_list(iirr, blk_size(iirr)) = m
    blk_pos(m) = blk_size(iirr)
end do

!==========================================================================
! STEP 4: Build full dipole matrix (cross-irrep pairs needed for intensities)
!==========================================================================
allocate(dip_mat(total_combinations, total_combinations, 3))
dip_mat = 0.d0

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP& SHARED(total_combinations, vec_combinations, N_modes, max_quanta, &
!$OMP&        modal_int, dipole_derivatives, second_dipole_derivatives, &
!$OMP&        dipole_vec_2d, final_index_dipole, dip_mat) &
!$OMP& PRIVATE(i, j, ii, jj, kk, pp, vm, vn, n_diff, &
!$OMP&         ovlp, prefix_prod, suffix_prod_loc, prod_except_i, prod_all, &
!$OMP&         step_prod, dip_local) &
!$OMP& SCHEDULE(dynamic, 64)
do i = 1, total_combinations
    do j = i, total_combinations
        vm(1:N_modes) = vec_combinations(i,1:N_modes)
        vn(1:N_modes) = vec_combinations(j,1:N_modes)

        ! Quick n_diff check (max 2 for dipole)
        n_diff = 0
        do ii = 1, N_modes
            if (vm(ii) /= vn(ii)) then
                n_diff = n_diff + 1
                if (n_diff > 2) exit
            end if
        end do
        if (n_diff > 2) cycle

        ! Overlaps and prefix/suffix products
        do ii = 1, N_modes
            ovlp(ii) = modal_int(ii, vm(ii), vn(ii), 0)
        end do
        prefix_prod(0) = 1.d0
        do ii = 1, N_modes
            prefix_prod(ii) = prefix_prod(ii-1) * ovlp(ii)
        end do
        prod_all = prefix_prod(N_modes)
        suffix_prod_loc(N_modes+1) = 1.d0
        do ii = N_modes, 1, -1
            suffix_prod_loc(ii) = suffix_prod_loc(ii+1) * ovlp(ii)
        end do
        do ii = 1, N_modes
            prod_except_i(ii) = prefix_prod(ii-1) * suffix_prod_loc(ii+1)
        end do

        dip_local = 0.d0

        ! One-mode dipole terms
        if (n_diff <= 1) then
            do ii = 1, N_modes
                if (abs(prod_except_i(ii)) < 1.d-30) cycle
                dip_local(1:3) = dip_local(1:3) + ( &
                    modal_int(ii,vm(ii),vn(ii),1) * dipole_derivatives(ii,1:3) &
                  + modal_int(ii,vm(ii),vn(ii),2) * second_dipole_derivatives(ii,ii,1:3) &
                    ) * prod_except_i(ii)
            end do
        end if

        ! Two-mode dipole terms (off-diagonal second derivatives)
        if (n_diff <= 2) then
            do kk = 1, N_modes*(N_modes+1)/2
                if (final_index_dipole(kk,1) == final_index_dipole(kk,2)) cycle
                if (abs(dipole_vec_2d(kk,1)) < 1.d-30 .and. &
                    abs(dipole_vec_2d(kk,2)) < 1.d-30 .and. &
                    abs(dipole_vec_2d(kk,3)) < 1.d-30) cycle
                ii = final_index_dipole(kk,1)
                jj = final_index_dipole(kk,2)
                if (abs(ovlp(ii)) < 1.d-300 .or. abs(ovlp(jj)) < 1.d-300) then
                    step_prod = 1.d0
                    do pp = 1, N_modes
                        if (pp /= ii .and. pp /= jj) &
                            step_prod = step_prod * ovlp(pp)
                    end do
                else
                    step_prod = prod_all / (ovlp(ii) * ovlp(jj))
                end if
                if (abs(step_prod) < 1.d-30) cycle
                step_prod = step_prod &
                    * modal_int(ii, vm(ii), vn(ii), 1) &
                    * modal_int(jj, vm(jj), vn(jj), 1)
                dip_local(1:3) = dip_local(1:3) + dipole_vec_2d(kk,1:3) * step_prod * 2.d0
            end do
        end if

        dip_mat(i,j,1:3) = dip_local(1:3)
        dip_mat(j,i,1:3) = dip_local(1:3)
    end do
end do
!$OMP END PARALLEL DO

!==========================================================================
! STEP 5: Allocate storage for all block eigensystems
!==========================================================================
allocate(blk_data(n_irreps))
do iirr = 1, n_irreps
    n_blk = blk_size(iirr)
    if (n_blk > 0) then
        allocate(blk_data(iirr)%ev(n_blk, n_blk))
        allocate(blk_data(iirr)%en(n_blk))
    else
        allocate(blk_data(iirr)%ev(0,0))
        allocate(blk_data(iirr)%en(0))
    end if
    blk_data(iirr)%ev = 0.d0
    blk_data(iirr)%en = 0.d0
end do

!==========================================================================
! STEP 6: For each irrep block build H and diagonalize
!==========================================================================
n_total_states = 0

do iirr = 1, n_irreps
    n_blk = blk_size(iirr)
    if (n_blk == 0) cycle

    write(*,'(A,A8,A,I8,A)') ' Building block for irrep ', &
        trim(irrep_names(iirr)), ', size = ', n_blk, ' ...'
    write(101,'(A,A8,A,I8,A)') ' Building block for irrep ', &
        trim(irrep_names(iirr)), ', size = ', n_blk, ' ...'

    !--- Build sparse pair list for this block (same-irrep only) ---
    n_sparse = 0
    do j = 1, n_blk
        do i = j, n_blk
            m     = blk_list(iirr, j)
            n_loc = blk_list(iirr, i)
            n_diff = 0
            do ii = 1, N_modes
                if (vec_combinations(m,ii) /= vec_combinations(n_loc,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff > 4) exit
                end if
            end do
            if (n_diff <= 4) n_sparse = n_sparse + 1
        end do
    end do

    allocate(sparse_m(n_sparse), sparse_n(n_sparse))
    allocate(sparse_ndiff(n_sparse))
    allocate(sparse_diff_modes(n_sparse, 4))
    sparse_diff_modes = 0
    idx = 0

    do j = 1, n_blk
        do i = j, n_blk
            m     = blk_list(iirr, j)
            n_loc = blk_list(iirr, i)
            n_diff = 0
            diff_modes = 0
            do ii = 1, N_modes
                if (vec_combinations(m,ii) /= vec_combinations(n_loc,ii)) then
                    n_diff = n_diff + 1
                    if (n_diff <= 4) diff_modes(n_diff) = ii
                    if (n_diff >  4) exit
                end if
            end do
            if (n_diff <= 4) then
                idx = idx + 1
                sparse_m(idx)          = j
                sparse_n(idx)          = i
                sparse_ndiff(idx)      = n_diff
                sparse_diff_modes(idx,1:4) = diff_modes(1:4)
            end if
        end do
    end do

    !--- Allocate and fill H_block ---
    allocate(H_block(n_blk, n_blk))
    H_block = 0.d0

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& SHARED(n_sparse, sparse_m, sparse_n, &
    !$OMP&        n_blk, blk_list, iirr, vec_combinations, N_modes, max_quanta, &
    !$OMP&        modal_int, Potential_3, Potential_4, HO_freq, &
    !$OMP&        Potential_3_vec, Potential_4_vec, &
    !$OMP&        check3, check4, total_3, total_4, &
    !$OMP&        final_index_3, count_index_3, n_unique_3, unique_modes_3, &
    !$OMP&        final_index_4, count_index_4, n_unique_4, unique_modes_4, &
    !$OMP&        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
    !$OMP&        quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
    !$OMP&        H_block) &
    !$OMP& PRIVATE(k, i, j, m, n_loc, ii, vm, vn, H_val) &
    !$OMP& SCHEDULE(dynamic, 64)
    do k = 1, n_sparse
        i     = sparse_m(k)
        j     = sparse_n(k)
        m     = blk_list(iirr, i)
        n_loc = blk_list(iirr, j)
        vm(1:N_modes) = vec_combinations(m,     1:N_modes)
        vn(1:N_modes) = vec_combinations(n_loc, 1:N_modes)

        call compute_H_element(m, n_loc, vm, vn, N_modes, max_quanta, &
            modal_int, Potential_3, Potential_4, HO_freq, &
            Potential_3_vec, Potential_4_vec, &
            check3, check4, total_3, total_4, &
            final_index_3, count_index_3, n_unique_3, unique_modes_3, &
            final_index_4, count_index_4, n_unique_4, unique_modes_4, &
            cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
            quartic_for_mode, n_quartic_for_mode, n_quartic_max, &
            H_val)

        H_block(i,j) = H_val
        H_block(j,i) = H_val
    end do
    !$OMP END PARALLEL DO

    deallocate(sparse_m, sparse_n, sparse_ndiff, sparse_diff_modes)

    !--- Diagonalize block ---
    allocate(evals_block(n_blk), evecs_block(n_blk, n_blk))
    evals_block = 0.d0; evecs_block = 0.d0
    call dsyevd_A(H_block, evals_block, evecs_block)
    deallocate(H_block)

    blk_data(iirr)%en(1:n_blk) = evals_block(1:n_blk)
    blk_data(iirr)%ev(1:n_blk, 1:n_blk) = evecs_block(1:n_blk, 1:n_blk)
    deallocate(evals_block, evecs_block)

    n_total_states = n_total_states + n_blk

    write(*,'(A,A8,A)') ' Block ', trim(irrep_names(iirr)), ' diagonalized.'
    write(101,'(A,A8,A)') ' Block ', trim(irrep_names(iirr)), ' diagonalized.'
end do

!==========================================================================
! STEP 7: Merge all eigenvalues and sort by energy
!==========================================================================
allocate(all_en(n_total_states))
allocate(all_irr(n_total_states))
allocate(all_st(n_total_states, 3))
allocate(all_co(n_total_states, 3))
allocate(blk_local_idx(n_total_states))

global_idx = 0
do iirr = 1, n_irreps
    n_blk = blk_size(iirr)
    if (n_blk == 0) cycle
    do i = 1, n_blk
        global_idx = global_idx + 1
        all_en(global_idx)  = blk_data(iirr)%en(i)
        all_irr(global_idx) = iirr
        blk_local_idx(global_idx) = i
        ! Top-3 coefficients, mapping block indices to global config indices
        call top_three_components(blk_data(iirr)%ev(:,i), n_blk, values, indices_3)
        all_co(global_idx, 1:3) = values(1:3)
        do k = 1, 3
            if (indices_3(k) >= 1 .and. indices_3(k) <= n_blk) then
                all_st(global_idx, k) = blk_list(iirr, indices_3(k))
            else
                all_st(global_idx, k) = 0
            end if
        end do
    end do
end do

! Sort by energy (insertion sort; n_total_states can be large but
! sorting once is cheap compared to H build)
allocate(sort_idx(n_total_states))
do i = 1, n_total_states
    sort_idx(i) = i
end do
do i = 1, n_total_states - 1
    do j = i + 1, n_total_states
        if (all_en(sort_idx(j)) < all_en(sort_idx(i))) then
            k = sort_idx(i)
            sort_idx(i) = sort_idx(j)
            sort_idx(j) = k
        end if
    end do
end do

zpe = all_en(sort_idx(1))

!==========================================================================
! STEP 8: Print  — same format as vibrational_ci + extra Sym column
!==========================================================================
write(*,'(A)') ' '
write(*,'(A)') '-------------------------------------------------------------------------------------------------------------------------------------------------'
write(*,'(A19,A22, A19, 6A12, 1A12)') 'E (cm-1)', 'E - ZPE', '  ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C', 'Sym Block'
write(*,'(A)') '-------------------------------------------------------------------------------------------------------------------------------------------------'
write(101,'(A)') ' '
write(101,'(A)') '-------------------------------------------------------------------------------------------------------------------------------------------------'
write(101,'(A19,A22, A19, 6A12, 1A12)') 'E (cm-1)', 'E - ZPE', '  ', 'Coeff A', 'Coeff B', 'Coeff C', 'State A', 'State B', 'State C', 'Sym Block'
write(101,'(A)') '-------------------------------------------------------------------------------------------------------------------------------------------------'

do i = 1, min(N_states_loc, n_total_states)
    idx = sort_idx(i)
    if (i <= 10) then
        write(*,'(F19.4, F22.4, A19, 3F12.6, 3I12, 1A12)') all_en(idx)/cm_to_hartree, (all_en(idx) - zpe)/cm_to_hartree, '      ||      ', all_co(idx,1), all_co(idx,2), all_co(idx,3), all_st(idx,1), all_st(idx,2), all_st(idx,3), trim(irrep_names(all_irr(idx)))
    end if
    write(101,'(F19.4, F22.4, A19, 3F12.6, 3I12, 1A12)')   all_en(idx)/cm_to_hartree, (all_en(idx) - zpe)/cm_to_hartree, '      ||      ', all_co(idx,1), all_co(idx,2), all_co(idx,3), all_st(idx,1), all_st(idx,2), all_st(idx,3), trim(irrep_names(all_irr(idx)))
end do

write(*,'(A)') '-------------------------------------------------------------------------------------------------------------------------------------------------'
write(101,'(A)') '-------------------------------------------------------------------------------------------------------------------------------------------------'

!==========================================================================
! STEP 9: Transition dipoles and intensities
! Normalization is global (not per block) as required.
! Ground state = sort_idx(1).  All transitions from GS to excited states.
!==========================================================================
write(*,*) ' -Calculating symmetry-adapted VCI transition dipoles'
write(101,*) ' -Calculating symmetry-adapted VCI transition dipoles'

! Count transitions below 4500 cm-1 from ZPE
number_to_print_int = 0
do i = 2, min(N_states_loc, n_total_states)
    if ((all_en(sort_idx(i)) - zpe)/cm_to_hartree < 4500.d0) &
        number_to_print_int = number_to_print_int + 1
end do
if (number_to_print_int < 1) number_to_print_int = 1

allocate(dip_final(number_to_print_int, 3))
allocate(intensity(number_to_print_int))
allocate(psi0(total_combinations))
allocate(psiI(total_combinations))
allocate(mu_psi0(total_combinations))
dip_final = 0.d0
intensity  = 0.d0

! Build ground-state wavefunction in full-space representation
! Ground state: sort_idx(1), irrep all_irr(sort_idx(1)), local index blk_local_idx(sort_idx(1))
psi0 = 0.d0
iirr     = all_irr(sort_idx(1))
local_gs = blk_local_idx(sort_idx(1))
n_blk    = blk_size(iirr)
do m = 1, n_blk
    psi0(blk_list(iirr, m)) = blk_data(iirr)%ev(m, local_gs)
end do

! For each Cartesian component compute mu|psi0> then project onto each excited state
do ii = 1, 3
    ! mu_psi0 = dip_mat(:,:,ii) * psi0
    call dgemv('N', total_combinations, total_combinations, &
               1.d0, dip_mat(1,1,ii), total_combinations, &
               psi0(1), 1, 0.d0, mu_psi0(1), 1)

    do i = 1, number_to_print_int
        idx      = sort_idx(i + 1)
        iirr     = all_irr(idx)
        local_i  = blk_local_idx(idx)
        n_blk    = blk_size(iirr)

        ! Build excited-state wavefunction in full-space representation
        psiI = 0.d0
        do m = 1, n_blk
            psiI(blk_list(iirr, m)) = blk_data(iirr)%ev(m, local_i)
        end do

        dip_final(i, ii) = ddot(total_combinations, psiI(1), 1, mu_psi0(1), 1)
    end do
end do

! Compute intensities with GLOBAL normalization
do i = 1, number_to_print_int
    intensity(i) = sum(dip_final(i,:)**2)
end do

! Write intensities file (same format as vibrational_ci, with extra sym column)
write(200,*) 'VCI'
write(200,*) number_to_print_int
do i = 1, number_to_print_int
    idx = sort_idx(i + 1)
    write(200,'(2F18.6,2X)') &
        (all_en(idx) - zpe)/cm_to_hartree, &
        100.d0 * intensity(i) / sqrt(sum(intensity(1:number_to_print_int)**2))!, &
       ! trim(irrep_names(all_irr(idx)))
end do

!==========================================================================
! STEP 10: Timing and cleanup
!==========================================================================
call cpu_time(end_time)
call system_clock(t_end)
elapsed_time = real(t_end - t_start) / real(count_rate_v)
write(*,'(A,F12.2,A)') ' VCI-SYM CPU time:     ', end_time - start_time, ' seconds'
write(*,'(A,F12.2,A)') ' VCI-SYM elapsed time: ', elapsed_time,          ' seconds'
write(101,'(A,F12.2,A)') ' VCI-SYM CPU time:     ', end_time - start_time, ' seconds'
write(101,'(A,F12.2,A)') ' VCI-SYM elapsed time: ', elapsed_time,          ' seconds'

! Deallocate all local allocatables
deallocate(modal_int)
deallocate(dipole_vec_2d, final_index_dipole)
deallocate(state_irrep_vec)
deallocate(blk_size, blk_list, blk_pos)
deallocate(dip_mat)
deallocate(all_en, all_irr, all_st, all_co, sort_idx, blk_local_idx)
deallocate(dip_final, intensity, psi0, psiI, mu_psi0)
do iirr = 1, n_irreps
    if (allocated(blk_data(iirr)%ev)) deallocate(blk_data(iirr)%ev)
    if (allocated(blk_data(iirr)%en)) deallocate(blk_data(iirr)%en)
end do
deallocate(blk_data)

end subroutine vibrational_ci_sym


subroutine top_three_components(vec, n, values, indices)
    ! This subroutine finds the three largest magnitude components of a vector
    ! and returns their values and indices.
    !
    ! Inputs:
    !   vec     : input vector of length n
    !   n       : actual size of the vector (vec may be larger than n)
    ! Outputs:
    !   values  : the three largest magnitude values (sorted descending)
    !   indices : the corresponding indices in the original vector
    implicit none
    real(8), intent(in)  :: vec(:)
    integer, intent(in)  :: n
    real(8), intent(out) :: values(3)
    integer, intent(out) :: indices(3)
    
    integer :: i, j, temp_idx
    real(8) :: temp_val
    real(8), allocatable :: vec_copy(:)
    integer, allocatable :: idx_copy(:)
    
    if (n < 3) then
        print *, "Warning: Vector size less than 3, only returning", n, "components"
        values = 0.0d0
        indices = 0
        do i = 1, n
            values(i) = vec(i)
            indices(i) = i
        end do
        return
    end if
    
    allocate(vec_copy(n))
    allocate(idx_copy(n))
    
    vec_copy = vec
    do i = 1, n
        idx_copy(i) = i
    end do
    
    do i = 1, 3
        do j = i+1, n
            if (abs(vec_copy(j)) > abs(vec_copy(i))) then
                temp_val = vec_copy(i)
                vec_copy(i) = vec_copy(j)
                vec_copy(j) = temp_val
                temp_idx = idx_copy(i)
                idx_copy(i) = idx_copy(j)
                idx_copy(j) = temp_idx
            end if
        end do
    end do
    
    do i = 1, 3
        values(i) = vec_copy(i)
        indices(i) = idx_copy(i)
    end do
    
    deallocate(vec_copy, idx_copy)
    
end subroutine top_three_components

end module vib_ci