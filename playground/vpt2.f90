module v_pt2
use combination  ! Provides count_unique_elements subroutine for processing unique mode indices
use jacobi_diagonalization  ! Provides diagonalization routines (dsyevr_A, dsyevd_A)
use omp_lib  ! Provides OpenMP parallelization functions
use symmetry_module
use read_input_file
contains

subroutine compute_H_element_vpt2(m_cfg, n_cfg, vm, vn, N_modes, max_quanta, &
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
           ! modal_int(ii, vm(ii), vn(ii), 5) & 
            + modal_int(ii, vm(ii), vn(ii), 4) * Potential_4(ii,ii,ii,ii) &
            + modal_int(ii, vm(ii), vn(ii), 3) * Potential_3(ii,ii,ii) &
           ! + modal_int(ii, vm(ii), vn(ii), 2) * HO_freq(ii) * 0.5d0 &
            ) * prod_except_i(ii)
    end do
end if
! The lines for modal_int(...,5) (kinetic energy) and modal_int(...,2)
! (harmonic potential, i.e., q^2) are commented out to keep consistency
! with the definition of the perturbation V used in VPT2 (see
! https://manual.q-chem.com/5.2/Ch11.S11.SS2.html). In that framework,
! V contains only the anharmonic terms (cubic, quartic, ...).
!
! The original implementation followed the Crystal paper, where the
! coupling term Vc was defined as the part of the potential where at
! least two mode indices differ (i.e., purely off-diagonal in the
! harmonic basis).
!
! Option A (current, commented lines):
!   Compute H_val = <m|V_anharm|n> only. In the VPT2 assembly, the
!   zeroth-order harmonic energy U_energies is added separately:
!       E_VPT2 = U_energies + <n|V_anharm|n> + sum_{m≠n} |<m|V_anharm|n>|^2 / (E_n^0 - E_m^0)
!   This is the standard VPT2 expression where V = anharmonic part.
!
! Option B (uncomment the two lines and remove U_energies addition):
!   H_val would then include the full harmonic Hamiltonian H0 (kinetic + harmonic potential)
!   as well as the anharmonic V. The harmonic part is diagonal in the HO basis:
!       <n|H0|n> = E_n^(0) = ω(n + 1/2)
!       <m|H0|n> = 0 for m ≠ n
!   Therefore the diagonal element becomes E_n^(0) + <n|V_anharm|n>,
!   identical to Option A. Off-diagonal elements remain purely anharmonic.
!
! The reason the harmonic Hamiltonian is diagonal is that, although q^2 and p^2
! individually have off-diagonal elements (connecting states differing by 2 quanta),
! when added together with the correct coefficients (1/2 p^2 + 1/2 ω^2 q^2), those
! off-diagonal contributions cancel exactly. Only the diagonal part survives.
! This cancellation is special to the harmonic oscillator and is why we can
! either include H0 explicitly or add its known eigenvalues separately.
!
! Thus both options produce the same VPT2 energies; the choice is a matter of: I want it like this :D
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

end subroutine compute_H_element_vpt2

! ============================================================================
! DIFFERENCE BETWEEN THIS IMPLEMENTATION AND STANDARD VPT2 (ORCA)
!
! Reference: Franke, Stanton, Douberly, J. Phys. Chem. A 2021, 125, 1301-1324
!            https://dx.doi.org/10.1021/acs.jpca.0c09526
!
! This code uses the full anharmonic potential as the perturbation:
!     V = V_cubic + V_quartic
! and applies ordinary second-order Rayleigh-Schrödinger perturbation theory.
! The energy of a state |n> is then:
!
!   E_n = E_n^(0) + <n|V|n> + sum_{m != n} |<m|V|n>|^2 / (E_n^(0) - E_m^(0))
!                                                                    (Eq. A)
!
! In this expression, the off-diagonal sum includes both cubic and quartic
! matrix elements. The diagonal term <n|V|n> includes the diagonal quartic
! (the cubic diagonal is zero for harmonic oscillator basis states).
!
! Standard VPT2, as used in ORCA and described in the above paper, introduces
! a formal perturbation parameter lambda and orders the potential as:
!     V = lambda * V_cubic + lambda^2 * V_quartic
! The energy is then truncated at total order lambda^2. This yields:
!
!   E_n = E_n^(0) + <n|V_quartic|n> 
!        + sum_{m != n} |<m|V_cubic|n>|^2 / (E_n^(0) - E_m^(0))   (Eq. B)
!
! The second term is the first-order correction from the quartic potential,
! while the sum is the second-order correction from the cubic potential.
! Off-diagonal quartic terms would contribute at order lambda^4 because the
! quartic operator is already lambda^2 and squaring it in the second-order
! sum gives lambda^4. Therefore they are omitted in standard VPT2.
!
! The present implementation (Eq. A) includes these off-diagonal quartic
! contributions, which is a legitimate higher-order perturbation treatment
! but is not the same as standard VPT2. This will produce different energies
! when compared to ORCA, especially for modes with sizable quartic coupling.
!
! =========================================================================
subroutine vpt2(Potential_3, Potential_4, N_modes, N_expansion, &
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
  real*8 :: H_val, denom
  integer :: vm(N_modes), vn(N_modes)
  real*8, external :: ddot

  ! Timing
  integer(8) :: t_start, t_end, count_rate, count_max
  real(8) :: elapsed_time, start_time, end_time

  integer(8) :: t_start_dav, t_end_dav, count_rate_dav, count_max_dav
  real(8) :: elapsed_time_dav, start_time_dav, end_time_dav

  real*8, allocatable :: U_energies(:), VPT2_energies(:) !zero-order (harmonic) energies and VPT2-corrected energies

  cm_to_hartree = 0.0000045563350d0
  call system_clock(count_rate=count_rate, count_max=count_max)
  call system_clock(t_start)
  call cpu_time(start_time)

  max_quanta_actual = N_quanta

  call classify_quartic_terms(total_4, final_index_4, n_unique_4, &
      check4, Potential_4_vec, is_semiquartic_ff)
  to_sparse_cut = merge(3, 4, is_semiquartic_ff)
  write(*,'(A,I2)') ' Sparse pair-list cutoff set to n_diff <= ', to_sparse_cut

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 0: Zero-order (harmonic) energies  U(i) = sum_k freq_k * (n_k + 1/2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(U_energies(total_combinations))
  write(*,'(A)') ' >>> Calculating zero-order (harmonic) energies of HO states'
  do j = 1, total_combinations
     sum_energy_HO = 0.d0
     do k = 1, N_modes
        sum_energy_HO = sum_energy_HO + HO_freq(k) * (real(vec_combinations(j, k), 8) + 0.5d0)
     end do
     U_energies(j) = sum_energy_HO
  end do

 ! write(*,'(A)') ' First 10 zero-order energies (cm-1):'
 ! do j = 1, min(10, total_combinations)
 !    write(*,*) U_energies(j) / cm_to_hartree
 ! end do

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
  ! STEP 4:  Precompute H_sparse  --  H_sparse(idx) = <vm|H_hat|vn>
  !          Since vec_combinations is a harmonic-oscillator product basis,
  !          H0 is exactly diagonal in this basis, so:
  !            m == n  ->  H_sparse = U(m) + <m|V|m>   (full diagonal energy)
  !            m /= n  ->  H_sparse = <m|V|n>           (pure anharmonic coupling)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(H_sparse(n_sparse))
        H_sparse = 0.d0
        write(*,'(A)') ' Precomputing sparse Hamiltonian matrix elements...'
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
           call compute_H_element_vpt2(m, n, vm, vn, N_modes, max_quanta_actual, &
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

  call cpu_time(end_time)
  call system_clock(t_end)
  elapsed_time = real(t_end - t_start) / real(count_rate)
  write(*,'(1A, 1F12.2, 1A)') " H build CPU time:     ", end_time - start_time, " seconds"
  write(*,'(1A, 1F12.2, 1A)') " H build elapsed time: ", elapsed_time, " seconds"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STEP 5: Assemble VPT2 energies
  !   E(i) = U(i) + <i|V|i> + sum_{j/=i} |<i|V|j>|^2 / (U(i) - U(j))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(VPT2_energies(total_combinations))
  VPT2_energies = 0.d0

  ! --- Diagonal contribution: H_sparse(i,i) already IS U(i) + <i|V|i> ---
  do i = 1, n_sparse
     if (sparse_m(i) == sparse_n(i)) then
        VPT2_energies(sparse_m(i)) = U_energies(sparse_m(i)) + H_sparse(i)
     end if
  end do
  H_threshold = 2.0*cm_to_hartree
  ! --- Second-order contribution from off-diagonal couplings ---   
  write(*,'(A)') ' Accumulating VPT2 second-order corrections ...'
  !$OMP PARALLEL DO PRIVATE(idx, m, n, H_val, denom) SCHEDULE(dynamic,256) DEFAULT(NONE) &
  !$OMP& SHARED(n_sparse, sparse_m, sparse_n, H_sparse, U_energies, VPT2_energies, H_threshold)
  do idx = 1, n_sparse
     m = sparse_m(idx)
     n = sparse_n(idx)
     if (m == n) cycle
     H_val = H_sparse(idx)
     denom = U_energies(m) - U_energies(n)
     if (abs(denom) < H_threshold) then
        ! Near-degenerate denominator: candidate Fermi/Darling-Dennison resonance.
        ! Flagged and skipped here rather than blown up; route these pairs to a
        ! variational (deperturbative) treatment instead of plain sum-over-states PT2.
        !$OMP CRITICAL
       ! write(*,'(A,I8,A,I8,A,ES12.4,A,ES12.4)') ' WARNING: near-degenerate pair (', m, ',', n, &
       !      ') |dE| in a.u. = ', abs(denom), '  Hij = ', H_val
        !$OMP END CRITICAL
        cycle
     end if
     !$OMP ATOMIC UPDATE
     VPT2_energies(m) = VPT2_energies(m) + H_val*H_val / denom
     !$OMP ATOMIC UPDATE
     VPT2_energies(n) = VPT2_energies(n) - H_val*H_val / denom
  end do
  !$OMP END PARALLEL DO

  VPT2_energies = VPT2_energies / cm_to_hartree
  U_energies    = U_energies    / cm_to_hartree
  

  !Here, the states are given in the order in which the vec_combinations is built, so itstart with 0, 0, 0... then 1 quanta per mode, so it is easier to print only the fundamentals.

  write(*,'(9A20)') ' States:',  'HO energy',          'VPT2 energy',         'VPT2 - HO',  'HO freq', 'VPT2 freq', 'VPT2 - HO', 'all in [cm-1]'
  do i = 1, N_modes + 1
    if(i == 1) then 
        write(*,'(1I20, 3F20.4, 1A20, 3F20.4)') i, U_energies(i), VPT2_energies(i), VPT2_energies(i) - U_energies(i), '0.0000', VPT2_energies(i) - VPT2_energies(1), VPT2_energies(i) - VPT2_energies(1) 
    else 
        write(*,'(1I20, 7F20.4)') i, U_energies(i), VPT2_energies(i), VPT2_energies(i) - U_energies(i), HO_freq(i-1)/cm_to_hartree, VPT2_energies(i) - VPT2_energies(1), VPT2_energies(i) - VPT2_energies(1) - HO_freq(i-1)/cm_to_hartree
    end if
    
  end do


  write(101,'(9A20)') ' States:',  'HO energy',          'VPT2 energy',         'VPT2 - HO',  'HO freq', 'VPT2 freq', 'VPT2 - HO', 'all in [cm-1]'
  do i = 1, total_combinations
    if(i == 1) then 
        write(101,'(1I20, 3F20.4, 1A20, 3F20.4)') i, U_energies(i), VPT2_energies(i), VPT2_energies(i) - U_energies(i), '0.0000', VPT2_energies(i) - VPT2_energies(1), VPT2_energies(i) - VPT2_energies(1) 
    else 
        write(101,'(1I20, 7F20.4)') i, U_energies(i), VPT2_energies(i), VPT2_energies(i) - U_energies(i), HO_freq(i-1)/cm_to_hartree, VPT2_energies(i) - VPT2_energies(1), VPT2_energies(i) - VPT2_energies(1) - HO_freq(i-1)/cm_to_hartree
    end if
    
  end do

  write(*,*) '-DONE VPT2'
  write(*,*)
end subroutine

end module