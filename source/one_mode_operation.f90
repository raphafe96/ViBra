module one_mode
use jacobi_diagonalization
use compute_integrals
use combination

contains

subroutine constant_one_mode(coefficients, Potential_3, Potential_4, N_modes, N_expansion, &
    HO_freq, new_coeff, store_integrals, total_energy, write_on_out, mode_excite, full_coef, &
    calc_intensity, dipole_derivatives, intensity, N_threads, &
    total_3, total_4, &
    Potential_3_vec, Potential_4_vec, &
    final_index_3, count_index_3, check3, &
    final_index_4, count_index_4, check4)
implicit none

integer,  intent(in) :: N_modes, N_expansion, write_on_out, calc_intensity, N_threads
integer,  intent(in) :: total_3, total_4
integer,  intent(in) :: mode_excite(N_expansion)
integer,  intent(in) :: final_index_3(total_3, 3), count_index_3(total_3, 3)
integer,  intent(in) :: check3(total_3)
integer,  intent(in) :: final_index_4(total_4, 4), count_index_4(total_4, 4)
integer,  intent(in) :: check4(total_4)
real*8,   intent(in) :: coefficients(N_modes, N_expansion)
real*8,   intent(in) :: Potential_3(N_modes, N_modes, N_modes)
real*8,   intent(in) :: Potential_4(N_modes, N_modes, N_modes, N_modes)
real*8,   intent(in) :: HO_freq(N_modes)
real*8,   intent(in) :: store_integrals(N_modes, N_expansion, N_expansion, 0:5)
real*8,   intent(in) :: Potential_3_vec(total_3)
real*8,   intent(in) :: Potential_4_vec(total_4)
real*8,   intent(in) :: dipole_derivatives(N_modes, 3)
real*8,   intent(inout) :: full_coef(N_modes, N_expansion, N_expansion)
real*8,   intent(out) :: new_coeff(N_modes, N_expansion)
real*8,   intent(out) :: total_energy
real*8,   intent(out) :: intensity

integer :: i, j, k, p, q, nu, mu, working_mode, for_output(N_expansion)
integer :: ii, jj, pp
logical :: is_vc_term
real*8  :: mid_integral, X_one_mode(N_modes, 0:4)
real*8  :: step(4), full_step(4)
real*8  :: Hamiltonian(N_modes, N_expansion, N_expansion)
real*8  :: Vc, u, full_term_val
real*8  :: dipole(N_modes, 3), step_multi_vec2(N_modes)
real*8  :: cm_to_hartree

! Overlap cache
real*8  :: overlap_cache(N_modes, 0:4)

! Vc accumulators
real*8  :: Vc_cubic, Vc_quartic

! For diagonalization
REAL*8  :: H(N_expansion, N_expansion)
REAL*8  :: eigenvalues(N_expansion)
REAL*8  :: eigenvectors(N_expansion, N_expansion)

cm_to_hartree = 0.0000045563350d0
Vc = 0.d0
X_one_mode = 0.d0

! Precompute overlap cache
overlap_cache = 0.d0
do k = 1, N_modes
    do p = 0, 4
        mid_integral = 0.d0
        do mu = 1, N_expansion
            do nu = 1, N_expansion
                mid_integral = mid_integral &
                    + coefficients(k, mu) * coefficients(k, nu) &
                    * store_integrals(k, mu, nu, p)
            end do
        end do
        overlap_cache(k, p) = mid_integral
    end do
end do

! No mixing contributions
do working_mode = 1, N_modes
    X_one_mode(working_mode, 2) = X_one_mode(working_mode, 2) + (HO_freq(working_mode)**1)/2.
end do    

! Cubic: X_one_mode and Vc in single pass
Vc_cubic = 0.d0
do ii = 1, total_3
    if (check3(ii) == 0) cycle

    ! Full contraction for Vc
    full_step = 1.0d0
    is_vc_term = .true.
    do pp = 1, 3
        if (final_index_3(ii,pp) > 0 .and. count_index_3(ii,pp) > 0) then
            full_step(pp) = overlap_cache(final_index_3(ii,pp), count_index_3(ii,pp))
            if (count_index_3(ii,pp) == 3) is_vc_term = .false.
        end if
    end do
    full_term_val = Potential_3_vec(ii) * full_step(1) * full_step(2) * full_step(3) * full_step(4)

    if (is_vc_term) Vc_cubic = Vc_cubic + full_term_val

    ! Per working_mode contributions
    do working_mode = 1, N_modes
        q = 0
        step = 1.0d0
        do pp = 1, 3
            if (count_index_3(ii,pp) > 0 .and. final_index_3(ii,pp) > 0 &
                .and. final_index_3(ii,pp) /= working_mode) then
                step(pp) = overlap_cache(final_index_3(ii,pp), count_index_3(ii,pp))
                if (count_index_3(ii,pp) == 3) step = 0.d0
            end if

            if (count_index_3(ii,pp) > 0 .and. final_index_3(ii,pp) == working_mode) then
                step(pp) = overlap_cache(final_index_3(ii,pp), 0)
                q = count_index_3(ii,pp)
            end if
        end do
        X_one_mode(working_mode, q) = X_one_mode(working_mode, q) &
            + Potential_3_vec(ii) * step(1) * step(2) * step(3) * step(4)
    end do
end do
Vc = Vc + Vc_cubic

! Quartic: X_one_mode and Vc in single pass
Vc_quartic = 0.d0
do ii = 1, total_4
    if (check4(ii) == 0) cycle

    ! Full contraction for Vc
    full_step = 1.0d0
    is_vc_term = .true.
    do pp = 1, 4
        if (final_index_4(ii,pp) > 0 .and. count_index_4(ii,pp) > 0) then
            full_step(pp) = overlap_cache(final_index_4(ii,pp), count_index_4(ii,pp))
            if (count_index_4(ii,pp) == 4) is_vc_term = .false.
        end if
    end do
    full_term_val = Potential_4_vec(ii) * full_step(1) * full_step(2) * full_step(3) * full_step(4)

    if (is_vc_term) Vc_quartic = Vc_quartic + full_term_val

    ! Per working_mode contributions
    do working_mode = 1, N_modes
        q = 0
        step = 1.0d0
        do pp = 1, 4
            if (count_index_4(ii,pp) > 0 .and. final_index_4(ii,pp) > 0 &
                .and. final_index_4(ii,pp) /= working_mode) then
                step(pp) = overlap_cache(final_index_4(ii,pp), count_index_4(ii,pp))
                if (count_index_4(ii,pp) == 4) step = 0.d0
            end if

            if (count_index_4(ii,pp) > 0 .and. final_index_4(ii,pp) == working_mode) then
                step(pp) = overlap_cache(final_index_4(ii,pp), 0)
                q = count_index_4(ii,pp)
            end if
        end do
        X_one_mode(working_mode, q) = X_one_mode(working_mode, q) &
            + Potential_4_vec(ii) * step(1) * step(2) * step(3) * step(4)
    end do
end do
Vc = Vc + Vc_quartic

! Build Hamiltonian
Hamiltonian = 0.d0
do k = 1, N_modes
    do i = 1, N_expansion
        do j = 1, N_expansion
            do p = 0, 4
                u = sqrt((HO_freq(k)**real(p)))
                Hamiltonian(k, i, j) = Hamiltonian(k, i, j) + X_one_mode(k, p)*store_integrals(k, i, j, p)
            end do
            u = HO_freq(k)
            Hamiltonian(k, i, j) = Hamiltonian(k, i, j) + store_integrals(k, i, j, 5)
        end do
    end do
end do

! Diagonalize
total_energy = 0.d0
new_coeff = 0.d0

do i = 1, N_expansion
    for_output(i) = i-1
end do

do k = 1, N_modes
    eigenvalues  = 0.d0
    eigenvectors = 0.d0
    H = 0.d0
    H(:,:) = Hamiltonian(k, :, :)

    call dsyevd_A(H, eigenvalues, eigenvectors)

    do i = 1, N_expansion
        H(i,:) = eigenvectors(i, :)
    end do

    if (write_on_out == 1) then
        write(101,'(1A30, 1I10)') 'Coefficient matrix for mode: ', k
        write(101,*) '---'
        write(101, '(100I10)') for_output(:)
        do i = 1, N_expansion
            write(101,'(100F10.5)') H(i, :)
            full_coef(k,i,:) = H(:, i)
        end do
        write(101,*)
    end if

    new_coeff(k, :) = H(:, mode_excite(k)+1)
    total_energy = total_energy + eigenvalues(mode_excite(k)+1)
end do

total_energy = (total_energy - Vc*real(N_modes-1))/cm_to_hartree

if (calc_intensity == 1) then
    dipole = 0.d0
    step_multi_vec2 = 1.d0

    do ii = 1, N_modes
        mid_integral = 0.d0
        do mu = 1, N_expansion
            mid_integral = mid_integral &
                + full_coef(ii, 1, mu)*new_coeff(ii,mu)*store_integrals(ii, mu, mu, 0)
        end do

        do jj = 1, N_modes
            if (jj /= ii) then
                step_multi_vec2(jj) = step_multi_vec2(jj)*mid_integral
            end if
        end do
    end do

    do ii = 1, N_modes
        do mu = 1, N_expansion
            do nu = 1, N_expansion
                dipole(ii, 1:3) = dipole(ii, 1:3) &
                    + full_coef(ii, 1, mu)*new_coeff(ii,nu) &
                    * store_integrals(ii, mu, nu, 1) &
                    * dipole_derivatives(ii,1:3) &
                    * step_multi_vec2(ii)
            end do
        end do
    end do

    intensity = ((sum(dipole(:, 1)))**2 + (sum(dipole(:, 2)))**2 + (sum(dipole(:, 3)))**2)*2000.
end if 

end subroutine constant_one_mode

end module one_mode