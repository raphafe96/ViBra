module combination
public
contains

    subroutine generate_combinations(N, M, vec, total_combinations, count_combo, numbering)
        implicit none
        integer, intent(in) :: N, M
        integer, dimension(N) :: combination
        integer :: total_combinations, numbering
        integer :: count_combo
        integer :: vec(total_combinations, N)
        
        count_combo = 0
        
        ! Start the recursive generation
        call recursive_generate(1, M, combination, N, count_combo, numbering, vec, total_combinations)
        
        write(101,*) count_combo, 'combinations for ', M, ' quantas' 
        
    end subroutine generate_combinations
    
    recursive subroutine recursive_generate(pos, remaining, combination, N, count_combo, numbering, vec, total_combinations)
        implicit none
        integer, intent(in) :: pos, remaining, N
        integer, intent(inout) :: combination(N), count_combo
        integer :: val, numbering, total_combinations
        integer :: vec(total_combinations, N)
        
        if (pos == N) then
            ! Last position: must take whatever is remaining
            combination(pos) = remaining
            count_combo = count_combo + 1
            call print_combination(combination, N, count_combo, numbering)
            vec(count_combo + numbering, :) = combination(:)
            return
        end if
        
        ! Try all possible values for the current position
        do val = remaining, 0, -1
            combination(pos) = val
            call recursive_generate(pos + 1, remaining - val, combination, N, count_combo, numbering, vec, total_combinations)
        end do
        
    end subroutine recursive_generate
    
    subroutine print_combination(arr, n, count, numbering)
        implicit none
        integer, intent(in) :: arr(:), n
        integer, intent(in) :: count
        integer :: i, numbering
        
        write(101, '(A,I6,A)', advance='no') ' state ', count+numbering, ': ('
        do i = 1, n
            if (i < n) then
                write(101, '(I0,A)', advance='no') arr(i), ','
            else
                write(101, '(I0,A)') arr(i), ')'
            end if
        end do
        
    end subroutine print_combination


    subroutine count_combinations(N, M, total) 
        implicit none
        integer, intent(in) :: N, M
        integer :: total
        integer :: i
        
        ! Calculate C(M+N-1, N-1) = C(M+N-1, M)
        total = 1
        do i = 1, min(N-1, M)
            total = total * (M + N - i) / i
        end do
        
    end subroutine


    subroutine count_unique_elements(arr, N, unique_values, counts, num_unique)
        implicit none
        integer, intent(in) :: arr(:)
        integer, intent(in) :: N
        integer, intent(out) :: unique_values(N)
        integer, intent(out) :: counts(N)
        integer, intent(out), optional :: num_unique

        integer :: i, j, local_num_unique
        logical :: found
        
        ! Inicializar
        unique_values = 0
        counts = 0
        local_num_unique = 0
        
        ! Para cada elemento do array
        do i = 1, N
            found = .false.
            
            ! Verificar se já está na lista de únicos
            do j = 1, local_num_unique
                if (unique_values(j) == arr(i)) then
                    ! Já existe, incrementar contagem
                    counts(j) = counts(j) + 1
                    found = .true.
                    exit
                end if
            end do
            
            ! Se não encontrado, adicionar como novo valor único
            if (.not. found) then
                local_num_unique = local_num_unique + 1
                unique_values(local_num_unique) = arr(i)
                counts(local_num_unique) = 1
            end if
        end do
        
        ! Se num_unique foi passado, atualizar seu valor
        if (present(num_unique)) then
            num_unique = local_num_unique
        end if
        
    end subroutine count_unique_elements


    !==========================================================================!
    ! Subroutine: build_inverted_index                                          !
    !                                                                           !
    ! Vectorizes the cubic and quartic potential tensors and builds the         !
    ! inverted index: for each mode, the list of potential terms that involve   !
    ! that mode. This index is precomputed once and passed to both the VSCF     !
    ! (one_mode) and VCI routines, avoiding redundant recomputation.            !
    !                                                                           !
    ! Outputs:                                                                  !
    !   Potential_3_vec  : cubic constants stored as 1D vector (i<=j<=k)        !
    !   Potential_4_vec  : quartic constants stored as 1D vector (i<=j<=k<=l)   !
    !   final_index_3/4  : unique mode indices for each potential term          !
    !   count_index_3/4  : multiplicity of each unique mode per term            !
    !   n_unique_3/4     : number of unique modes per term                      !
    !   unique_modes_3/4 : the unique mode indices per term                     !
    !   check3/4         : 1 if term is nonzero, 0 otherwise                    !
    !   cubic_for_mode   : inverted index: cubic terms involving each mode      !
    !   n_cubic_for_mode : number of cubic terms for each mode                  !
    !   quartic_for_mode : inverted index: quartic terms involving each mode    !
    !   n_quartic_for_mode: number of quartic terms for each mode               !
    !   n_cubic_max      : max entries per mode in cubic inverted index         !
    !   n_quartic_max    : max entries per mode in quartic inverted index       !
    !==========================================================================!
    subroutine build_inverted_index(N_modes, Potential_3, Potential_4, &
        total_3, total_4, &
        Potential_3_vec, Potential_4_vec, &
        final_index_3, count_index_3, n_unique_3, unique_modes_3, check3, &
        final_index_4, count_index_4, n_unique_4, unique_modes_4, check4, &
        cubic_for_mode, n_cubic_for_mode, n_cubic_max, &
        quartic_for_mode, n_quartic_for_mode, n_quartic_max)
    implicit none

    integer,  intent(in)  :: N_modes, total_3, total_4
    real*8,   intent(in)  :: Potential_3(N_modes, N_modes, N_modes)
    real*8,   intent(in)  :: Potential_4(N_modes, N_modes, N_modes, N_modes)

    real*8,   intent(out) :: Potential_3_vec(total_3)
    real*8,   intent(out) :: Potential_4_vec(total_4)
    integer,  intent(out) :: final_index_3(total_3, 3)
    integer,  intent(out) :: count_index_3(total_3, 3)
    integer,  intent(out) :: n_unique_3(total_3)
    integer,  intent(out) :: unique_modes_3(total_3, 3)
    integer,  intent(out) :: check3(total_3)
    integer,  intent(out) :: final_index_4(total_4, 4)
    integer,  intent(out) :: count_index_4(total_4, 4)
    integer,  intent(out) :: n_unique_4(total_4)
    integer,  intent(out) :: unique_modes_4(total_4, 4)
    integer,  intent(out) :: check4(total_4)
    integer,  intent(out) :: n_cubic_max, n_quartic_max
    integer,  allocatable, intent(out) :: cubic_for_mode(:,:)
    integer,  allocatable, intent(out) :: n_cubic_for_mode(:)
    integer,  allocatable, intent(out) :: quartic_for_mode(:,:)
    integer,  allocatable, intent(out) :: n_quartic_for_mode(:)

    integer :: i, j, k, l, ii, jj, pp, mode_idx
    integer :: array_index(4), final_index(4), count_index(4)

    ! --- Initialize ---
    check3 = 0
    check4 = 0
    Potential_3_vec = 0.d0
    Potential_4_vec = 0.d0
    final_index_3   = -1
    count_index_3   = -1
    final_index_4   = -1
    count_index_4   = -1
    n_unique_3      = 0
    unique_modes_3  = 0
    n_unique_4      = 0
    unique_modes_4  = 0

    ! --- Vectorize potentials and extract unique-mode lists ---
    ii = 0
    jj = 0
    do i = 1, N_modes
        do j = i, N_modes
            do k = j, N_modes
                ii = ii + 1
                Potential_3_vec(ii) = Potential_3(i,j,k)
                if (abs(Potential_3(i,j,k)) > 1.d-26) check3(ii) = 1

                array_index    = -1
                array_index(1) = i
                array_index(2) = j
                array_index(3) = k
                final_index    = -1
                count_index    = -1
                call count_unique_elements(array_index(1:3), 3, final_index, count_index)
                final_index_3(ii, 1:3) = final_index(1:3)
                count_index_3(ii, 1:3) = count_index(1:3)
                n_unique_3(ii) = 0
                unique_modes_3(ii, :) = 0
                do pp = 1, 3
                    if (final_index(pp) > 0) then
                        n_unique_3(ii) = n_unique_3(ii) + 1
                        unique_modes_3(ii, n_unique_3(ii)) = final_index(pp)
                    end if
                end do

                do l = k, N_modes
                    jj = jj + 1
                    Potential_4_vec(jj) = Potential_4(i,j,k,l)
                    if (abs(Potential_4(i,j,k,l)) > 1.d-26) check4(jj) = 1

                    array_index(1) = i
                    array_index(2) = j
                    array_index(3) = k
                    array_index(4) = l
                    final_index    = -1
                    count_index    = -1
                    call count_unique_elements(array_index(1:4), 4, final_index, count_index)
                    final_index_4(jj, 1:4) = final_index(1:4)
                    count_index_4(jj, 1:4) = count_index(1:4)
                    n_unique_4(jj) = 0
                    unique_modes_4(jj, :) = 0
                    do pp = 1, 4
                        if (final_index(pp) > 0) then
                            n_unique_4(jj) = n_unique_4(jj) + 1
                            unique_modes_4(jj, n_unique_4(jj)) = final_index(pp)
                        end if
                    end do
                end do ! l
            end do ! k
        end do ! j
    end do ! i

    ! --- Count how many terms involve each mode (first pass) ---
    allocate(n_cubic_for_mode(N_modes), n_quartic_for_mode(N_modes))
    n_cubic_for_mode  = 0
    n_quartic_for_mode = 0

    do jj = 1, total_3
        if (check3(jj) == 0) cycle
        do pp = 1, n_unique_3(jj)
            n_cubic_for_mode(unique_modes_3(jj,pp)) = &
                n_cubic_for_mode(unique_modes_3(jj,pp)) + 1
        end do
    end do
    do jj = 1, total_4
        if (check4(jj) == 0) cycle
        do pp = 1, n_unique_4(jj)
            n_quartic_for_mode(unique_modes_4(jj,pp)) = &
                n_quartic_for_mode(unique_modes_4(jj,pp)) + 1
        end do
    end do

    n_cubic_max   = max(1, maxval(n_cubic_for_mode))
    n_quartic_max = max(1, maxval(n_quartic_for_mode))

    ! --- Fill the inverted index arrays ---
    allocate(cubic_for_mode(n_cubic_max, N_modes))
    allocate(quartic_for_mode(n_quartic_max, N_modes))
    cubic_for_mode  = 0
    quartic_for_mode = 0
    n_cubic_for_mode  = 0
    n_quartic_for_mode = 0

    do jj = 1, total_3
        if (check3(jj) == 0) cycle
        do pp = 1, n_unique_3(jj)
            mode_idx = unique_modes_3(jj,pp)
            n_cubic_for_mode(mode_idx) = n_cubic_for_mode(mode_idx) + 1
            cubic_for_mode(n_cubic_for_mode(mode_idx), mode_idx) = jj
        end do
    end do
    do jj = 1, total_4
        if (check4(jj) == 0) cycle
        do pp = 1, n_unique_4(jj)
            mode_idx = unique_modes_4(jj,pp)
            n_quartic_for_mode(mode_idx) = n_quartic_for_mode(mode_idx) + 1
            quartic_for_mode(n_quartic_for_mode(mode_idx), mode_idx) = jj
        end do
    end do

    end subroutine build_inverted_index


    subroutine apply_degeneracy_factors(N_modes, Potential_3, Potential_4)
        implicit none
        
        ! Input/Output
        integer, intent(in) :: N_modes
        real(kind=8), intent(inout) :: Potential_3(N_modes, N_modes, N_modes)
        real(kind=8), intent(inout) :: Potential_4(N_modes, N_modes, N_modes, N_modes)
        
        ! Local variables
        integer :: i, j, k, l, i_perm, j_perm, k_perm, l_perm, p1, p2, p3, p4
        integer :: N_unique, perm_factor_3, perm_factor_4
        integer :: array_index(4), unique_values(4), counts(4)
        real(kind=8) :: temp_4(N_modes, N_modes, N_modes, N_modes), temp_3(N_modes, N_modes, N_modes)
        logical :: found_value

        !Orca por exemplo fornce termos 1322 mas ao organizar apenas a triangular o termo ficaria 1223, que nao eh fornecido pelo ORCA e portanto o termo seria cancelado. O que deve ser feito primeiro eh simetrizar por completo todas as components.
        !Actually I do not know how the constants are being stored... when I symmetrize completely it simply does not work for the quartic. I would rather use the constants the way they are from ORCA rathar than symmetrizing
        !I will keep this part for future testings...
        !Ok, new update, it seems I was testing for amonia, with degenerated modes, so the problem was not symmetrization????
        !NEw update :::::
        !It works :D

        do i = 1, N_modes
            do j = 1, N_modes
                do k = 1, N_modes
                    do l = 1, N_modes
                        
                        ! Skip if already processed
                        if (Potential_4(i,j,k,l) /= 0.0d0) cycle
                        
                        ! Check if all indices are the same - no permutation needed
                        if (i == j .and. j == k .and. k == l) cycle
                        
                        ! Store indices in array for permutation
                        array_index(1) = i
                        array_index(2) = j
                        array_index(3) = k
                        array_index(4) = l
                        
                        ! Try to find a non-zero value in permutations
                        found_value = .false.
                        
                        ! Generate all unique permutations
                        do p1 = 1, 4
                            do p2 = 1, 4
                                if (p2 == p1) cycle
                                do p3 = 1, 4
                                    if (p3 == p1 .or. p3 == p2) cycle
                                    p4 = 10 - p1 - p2 - p3  ! Since 1+2+3+4=10
                                    
                                    ! Get the permuted indices
                                    i_perm = array_index(p1)
                                    j_perm = array_index(p2)
                                    k_perm = array_index(p3)
                                    l_perm = array_index(p4)
                                    
                                    ! Check if this permutation has a non-zero value
                                    if (Potential_4(i_perm, j_perm, k_perm, l_perm) /= 0.0d0) then
                                        Potential_4(i,j,k,l) = Potential_4(i_perm, j_perm, k_perm, l_perm)
                                        found_value = .true.
                                        exit
                                    end if
                                end do
                                if (found_value) exit
                            end do
                            if (found_value) exit
                        end do
                        
                    end do
                end do
            end do
        end do

        !Now in case of potential ijk is not given entirely
        do i = 1, N_modes
            do j = 1, N_modes
                do k = 1, N_modes
                    
                    ! Skip if already processed
                    if (Potential_3(i,j,k) /= 0.0d0) cycle
                    
                    ! Check if all indices are the same - no permutation needed
                    if (i == j .and. j == k) cycle
                    
                    ! Store indices in array for permutation
                    array_index(1) = i
                    array_index(2) = j
                    array_index(3) = k
                    
                    ! Try to find a non-zero value in permutations
                    found_value = .false.
                    
                    ! Generate all unique permutations of 3 indices
                    do p1 = 1, 3
                        do p2 = 1, 3
                            if (p2 == p1) cycle
                            do p3 = 1, 3
                                if (p3 == p1 .or. p3 == p2) cycle
                                
                                ! Get the permuted indices
                                i_perm = array_index(p1)
                                j_perm = array_index(p2)
                                k_perm = array_index(p3)
                                
                                ! Check if this permutation has a non-zero value
                                if (Potential_3(i_perm, j_perm, k_perm) /= 0.0d0) then
                                    Potential_3(i,j,k) = Potential_3(i_perm, j_perm, k_perm)
                                    found_value = .true.
                                    exit
                                end if
                            end do
                            if (found_value) exit
                        end do
                        if (found_value) exit
                    end do
                    
                end do
            end do
        end do

        ! SALVAR os valores originais de Potential_4 ANTES de zerar qualquer coisa
        temp_4 = Potential_4
        temp_3 = Potential_3
        
        ! ZERAR todos os componentes de Potential_4 primeiro 
        Potential_4 = 0.0d0
        Potential_3 = 0.0d0
        
        ! Process 3-index terms (Potential_3)
        ! Primeiro, aplicar os fatores nos valores originais
        do i = 1, N_modes
            do j = i, N_modes
                do k = j, N_modes
                    
                    ! Analisar degenerescência para 3 índices
                    array_index(1) = i
                    array_index(2) = j
                    array_index(3) = k
                    
                    ! Contar elementos únicos
                    call count_unique_elements(array_index(1:3), 3, unique_values, counts, N_unique)
                    
                    ! Aplicar fatores de degenerescência para Potential_3
                    select case (N_unique)
                        case (1)
                            perm_factor_3 = 1
                        case (2)
                            perm_factor_3 = 3
                        case (3)
                            perm_factor_3 = 6
                        case default
                            perm_factor_3 = 1
                    end select
                    
                    Potential_3(i,j,k) = temp_3(i,j,k) * dble(perm_factor_3)
                    
                end do
            end do
        end do
        
        ! Process 4-index terms (Potential_4)
        ! Primeiro, restaurar os valores salvos e aplicar os fatores
        do i = 1, N_modes
            do j = i, N_modes
                do k = j, N_modes
                    do l = k, N_modes
                        
                        ! Analisar degenerescência para 4 índices
                        array_index(1) = i
                        array_index(2) = j
                        array_index(3) = k
                        array_index(4) = l
                        
                        ! Contar elementos únicos e suas frequências
                        call count_unique_elements(array_index, 4, unique_values, counts, N_unique)
                        
                        ! Aplicar fatores de degenerescência para Potential_4 baseado no padrão
                        select case (N_unique)
                            case (1)
                                perm_factor_4 = 1
                            case (2)
                                if (maxval(counts(:)) == 3) then
                                    perm_factor_4 = 4  ! 4!/(3!1!) = 4
                                else if (maxval(counts(:)) == 2) then
                                    perm_factor_4 = 6  ! 4!/(2!2!) = 6
                                else
                                    perm_factor_4 = 1
                                end if
                            case (3)
                                perm_factor_4 = 12  ! 4!/(2!1!1!) = 12
                            case (4)
                                perm_factor_4 = 24  ! 4! = 24
                            case default
                                perm_factor_4 = 1
                        end select
                        
                        ! Aplicar fator
                        Potential_4(i,j,k,l) = temp_4(i,j,k,l) * dble(perm_factor_4)
                        
                    end do
                end do
            end do
        end do

        ! (Os componentes fora já estão zero porque inicializamos Potential_4 = 0)
        
    end subroutine apply_degeneracy_factors


    !=============================================================================!
! Subroutine: classify_quartic_terms                                        !
!                                                                             !
! Classifies every nonzero quartic force-constant term by how many DISTINCT !
! mode indices it touches (n_unique_4), and reports how many fall into each !
! class (1/2/3/4 distinct modes). A term with n_unique_4 == 4 is a genuine  !
! quartic (Phi_ijkl, all indices different) -- the ONLY case that requires  !
! n_diff <= 4 in the sparse-pair-list cutoff. Everything with n_unique_4 <=3!
! is "semi-quartic" and only ever needs n_diff <= 3.                        !
!                                                                             !
!=============================================================================!
subroutine classify_quartic_terms(total_4, final_index_4, n_unique_4, &
    check4, Potential_4_vec, is_semiquartic)
    implicit none

    integer, intent(in)  :: total_4
    integer, intent(in)  :: final_index_4(total_4, 4)
    integer, intent(in)  :: n_unique_4(total_4)
    integer, intent(in)  :: check4(total_4)
    real*8,  intent(in)  :: Potential_4_vec(total_4)
    logical, intent(out) :: is_semiquartic

    integer :: jj, count_1, count_2, count_3, count_4, n_active
    real*8, parameter :: report_tol = 1.0d-16

    count_1 = 0
    count_2 = 0
    count_3 = 0
    count_4 = 0
    n_active = 0

    write(*,'(A)') ' '
    write(*,'(A)') ' >>> QUARTIC FORCE FIELD CLASSIFICATION (by distinct mode count)'

    do jj = 1, total_4
        if (check4(jj) == 0) cycle
        if (abs(Potential_4_vec(jj)) < report_tol) cycle
        n_active = n_active + 1

        select case (n_unique_4(jj))
        case (1)
            count_1 = count_1 + 1
        case (2)
            count_2 = count_2 + 1
        case (3)
            count_3 = count_3 + 1
        case (4)
            count_4 = count_4 + 1
            ! Print every genuine quartic term found -- there should be none
            ! for a semi-quartic-only force field, so this list should be empty.
            write(*,'(A,I8,A,4I5,A,ES14.4)') '  TRUE QUARTIC  term #', jj, &
                '  modes (i,j,k,l) = ', final_index_4(jj,1:4), &
                '  Phi = ', Potential_4_vec(jj)
        end select
    end do

    write(*,'(A)') ' ------------------------------------------------------------'
    write(*,'(A,I8)') '  Active (nonzero) quartic terms        : ', n_active
    write(*,'(A,I8)') '    n_unique = 1 (iiii)                 : ', count_1
    write(*,'(A,I8)') '    n_unique = 2 (iijj / iiik)          : ', count_2
    write(*,'(A,I8)') '    n_unique = 3 (iijk)                 : ', count_3
    write(*,'(A,I8)') '    n_unique = 4 (ijkl, TRUE quartic)   : ', count_4
    write(*,'(A)') ' ------------------------------------------------------------'

    is_semiquartic = (count_4 == 0)

    if (is_semiquartic) then
        write(*,'(A)') '  >>> No true 4-distinct-mode quartic terms found.'
        write(*,'(A)') '  >>> Safe to build the sparse pair list with n_diff <= 3.'
    else
        write(*,'(A,I6,A)') '  >>> WARNING: ', count_4, &
            ' true quartic term(s) found (n_unique_4 == 4).'
        write(*,'(A)') '  >>> This is NOT a semi-quartic-only force field.'
        write(*,'(A)') '  >>> The sparse pair list MUST use n_diff <= 4, or these'
        write(*,'(A)') '  >>> couplings will be silently dropped.'
    end if
    write(*,'(A)') ' ============================================================'
    write(*,'(A)') ' '

    end subroutine classify_quartic_terms

end module combination