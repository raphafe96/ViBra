module jacobi_diagonalization 
use read_input_file
  implicit none
  public dsyevd_A
  public dsyevr_A

contains

  !===========================================================
  ! Diagonalize a symmetric matrix using LAPACK's DSYEVD
  ! Input:  A - symmetric matrix (NxN), will be destroyed
  ! Output: eigvals - eigenvalues (N)
  !         eigvecs - eigenvectors (NxN) in columns
  !===========================================================
  subroutine dsyevd_A(A, eigvals, eigvecs)
    real(8), intent(inout) :: A(:,:)
    real(8), intent(out)   :: eigvals(:), eigvecs(:,:)
    
    integer :: N, info, lwork, liwork, mkl_threads, mkl_dynamic
    real(8), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    
    N = size(A, 1)
    
    ! Check if matrix is square
    if (size(A, 2) /= N) then
      print *, "Error: Matrix must be square"
      return
    end if
    
    ! Copy A to eigvecs (LAPACK will overwrite with eigenvectors)
    eigvecs = A
    
    ! Query optimal workspace size
    allocate(work(1))
    allocate(iwork(1))

    call dsyevd('V', 'L', N, eigvecs, N, eigvals, work, -1, iwork, -1, info)
    
    if (info /= 0) then
      print *, "Error: DSYEVD workspace query failed, info =", info
      deallocate(work, iwork)
      return
    end if
    
    lwork = int(work(1))
    liwork = iwork(1)
    
    deallocate(work, iwork)
    
    ! Allocate optimal workspace
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Solve eigenvalue problem
    call dsyevd('V', 'L', N, eigvecs, N, eigvals, work, lwork, iwork, liwork, info)
    
    call normalize_eigenvectors(eigvecs, size(eigvals), size(eigvals))

    if (info /= 0) then
      print *, "Error: DSYEVD diagonalization failed, info =", info
      if (info < 0) then
        print *, "The ", -info, "-th argument had an illegal value"
      else
        print *, "The algorithm failed to converge"
      end if
      deallocate(work, iwork)
      return
    end if
    
    deallocate(work, iwork)
    
  end subroutine 



subroutine dsyevr_A(A, nev, eigvals, eigvecs, jobz)
    real(8), intent(inout) :: A(:,:)
    integer, intent(in)    :: nev
    real(8), intent(out)   :: eigvals(:)
    real(8), intent(out)   :: eigvecs(:,:)
    character :: jobz      ! Compute eigenvalues and eigenvectors

    integer :: N, info, lwork, liwork, m
    real(8), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    integer, allocatable :: isuppz(:)
    
    ! Local variables for DSYEVR
    character :: range = 'I'     ! Get eigenvalues with indices from IL to IU
    character :: uplo = 'L'      ! Upper triangular part of A is stored
    integer :: il, iu            ! IL = 1, IU = nev (get first nev eigenvalues)
    real(8) :: vl, vu            ! Not referenced when range='I'
    real(8) :: abstol            ! Absolute error tolerance
    integer :: ldz               ! Leading dimension of z
    integer, allocatable :: ifail(:)  ! Not used with DSYEVR

    
    ! External function declaration
    real(8) :: dlamch
    external dlamch
    
    N = size(A, 1)
    
    ! Check if matrix is square
    if (size(A, 2) /= N) then
        print *, "Error: Matrix must be square"
        return
    end if
    
    ! Check if nev is valid
    if (nev < 1 .or. nev > N) then
        print *, "Error: nev must be between 1 and ", N
        return
    end if
    
    ! Set parameters for DSYEVR
    il = 1
    iu = nev
    ldz = N
    
    ! Set absolute tolerance (use LAPACK's default)
    abstol = 2.0d0 * dlamch('S')
    
    ! Allocate support array
    allocate(isuppz(2*N))
    
    ! Query optimal workspace size
    allocate(work(1))
    allocate(iwork(1))
    
    call dsyevr(jobz, range, uplo, N, A, N, vl, vu, il, iu, abstol, &
                m, eigvals, eigvecs, ldz, isuppz, work, -1, iwork, -1, info)
    
    if (info /= 0) then
        print *, "Error: DSYEVR workspace query failed, info =", info
        deallocate(work, iwork, isuppz)
        return
    end if
    
    lwork = int(work(1))
    liwork = iwork(1)
    
    deallocate(work, iwork)
    
    ! Allocate optimal workspace
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    ! Solve eigenvalue problem for first nev eigenvalues/vectors
    call dsyevr(jobz, range, uplo, N, A, N, vl, vu, il, iu, abstol, &
                m, eigvals, eigvecs, ldz, isuppz, work, lwork, iwork, liwork, info)

    call normalize_eigenvectors(eigvecs, nev, size(eigvals))
    
    if (info /= 0) then
        print *, "Error: DSYEVR diagonalization failed, info =", info
        if (info < 0) then
            print *, "The ", -info, "-th argument had an illegal value"
        else
            print *, "The algorithm failed to converge"
        end if
        deallocate(work, iwork, isuppz)
        return
    end if
    
    ! Check if we got the requested number of eigenvalues
    if (m /= nev) then
        print *, "Warning: DSYEVR returned", m, "eigenvalues, requested", nev
    end if
    
    deallocate(work, iwork, isuppz)
    
end subroutine

subroutine normalize_eigenvectors(eigvecs, nev, n)
    implicit none
    real(8), intent(inout) :: eigvecs(:,:)  ! Eigenvectors in columns
    integer, intent(in)    :: nev            ! Number of vectors
    integer, intent(in)    :: n              ! Dimension
    
    integer :: i
    real(8) :: norm
    
    do i = 1, nev
        norm = sqrt(dot_product(eigvecs(:,i), eigvecs(:,i)))
        if (norm > 1.0d-10) then  ! Avoid division by zero
            eigvecs(:,i) = eigvecs(:,i) / norm
        end if
    end do
    
end subroutine normalize_eigenvectors

!=============================================================================!
! Subroutine: jacobi_davidson_eigensolver                                     !
!                                                                             !
! Purpose:                                                                    !
!   Solve H * x = lambda * x for the lowest Nfirst eigenvalues using the      !
!   block Jacobi‑Davidson method with PCG inner solver.                       !
!                                                                             !
!   The correction equation is solved with a few steps of preconditioned      !
!   conjugate gradient, using the diagonal (H_diag - lambda) as               !
!   preconditioner.  This avoids the stagnation.                              !
!=============================================================================!
subroutine jacobi_davidson_eigensolver(N_dim, Nfirst, &
    n_sparse, sparse_m, sparse_n, H_sparse, &
    tol, max_iter, eigenvalues, eigenvectors)
  use omp_lib
  implicit none
  integer, intent(in) :: N_dim, Nfirst, n_sparse, max_iter
  integer(kind=4), intent(in) :: sparse_m(n_sparse), sparse_n(n_sparse)
  real*8,  intent(in) :: H_sparse(n_sparse)
  real*8,  intent(in) :: tol
  real*8,  intent(out) :: eigenvalues(Nfirst)
  real*8,  intent(out) :: eigenvectors(N_dim, Nfirst)


  ! --- local variables ---
  integer :: i, j, info, iter
  real*8, allocatable :: H_diag(:)
  real*8, allocatable :: V(:,:), W(:,:)
  real*8, allocatable :: H_sub(:,:), theta(:)
  real*8, allocatable :: r(:), t(:), u(:)
  real*8, allocatable :: Ht(:), Hu(:), p(:), Ap(:), z(:)
  real*8 :: alpha, norm_r, lambda
  integer, allocatable :: ind(:)
  logical :: converged(Nfirst)
  integer :: k_current, k_old, nb
  real*8, allocatable :: work(:)
  integer :: lwork
  real*8, allocatable :: dummy(:,:), dummy_eig(:)
  real*8, allocatable :: res_norms(:)
  real*8 :: max_res
  integer :: max_res_idx
  integer :: buffer

  ! inner PCG parameters
  integer, parameter :: max_inner_iter = 20    ! maximum CG iterations per correction
  real*8, parameter  :: inner_tol = 1.0d-5    ! relative tolerance for inner solver

  ! subspace dimensions for restarting
  integer :: max_basis, restart_size


  call read_intg('DAVBUF', buffer, 4000)
  max_basis    = max(buffer, 20 * Nfirst)        ! do not let the subspace exceed this size
  restart_size = min(max_basis, Nfirst * 6)    ! keep this many vectors after a restart



  !----------------------------------------------------------------------
  ! 1.  Extract diagonal
  !     The diagonal of H is used as the preconditioner and also for the
  !     initial subspace. We pull it out of the sparse storage.
  !----------------------------------------------------------------------
  allocate(H_diag(N_dim))
  H_diag = 0.0d0
  do i = 1, n_sparse
     if (sparse_m(i) == sparse_n(i)) H_diag(sparse_m(i)) = H_sparse(i)
  end do

  !----------------------------------------------------------------------
  ! 2.  Initial subspace
  !     We build an orthonormal basis V from unit vectors placed at the
  !     configurations with the smallest diagonal elements.  The idea
  !     is that the lowest eigenstates have large components on these
  !     configurations, so they provide a good starting guess.
  !----------------------------------------------------------------------
  allocate(V(N_dim, max_basis), W(N_dim, max_basis))
  V = 0.0d0; W = 0.0d0

  write(*,'(A,I6)') ' Jacobi-Davidson (PCG) solver: building initial subspace...'
  allocate(ind(N_dim))
  do i = 1, N_dim
     ind(i) = i
  end do
  call qsort_diagB(H_diag, ind, N_dim, Nfirst)          ! index array sorted by increasing H_diag

  k_current = 0
  do j = 1, Nfirst
     k_current = k_current + 1
     V(ind(j), k_current) = 1.0d0              ! unit vector for j‑th smallest diagonal
  end do

  call mgs(V, N_dim, k_current)               ! orthonormalise the initial basis
  call matvec_block(N_dim, k_current, V(:,1:k_current), W(:,1:k_current), &
                    n_sparse, sparse_m, sparse_n, H_sparse)   ! W = H * V

  iter = 0
  converged = .false.
  write(*,'(A)') ' Iter   Subspace   Max residual and its index'

  ! workspace query for the dense diagonalisation (DSYEV)
  allocate(dummy(max_basis, max_basis), dummy_eig(max_basis))
  allocate(work(1))
  lwork = -1
  call dsyev('V', 'U', max_basis, dummy, max_basis, dummy_eig, work, lwork, info)
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))
  deallocate(dummy, dummy_eig)

  ! allocate working vectors for the correction equations
  allocate(r(N_dim), t(N_dim), u(N_dim), Ht(N_dim), Hu(N_dim), &
           p(N_dim), Ap(N_dim), z(N_dim))

  !----------------------------------------------------------------------
  ! Main iteration loop – runs until all eigenvalues are converged or
  ! the maximum iteration count is reached.
  !----------------------------------------------------------------------
  do while (iter < max_iter .and. .not. all(converged))
     iter = iter + 1

     ! Rayleigh–Ritz step:
     !   - Build the projected matrix H_sub = V^T * W
     !   - Diagonalise it to obtain Ritz values theta and Ritz vectors y.
     !     After dsyev, H_sub holds the eigenvectors of the small matrix.
     allocate(H_sub(k_current, k_current), theta(k_current))
     call dgemm('T', 'N', k_current, k_current, N_dim, 1.0d0, V, N_dim, W, N_dim, 0.0d0, H_sub, k_current)
     call dsyev('V', 'U', k_current, H_sub, k_current, theta, work, lwork, info)

     k_old = k_current     ! size before adding new vectors

     ! Loop over the requested eigenvalues and compute corrections
     do j = 1, Nfirst
        lambda = theta(j)                   ! current Ritz value
        ! Form the Ritz vector u = V * y_j
        u = 0.0d0
        call dgemv('N', N_dim, k_current, 1.0d0, V, N_dim, H_sub(1,j), 1, 0.0d0, u, 1)
        ! Compute residual r = W*y_j - lambda * u   (and keep Hu = H*u = W*y_j)
        r = 0.0d0
        call dgemv('N', N_dim, k_current, 1.0d0, W, N_dim, H_sub(1,j), 1, 0.0d0, r, 1)
        Hu = r
        r = r - lambda * u
        norm_r = sqrt(dot_product(r,r))
        if (norm_r < tol) then
           converged(j) = .true.
           cycle
        end if
        converged(j) = .false.

        !---------------------------------------------------------------
        ! PCG inner solver – approximately solve the correction equation
        !   (I - u u')(H - lambda I)(I - u u') t = -r
        ! using the diagonal (H_diag - lambda) as preconditioner.
        ! The PCG iteration is carried out in the subspace orthogonal to u.
        !---------------------------------------------------------------
        t = 0.0d0
        block
          real*8 :: res(N_dim), res_new(N_dim), q(N_dim)
          real*8 :: rho, rho_old, beta, pap, alpha_cg, up
          integer :: cg_iter
          res = -r                          ! initial residual (b - A*0)
          do cg_iter = 1, max_inner_iter
             ! apply preconditioner: z = M^{-1} res, M = diag(H_diag - lambda)
             do i = 1, N_dim
                if (abs(H_diag(i)-lambda) > 1.d-12) then
                   z(i) = res(i) / (H_diag(i) - lambda)
                else
                   z(i) = 0.0d0
                end if
             end do
             rho = dot_product(res, z)
             if (cg_iter == 1) then
                p = z                           ! initial search direction
             else
                beta = rho / rho_old
                p = z + beta * p               ! updated search direction
             end if

             ! compute q = A * p, where A = (I - u u')(H - lambda I)(I - u u')
             call matvec_block(N_dim, 1, p, Ap, n_sparse, sparse_m, sparse_n, H_sparse)
             up = dot_product(u, p)
             q = p - up * u                    ! (I - u u') p
             q = (Ap - up * Hu) - lambda * q   ! (H - lambda I) applied
             up = dot_product(u, q)
             q = q - up * u                    ! (I - u u') again

             pap = dot_product(p, q)
             alpha_cg = rho / pap
             t = t + alpha_cg * p              ! update correction vector
             res_new = res - alpha_cg * q      ! update residual
             if (sqrt(dot_product(res_new, res_new)) < inner_tol * norm_r) then
                res = res_new
                exit                            ! PCG converged for this root
             end if
             rho_old = rho
             res = res_new
          end do
        end block

        ! Add the correction vector to the subspace (if it is not numerically zero)
        if (mgs_add(V, N_dim, k_current, t)) then
           if (k_current < max_basis) then
              k_current = k_current + 1
              V(:,k_current) = t
           end if
        end if
     end do

     ! Update W for the newly added basis vectors
     nb = k_current - k_old
     if (nb > 0) then
        call matvec_block(N_dim, nb, V(:, k_old+1:k_current), W(:, k_old+1:k_current), &
                          n_sparse, sparse_m, sparse_n, H_sparse)
     end if

          ! ---- Convergence monitoring 
     allocate(res_norms(Nfirst))
     do j = 1, Nfirst
        u = 0.0d0
        call dgemv('N', N_dim, k_current, 1.0d0, V, N_dim, H_sub(1,j), 1, 0.0d0, u, 1)
        r = 0.0d0
        call dgemv('N', N_dim, k_current, 1.0d0, W, N_dim, H_sub(1,j), 1, 0.0d0, r, 1)
        r = r - theta(j) * u
        res_norms(j) = sqrt(dot_product(r,r))
     end do

    BLOCK
    integer            :: nconv, worst_idx(5), best_idx(5), pos, k, j
    real*8             :: worst_val(5), best_val(5), val

    nconv = count(res_norms(1:Nfirst) < tol)

    ! Print header at first iteration and every 20 iterations
    if (iter == 1 .or. mod(iter, 20) == 0) then
        write(*,'(A)')    '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
        write(*,'(A)') ' Iter  Basis  Conv  Worst 5 (root:res)                                                                  |   Best 5 among unconverged'
        write(*,'(A)')    '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    end if

    ! Begin the output line
    write(*,'(I6,1X,I5,1X,I5,2X,A)',advance='no') iter, k_current, nconv, '['

    ! ---- Find the 5 largest (worst) residuals ----
    worst_val = -1.0d0
    worst_idx = 0
    do j = 1, Nfirst
        val = res_norms(j)
        if (val > worst_val(5)) then
            pos = 5
            do k = 4, 1, -1
                if (val > worst_val(k)) then
                    worst_val(k+1) = worst_val(k)
                    worst_idx(k+1) = worst_idx(k)
                    pos = k
                else
                    exit
                end if
            end do
            worst_val(pos) = val
            worst_idx(pos) = j
        end if
    end do

    ! Print worst pairs – fixed width (I4,":",ES9.2) gives exactly 14 characters
    do k = 1, min(5, Nfirst)
        if (k == 1) then
            write(*,'(I4,":",ES9.2)',advance='no') worst_idx(k), worst_val(k)
        else
            write(*,'(", ",I4,":",ES9.2)',advance='no') worst_idx(k), worst_val(k)
        end if
    end do
    write(*,'(A)',advance='no') ' ]'

    ! Separator
    write(*,'(A)',advance='no') '   |   '

    ! ---- Find and print the 5 smallest unconverged (best) residuals ----
    if (nconv < Nfirst) then
        write(*,'(A)',advance='no') '['
        best_val = huge(1.0d0)
        best_idx = 0
        do j = 1, Nfirst
            if (res_norms(j) >= tol) then
                val = res_norms(j)
                if (val < best_val(5)) then
                    pos = 5
                    do k = 4, 1, -1
                        if (val < best_val(k)) then
                            best_val(k+1) = best_val(k)
                            best_idx(k+1) = best_idx(k)
                            pos = k
                        else
                            exit
                        end if
                    end do
                    best_val(pos) = val
                    best_idx(pos) = j
                end if
            end if
        end do

        do k = 1, min(5, Nfirst - nconv)
            if (k == 1) then
                write(*,'(I4,":",ES9.2)',advance='no') best_idx(k), best_val(k)
            else
                write(*,'(", ",I4,":",ES9.2)',advance='no') best_idx(k), best_val(k)
            end if
        end do
        write(*,'(A)') ' ]'    ! closes the line
    else
        write(*,'(A)') 'All converged!'
    end if
    

    deallocate(res_norms)
   END BLOCK

     !------------------------------------------------------------------
     ! Restart if subspace is full and not yet converged.
     ! Keep only the restart_size best Ritz vectors, then re‑orthonormalise.
     !
     ! H_sub (already diagonalised above via dsyev) has ORTHONORMAL
     ! columns. V's columns are already orthonormal too (invariant kept
     ! by mgs/mgs_add throughout the run). So V_new = V * H_sub(:,1:keep)
     ! is orthonormal again *by construction* -- and because H is linear,
     ! applying the SAME projection to W gives W_new = H * V_new exactly,
     ! with NO new sparse matrix-vector product required. This replaces
     ! the previous fresh matvec_block(N_dim, k_current, ...) restart
     ! call, which is what caused the large transient memory spike (it
     ! forced matvec_block's per-thread duplicated buffer up to
     ! N_dim * restart_size * nthreads).
     !
     ! A full mgs() re-orthogonalisation of V here (as before) is
     ! deliberately NOT used: mgs mixes columns, and that mixing is not
     ! reflected in W, which would silently decouple W from H*V. Instead
     ! we only rescale each column by its own norm (round-off cleanup),
     ! which preserves w = H*v exactly since it is a per-column scalar
     ! operation, not a cross-column rotation.
     !------------------------------------------------------------------
     if (k_current >= max_basis .and. .not. all(converged)) then
        block
          real*8, allocatable :: V_new(:,:), W_new(:,:)
          integer :: keep, col
          real*8 :: colnorm
          keep = restart_size
          allocate(V_new(N_dim, keep), W_new(N_dim, keep))

          call dgemm('N', 'N', N_dim, keep, k_current, 1.0d0, V, N_dim, &
                     H_sub, k_current, 0.0d0, V_new, N_dim)
          call dgemm('N', 'N', N_dim, keep, k_current, 1.0d0, W, N_dim, &
                     H_sub, k_current, 0.0d0, W_new, N_dim)

          V(:, 1:keep) = V_new
          W(:, 1:keep) = W_new
          V(:, keep+1:max_basis) = 0.0d0
          W(:, keep+1:max_basis) = 0.0d0
          k_current = keep

          ! Defensive round-off cleanup only (see note above): rescale,
          ! do not re-orthogonalise across columns.
          do col = 1, k_current
             colnorm = sqrt(dot_product(V(:,col), V(:,col)))
             if (colnorm > 1.d-12) then
                V(:,col) = V(:,col) / colnorm
                W(:,col) = W(:,col) / colnorm
             end if
          end do

          deallocate(V_new, W_new)
        end block
        write(*,'(A,I4,A,I4)') ' *** Restarting: subspace reduced to ', restart_size, &
                               ' vectors at iteration ', iter
     end if

     deallocate(H_sub, theta)
  end do

  ! Final extraction of eigenvalues and eigenvectors from the converged subspace
  allocate(theta(k_current), H_sub(k_current, k_current))
  call dgemm('T', 'N', k_current, k_current, N_dim, 1.0d0, V, N_dim, W, N_dim, 0.0d0, H_sub, k_current)
  call dsyev('V', 'U', k_current, H_sub, k_current, theta, work, lwork, info)
  eigenvalues(1:Nfirst) = theta(1:Nfirst)
  call dgemm('N', 'N', N_dim, Nfirst, k_current, 1.0d0, V, N_dim, H_sub(1,1), k_current, 0.0d0, eigenvectors, N_dim)

  deallocate(work)
  deallocate(H_diag, V, W, theta, H_sub, ind)
  deallocate(r, t, u, Ht, Hu, p, Ap, z)
  write(*,'(A)')    '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  write(*,'(A,I4,A)') ' Jacobi-Davidson (PCG) converged after ', iter, ' iterations.'
  
end subroutine jacobi_davidson_eigensolver

!=============================================================================!
! Auxiliary routines with tightened rejection threshold (1.d-8)               !
!=============================================================================!

! Simple insertion sort to order the diagonal elements.  It returns indices
! in ascending order of the corresponding H_diag values. OLD OLD OLD Deprecate, I will keep just for... why not.
subroutine qsort_diag(arr, ind, n)
  implicit none
  integer, intent(in) :: n
  real*8,  intent(in) :: arr(n)
  integer, intent(inout) :: ind(n)
  integer :: i, j, tmp
  real*8 :: tmp_val
  do i = 2, n
     tmp_val = arr(ind(i))
     tmp = ind(i)
     j = i - 1
     do while (j >= 1 .and. arr(ind(j)) > tmp_val)
        ind(j+1) = ind(j)
        j = j - 1
     end do
     ind(j+1) = tmp
  end do
end subroutine qsort_diag


subroutine qsort_diagB(arr, ind, n, k)
  implicit none
  integer, intent(in)  :: n, k  
  real*8,  intent(in)  :: arr(n)     
  integer, intent(out) :: ind(n)     
                                 

  real*8  :: temp(n)                 
  integer :: i, idx
  ind = 0

  if (k > n) then
     print *, "Error: k > n in qsort_diag"
     stop
  end if
  temp = arr

  do i = 1, k
     idx = minloc(temp, dim=1)       
     ind(i) = idx                    
     temp(idx) = huge(1.0d0)         
  end do

end subroutine qsort_diagB

! Modified Gram–Schmidt orthonormalisation for a set of vectors.
! Vectors whose norm drops below 1.d-8 are considered numerically zero
! and are not normalised; this stops round‑off noise from entering the basis.
subroutine mgs(V, m, n)
  implicit none
  integer, intent(in) :: m, n
  real*8, intent(inout) :: V(m, n)
  integer :: i, j
  real*8 :: r
  do i = 1, n
     do j = 1, i-1
        r = dot_product(V(:,i), V(:,j))
        V(:,i) = V(:,i) - r * V(:,j)
     end do
     r = sqrt(dot_product(V(:,i), V(:,i)))
     if (r > 1.d-8) V(:,i) = V(:,i) / r
  end do
end subroutine mgs

!----------------------------------------------------------------------
! Function: mgs_add  (improved – iterative MGS for better orthogonality)
!   Orthogonalise u against V(1:n) using a twice‑is‑enough MGS.
!   If the norm after the first pass drops by a factor > 0.5,
!   a second pass is performed.  The final vector is normalised if its
!   norm is > 1.d-8.
!
!   Returns .true. if the vector survived (norm > 1.d-8) and was added,
!   .false. if it became too small after orthogonalisation and should be
!   discarded to avoid polluting the subspace with numerical noise.
!----------------------------------------------------------------------
logical function mgs_add(V, m, n, u)
  implicit none
  integer, intent(in)    :: m, n
  real*8,  intent(inout) :: V(m, n)
  real*8,  intent(inout) :: u(m)
  integer :: i, pass
  real*8  :: r, norm_old, norm_new

  norm_old = sqrt(dot_product(u, u))
  if (norm_old < 1.d-12) then
     mgs_add = .false.        ! vector already too small – reject
     return
  end if

  do pass = 1, 2
     do i = 1, n
        r = dot_product(u, V(:,i))
        u = u - r * V(:,i)
     end do
     norm_new = sqrt(dot_product(u, u))
     if (pass == 1 .and. norm_new > 0.5d0 * norm_old) exit   ! first pass was enough
     norm_old = norm_new
  end do

  if (norm_new > 1.d-8) then
     u = u / norm_new
     mgs_add = .true.
  else
     u = 0.0d0
     mgs_add = .false.
  end if
end function mgs_add

!=============================================================================!
! Subroutine: matvec_block                                                    !
!                                                                             !
! Sparse matrix-times-dense-block product: W_out = H_sparse * V_in.          !
! The sparse matrix is stored in symmetric upper-triangular form (m <= n).   !
!                                                                             !
! Two parallelisation strategies, chosen by the block width nb:              !
!                                                                             !
!  nb > 1  (initial subspace build, restart, batch of new correction        !
!           vectors): parallelise over the OUTPUT COLUMNS. Each thread owns !
!           whole columns of W_out exclusively, so there is no race and NO  !
!           per-thread duplication of the (N_dim, nb) block. Memory is      !
!           O(N_dim * nb), independent of the thread count -- this is what  !
!           avoids the large transient spike previously caused by calling   !
!           the entries-parallel path with nb = restart_size at restart.    !
!                                                                             !
!  nb == 1 (the PCG inner-loop matvec, called very frequently): there is    !
!           only one column, so column-parallelism gives no speedup.        !
!           Parallelise over the SPARSE ENTRIES instead, as before, but     !
!           with a private accumulator sized N_dim (not N_dim*nb*nthreads   !
!           -- nb is always 1 on this path, so this is unchanged from the   !
!           original cost for this call).                                  !
!=============================================================================!
subroutine matvec_block(N_dim, nb, V_in, W_out, n_sparse, sparse_m, sparse_n, H_sparse)
  use omp_lib
  implicit none
  integer, intent(in)  :: N_dim, nb, n_sparse
  real*8,  intent(in)  :: V_in(N_dim, nb)
  real*8,  intent(out) :: W_out(N_dim, nb)
  integer(kind=4), intent(in)  :: sparse_m(n_sparse), sparse_n(n_sparse)
  real*8,  intent(in)  :: H_sparse(n_sparse)

  integer :: idx, m, n, col, tid, nthreads
  real*8, allocatable :: W_priv(:,:)

  if (nb > 1) then
     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP& SHARED(n_sparse, sparse_m, sparse_n, H_sparse, V_in, W_out, nb) &
     !$OMP& PRIVATE(col, idx, m, n) &
     !$OMP& SCHEDULE(dynamic)
     do col = 1, nb
        W_out(:,col) = 0.0d0
        do idx = 1, n_sparse
           m = sparse_m(idx)
           n = sparse_n(idx)
           W_out(m,col) = W_out(m,col) + H_sparse(idx) * V_in(n,col)
           if (m /= n) W_out(n,col) = W_out(n,col) + H_sparse(idx) * V_in(m,col)
        end do
     end do
     !$OMP END PARALLEL DO
  else
     nthreads = omp_get_max_threads()
     allocate(W_priv(N_dim, 0:nthreads-1))
     W_priv = 0.0d0

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP& SHARED(n_sparse, sparse_m, sparse_n, H_sparse, V_in, W_priv) &
     !$OMP& PRIVATE(idx, m, n, tid) &
     !$OMP& SCHEDULE(dynamic, 512)
     do idx = 1, n_sparse
        tid = omp_get_thread_num()
        m = sparse_m(idx)
        n = sparse_n(idx)
        W_priv(m, tid) = W_priv(m, tid) + H_sparse(idx) * V_in(n, 1)
        if (m /= n) W_priv(n, tid) = W_priv(n, tid) + H_sparse(idx) * V_in(m, 1)
     end do
     !$OMP END PARALLEL DO

     W_out(:,1) = 0.0d0
     do tid = 0, nthreads-1
        W_out(:,1) = W_out(:,1) + W_priv(:,tid)
     end do
     deallocate(W_priv)
  end if
end subroutine matvec_block

end module jacobi_diagonalization