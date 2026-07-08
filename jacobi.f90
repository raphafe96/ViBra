module jacobi_diagonalization 
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



end module jacobi_diagonalization