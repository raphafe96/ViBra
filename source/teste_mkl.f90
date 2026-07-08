program test_dsyevd_parallel
    ! Include MKL module for Fortran
    use mkl_service  ! This provides mkl_get_max_threads, mkl_set_num_threads, etc.
    !use mkl_version  ! For version info
    implicit none
    
    integer, parameter :: n = 5000
    real(8), allocatable :: a(:,:), eigenvectors(:,:), eigenvalues(:), work(:)
    integer, allocatable :: iwork(:)
    integer :: lwork, liwork, info, i, j
    real(8) :: start_time, end_time
    character(200) :: mkl_version_string
    integer :: mkl_threads, mkl_dynamic
    
    ! Get MKL version
    call mkl_get_version_string(mkl_version_string)
    write(*,*) '========================================'
    write(*,*) 'MKL Version: ', trim(mkl_version_string)
    write(*,*) '========================================'
    
    ! Check thread settings using mkl_service functions
    mkl_threads = mkl_get_max_threads()
    mkl_dynamic = mkl_get_dynamic()
    
    write(*,*) 'Current MKL settings:'
    write(*,*) 'MKL max threads:', mkl_threads
    write(*,*) 'MKL dynamic:', mkl_dynamic
    write(*,*) ''
    
    ! Force thread settings for testing
    call mkl_set_dynamic(0)        ! Disable dynamic adjustment (0 = FALSE)
    call mkl_set_num_threads(8)    ! Set to desired number
    
    ! Verify settings were applied
    mkl_threads = mkl_get_max_threads()
    mkl_dynamic = mkl_get_dynamic()
    
    write(*,*) 'After forcing settings:'
    write(*,*) 'MKL max threads:', mkl_threads
    write(*,*) 'MKL dynamic:', mkl_dynamic
    write(*,*) ''
    
    ! Allocate matrices
    allocate(a(n,n), eigenvectors(n,n), eigenvalues(n))
    
    ! Create a symmetric matrix
    write(*,*) 'Creating symmetric matrix of size', n, 'x', n
    call random_number(a)
    do i = 1, n
        a(i,i) = a(i,i) * 1.0d0  ! Increase diagonal dominance
        do j = i+1, n
            a(j,i) = a(i,j)  ! Make symmetric
        end do
    end do
    eigenvectors = a  ! Save copy for dsyevd
    
    ! Query optimal workspace
    allocate(work(1), iwork(1))
    lwork = -1
    liwork = -1
    call dsyevd('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, iwork, liwork, info)
    
    if (info /= 0) then
        write(*,*) 'Warning: DSYEVD workspace query returned info =', info
    end if
    
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)
    
    ! Allocate workspace
    allocate(work(lwork), iwork(liwork))
    write(*,*) 'Workspace allocated: lwork =', lwork, 'liwork =', liwork
    write(*,*) ''
    
    ! Run dsyevd and time it
    write(*,*) 'Running DSYEVD (divide-and-conquer eigensolver)...'
    write(*,*) 'This should use multiple threads if MKL is configured correctly'
    write(*,*) 'Check MKL_VERBOSE output if enabled'
    write(*,*) ''
    
    call cpu_time(start_time)
    call dsyevd('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, iwork, liwork, info)
    call cpu_time(end_time)
    
    ! Report results
    write(*,*) '========================================'
    write(*,*) 'DSYEVD completed'
    write(*,*) 'Info:', info
    
    if (info == 0) then
        write(*,*) 'Time:', end_time - start_time, 'seconds'
        write(*,*) 'First few eigenvalues:', eigenvalues(1:min(5,n))
        write(*,*) 'Last few eigenvalues:', eigenvalues(max(1,n-4):n)
    else
        write(*,*) 'DSYEVD failed with info =', info
    end if
    write(*,*) '========================================'
    
    ! Clean up
    deallocate(a, eigenvectors, eigenvalues, work, iwork)
    
end program test_dsyevd_parallel
