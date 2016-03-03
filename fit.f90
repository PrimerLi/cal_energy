subroutine fit_cc_real(istart, iend, a, b, sigma_cc_real, xiom)
    implicit none
    integer, intent(in) :: istart, iend
    real(kind = 8), intent(out) :: a, b
    real(kind = 8), intent(in) :: sigma_cc_real(:)
    real(kind = 8), intent(in) :: xiom(:)
    real(kind = 8) :: s
    integer, parameter :: dp = 8
    real(kind = 8), dimension(:, :), allocatable :: matrix
    real(kind = 8), dimension(:, :), allocatable :: matrix_inverse
    real(kind = 8) :: determinant
    integer :: Norb
    integer :: i, j
    real(kind = 8), dimension(:), allocatable :: vector

    norb = 2
    allocate(matrix(norb, norb))
    allocate(matrix_inverse(norb, norb))
    allocate(vector(norb))
    matrix(1,1) = iend - istart + 1
    s = 0.d0
    do i = istart, iend
        s = s + 1.d0/(xiom(i)**2)
    enddo
    matrix(1,2) = s
    matrix(2,1) = s
    s = 0.d0
    do i = istart, iend
        s = s + 1.d0/(xiom(i)**4)
    enddo
    matrix(2,2) = s
    determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
    if (abs(determinant) < 1.0e-10) then
        write(*, *) "singular matrix in fit function. "
	write(*, *) matrix(1,1), matrix(1,2)
   	write(*, *) matrix(2,1), matrix(2,2)
        deallocate(matrix, matrix_inverse)
        deallocate(vector)
        stop
    endif

    matrix_inverse(1,1) = matrix(2,2)/determinant
    matrix_inverse(1,2) = -matrix(1,2)/determinant
    matrix_inverse(2,1) = -matrix(2,1)/determinant
    matrix_inverse(2,2) = matrix(1,1)/determinant

    s = 0.d0
    do i = istart, iend
        s = s + sigma_cc_real(i)
    enddo
    vector(1) = s
    s = 0.d0
    do i = istart, iend
        s = s + sigma_cc_real(i)/(xiom(i)**2)
    enddo
    vector(2) = s
    
    s = 0.d0
    do i = 1, norb
        s = s + matrix_inverse(1, i)*vector(i)
    enddo
    a = s
    s = 0.d0
    do i = 1, norb
        s = s + matrix_inverse(2,i)*vector(i)
    enddo
    b = s

    deallocate(matrix)
    deallocate(matrix_inverse)
    deallocate(vector)
end subroutine fit_cc_real

subroutine fit_cc_imag(istart, iend, b, sigma_cc_imag, xiom)
    implicit none
    integer, intent(in) :: istart, iend
    real(kind = 8), intent(out) :: b
    real(kind = 8), intent(in) :: sigma_cc_imag(:)
    real(kind = 8), intent(in) :: xiom(:)
    real(kind = 8) :: s
    real(kind = 8) :: left, right
    integer :: i

    s = 0.d0
    do i = istart, iend
        s = s + 1.d0/xiom(i)**2
    enddo
    left = s
    s = 0.d0
    do i = istart, iend
        s = s + sigma_cc_imag(i)/xiom(i)
    enddo
    right = s
    b = right/left
end subroutine fit_cc_imag
