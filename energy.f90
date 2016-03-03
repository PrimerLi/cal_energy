program main
implicit none
    interface 
        subroutine inverse(matrix, matrixInverse)
        implicit none
            complex(kind = 8), intent(in) :: matrix(:, :)
            complex(kind = 8), intent(out) :: matrixInverse(:, :)
        end subroutine inverse
        subroutine matrixProduct(matrixA, matrixB, matrixResult)
        implicit none
            complex(kind = 8), intent(in) :: matrixA(:, :), matrixB(:, :)
            complex(kind = 8), intent(out) :: matrixResult(:, :)
        end subroutine matrixProduct
        subroutine fit_cc_real(istart, iend, a, b, sigma_cc_real, xiom)
        implicit none
            integer, intent(in) :: istart, iend
            real(kind = 8), intent(out) :: a, b
            real(kind = 8), intent(in) :: sigma_cc_real(:), xiom(:)
        end subroutine fit_cc_real
        subroutine fit_cc_imag(istart, iend, b, sigma_cc_real, xiom)
        implicit none
            integer, intent(in) :: istart, iend
            real(kind = 8), intent(out) :: b
            real(kind = 8), intent(in) :: sigma_cc_real(:), xiom(:)
        end subroutine fit_cc_imag
        real(kind = 8) function cNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: mu, eps_f, beta, V
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function cNumber
        real(kind = 8) function fNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: mu, eps_f, beta, V
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function fNumber
        real(kind = 8) function kinetic_energy(beta, V, xiom, Giom, Niom, eps_f,&
                &chem, Sigma)
        implicit none
            real(kind = 8), intent(in) :: beta, V, eps_f, chem
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Giom(:, :, :), Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function kinetic_energy
        real(kind = 8) function potential_energy(beta, V, xiom, Giom, Niom,&
                &eps_f, chem, Sigma, cc_real_a, ff_real_a, Giom_0)
        implicit none
            real(kind = 8), intent(in) :: beta, V, eps_f, chem
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Giom(:, :, :), Sigma(:, :, :)
            complex(kind = 8), intent(in) :: Giom_0(:, :, :)
            integer, intent(in) :: Niom
            real(kind = 8), intent(in) :: cc_real_a, ff_real_a
        end function potential_energy
    end interface
    integer :: Niom
    integer :: matrixDimension
    integer :: i, j, nw
    integer :: row, col
    complex(kind = 8), dimension(:, :, :), allocatable :: Giom_0, Giom
    complex(kind = 8), dimension(:, :, :), allocatable :: Sigma
    complex(kind = 8), dimension(:, :, :), allocatable :: Sigma_temp, Giom_temp
    real(kind = 8), dimension(:), allocatable :: xiom
    real(kind = 8) :: real, imag
    complex(kind = 8), dimension(:, :), allocatable :: G_mat, G_mat_inv
    complex(kind = 8), dimension(:, :), allocatable :: S_mat, S_mat_inv
    complex(kind = 8), dimension(:, :), allocatable :: Z_mat, Z_mat_inv
    real(kind = 8) :: beta
    complex(kind = 8) :: s
    real(kind = 8) :: mu, chem, eps_f, temp
    integer :: Nuse_fit
    real(kind = 8) :: ratio 
    integer, parameter :: dp = 8
    real(kind = 8), dimension(:), allocatable :: sigma_cc_real, sigma_cc_imag
    real(kind = 8), dimension(:), allocatable :: sigma_ff_real, sigma_ff_imag
    integer :: istart, iend
    real(kind = 8) :: a, b
    real(kind = 8) :: cc_real_a, cc_real_b, ff_real_a, ff_real_b
    real(kind = 8) :: cc_imag_b, ff_imag_b

    complex(kind = 8) :: f
    complex(kind = 8) :: alpha_n, beta_n, gamma_n
    real(kind = 8) :: Gaussian
    real(kind = 8) :: V
    real(kind = 8) :: kinetic
    real(kind = 8) :: potential
    real(kind = 8) :: total_energy
    real(kind = 8) :: fill_c, fill_f, fill_total

    V = 0.955_dp

    ratio = 0.3_dp

    open(unit = 12, file = "mu_fill", action = "read")
    read(12, *) chem
    read(12, *) temp
    read(12, *) eps_f
    close(12)
    mu = chem

    open(unit = 20, file = "betaInfo", action = "read")
    read(20, *) beta
    close(20)

    open(unit = 12, file = "Niom", action = "read")
    read(12, *) Niom
    close(12)

    Nuse_fit = int(Niom*ratio)
    allocate(sigma_cc_real(Nuse_fit), sigma_cc_imag(Nuse_fit))
    allocate(sigma_ff_real(Nuse_fit), sigma_ff_imag(Nuse_fit))
    
    matrixDimension = 2
    allocate(Giom_0(Niom, matrixDimension, matrixDimension))
    allocate(Giom(Niom, matrixDimension, matrixDimension))
    allocate(Sigma(Niom, matrixDimension, matrixDimension))
    allocate(Sigma_temp(-Niom+1:Niom, matrixDimension, matrixDimension))
    allocate(Giom_temp(-Niom+1:Niom, matrixDimension, matrixDimension))
    allocate(xiom(Niom))
    allocate(G_mat(matrixDimension, matrixDimension))
    allocate(G_mat_inv(matrixDimension, matrixDimension))
    allocate(S_mat(matrixDimension, matrixDimension))
    allocate(S_mat_inv(matrixDimension, matrixDimension))
    allocate(Z_mat(matrixDimension, matrixDimension))
    allocate(Z_mat_inv(matrixDimension, matrixDimension))

    open(unit = 14, file = "g_iom_0_new", action = "read")
    do i = 1, matrixDimension
    do j = 1, matrixDimension
        read(14, *) row, col
        do nw = 1, Niom
            read(14, *) xiom(nw), Giom_0(nw, i, j)
        enddo
    enddo
    enddo
    close(14)

    open(unit = 16, file = "Self_energy", action = "read")
    do i = 1, matrixDimension
    do j = 1, matrixDimension
        read(16, *) row, col
        do nw = 1, Niom
            read(16, *) xiom(nw), real, imag
            Sigma(nw, i, j) = dcmplx(real, imag)
        enddo
    enddo
    enddo
    close(16)

    do nw = 1, Nuse_fit
        sigma_cc_real(nw) = dble(Sigma(nw,1,1))
        sigma_cc_imag(nw) = aimag(Sigma(nw,1,1))
        sigma_ff_real(nw) = dble(Sigma(nw,2,2))
        sigma_ff_imag(nw) = aimag(Sigma(nw,2,2))
    enddo
    istart = int(Nuse_fit*0.9)
    iend = Nuse_fit
    call fit_cc_real(istart, iend, a, b, sigma_cc_real, xiom)
    cc_real_a = a
    cc_real_b = b
    call fit_cc_imag(istart, iend, b, sigma_cc_imag, xiom)
    cc_imag_b = b
    call fit_cc_real(istart, iend, a, b, sigma_ff_real, xiom)
    ff_real_a = a
    ff_real_b = b
    call fit_cc_imag(istart, iend, b, sigma_ff_imag, xiom)
    ff_imag_b = b

  if (.false.) then
    open(unit = 20, file = "sigma_cc_real", action = "write")
    do nw = Niom/10, Niom
        write(20, *) xiom(nw), dble(Sigma(nw,1,1))
    enddo
    close(20)
    open(unit = 22, file = "fit_cc_real", action = "write")
    do nw = Niom/10, Niom
        write(22, *) xiom(nw), cc_real_a + cc_real_b/xiom(nw)**2
    enddo
    close(22)
    open(unit = 12, file = "sigma_cc_imag", action = "write")
    do nw = Niom/10, Niom
        write(12, *) xiom(nw), aimag(Sigma(nw,1,1))
    enddo
    close(12)
    open(unit = 14, file = "fit_cc_imag", action = "write")
    do nw = Niom/10, Niom
        write(14, *) xiom(nw), cc_imag_b/xiom(nw)
    enddo
    close(14)
  endif

    do nw = 1, Niom
        do i = 1, matrixDimension
        do j = 1, matrixDimension
            G_mat(i, j) = Giom_0(nw, i, j)
            S_mat(i, j) = Sigma(nw, i, j)
        enddo
        enddo
        call inverse(G_mat, G_mat_inv)
        Z_mat = G_mat_inv - S_mat
        call inverse(Z_mat, Z_mat_inv)
        Giom(nw, :, :) = Z_mat_inv(:, :)
    enddo

    do nw = 1, Niom
        Sigma_temp(nw, :, :) = Sigma(nw, :, :)
        Sigma_temp(1-nw, :, :) = conjg(Sigma(nw, :, :))
        Giom_temp(nw, :, :) = Giom(nw, :, :)
        Giom_temp(1-nw, :, :) = conjg(Giom(nw, :, :))
    enddo

  if (.false.) then
    open(unit = 18, file = "Giom_output", action = "write")
    do i = 1, matrixDimension
    do j = 1, matrixDimension
        write(18, *) i, j
        do nw = 1, Niom
            write(18, *) xiom(nw), Giom(nw, i, j)
        enddo
    enddo
    enddo
    close(18)
  endif

  if (.true.) then
    fill_c = 2*cNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
    fill_f = 2*fNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
    fill_total = fill_c + fill_f
    write(*, *) "c electron occupation number = ", fill_c
    write(*, *) "f electron occupation number = ", fill_f
    write(*, *) "Total filling number = ", fill_total
    potential = potential_energy(beta, V, xiom, Giom, Niom, eps_f, &
            &chem, Sigma, cc_real_a, ff_real_a, Giom_0)
    kinetic = 2.d0*kinetic_energy(beta, V, xiom, Giom, Niom, eps_f, &
            &chem, Sigma)
    total_energy = potential + kinetic
    write(*, *) "Kinetic energy = ", kinetic
    write(*, *) "Potential energy = ", potential
    write(*, *) "Total energy = ", total_energy
  endif

    deallocate(Giom, Giom_0, Sigma)
    deallocate(xiom)
    deallocate(G_mat, G_mat_inv)
    deallocate(S_mat, S_mat_inv)
    deallocate(Z_mat, Z_mat_inv)
    deallocate(sigma_cc_real, sigma_cc_imag)
    deallocate(sigma_ff_real, sigma_ff_imag)
    deallocate(Sigma_temp, Giom_temp)
end program main

subroutine inverse(matrix, matrixInverse)
implicit none
    complex(kind = 8), intent(in) :: matrix(:, :)
    complex(kind = 8), intent(out) :: matrixInverse(:, :)
    integer :: i, j, matrixDimension
    complex(kind = 8) :: determinant

    matrixDimension = size(matrix, 1)
    if (matrixDimension .ne. size(matrix, 2)) then
        write(*, *) "Only square matrix is accepted. "
        stop
    endif
    if (matrixDimension .ne. 2) then
        write(*, *) "Only square matrix of dimesion 2 is accepted. "
        stop
    endif

    determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
    if (dble(determinant) == 0 .and. aimag(determinant) == 0) then
        write(*, *) "The matrix is singular. "
        stop
    endif

    matrixInverse(1,1) = matrix(1,1)/determinant
    matrixInverse(1,2) = -matrix(1,2)/determinant
    matrixInverse(2,1) = -matrix(2,1)/determinant
    matrixInverse(2,2) = matrix(1,1)/determinant
end subroutine inverse

subroutine matrixProduct(matrixA, matrixB, matrixResult)
    implicit none
    complex(kind = 8), intent(in) :: matrixA(:, :), matrixB(:, :)
    complex(kind = 8), intent(out) :: matrixResult(:, :)
    complex(kind = 8) :: s
    integer :: i, j, k, matrixDimension

    matrixDimension = size(matrixA, 1)
    if (size(matrixA, 1) .ne. size(matrixA, 2) .or. &
    size(matrixA, 2) .ne. size(matrixB, 1)) then
        write(*, *) "Matrix dimension wrong. "
        stop
    endif

    do i = 1, matrixDimension
    do j = 1, matrixDimension
        s = dcmplx(0, 0)
        do k = 1, matrixDimension
            s = s + matrixA(i, k)*matrixB(k, j)
        enddo
        matrixResult(i, j) = s
    enddo
    enddo
end subroutine matrixProduct

real(kind = 8) function Gaussian(energy)
implicit none
    real(kind = 8), intent(in) :: energy
    real(kind = 8), parameter :: pi = acos(-1.d0)
    Gaussian = 1.d0/sqrt(pi)*exp(-energy**2)
end function Gaussian

complex(kind = 8) function f(gamma_n)
implicit none
    complex(kind = 8), intent(in) :: gamma_n
    real(kind = 8) :: energyUpper, energyLower
    integer :: Nenergy, ienergy
    real(kind = 8) :: denergy
    complex(kind = 8) :: s
    real(kind = 8) :: energy
    integer, parameter :: dp = 8
    real(kind = 8) :: Gaussian
    
    energyUpper = 20.d0
    energyLower = -energyUpper
    Nenergy = 1000
    denergy = (energyUpper - energyLower)/dble(Nenergy)
    s = dcmplx(0, 0)
    do ienergy = 1, Nenergy
        energy = energyLower + denergy*(ienergy-0.5_dp)
        s = s + (energy)/(gamma_n - energy)*Gaussian(energy)*denergy
    enddo
    f = s
end function f

!Function scc is the summation result for bare Green function Gcc
real(kind = 8) function scc(eps, eps_f, Delta, beta, Zplus, Zminus)
implicit none
    real(kind = 8), intent(in) :: eps, eps_f, Delta, beta
    real(kind = 8), intent(in) :: Zplus, Zminus
    integer, parameter :: dp = 8
    real(kind = 8) :: first, second
    first = 0.5_dp*(eps - eps_f +&
            &sqrt(Delta))/(sqrt(Delta))*(1.d0/(exp(beta*Zplus) + 1))
    second = 0.5_dp*(eps - eps_f - sqrt(Delta))/(-sqrt(Delta))*(1.d0/(exp(beta*Zminus) + 1))
    scc = first + second
end function scc

!Function cintegrand is to be used to calculate the kinetic energy contributed by
!c electrons. Integration over eps*cintegrand*rho gives the c electron kinetic energy. 
!Integration over cintegrand*rho gives the c electron occupation number.
!cintegrand = 1/beta*sum(1/(gamma_n - epsilon)). High frequency conditioning is
!added. 
real(kind = 8) function cintegrand(eps, mu, eps_f, beta, V, xiom, Sigma, Niom)
implicit none
    real(kind = 8), intent(in) :: mu
    real(kind = 8), intent(in) :: eps, eps_f, beta, V
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Sigma(:, :, :)
    integer, intent(in) :: Niom

    integer :: nw
    complex(kind = 8) :: alpha_n, beta_n, gamma_n
    complex(kind = 8) :: s
    real(kind = 8) :: scc
    real(kind = 8) :: tempResult
    real(kind = 8) :: Delta, Zplus, Zminus
    integer, parameter :: dp = 8

    tempResult = 0
    s = dcmplx(0, 0)
    do nw = 1, Niom
        alpha_n = dcmplx(0, xiom(nw)) - eps_f + mu
        beta_n = dcmplx(0, xiom(nw)) + mu
        gamma_n = beta_n - V**2/alpha_n
        s = s + 1.d0/(gamma_n - eps)
    enddo
    s = s/beta

    Delta = (eps - eps_f)**2 + 4*V**2
    Zplus = -mu + 0.5_dp*(eps + eps_f + sqrt(Delta))
    Zminus = -mu + 0.5_dp*(eps + eps_f - sqrt(Delta))
    tempResult = scc(eps, eps_f, Delta, beta, Zplus, Zminus) - 2.d0*dble(s)

    s = dcmplx(0, 0)
    do nw = 1, Niom
        alpha_n = dcmplx(0, xiom(nw)) - eps_f + mu - Sigma(nw, 2,2)
        beta_n = dcmplx(0, xiom(nw)) + mu - Sigma(nw, 1,1)
        gamma_n = beta_n - (V + Sigma(nw, 1,2))*(V + Sigma(nw,2,1))/alpha_n
        s = s + 1.d0/(gamma_n - eps)
    enddo
    s = s/beta
    tempResult = tempResult + 2.d0*dble(s)
    cintegrand = tempResult
end function cintegrand

!Function sff is the summation result for bare Green function Gff
real(kind = 8) function sff(eps, eps_f, Delta, beta, Zplus, Zminus)
implicit none
    real(kind = 8), intent(in) :: eps, eps_f, Delta, beta
    real(kind = 8), intent(in) :: Zplus, Zminus
    integer, parameter :: dp = 8
    real(kind = 8) :: first, second
    first = 0.5_dp*(-eps + eps_f +&
            &sqrt(Delta))/(sqrt(Delta))*(1.d0/(exp(beta*Zplus) + 1))
    second = 0.5_dp*(-eps + eps_f - sqrt(Delta))/(-sqrt(Delta))*(1.d0/(exp(beta*Zminus) + 1))
    sff = first + second
end function sff

!Integration over eps_f*fintegrand*rho gives the f electron on-site energy.
!Integration over fintengrand*rho gives the f electron occupation number.
!fintegrand = 1/beta*sum(1/alpha_n*(beta_n - epsilon)/(gamma_n - epsilon)). 
!High frequency conditioning is added. 
real(kind = 8) function fintegrand(eps, mu, eps_f, beta, V, xiom, Sigma, Niom)
implicit none
    real(kind = 8), intent(in) :: eps, mu, eps_f, beta, V
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Sigma(:, :, :)
    integer, intent(in) :: Niom

    real(kind = 8) :: sff
    real(kind = 8) :: Delta, Zplus, Zminus
    integer :: nw
    complex(kind = 8) :: s
    real(kind = 8) :: tempResult
    complex(kind = 8) :: alpha_n, beta_n, gamma_n
    integer, parameter :: dp = 8
    
    tempResult = 0
    s = dcmplx(0, 0)
    do nw = 1, Niom
        alpha_n = dcmplx(0, xiom(nw)) + mu - eps_f
        beta_n = dcmplx(0, xiom(nw)) + mu
        gamma_n = beta_n - V**2/alpha_n
        s = s + 1.d0/alpha_n*(beta_n - eps)/(gamma_n - eps)
    enddo
    s = s/beta

    Delta = (eps - eps_f)**2 + 4*V**2
    Zplus = -mu + 0.5_dp*(eps + eps_f + sqrt(Delta))
    Zminus = -mu + 0.5_dp*(eps + eps_f - sqrt(Delta))
    tempResult = sff(eps, eps_f, Delta, beta, Zplus, Zminus) - 2.d0*dble(s)

    s = dcmplx(0, 0)
    do nw = 1, Niom
        alpha_n = dcmplx(0, xiom(nw)) - eps_f + mu - Sigma(nw, 2,2)
        beta_n = dcmplx(0, xiom(nw)) + mu - Sigma(nw, 1,1)
        gamma_n = beta_n - (V + Sigma(nw, 1,2))*(V + Sigma(nw, 2,1))/alpha_n
        s = s + 1.d0/alpha_n*(beta_n - eps)/(gamma_n - eps)
    enddo
    s = s/beta
    tempResult = tempResult + 2.d0*dble(s)
    fintegrand = tempResult
end function fintegrand

real(kind = 8) function cKineticEnergy(mu, eps_f, beta, V, xiom, Sigma, Niom)
implicit none
    interface
        real(kind = 8) function cintegrand(eps, mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: eps, mu, eps_f, V, beta
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function cintegrand
    end interface
    real(kind = 8), intent(in) :: mu, eps_f, beta, V
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Sigma(:, :, :)
    integer, intent(in) :: Niom
    real(kind = 8) :: Gaussian
    real(kind = 8) :: energyUpper, energyLower, denergy
    integer :: Nenergy, ienergy
    real(kind = 8) :: eps
    real(kind = 8) :: s
    integer, parameter :: dp = 8
    
    s = 0
    Nenergy = 1000
    energyUpper = 20
    energyLower = -energyUpper
    denergy = (energyUpper - energyLower)/dble(Nenergy)
    do ienergy = 1, Nenergy
        eps = energyLower + denergy*(ienergy + 0.5_dp)
        s = s + eps*Gaussian(eps)*cintegrand(eps, mu, eps_f, beta, V, xiom, Sigma,&
                &Niom)*denergy
    enddo
    cKineticEnergy = s
end function cKineticEnergy

real(kind = 8) function cNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
implicit none
    interface
        real(kind = 8) function cintegrand(eps, mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: eps, mu, eps_f, V, beta
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function cintegrand
    end interface
    real(kind = 8), intent(in) :: mu, eps_f, beta, V
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Sigma(:, :, :)
    integer, intent(in) :: Niom
    real(kind = 8) :: Gaussian
    real(kind = 8) :: energyUpper, energyLower, denergy
    integer :: Nenergy, ienergy
    real(kind = 8) :: eps
    real(kind = 8) :: s
    integer, parameter :: dp = 8
    
    s = 0
    Nenergy = 1000
    energyUpper = 20
    energyLower = -energyUpper
    denergy = (energyUpper - energyLower)/dble(Nenergy)
    do ienergy = 1, Nenergy
        eps = energyLower + denergy*(ienergy + 0.5_dp)
        s = s + Gaussian(eps)*cintegrand(eps, mu, eps_f, beta, V, xiom, Sigma,&
                &Niom)*denergy
    enddo
    cNumber = s
end function cNumber

real(kind = 8) function fKineticEnergy(mu, eps_f, beta, V, xiom, Sigma, Niom)
implicit none
    interface
        real(kind = 8) function fintegrand(eps, mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: eps, mu, eps_f, beta, V
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function fintegrand
    end interface
    real(kind = 8), intent(in) :: mu, eps_f, beta, V
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Sigma(:, :, :)
    integer, intent(in) :: Niom
    real(kind = 8) :: Gaussian
    real(kind = 8) :: energyUpper, energyLower, denergy
    integer :: Nenergy, ienergy
    real(kind = 8) :: eps
    real(kind = 8) :: s
    integer, parameter :: dp = 8
    
    s = 0
    Nenergy = 1000
    energyUpper = 20
    energyLower = -energyUpper
    denergy = (energyUpper - energyLower)/dble(Nenergy)
    do ienergy = 1, Nenergy
        eps = energyLower + denergy*(ienergy + 0.5_dp)
        s = s + Gaussian(eps)*fintegrand(eps, mu, eps_f, beta, V, xiom, Sigma,&
                &Niom)*denergy
    enddo
    fKineticEnergy = s*eps_f
end function fKineticEnergy

real(kind = 8) function fNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
implicit none
    interface
        real(kind = 8) function fintegrand(eps, mu, eps_f, beta, V, &
                &xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: eps, mu, eps_f, beta, V
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function fintegrand
    end interface
    real(kind = 8), intent(in) :: mu, eps_f, beta, V
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Sigma(:, :, :)
    integer, intent(in) :: Niom
    real(kind = 8) :: Gaussian
    real(kind = 8) :: energyUpper, energyLower, denergy
    integer :: Nenergy, ienergy
    real(kind = 8) :: eps
    real(kind = 8) :: s
    integer, parameter :: dp = 8
    
    s = 0
    Nenergy = 1000
    energyUpper = 20
    energyLower = -energyUpper
    denergy = (energyUpper - energyLower)/dble(Nenergy)
    do ienergy = 1, Nenergy
        eps = energyLower + denergy*(ienergy + 0.5_dp)
        s = s + Gaussian(eps)*fintegrand(eps, mu, eps_f, beta, V, xiom, Sigma,&
                &Niom)*denergy
    enddo
    fNumber = s
end function fNumber

real(kind = 8) function kinetic_energy(beta, V, xiom, Giom, Niom, eps_f, chem, Sigma)
implicit none
    interface
        real(kind = 8) function cKineticEnergy(mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: mu, eps_f, beta, V
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function cKineticEnergy
        real(kind = 8) function fKineticEnergy(mu, eps_f, beta, V, xiom, Sigma, Niom)
        implicit none
            real(kind = 8), intent(in) :: mu, eps_f, beta, V
            real(kind = 8), intent(in) :: xiom(:)
            complex(kind = 8), intent(in) :: Sigma(:, :, :)
            integer, intent(in) :: Niom
        end function fKineticEnergy
    end interface
    real(kind = 8), intent(in) :: beta, V, eps_f, chem
    integer, intent(in) :: Niom
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Giom(:, :, :), Sigma(:, :, :)

    complex(kind = 8) :: f
    complex(kind = 8) :: s
    integer :: nw
    complex(kind = 8) :: alpha_n, beta_n, gamma_n
    !real(kind = 8) :: cKineticEnergy, fKineticEnergy

    s = dcmplx(0, 0)
    do nw = 1, Niom 
        !alpha_n = dcmplx(0, xiom(nw)) - (eps_f - chem) + Sigma(nw,2,2)
        !beta_n = dcmplx(0, xiom(nw)) + chem - Sigma(nw, 1, 1)
        !gamma_n = beta_n - (V + Sigma(nw,1,2))*(V + Sigma(nw, 2,1))/alpha_n
        s = s + V*Giom(nw, 2,1) + V*Giom(nw, 1,2)
    enddo
    s = s/dcmplx(beta, 0)
    kinetic_energy = 2.d0*dble(s) + cKineticEnergy(chem, eps_f, beta, V, xiom,&
            &Sigma, Niom) + fKineticEnergy(chem, eps_f, beta, V, xiom, Sigma, Niom)
end function kinetic_energy

real(kind = 8) function potential_energy(beta, V, xiom, Giom, Niom, eps_f,&
        &chem, Sigma, cc_real_a, ff_real_a, Giom_0)
implicit none
    interface 
    real(kind = 8) function cNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
    implicit none
        real(kind = 8), intent(in) :: mu, eps_f, beta, V
        real(kind = 8), intent(in) :: xiom(:)
        complex(kind = 8), intent(in) :: Sigma(:, :, :)
        integer, intent(in) :: Niom
    end function cNumber
    real(kind = 8) function fNumber(mu, eps_f, beta, V, xiom, Sigma, Niom)
    implicit none
        real(kind = 8), intent(in) :: mu, eps_f, beta, V
        real(kind = 8), intent(in) :: xiom(:)
        complex(kind = 8), intent(in) :: Sigma(:, :, :)
        integer, intent(in) :: Niom
    end function fNumber
    end interface 
    real(kind = 8), intent(in) :: beta, V, eps_f, chem
    real(kind = 8), intent(in) :: xiom(:)
    complex(kind = 8), intent(in) :: Giom(:, :, :), Sigma(:, :, :)
    complex(kind = 8), intent(in) :: Giom_0(:, :, :)
    integer, intent(in) :: Niom
    real(kind = 8), intent(in) :: cc_real_a, ff_real_a

    complex(kind = 8) :: s
    integer :: nw
    integer :: i, j, row, col, matrixDimension
    real(kind = 8) :: temp_result
    real(kind = 8) :: scc, sff
    real(kind = 8) :: Gaussian
    real(kind = 8) :: Delta, Zplus, Zminus
    real(kind = 8) :: integral
    real(kind = 8) :: energy, denergy
    integer :: Nenergy, ienergy
    real(kind = 8) :: energyLower, energyUpper
    integer, parameter :: dp = 8

    energyLower = -20.d0
    energyUpper = -energyLower
    Nenergy = 1000
    denergy = (energyUpper - energyLower)/dble(Nenergy)

    Delta = (chem - eps_f)**2 + 4*V**2
    Zplus = -chem + 0.5*(chem + eps_f + sqrt(Delta))
    Zminus = -chem + 0.5*(chem + eps_f - sqrt(Delta))
    
    s = dcmplx(0, 0)
    do nw = 1, Niom
        s = s + Sigma(nw, 1,1)*Giom(nw,1,1) + Sigma(nw,1,2)*Giom(nw,2,1)&
        &+Sigma(nw,2,1)*Giom(nw,1,2) + Sigma(nw,2,2)*Giom(nw,2,2)
    enddo
    s = s/beta
    temp_result = 2.d0*dble(s)
    
    integral = 0.d0
    do ienergy = 1, Nenergy
        energy = energyLower + denergy*(ienergy - 0.5_dp)
        integral = integral + Gaussian(energy)*scc(energy, eps_f, Delta, beta,&
                &Zplus, Zminus)*denergy
    enddo
    temp_result = temp_result + cc_real_a*integral

    integral = 0.d0
    do ienergy = 1, Nenergy
        energy = energyLower + denergy*(ienergy - 0.5_dp)
        integral = integral + Gaussian(energy)*sff(energy, eps_f, Delta, beta,&
                &Zplus, Zminus)*denergy
    enddo
    temp_result = temp_result + ff_real_a*integral

    s = dcmplx(0, 0)
    do nw = 1, Niom
        s = s + Giom_0(nw, 1,1)
    enddo
    s = s/beta
    temp_result = temp_result - 2.d0*dble(s)*cc_real_a

    s = dcmplx(0, 0)
    do nw = 1, Niom
        s = s + Giom_0(nw,2,2)
    enddo
    s = s/beta
    temp_result = temp_result - 2.d0*dble(s)*ff_real_a

    potential_energy = temp_result
end function potential_energy
