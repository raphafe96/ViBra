Module compute_integrals

contains

subroutine sobrepos(i, j, coeff_i, coeff_j, integral, freq, power)
implicit none

integer :: i, j, power
real*8 :: coeff_i, coeff_j, integral, freq
real*8 :: omega, temp

integral = 0.d0
omega = freq

if (power .lt. 5) then
    omega = 1.0
end if

! Overlap (power = 0)
if (power == 0) then
    if (i == j) then
        integral = coeff_i * coeff_j
    end if
    return
end if

! Position operator Q (power = 1)
if (power == 1) then
    if (j == i + 1) then
        temp = sqrt(dble(i + 1) / (2.0d0 * omega))
    else if (j == i - 1) then
        temp = sqrt(dble(i) / (2.0d0 * omega))
    else
        temp = 0.0d0
    end if
    integral = coeff_i * coeff_j * temp
    return
end if

! Q² operator (power = 2)
if (power == 2) then
    if (i == j) then
        temp = (dble(i) + 0.5d0) / omega
    else if (j == i + 2) then
        temp = sqrt(dble((i + 1) * (i + 2))) / (2.0d0 * omega)
    else if (j == i - 2) then
        if (i >= 2) then
            temp = sqrt(dble(i * (i - 1))) / (2.0d0 * omega)
        else
            temp = 0.0d0
        end if
    else
        temp = 0.0d0
    end if
    integral = coeff_i * coeff_j * temp
    return
end if

! Q³ operator (power = 3)
if (power == 3) then
    if (j == i + 1) then
        temp = (3.0d0 * dble(i + 1) * sqrt(dble(i + 1))) / (2.0d0 * omega * sqrt(2.0d0 * omega))
    else if (j == i - 1) then
        temp = (3.0d0 * dble(i) * sqrt(dble(i))) / (2.0d0 * omega * sqrt(2.0d0 * omega))
    else if (j == i + 3) then
        temp = sqrt(dble((i + 1) * (i + 2) * (i + 3))) / (2.0d0 * omega * sqrt(2.0d0 * omega))
    else if (j == i - 3) then
        if (i >= 3) then
            temp = sqrt(dble(i * (i - 1) * (i - 2))) / (2.0d0 * omega * sqrt(2.0d0 * omega))
        else
            temp = 0.0d0
        end if
    else
        temp = 0.0d0
    end if
    integral = coeff_i * coeff_j * temp
    return
end if

! Q⁴ operator (power = 4)
if (power == 4) then
    if (i == j) then
        temp = (6.0d0 * dble(i) * dble(i + 1) + 3.0d0) / (4.0d0 * omega * omega)
    else if (j == i + 2) then
        temp = sqrt(dble((i + 1) * (i + 2))) * dble(2*i + 3) / (2.0d0 * omega * omega)
    else if (j == i - 2) then
        if (i >= 2) then
            temp = sqrt(dble(i * (i - 1))) * dble(2*i - 1) / (2.0d0 * omega * omega)
        else
            temp = 0.0d0
        end if
    else if (j == i + 4) then
        temp = sqrt(dble((i + 1) * (i + 2) * (i + 3) * (i + 4))) / (4.0d0 * omega * omega)
    else if (j == i - 4) then
        if (i >= 4) then
            temp = sqrt(dble(i * (i - 1) * (i - 2) * (i - 3))) / (4.0d0 * omega * omega)
        else
            temp = 0.0d0
        end if
    else
        temp = 0.0d0
    end if
    integral = coeff_i * coeff_j * temp
    return
end if

! KINETIC ENERGY OPERATOR T = -½ ∂²/∂Q² (power = 5)
if (power == 5) then
    if (i == j) then
        ! Diagonal: T_nn = ω(n+1/2)/2
        temp = omega * (dble(i) + 0.5d0) / 2.0d0
    else if (j == i + 2) then
        ! T_{n,n+2} = -ω/4 * √[(n+1)(n+2)]
        temp = -omega * sqrt(dble((i + 1) * (i + 2))) / 4.0d0
    else if (j == i - 2) then
        ! T_{n,n-2} = -ω/4 * √[n(n-1)]
        if (i >= 2) then
            temp = -omega * sqrt(dble(i * (i - 1))) / 4.0d0
        else
            temp = 0.0d0
        end if
    else
        temp = 0.0d0
    end if
    integral = coeff_i * coeff_j * temp
    return
end if

! If ∂²/∂Q² itself (not -½ of it), use power = 6
if (power == 6) then
    if (i == j) then
        ! ⟨n|∂²/∂Q²|n⟩ = -ω(n+1/2)
        temp = -omega * (dble(i) + 0.5d0)
    else if (j == i + 2) then
        ! ⟨n+2|∂²/∂Q²|n⟩ = ω√[(n+1)(n+2)]/2
        temp = omega * sqrt(dble((i + 1) * (i + 2))) / 2.0d0
    else if (j == i - 2) then
        if (i >= 2) then
            temp = omega * sqrt(dble(i * (i - 1))) / 2.0d0
        else
            temp = 0.0d0
        end if
    else
        temp = 0.0d0
    end if
    integral = coeff_i * coeff_j * temp
    return
end if

end subroutine

end module