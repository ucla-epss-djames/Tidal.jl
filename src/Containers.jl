module Containers

export m, n, aggregate_matrix

using LinearAlgebra

const m = 6
const n = 3

function aggregate_matrix(l::Integer, layers::Integer, data::Matrix{Complex},
                          flag::Bool)

    B = zeros(Complex, m, n)
    Bi = zeros(Complex, m, n, layers)

    r = real(data[:,1])
    μ = data[:,2]
    a = real(data[:,3])
    ρ = real(data[:,4])

    if flag
        B = core_GD_matrix(l, r[1], a[1], ρ[1])
    else
        B = core_HH_matrix(l, r[1], μ[1], a[1], ρ[1])
    end

    Bi[:,:,1] = B

    Yi1 = zeros(Complex, m, m)
    Yi = zeros(Complex, m, m)

    for i in 2:layers

        Bi1 = B

        Yi1 = inverse_matrix!(Yi1, l, r[i-1], μ[i], a[i-1], ρ[i])
        Yi = propagator_matrix!(Yi, l, r[i], μ[i], a[i], ρ[i])

        Yi1 = Yi * Yi1
        B = Yi1 * Bi1

        Bi[:,:,i] = B

    end

    return B, Bi

end

function core_GD_matrix(l::Integer, r::Real, a::Real, ρ::Real)

    B = zeros(Complex, m, n)

    Ac = a / r

    B[1,1] = -r^l / a
    B[5,1] = r^l
    B[6,1] = 2*(l - 1)*r^(l - 1)

    B[2,2] = 1

    B[1,3] = 1
    B[3,3] = ρ * r * Ac
    B[6,3] = 3 * Ac

    return B

end

function core_HH_matrix(l::Integer, r::Real, μ::Complex, a::Real, ρ::Real)

    B = zeros(Complex, m, n)

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3
    lm1 = l - 1
    lm2 = l - 2

    B[1,1] = l * r^lp1 / (2*(2*l + 3))
    B[2,1] = lp3 * r^lp1 / (2*lp1 * (2*l + 3))
    B[3,1] = (l*ρ*a*r + 2*μ*(l^2 - lp3))*r^l / (2*(2*l + 3))
    B[4,1] = l*lp2*μ * r^l / (lp1 * (2*l + 3))
    B[6,1] = 2*π * ρ*l * r^lp1 / (2*l + 3)

    B[1,2] = r^lm1
    B[2,2] = r^lm1 / l
    B[3,2] = (r*ρ*a + 2*μ*lm1) * r^lm2
    B[4,2] = 2*μ*lm1 * r^lm2 / l
    B[6,2] = 4*π * ρ * r^lm1

    B[3,3] = -ρ * r^l
    B[5,3] = -r^l
    B[6,3] = -r^lm1 * (2*l + 1)

    return B

end

function propagator_matrix!(Y::Matrix{Complex}, l::Integer, r::Real,
                            μ::Complex, a::Real, ρ::Real)

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3

    Y[:,:] .= 0

    Y[1:m,1:n] = core_HH_matrix(l, r, μ, a, ρ)

    Y[1,4] = lp1 * r^-l / (2*(2*l - 1))
    Y[2,4] = (2 - l)*r^-l / (2*l * (2*l - 1))
    Y[3,4] = (r*ρ*a*lp1 - 2*μ*(l^2 + 3*l - 1)) / (2*r^lp1 * (2*l - 1))
    Y[4,4] = μ*(l^2 - 1) / (l*r^lp1 * (2*l - 1))
    Y[6,4] = 2*π * ρ*lp1 / (r^l * (2*l - 1))

    Y[1,5] = r^-lp2
    Y[2,5] = -r^-lp2 / lp1
    Y[3,5] = (r*ρ*a - 2*μ * lp2) / r^lp3
    Y[4,5] = 2*lp2 * μ / (lp1 * r^lp3)
    Y[6,5] = 4*π * ρ / r^lp2

    Y[3,6] = -ρ / r^lp1
    Y[5,6] = -1 / r^lp1

    return Y

end

function inverse_matrix!(Y::Matrix, l::Integer, r::Real, μ::Complex, a::Real, ρ::Real)


    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3

    Y[:,:] .= 0

    Y[1,1] = ρ*a*r / μ - 2*lp2
    Y[2,1] = -ρ*a*r / μ + 2*(l^2 + 3*l - 1) / lp1
    Y[3,1] = 4*π * ρ
    Y[4,1] = ρ*a*r / μ + 2*(l - 1)
    Y[5,1] = -ρ*a*r / μ - 2*(l^2 - lp3) / l
    Y[6,1] = 4*π * ρ * r

    Y[1,2] = 2*l * lp2
    Y[2,2] = -2*(l^2 - 1)
    Y[4,2] = 2*(l^2 - 1)
    Y[5,2] = -2*l * lp2

    Y[1,3] = -r / μ
    Y[2,3] = r / μ
    Y[4,3] = -r / μ
    Y[5,3] = r / μ

    Y[1,4] = l*r / μ
    Y[2,4] = r*(2 - l) / μ
    Y[4,4] = -lp1*r / μ
    Y[5,4] = lp3*r / μ

    Y[1,5] = ρ*r / μ
    Y[2,5] = -ρ*r / μ
    Y[4,5] = ρ*r / μ
    Y[5,5] = -ρ*r / μ
    Y[6,5] = 2*l + 1

    Y[3,6] = -1
    Y[6,6] = -r

    V = zeros(Complex, 6)

    V[1] = lp1 / r^lp1
    V[2] = l*lp1 / (2*r^(l - 1) * (2*l - 1))
    V[3] = 1 / r^(l - 1)
    V[4] = l * r^l
    V[5] = l*lp1 * r^lp2 / (2*(2*l + 3))
    V[6] = -r^lp1

    V /= 2*l + 1

    D = Diagonal(V)

    Y = D * Y

    return Y

end

end # module
