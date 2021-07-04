module Containers

export m, n, aggregate_matrix!

using LinearAlgebra

const m = 6
const n = 3

function aggregate_matrix!(l::Integer, layers::Integer, data::Matrix{Comlpex}, B::Matrix{Complex}, Bi::Matrix{Complex})

    r = data[:,1]
    μ = data[:,2]
    a = data[:,3]
    ρ = data[:,4]

    B = core_matrix(l, r[1], μ[1], a[1], ρ[1])
    Bi[:,:,1] = B

    Yi1 = zeros(Complex, m, m)
    Yi = zeros(Complex, m, m)

    for i in 2,layers

        Bi1 = B

        inverse_matrix!(l, r[i-1], μ[i], a[i-1], ρ[i], Yi1)
        propagtor_matrix!(l, r[i], μ[i], a[i], ρ[i], Yi)

        Yi1 = Yi * Yi1
        B = Yi1 * Bi1

        Bi[:,:,i] = B

    end

end


function core_matrix(l::Integer, r::Real, μ::Complex, a::Real, ρ::Real)

    B = zeros(Complex, m, n)

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3
    lm1 = l - 1
    lm2 = l - 2

    B[1,1] = l * r^lp1 / (2(2l + 3))
    B[2,1] = lp3 * r^lp1 / (2lp1 * (2l + 3))
    B[3,1] = (l*ρ*a*r + 2μ*(l^2 - lp3))*r^l / (2(2l + 3))
    B[4,1] = l*lp2*μ * r^l / (lp1 * (2l + 3))
    B[6,1] = 2π * ρ*l * r^lp1 / (2l + 3)

    B[1,2] = r^lm1
    B[2,2] = r^lm1 / l
    B[3,2] = (r*ρ*a + 2μ*lm1) * r^lm2
    B[4,2] = 2μ*lm1 * r^lm2 / l
    B[6,2] = 4π * ρ * r^lm1

    B[3,3] = -ρ * r^l
    B[5,3] = -r^l
    B[6,3] = -r^lm1 * (2l + 1)

    return B

end

function propagator_matrix!(l::Integer, r::Real, μ::Complex, a::Real, ρ::Real, Y::Matrix{Complex})

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3

    Y[1:m,1:n] = core_matrix(l, r, μ, a, ρ)

    Y[1,4] = lp1 * r^-l / (2(2l - 1))
    Y[2,4] = (2 - l)*r^-l / (2l * (2l - 1))
    Y[3,4] = (r*ρ*a*lp1 - 2μ*(l^2 + 3l - 1)) / (2r^lp1 * (2l - 1))
    Y[4,4] = μ*(l^2 - 1) / (l*r^lp1 * (2l - 1))
    Y[6,4] = 2π * ρ*lp1 / (r^l * (2l - 1))

    Y[1,5] = r^-lp2
    Y[2,5] = -r^-lp2 / lp1
    Y[3,5] = (r*ρ*a - 2μ * lp2) / r^lp3
    Y[4,5] = 2lp2 * μ / (lp1 * r^lp3)
    Y[6,5] = 4π * ρ / r^lp2

    Y[3,6] = -ρ / r^lp1
    Y[5,6] = -1 / r^lp1

end

function inverse_matrix!(l::Integer, r::Real, μ::Complex, a::Real, ρ::Real, Y::Matrix{Complex})

    lp1 = l + 1
    lp2 = l + 2
    lp3 = l + 3

    Y[:,:] .= 0

    Y[1,1] = ρ*a*r / μ - 2lp2
    Y[2,1] = -ρ*a*r / μ + 2(l^2 + 3l - 1) / lp1
    Y[3,1] = 4π * ρ
    Y[4,1] = ρ*a*r / μ + 2(l - 1)
    Y[5,1] = -ρ*a*r / μ - 2(l^2 - lp3) / l
    Y[6,1] = 4π * ρ * r

    Y[1,2] = 2l * lp2
    Y[2,2] = -2(l^2 - 1)
    Y[4,2] = 2(l^2 - 1)
    Y[5,2] = -2l * lp2

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
    Y[6,5] = 2l + 1

    Y[3,6] = -1
    Y[6,6] = -r

    V = zeros(Complex, 6)

    V[1] = lp1 / r^lp1
    V[2] = l*lp1 / (2r^(l - 1) * (2l - 1))
    V[3] = 1 / r^(l - 1)
    V[4] = l * r^l
    V[5] = l*lp1 * r^lp2 / (2(2l + 3))
    V[6] = -r^lp1

    V /= 2l + 1

    D = Diagonal(V)

    Y = D * Y

end

end # module
