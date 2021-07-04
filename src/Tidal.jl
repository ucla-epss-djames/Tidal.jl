module Tidal

include("Containers.jl")

using .Containers

export propagtor_method, normalize!

const G = 6.67408e-11

# Functions to add here
# propagator_method
function propagtor_method(l::Integer, layers::Integer, data::Matrix{Complex})

    B = zeros(Complex, m, n)
    Bi = zeros(Complex, m, n, layers)

    aggregate_matrix!(l, layers, data, B, Bi)

    bound = solve_vector()

    tidal = solve_layer()

end

# solve_vector
function solve_vector(l::Integer, r::Real, B::Matrix{Complex})

    b = zeros(Complex, n)
    b[3] = (-2l - 1) / r

    M = zeros(Complex, n, n)

    M[1,:] = B[3,:]
    M[2,:] = B[4,:]
    M[3,:] = B[6,:]

    c = M \ b

    return c

end

# solve_layer
function solve_layer(l::Integer, layers::Integer, data::Matrix{Complex}, c::Vector{Complex}, Bi::Array{Matrix})

    r = data[:,1]
    a = data[:,3]
    y = zeros(Complex, m)
    tidal = zeros(Complex, layers, 3)

    for i in 1, layers

        y = Bi[:,:,i] * c

        kl = -r[i]^l - y[5]
        hl = a[i] * y[1]
        ll = a[i] * y[2]

        tidal[i,:] = [kl, hl, ll]

    end

end

# normalize
function normalize!(mass::Real, r::Real, data::Matrix{Complex})

    a = mass / r^2 * G
    ρ = mass / r^3
    p = ρ * a * r

    data[:,1] /= r
    data[:,2] /= p
    data[:,3] /= a
    data[:,4] /= ρ

end


end # module
