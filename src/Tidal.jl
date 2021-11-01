module Tidal

include("Containers.jl")

using .Containers

using PhysicalConstants.CODATA2018: G

export propagator_method, normalize!

# Functions to add here
# propagator_method
function propagator_method(l::Integer, layers::Integer, data::Matrix{Complex},
                          flag::Bool)

    B, Bi = aggregate_matrix(l, layers, data, flag)

    bound = solve_vector(l, real(data[layers,1]), B)

    tidal = solve_layer(l, layers, data, bound, Bi)

end

# solve_vector
function solve_vector(l::Integer, r::Real, B::Matrix)

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
function solve_layer(l::Integer, layers::Integer, data::Matrix{Complex}, c::Vector, Bi::Array{Complex})

    r = data[:,1]
    a = data[:,3]
    y = zeros(Complex, m)
    tidal = zeros(Complex, layers, 3)

    for i in 1:layers

        y = Bi[:,:,i] * c

        kl = -r[i]^l - y[5]
        hl = a[i] * y[1]
        ll = a[i] * y[2]

        tidal[i,:] = [kl, hl, ll]

    end

    return tidal

end

function use_GD_core(r::Real, μ::Real)

    if r > 0.0 && μ == 0
        flag = true
    else
        flag = false
    end

    return flag

end


# normalize
function normalize!(mass::Real, r::Real, data::Matrix{Complex})

    a = mass / r^2 * G.val
    ρ = mass / r^3
    p = ρ * a * r

    data[:,1] /= r
    data[:,2] /= p
    data[:,3] /= a
    data[:,4] /= ρ

end


end # module
