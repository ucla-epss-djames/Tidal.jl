module Tidal
# Author: David James, davidabraham@ucla.edu
# Background: The following is the tidal love number calculations.
# Refer to Henning & Hurford '14 (HH'14) for a full description of
# the method. Core matrices come from HH'14 or the "Global Dynamics"
# book by Sabadini (GD). Furthermore, corrections taken from Spada
# et al. '92 and Vermeersen et al. '96.

include("containers.jl")

using PhysicalConstants.CODATA2018: G

export propagator_method, normalize!

"""
    propagator_method(l::Integer, layers::Integer,
                      data::Matrix{Complex}, flag::Bool)

Performs the propagator method from HH'14 to produce the tidal love
numbers of harmonic degree `l` for a planetary body.

The body is described in the `data` matrix where `layers` is the
number of rows and each column is the radius, complex shear modulus,
gravitational acceleration, and density, respectively. All units in
the matrix should be SI and normalized by `normalize!`.

Finally, the `flag` designates which to use. `True` uses the GD core
and `False` uses the HH core. Refer to `core_GD_matrix`,
`core_HH_matrix`, and/or `use_GD_core` for more information.
"""
function propagator_method(l::Integer, layers::Integer,
                           data::Matrix{Complex}, flag::Bool)

    B, Bi = aggregate_matrix(l, layers, data, flag)

    bound = solve_vector(l, real(data[layers,1]), B)

    tidal = solve_layer(l, layers, data, bound, Bi)

end

"""
    solve_vector(l::Integer, r::Real, B::Matrix)

Solves out the bound vector given a harmonic degree `l`, the final
radiual value `r` of the `data` matrix, and the uppermost aggregate
matrix `B` outputted from `aggregate_matrix`.

Refer to A6-A7 of HH'14.
"""
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

"""
    solve_layer(l::Integer, layers::Integer, data::Matrix{Complex},
                c::Vector, Bi::Array{Complex})

Solves out the tidal values given a harmonic degree `l`, rows of the
`data` matrix known as `layers`, the structure of the planet `data`,
the solved vector `c` from `solve_vector`, and the aggregate matrices
`Bi` from `aggregate_matrix`.

Refer to HH'14 A8-A11. Furthermore, refer to `propagator_method` for
information on the `data` matrix.
"""
function solve_layer(l::Integer, layers::Integer,
                     data::Matrix{Complex}, c::Vector,
                     Bi::Array{Complex})

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

"""
    use_GD_core(r::Real, μ::Real)

If the core is fluid then GD core will be used, otherwise the HH core
will be used for the propagation method given the initial radial and
shear modulus values from `data`.

NOTE: more research needs to be done here for there isn't complete
understanding why one core is better than the other.
"""
function use_GD_core(r::Real, μ::Real)

    if r > 0.0 && μ == 0
        flag = true
    else
        flag = false
    end

    return flag

end

"""
    normalize!(mass::Real, r::Real, data::Matrix{Complex})

Calculates the normalization units followed by normalizing the `data`
matrix given the `mass` and radius `r` of the planet.

NOTE: This should be done first before inserting the `data` matrix
into the `propagator_method`.
"""
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
