using PhysicalConstants.CODATA2014: G
using PlanetEvolution
using Planets
using Tidal

using Plots
using Printf

const Q_c = 1e4
const R_e = 6371e3
const M_e = 5.972e24

# Volume of a Sphere
V(r::Real) = 4/3 * π * r^3

function fluid_earth_param()
    name = "Fluid Super Earth"
    layers = 100
    R = 2.68 * R_e
    M = 6.55 * M_e
    GM = G.val * M
    omega, Cp, alpha, k, kappa = 0.0, 0.0, 0.0, 0.0, 0.0
    T0, P0, a, b, B, nabla = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    Teq, Tef, T1, eta, A, Ra = 0.0, 0.0, 0.0, 0.0, 0, 0
    mu_f = (0, 0, 0)
    model = 1
    plnt = Planet(name, layers, R, M, GM, omega, Cp, alpha, k, kappa, T0, P0, a,
                  b, B, nabla, Teq, Tef, T1, eta, A, Ra, mu_f, model)
    mn = Moon("No Moon", 0, 0, 0, 0, false)

    return plnt, mn
end

function visco_earth_param()

    name = "Viscoelastic Super Earth"
    layers = 100
    R = 2.68 * R_e
    M = 6.55 * M_e
    GM = G.val * M
    omega = 2π / (17.24 * 3600)
    Cp, alpha, k, kappa = 0.0, 0.0, 0.0, 0.0
    T0, P0, a, b, B, nabla = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    Teq, Tef, T1, eta, A, Ra = 0.0, 0.0, 0.0, 0.0, 0, 0
    mu_f = (0, 0, 0)
    model = 1
    plnt = Planet(name, layers, R, M, GM, omega, Cp, alpha, k, kappa, T0, P0, a,
                  b, B, nabla, Teq, Tef, eta, T1, A, Ra, mu_f, model)
    return plnt
end

function fig_4(save_data=false)

    function k2_remus(R_p, R_c, ρ_m, ρ_c, g_c, μ_c)
        α = eq13_α(R_p, R_c, ρ_m, ρ_c)
        β = eq13_β(R_p, R_c, α)
        ϵ = eq17_ϵ(R_c, ρ_m, ρ_c, g_c, μ_c, α, β)
        k2 = eq18_k2_p(α, β, ϵ)
        return k2
    end

    plnt, mn = fluid_earth_param()
    layers_c = floor(Int, .1*plnt.layers)
    layers_m = plnt.layers - layers_c
    rho_ratio = [0.73, 0.01]
    rho = plnt.M / V(plnt.R)

    d = Dict()
    supe = zeros(plnt.layers, 4)
    println("FLUID STRUCTURE")
    it = 1
    for i in rho_ratio
        rho_m = i * rho
        data = zeros(plnt.layers, 9)
        for j in 1:plnt.layers
            R_c = plnt.R * (j / plnt.layers)
            supe[1:layers_c,1] = range(0.1, R_c, length=layers_c)
            dr = supe[2,1] - supe[1,1]
            supe[layers_c+1:plnt.layers,1] = range(R_c + dr, plnt.R,
                                                   length=layers_m)
            mass_m = rho_m * (V(plnt.R) - V(R_c))
            rho_c = (plnt.M - mass_m) / V(R_c)

            supe[1:layers_c,2] .= rho_c
            supe[layers_c+1:end,2] .= rho_m
            sd, mass = planet_structure(plnt, mn, supe)

            cmu = sd[layers_c,2]
            a = real(sd[layers_c,3])

            normalize!(mass, real(sd[end,1]), sd)
            tidal = propagator_method(2, plnt.layers, sd, false)
            k2r = k2_remus(plnt.R, R_c, rho_m, rho_c, a, cmu)
            data[j,1] = real(sd[layers_c,1])
            data[j,2] = real(tidal[plnt.layers,1])
            data[j,3] = imag(tidal[plnt.layers,1])
            data[j,4] = rho_c
            data[j,5] = real(cmu)
            data[j,6] = imag(cmu)
            data[j,7] = a
            data[j,8] = real(k2r)
            data[j,9] = imag(k2r)
        end
        d[it] = (type="fluid", rho_m=rho_m, data=data)
        it += 1
    end

    println("VISCOELASTIC STRUCTURE")
    plnt = visco_earth_param()
    for i in rho_ratio
        rho_m = i * rho
        data = zeros(plnt.layers, 9)
        for j in 1:plnt.layers
            R_c = plnt.R * (j / plnt.layers)
            supe[1:layers_c,1] = range(0.1, R_c, length=layers_c)
            dr = supe[2,1] - supe[1,1]
            supe[layers_c+1:plnt.layers,1] = range(R_c + dr, plnt.R,
                                                   length=layers_m)

            mass_m = rho_m * (V(plnt.R) - V(R_c))
            rho_c = (plnt.M - mass_m) / V(R_c)

            supe[1:layers_c,2] .= rho_c
            supe[layers_c+1:end,2] .= rho_m

            a = planet_g(R_c, (plnt.M - mass_m)*G.val)
            mu = planet_mu(R_c, a, rho_c)
            tau_m = Q_c / plnt.ω
            eta = mu * tau_m

            supe[:,3] .= mu
            supe[1:layers_c,4] .= eta
            supe[layers_c+1:end,4] .= 0

            sd, mass = planet_structure(plnt, mn, supe)

            cmu = sd[layers_c,2]
            a = real(sd[layers_c,3])

            normalize!(mass, real(sd[end,1]), sd)
            tidal = propagator_method(2, plnt.layers, sd, false)
            k2r = k2_remus(plnt.R, R_c, rho_m, rho_c, a, cmu)
            data[j,1] = real(sd[layers_c,1])
            data[j,2] = real(tidal[plnt.layers,1])
            data[j,3] = imag(tidal[plnt.layers,1])
            data[j,4] = rho_c
            data[j,5] = real(cmu)
            data[j,6] = imag(cmu)
            data[j,7] = a
            data[j,8] = real(k2r)
            data[j,9] = imag(k2r)
        end
        d[it] = (type="viscoelastic", rho_m=rho_m, data=data)
        it += 1
    end

    # SAVING DATA
    if(save_data == true)
        open("stixrude+2021_data.txt", "w") do io
            for i in 1:4
                @printf(io, "%s rho_m = %12.6f\n", d[i].type, d[i].rho_m)
                for j in 1:plnt.layers
                    @printf(io, "%6.4f %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
                            d[i].data[j,1], d[i].data[j,2], d[i].data[j,3],
                            d[i].data[j,4], d[i].data[j,5], d[i].data[j,6],
                            d[i].data[j,7], d[i].data[j,8], d[i].data[j,9])
                end
            end
        end
    end

    axlabel = ["k₂", "Tidal Quality Q", "Core Radius c/R"]
    colors = [:red, :green, :blue]

    plt_1 = plot(yaxis=(axlabel[1], (0.0, 1.4)), legend=:none)
    it = 1
    for i in [1, 3]
        c = d[i].data[:,1]
        k2 = d[i].data[:,2]
        k2r = d[i].data[:,8]
        plt_1 = plot!(c, k2r, color=colors[it])
        plt_1 = scatter!(c, k2, markers=:+, color=colors[it])
        it += 1
    end

    plt_2 = plot(yaxis=(axlabel[1], (0.0, 1.4)), legend=:none)
    it = 1
    for i in [2, 4]
        c = d[i].data[:,1]
        k2 = d[i].data[:,2]
        k2r = d[i].data[:,8]
        plt_2 = plot!(c, k2r, color=colors[it])
        plt_2 = scatter!(c, k2, markers=:+, color=colors[it])
        it += 1
    end

    plt_3 = plot(xaxis=(axlabel[3], 0:0.2:1),
                 yaxis=(axlabel[2], (1e4, 1e13), :log10), legend=:none,
                 minorgrid=true)
    it = 3
    for i in [3, 4]
        c = d[i].data[:,1]
        k2r = d[i].data[:,2]
        k2i = d[i].data[:,3]
        Q = abs.(k2r ./ k2i)
        k2r = d[i].data[:,8]
        k2i = d[i].data[:,9]
        Qr = abs.(k2r ./ k2i)
        plt_3 = plot!(c, Qr, color=colors[it])
        plt_3 = scatter!(c, Q, markers=:+, color=colors[it])
        it -= 2
    end

    plt = plot(plt_1, plt_2, plt_3, layout=(3, 1), size=(1080, 1350), dpi=600,
               left_margin=5Plots.mm)
end
