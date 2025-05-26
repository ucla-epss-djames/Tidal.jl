using PhysicalConstants.CODATA2014: G
using PlanetEvolution
using Planets
using Tidal

using Plots
using Printf

# Volume of a Sphere
V(r::Real) = 4/3 * π * r^3

# helled+guillot2018
function satn_param()
    name = "2-layer Saturn"
    layers = 500
    R = 58232e3
    M = 568.319e24
    GM = M * G.val
    omega = 2π / (10.656 * 3600)
    Cp, alpha, k, kappa = 0.0, 0.0, 0.0, 0.0
    T0, P0, a, b, B, nabla = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    Teq, Tef, T1, eta, A, Ra = 0.0, 0.0, 0.0, 0.0, 0, 0
    mu_f = (0, 0, 0)
    model = 1
    plnt = Planet(name, layers, R, M, GM, omega, Cp, alpha, k, kappa, T0, P0, a,
                  b, B, nabla, Teq, Tef, eta, T1, A, Ra, mu_f, model)
    mn = Moon("No Moon", 0, 0, 0, 0, false)

    return (plnt, mn)
end

function fig_A2(save_data=false)

    function k2_remus(R_p, R_c, ρ_m, ρ_c, g_c, μ_c)
        α = eq13_α(R_p, R_c, ρ_m, ρ_c)
        β = eq13_β(R_p, R_c, α)
        ϵ = eq17_ϵ(R_c, ρ_m, ρ_c, g_c, μ_c, α, β)
        k2 = eq18_k2_p(α, β, ϵ)
        return k2
    end

    plnt, mn = satn_param()
    layers_c = floor(Int, .2*plnt.layers)
    layers_m = plnt.layers - layers_c
    rho = plnt.M / V(plnt.R)
    rho_c = 14500
    mu = 1000e9
    eta = 1e15

    data = zeros(111, 9)
    it = 1
    for i in 50:160 # 50:160 # i - core space
        if(i > 100) rho_c = 5500 end
        satn = zeros(plnt.layers, 4)
        R_c = plnt.R * (i / plnt.layers)
        satn[1:layers_c,1] = range(0.1, R_c, length=layers_c)
        dr = satn[2,1] - satn[1,1]
        satn[layers_c+1:plnt.layers,1] = range(R_c + dr, plnt.R,
                                               length=layers_m)
        mass_c = rho_c * V(R_c)
        mass_m = plnt.M - mass_c
        rho_m = mass_m / (V(plnt.R) - V(R_c))
        satn[1:layers_c,2] .= rho_c
        satn[layers_c+1:end,2] .= rho_m
        satn[:,3] .= mu
        satn[1:layers_c,4] .= eta
        sd, mass = planet_structure(plnt, mn, satn)
        println(mass)

        cmu = sd[layers_c,2]
        g_c = real(sd[layers_c,3])

        normalize!(mass, real(sd[end,1]), sd)
        tidal = propagator_method(2, plnt.layers, sd, false)
        k2r = k2_remus(plnt.R, R_c, rho_m, rho_c, g_c, cmu)
        data[it,1] = R_c
        data[it,2] = real(tidal[plnt.layers,1])
        data[it,3] = imag(tidal[plnt.layers,1])
        data[it,4] = rho_c
        data[it,5] = real(cmu)
        data[it,6] = imag(cmu)
        data[it,7] = g_c
        data[it,8] = real(k2r)
        data[it,9] = imag(k2r)
        it += 1
    end

    # SAVING DATA
    if(save_data == true)
        open("lainey+2017_data.txt", "w") do io
            for i in 1:111
                @printf(io, "%12.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
                        data[i,1], data[i,2], data[i,3], data[i,4], data[i,5],
                        data[i,6], data[i,7], data[i,8], data[i,9])
            end
        end
    end

    plt_1 = plot(xaxis=("Core Radius [x 1000 km]", (7.5, 16.5), 8:2:16),
                 yaxis=("k₂", (0.9, 1.4), 1.0:0.1:1.4), minorgrid=true,
                 legend=:none)
    c = data[:,1] ./ 1000000
    k2 = data[:,8]
    plt_1 = scatter!(c, k2, markershape=:square, color=:orange)
    k2 = data[:,2]
    plt_1 = scatter!(c, k2, markershape=:+, color=:black)

    plt_2 = plot(xaxis=("Core Radius [x 1000 km]", (7.5, 16.5), 8:2:16),
                 yaxis=("Q", :log10), minorgrid=true, legend=:none)
    Q = data[:,8] ./ abs.(data[:,9])
    plt_2 = scatter!(c, Q, markershape=:square, color=:orange)
    Q = data[:,2] ./ abs.(data[:,3])
    plt_2 = scatter!(c, Q, markershape=:+, color=:black)

    plt = plot(plt_1, plt_2, layout=(1, 2), framestyle=:box, size=(1080, 566), dpi=600,
               left_margin=5Plots.mm, bottom_margin=5Plots.mm)

end
