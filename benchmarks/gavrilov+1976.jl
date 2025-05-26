# PlanetEvolution has a wrapper for Tidal
using PlanetEvolution
using Planets

using Plots
using Printf

# profiles
line_rho(ρ::Real, s::Real, s1::Real) = 4*ρ * (1 - (s/s1))
quad_rho(ρ::Real, s::Real, s1::Real) = 2.5*ρ * (1 - (s/s1)^2)
poly_rho(ρ::Real, s::Real, s1::Real) = π^2/3*ρ * sin(π*(s/s1)) / (π*(s/s1))

function jupiter_param()
    name = "Jupiter"
    layers = 1000
    R = 69.911e6
    M = 1.8987e27
    GM = 0.0
    omega, Cp, alpha, k, kappa = 0.0, 0.0, 0.0, 0.0, 0.0
    T0, P0, a, b, B,nabla = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    Teq, Tef, T1, eta, A, Ra = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    mu_f = (0, 0, 0)
    model = 1
    plnt = Planet(name, layers, R, M, GM, omega, Cp, alpha, k, kappa, T0, P0, a, b,
                  B,
                  nabla, Teq, Tef, T1, eta, A, Ra, mu_f, model)
    mn = Moon("No Moon", 0, 0, 0, 0, false)

    return plnt, mn
end

# set to true to save the file
function fig_1(save_data=false)

    plnt, mn = jupiter_param()
    ri = 0.1
    ρ_j = 1326
    jupt = zeros(plnt.layers, 4)
    ln = [2, 3, 4, 5, 6]
    n = 5
    rx = range(ri, plnt.R, length=plnt.layers)
    jupt[:,1] = rx

    # linear calc
    res_l = zeros(plnt.layers, n)
    rho_l = line_rho.(ρ_j, real.(jupt[:,1]), plnt.R)
    jupt[:,2] = rho_l
    for i in 1:n
        tidal = tidal_resp(plnt, mn, jupt, false, l=ln[i])
        res_l[:,i] = real.(tidal[:,1])
    end

    # quad calc
    res_q = zeros(plnt.layers, n)
    rho_q = quad_rho.(ρ_j, real.(jupt[:,1]), plnt.R)
    jupt[:,2] = rho_q
    for i in 1:n
        tidal = tidal_resp(plnt, mn, jupt, false, l=ln[i])
        res_q[:,i] = real.(tidal[:,1])
    end

    # poly calc
    res_p = zeros(plnt.layers, n)
    rho_p = poly_rho.(ρ_j, real.(jupt[:,1]), plnt.R)
    jupt[:,2] = rho_p
    for i in 1:n
        tidal = tidal_resp(plnt, mn, jupt, false, l=ln[i])
        res_p[:,i] = real.(tidal[:,1])
    end

    # SAVING DATA
    if(save_data == true)
        open("gavrilov+1976_data.txt", "w") do io
            for i in 1:plnt.layers
                @printf(io, "%12.6f %12.6f%12.6f%12.6f%12.6f%12.6f %12.6f%12.6f%12.6f%12.6f%12.6f %12.6f%12.6f%12.6f%12.6f%12.6f\n",
                        jupt[i,1]/plnt.R,
                        res_l[i,:]..., res_q[i,:]..., res_p[i,:]...)
            end
        end
    end

    plt_1 = plot(xticks=([0:0.2:1;], ["", "", "", "", "", ""]),
                 yaxis=("kₙ(x)", 0:0.1:0.7), legend=:none)
    for i in 1:n
        plt_1 = plot!(jupt[:,1] ./ plnt.R, res_l[:,i], color=:black)
        plt_1 = plot!(jupt[:,1] ./ plnt.R, res_q[:,i], color=:black,
                      linestyle=:dash)
        plt_1 = plot!(jupt[:,1] ./ plnt.R, res_p[:,i], color=:black,
                      linestyle=:dot)
    end
    plt_2 = plot(xaxis=("x", 0:0.2:1.0), yaxis=("Density, [g/cm³]", 1:5))
    plt_2 = plot!(jupt[:,1]./plnt.R, rho_l./1000, color=:black, label="linear")
    plt_2 = plot!(jupt[:,1]./plnt.R, rho_q./1000, color=:black,
                  linestyle=:dash, label="quadratic")
    plt_2 = plot!(jupt[:,1]./plnt.R, rho_p./1000, color=:black,
                  linestyle=:dot, label="polytrope")
    l = @layout [a; b{0.3h}]
    plt = plot(plt_1, plt_2, layout=l, size=(1080, 1350), left_margin=3Plots.mm,
               dpi=600)


end
