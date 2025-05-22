# analytical forms of Q and k2 for a 2-layer planet
# equations come from Remus+2015
# variable names and definitions should follow the paper
# c-subscript is core
# p-subscript is planet
# ρ_0 is the density of the envelope
# these equations are benchmark equations for the tidal PMM code
# benchmarks can be found for lainey+2017-figA2, stixrude+2021-fig4

export eq13_α, eq13_β, eq20b_G_p, eq17_ϵ
export eq18_k2_p, eq18_k2_c, eq1_Q_p, eq1_Q_c, eq20a_Q_p, eq20a_Q_p

eq13_α(R_p, R_c, ρ_0, ρ_c) = 1 + 5/2*(ρ_c/ρ_0 - 1)*(R_c/R_p)^3
eq13_β(R_p, R_c, α) = 3/5*(R_c/R_p)^2 * (α - 1)
eq20b_G_p(R_p, R_c, α) = (α + 3/2) / (α + 3/2*(R_c/R_p)^5)

function eq17_ϵ(R_c, ρ_0, ρ_c, g_c, μ_c, α, β)
    μ_f = 19*μ_c / (2*ρ_c*g_c*R_c)
    num =  μ_f + ρ_0/ρ_c * (1 - ρ_0/ρ_c) * (β + 3/2) + (1 - ρ_0/ρ_c)
    dem = (α + 3/2) * (1 - ρ_0/ρ_c) * ρ_0/ρ_c
    ϵ = num / dem
    return ϵ
end

eq18_k2_p(α, β, ϵ) = 3/2 * (ϵ + 2/3*β) / (α*ϵ - β)
eq18_k2_c(ρ_0, ρ_c, α, β, ϵ) = 3/2 * (ϵ - 1 + ρ_c/ρ_0) / (α*ϵ - β)

eq1_Q_p(R_c, R_p, k2_p, k2_c, Q_c) = Q_c * k2_p/k2_c * (R_p/R_c)^5
eq1_Q_c(R_c, R_p, k2_p, k2_c, Q_p) = Q_p * k2_c/k2_p * (R_c/R_p)^5
eq20a_Q_p(R_c, R_p, G_p, k2_p, k2_c, Q_c) = eq1_Q_p(R_c, R_p, k2_p, k2_c, Q_c)/G_p
eq20a_Q_c(R_c, R_p, G_p, k2_p, k2_c, Q_p) = eq1_Q_c(R_c, R_p, k2_p, k2_c, Q_p)*G_p
