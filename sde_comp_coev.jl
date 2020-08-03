#
# using numerical solutions to estimate
# Cov(α,|β|), Cov(α,γ) Cov(α,ℭ) and
# Cor(α,|β|), Cor(α,γ) Cor(α,ℭ)
#


using Parameters,
      Statistics,
      Random,
      LinearAlgebra,
      Distributions,
      DifferentialEquations,
      StatsBase,
      StatsPlots,
      Plots,
      DataFrames,
      CSV,
      Optim

#include("/home/bb/Gits/branching.brownian.motion.and.spde/sde_functions.jl")

# background parameters
S = 100
w = fill(0.1, S)  # niche breadths
U = fill(1.0, S)  # total niche use
Ω = sum(U) # niche use scaling
η = fill(1.0, S)  # segregation variances
μ = fill(1e-5, S) # mutation rates
V = fill(1.0, S)  # magnitudes of drift
R = fill(1.0, S) # innate rate of growth
θ = fill(0.0, S)  # phenotypic optima

# set the timespan
# T₁ corresponds to the 'burn-in'
# T₂ corresponds the region we estimate the covariances over
T₁ = 5e2
T₂ = 5e2
tspan = (0.0, T₁ + T₂)

# containers for covariances and correlations
Cov_αβ_mean = zeros(0)
Cov_αabsβ_mean = zeros(0)
Cov_αγ_mean = zeros(0)
Cov_αℭ_mean = zeros(0)
Cor_αβ_mean = zeros(0)
Cor_αabsβ_mean = zeros(0)
Cor_αγ_mean = zeros(0)
Cor_αℭ_mean = zeros(0)

# containers for variance in mean traits among species
# and mean of the variance within species
Vₓ_mean = zeros(0)
σ²_mean = zeros(0)

# counts number of points accepted
num = 0

# accumulate 1000 points
while num < 1000

    # random draws for a and c
    adraw = exp.(rand(Normal(-10,6),1))[1]
    cdraw = exp.(rand(Normal(-20,4),1))[1]

    # continue to draw a until it is less than 1e-2
    # this prevents small N
    while maximum(adraw)>1e-2
        adraw = exp.(rand(Normal(-10,6),1))[1]
    end

    # continue to draw c until it is less than 1e-4
    # this prevents astronomical Vₓ and small N
    while maximum(cdraw)>1e-4
        cdraw = exp.(rand(Normal(-20,4),1))[1]
    end

    a = fill(adraw,S) # strength of abiotic selection
    c = fill(cdraw,S) # strengths of competition

    pars = ModelParameters(
        S = S,
        w = w,
        U = U,
        η = η,
        c = c,
        a = a,
        μ = μ,
        V = V,
        R = R,
        θ = θ,
        Ω = Ω,
    )

    # initial condition
    u₀ = cat(θ, fill(10.0, S), fill(1000.0, S), dims = 1)

    # numerically solve SDE
    prob = SDEProblem(f, g, u₀, tspan, pars)
    dumm = false
    sol = try solve(prob,maxiters=1e8)
    catch y
        if isa(y,DomainError)
            dumm = true
        end
    end

    if !dumm

        index_T₁ = findmin(abs.(sol.t .- T₁))[2]
        index_T₂ = length(sol.t)

        Cov_αβ = zeros(index_T₂ - index_T₁ + 1)
        Cov_αabsβ = zeros(index_T₂ - index_T₁ + 1)
        Cov_αγ = zeros(index_T₂ - index_T₁ + 1)
        Cov_αℭ = zeros(index_T₂ - index_T₁ + 1)

        Cor_αβ = zeros(index_T₂ - index_T₁ + 1)
        Cor_αabsβ = zeros(index_T₂ - index_T₁ + 1)
        Cor_αγ = zeros(index_T₂ - index_T₁ + 1)
        Cor_αℭ = zeros(index_T₂ - index_T₁ + 1)

        Vₓ = zeros(index_T₂ - index_T₁ + 1)
        σ² = zeros(index_T₂ - index_T₁ + 1)

        for j = index_T₁:index_T₂
            t = sol.t[j]

            k = j - index_T₁ + 1

            α = alpha(sol(t), pars)
            β = beta(sol(t), pars)
            γ = gamma(sol(t), pars)
            ℭ = coevolution(sol(t), pars)

            ᾱ = mean(α)
            β̄ = mean(β)
            γ̄ = mean(γ)
            ℭ̄ = mean(ℭ)

            Cov_αβ[k] = cov(vec(α), vec(β))
            Cov_αabsβ[k] = cov(vec(α), vec(abs.(β)))
            Cov_αγ[k] = cov(vec(α), vec(γ))
            Cov_αℭ[k] = cov(vec(α), vec(ℭ))

            Cor_αβ[k] = cor(vec(α), vec(β))
            Cor_αabsβ[k] = cor(vec(α), vec(abs.(β)))
            Cor_αγ[k] = cor(vec(α), vec(γ))
            Cor_αℭ[k] = cor(vec(α), vec(ℭ))

            Vₓ[k] = var(sol(t)[1:S])
            σ²[k] = mean(sol(t)[(S+1):(2*S)] .+ pars.η)

        end

        # append results to containers
        append!(Cov_αβ_mean,   mean(Cov_αβ))
        append!(Cov_αabsβ_mean,mean(Cov_αabsβ))
        append!(Cov_αγ_mean,   mean(Cov_αγ))
        append!(Cov_αℭ_mean,   mean(Cov_αℭ))

        append!(Cor_αβ_mean,   mean(Cor_αβ))
        append!(Cor_αabsβ_mean,mean(Cor_αabsβ))
        append!(Cor_αγ_mean,   mean(Cor_αγ))
        append!(Cor_αℭ_mean,   mean(Cor_αℭ))

        append!(Vₓ_mean,mean(Vₓ))
        append!(σ²_mean,mean(σ²))

        global num += 1

    end

end

# check it out
scatter(Vₓ_mean./σ²_mean,Cor_αβ_mean)
scatter(Vₓ_mean./σ²_mean,Cor_αabsβ_mean)
scatter(Vₓ_mean./σ²_mean,Cor_αγ_mean)
scatter(Vₓ_mean./σ²_mean,Cor_αℭ_mean)

# build dataframe
df = DataFrame( Cαβ = Cor_αβ_mean, Cαabsβ = Cor_αabsβ_mean,
	Cαγ = Cor_αγ_mean, Cαℭ = Cor_αℭ_mean, V = Vₓ_mean, σ = σ²_mean)

# export to csv for ggplot
CSV.write("/home/bb/Gits/branching.brownian.motion.and.spde/corrs.csv", df)
