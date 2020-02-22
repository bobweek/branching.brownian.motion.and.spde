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

include("/home/bob/Research/Branching Brownian Motion/bbm_functions_structs.jl")

#
# shows moments are finite across a range of c
#

# number of generations to halt at
T = 200

# initial trait values
x₀ = rand(Normal(0.0, 1.0), n)

# parameter values
r = 2.0
a = 0.01
θ = 0.0
crange = [1e-5,1e-4,1e-3,1e-2]
μ = 1.0
V = 20.0

# how many times to repeat for each c?
reps = 10

# containers to store means and variances for each run
Means = zeros(length(crange), reps)
Varis = zeros(length(crange), reps)

for c in crange
    for k in 1:reps

        # equilibrium abundance
        # we use this as the initial abundance
        N̂ = Int64(floor((r - 0.5 * √(μ * a)) / c))

        # initial trait values
        x₀ = rand(Normal(0.0,1.0),N̂)

        ##
        ## VERY IMPORTANT REQUIREMENT   -->  V >= exp(r)
        ##
        ## this inequality must be satisfied to use negative binomial sampling
        ##
        ##
        ## TWO MORE IMPORTANT REQUIREMENTS --> 2*r > √(μ*a) && c > r - √(μ*a)/2
        ##
        ## these inequalities must be satisfied for positive equilibrium abundance
        ##

        # set up initial population
        X = population(
            x = x₀,
            x̄ = mean(x₀),
            σ² = var(x₀),
            N = N̂,
            r = r,
            a = a,
            θ = θ,
            c = c,
            μ = μ,
            V = V,
        )

        # set up history of population
        Xₕ = fill(X, T)

        # simulate
        for i = 2:T
            if Xₕ[i-1].N > 0

                Xₕ[i] = update(Xₕ[i-1])

            else

                Xₕ[i] = Xₕ[i-1]

            end

        end

        # store
        j = findall(μrange .== μ)[1]

        global Means[j, k] = Xₕ[T].x̄
        global Varis[j, k] = Xₕ[T].σ²

    end

end

# compute mean and variance of means and variances 😱
x̄_mean = zeros(length(μrange))
x̄_vari = zeros(length(μrange))
σ²_mean = zeros(length(μrange))
σ²_vari = zeros(length(μrange))

for i in 1:length(μrange)

    x̄_mean[i] = mean(Means[i,:])
    x̄_vari[i] =  var(Means[i,:])
    σ²_mean[i] = mean(Varis[i,:])
    σ²_vari[i] =  var(Varis[i,:])

end

boxplot(transpose(Means),ylims=(-1,1))

boxplot(transpose(Varis))
