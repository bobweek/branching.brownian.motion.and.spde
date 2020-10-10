#
# Generates sample paths to compare with ind-based simulation
#

using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	DifferentialEquations, StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bb/Gits/branching.brownian.motion.and.spde/sde_functions.jl")

j=1
k=1
# strengths of competition
cc = [1e-7,5e-6]

# for S = 1000
cc = [1e-8,5e-7]


# durations
TT = [1e3,5e3]

tspan = (0.0,TT[j])

# background parameters
S = 1000
w = fill(0.1, S)  # niche breadths
U = fill(1.0, S)  # total niche use
c = fill(cc[k],S)      # strengths of competition
Ω = sum(U) # niche use scaling
η = fill(1.0, S)  # segregation variances
μ = fill(1e-7, S) # mutation rates
V = fill(5.0, S)  # magnitudes of drift
R = fill(1.0, S) # innate rate of growth
a = fill(1e-2,S)       # strengths of abiotic selection
θ = fill(0.0, S)  # phenotypic optima

pars = ModelParameters(S=S, w=w, U=U, η=η, c=c, a=a, μ=μ, V=V, R=R, θ=θ, Ω=Ω)

#
# find deterministic equilibrium
#

# approximation for equilibrium abundance
C = c.*.√( U.^2 ./ .√(4*π.*w) )
N₀ = Int64.(floor.( (R.-0.5*.√(μ.*a))./C ) )

# initial condition
u₀ = cat(θ,fill(10.0, S),fill(1000.0, S),dims=1)

# numerically solve SDE
prob = SDEProblem(f,g,u₀,tspan,pars)
sol = solve(prob)

# extract state variables
x̄_sol = copy(.-sol[(0*S+1):(1*S),:])
G_sol = copy(sol[(1*S+1):(2*S),:])
N_sol = copy(sol[(2*S+1):(3*S),:])

# build dataframes
spp = string("Species", 1)
df = DataFrame(spp = spp, x = x̄_sol[1,:], G = G_sol[1,:], N = N_sol[1,:], time = sol.t)
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i,:], G = G_sol[i,:],
	N = N_sol[i,:], time = sol.t)
	global df = append!(df,df_dummy)
end

# weak competition
#CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_wc.csv", df)

# moderate competition
#CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_mc.csv", df)

# for S = 1000

LT = length(sol)

# extract state variables
x̄_sol = copy(.-sol[(0*S+1):(1*S),LT])
G_sol = copy(sol[(1*S+1):(2*S),LT])
N_sol = copy(sol[(2*S+1):(3*S),LT])

# build dataframes
spp = string("Species", 1)
df = DataFrame(spp = spp, x = x̄_sol[1], G = G_sol[1], N = N_sol[1])
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i], G = G_sol[i],
	N = N_sol[i])
	global df = append!(df,df_dummy)
end


# weak competition
CSV.write("/home/bb/Gits/branching.brownian.motion.and.spde/sample_path_wc1e3.csv", df)

# moderate competition
#CSV.write("/home/bb/Gits/branching.brownian.motion.and.spde/sample_path_mc1e3.csv", df)
