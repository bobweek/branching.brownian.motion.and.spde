using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	DifferentialEquations, StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bob/Research/Branching Brownian Motion/sde_functions.jl")

# weak competition
wc = 1e-5

# moderate competition
mc = 5e-5

# defines 'burn-in' for locating deterministic equilibrium
# assumes convergence when δ < ε
S = 50
T = 1.0e4
tspan = (0.0,T)
ε = 1.0e-3
ΔT = 100.0 # must be less than T

# pick some parameters
w = fill(0.1,S)        # niche breadths
U = fill(1.0,S)       # total niche use
Ω = sum(U)							# niche use scaling
η = fill(1.0,S)       # segregation variances
c = fill(wc,S)/S      # strengths of competition
a = fill(1e-4,S)       # strengths of abiotic selection
μ = fill(1e-3,S)      # mutation rates
V = fill(1.0,S)       # magnitudes of drift
R = fill(0.01,S)    	# innate rate of growth
θ = rand(Normal(0.0,1.0),S)             # phenotypic optima
pars = cat(S, w, U, η, c, a, μ, V, R, θ, dims=1)

pars = ModelParameters(S=S, w=w, U=U, η=η, c=c, a=a, μ=μ, V=V, R=R, θ=θ, Ω=Ω)

#
# find deterministic equilibrium
#

# initial condition
u₀ = cat(θ,fill(10.0, S),fill(1000.0, S),dims=1)

# numerically solve SDE
tspans = (0.0,1e3)
probhat = SDEProblem(f,g,u₀,tspans,pars)
solhat = solve(probhat,maxiters=1e8)

plot(solhat)

# state at time T
endT = maximum(solhat.t)
x̃ = sol(endT)

# distritbution of mean traits
histogram(x̃[(0*S+1):(1*S)])

# distritbution of trait variances
histogram(log.(x̃[(1*S+1):(2*S)]))

# distritbution of abundance
histogram(log.(x̃[(2*S+1):(3*S)]))

# extract state variables
x̄_sol = copy(solhat[(0*S+1):(1*S),:])
G_sol = copy(solhat[(1*S+1):(2*S),:])
N_sol = copy(solhat[(2*S+1):(3*S),:])

# build dataframes
spp = string("Species", 1)
df = DataFrame(spp = spp, x = x̄_sol[1,:], G = G_sol[1,:],
	N = N_sol[1,:], time = solhat.t)
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i,:], G = G_sol[i,:],
		N = N_sol[i,:], time = solhat.t)
	global df = append!(df,df_dummy)
end

# export to csv for ggplot
CSV.write("/home/bob/Research/WNEE/sample_path_wc.csv", df)

#
# moderate competition
#

c = fill(mc,S)/S      # strengths of competition
pars = cat(S, w, U, η, c, a, μ, v, r, θ, dims=1)
probhat = SDEProblem(f,g,u₀,tspans,pars)
solhat = solve(probhat,maxiters=1e8)

# extract state variables
x̄_sol = copy(solhat[(0*S+1):(1*S),:])
G_sol = copy(solhat[(1*S+1):(2*S),:])
N_sol = copy(solhat[(2*S+1):(3*S),:])

# build dataframes
spp = string("Species", 1)
df = DataFrame(spp = spp, x = x̄_sol[1,:], G = G_sol[1,:],
	N = N_sol[1,:], time = solhat.t)
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i,:], G = G_sol[i,:],
		N = N_sol[i,:], time = solhat.t)
	global df = append!(df,df_dummy)
end

# export to csv for ggplot
CSV.write("/home/bob/Research/WNEE/sample_path_mc.csv", df)
