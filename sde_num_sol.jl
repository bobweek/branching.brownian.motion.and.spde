using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	DifferentialEquations, StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bob/Research/Branching Brownian Motion/sde_functions.jl")

j=1
k=2
# strengths of competition
cc = [1e-5,1e-3]

# durations
TT = [1e3,5e3]

dumm = true

while dumm

	try for j in 1:2

		for k in 1:2

			tspan = (0.0,TT[k])

			# background parameters
			S = 100
			w = fill(0.1, S)  # niche breadths
			U = fill(1.0, S)  # total niche use
			c = fill(cc[i],S)/S      # strengths of competition
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
			sol = solve(prob,maxiters=1e8)

			plot(sol)

			# state at time T
			endT = maximum(sol.t)
			x̃ = sol(endT)

			# extract state variables
			x̄_sol = copy(sol[(0*S+1):(1*S),:])
			G_sol = copy(sol[(1*S+1):(2*S),:])
			N_sol = copy(sol[(2*S+1):(3*S),:])

			# build dataframes
			spp = string("Species", 1)
			global df = DataFrame(spp = spp, x = x̄_sol[1,:], G = G_sol[1,:],
			N = N_sol[1,:], time = sol.t)
			for i in 2:S
				spp = string("Species", i)
				df_dummy = DataFrame(spp = spp, x = x̄_sol[i,:], G = G_sol[i,:],
				N = N_sol[i,:], time = sol.t)
				global df = append!(df,df_dummy)
			end

			df.spp

			if j==1 & k==1
				# weak competition
				CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_wc.csv", df)
			end

			if j==2 & k==1
				# strong competition
				CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_mc.csv", df)
			end

			if j==1 & k==2
				# weak competition long duration
				CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_wc_l.csv", df)
			end

			if j==2 & k==2
				# strong competition long duration
				CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_mc_l.csv", df)
			end

		end

	end
	catch y
		if isa(y,DomainError)
			dumm = true
		end
	end

end

tspan = (0.0,TT[k])

# background parameters
S = 100
w = fill(0.1, S)  # niche breadths
U = fill(1.0, S)  # total niche use
c = fill(cc[i],S)/S      # strengths of competition
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
sol = solve(prob,maxiters=1e8)

plot(sol)

# state at time T
endT = maximum(sol.t)
x̃ = sol(endT)

# extract state variables
x̄_sol = copy(sol[(0*S+1):(1*S),:])
G_sol = copy(sol[(1*S+1):(2*S),:])
N_sol = copy(sol[(2*S+1):(3*S),:])

# build dataframes
spp = string("Species", 1)
global df = DataFrame(spp = spp, x = x̄_sol[1,:], G = G_sol[1,:],
N = N_sol[1,:], time = sol.t)
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i,:], G = G_sol[i,:],
	N = N_sol[i,:], time = sol.t)
	global df = append!(df,df_dummy)
end

df.spp
