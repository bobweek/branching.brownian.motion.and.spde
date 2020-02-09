using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	DifferentialEquations, StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bob/Research/White\ Noise/bbm_functions_structs.jl")

#
# a branching Brownian motion
#

# number of generations to halt at
T = 1000

# initial number of individuals
N₀ = 100

# initial trait values
x₀=rand(Normal(0.0,1.0),N₀)

# set up initial population
X = population(x=x₀, x̄=mean(x₀), σ²=var(x₀), N=N₀,
	r=1.0, a=0.01, θ=0.0, c=0.01, μ=2.0, V=exp(1.0))

##
## VERY IMPORTANT REQUIREMENT  -->   V >= exp(r)
##

# always a good idea to inspect a single iteration
update(X)

# set up history of population
Xₕ = fill(X,T)

# set up containers for paths of N, x\bar and σ²
Nₕ = fill(X.N,T)
x̄ₕ = fill(X.x̄,T)
σ²ₕ = fill(X.σ²,T)

# simulate
for i in 2:T

	if Xₕ[i-1].N > 0

		Xₕ[i] = update(Xₕ[i-1])

	else

		Xₕ[i] = Xₕ[i-1]

	end

	Nₕ[i] = Xₕ[i].N
	x̄ₕ[i] = Xₕ[i].x̄
	σ²ₕ[i] = Xₕ[i].σ²

end

# plot some results

plot(1:T,Nₕ)

plot(1:T,x̄ₕ)

plot(1:T,σ²ₕ)

# rescale time
Tₛ = fill(1/N₀,T)
for i in 2:T
	Tₛ[i] = Tₛ[i-1]+1/Nₕ[i]
end

# plotting the rescaled process
plot(Tₛ,x̄ₕ./.√Nₕ,ylims=(-1.0,1.0))
plot(Tₛ,σ²ₕ./.√Nₕ)
