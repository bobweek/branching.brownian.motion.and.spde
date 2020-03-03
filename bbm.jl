using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bob/Research/Branching Brownian Motion/bbm_functions_structs.jl")

#
# a branching Brownian motion
#

# parameter values
S = 1
R = fill(1.0,S)
a = fill(0.01,S)
θ = fill(0.0,S)
c = fill(1.0e-4,S)
w = fill(0.1, S)  # niche breadths
U = fill(1.0, S)  # total niche use
η = fill(1.0, S)  # segregation variances
μ = fill(1.0e-5,S)
V = fill(5.0,S)

# equilibrium abundance an the absence of interspecific interactions
# we use this as the initial abundance

C = c.*.√( U.^2 ./ .√(4*π.*w) )

N₀ = Int64.(floor.( (R.-0.5*.√(μ.*a))./C ) )

# initial breeding values
g₀ = rand.(Normal(0.0,1.0),N₀)

# initial trait values
x₀ = fill(zeros(0),S)
for i in 1:S
	ηₘ = √η[i]*Matrix(I, N₀[i], N₀[i])
	x₀[i] = vec(rand(MvNormal(g₀[i],ηₘ),1))
end

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
var.(x₀)
# set up initial population
X = community(S=S, x=x₀, g=g₀, N=N₀, n=N₀, x̄=mean.(x₀), σ²=var.(x₀),
	G=var.(g₀), R=R, a=a, θ=θ, c=c, w=w, U=U, η=η, μ=μ, V=V)

# always a good idea to inspect a single iteration
rescaled_lower(X)

# number of generations to halt at
T = 1000

# set up history of population
Xₕ = fill(X,T)

# simulate
for i in 2:T

	if prod( log.( 1 .+ Xₕ[i-1].N ) ) > 0

		Xₕ[i] = rescaled_lower(Xₕ[i-1])

	else

		Xₕ[i] = Xₕ[i-1]

	end

end

# set up containers for paths of N, x̄ and σ²
Nₕ = zeros(S,T)
x̄ₕ = zeros(S,T)
σ²ₕ= zeros(S,T)
Gₕ = zeros(S,T)

# fill them in
for i in 1:S
	for j in 1:T
		Nₕ[i,j] =Xₕ[j].N[i]
		x̄ₕ[i,j] =Xₕ[j].x̄[i]
		σ²ₕ[i,j]=Xₕ[j].σ²[i]
		Gₕ[i,j] =Xₕ[j].G[i]
	end
end

plot((1:T)./N₀[1],Nₕ[1,:])
plot((1:T)./N₀[1],x̄ₕ[1,:])
plot((1:T)./N₀[1],σ²ₕ[1,:])

plot(1:T,x̄ₕ,title="Mean Trait: N ~ 195000",ylabel="Value",xlabel="Iteration",ylim=(-1,1))
png("mean2.png")
plot(1:T,σ²ₕ,title="Trait Variance: N ~ 195000",ylabel="Value",xlabel="Iteration",ylim=(0,15))
png("var2.png")

# rescale time
Tₛ = fill(1/N₀,T)
for i in 2:T
	Tₛ[i] = Tₛ[i-1]+1/n
end

# plotting the rescaled process

plot(Tₛ,Nₕ)
plot(Tₛ,x̄ₕ)
plot(Tₛ,σ²ₕ)
