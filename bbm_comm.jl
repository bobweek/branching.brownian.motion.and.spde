using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bob/Research/Branching Brownian Motion/bbm_functions_structs.jl")

########################################################
#                                                      #
#  An individual-based model for the entire community  #
#                                                      #
#  Used to compare with diffusion approximation        #
#                                                      #
########################################################

# parameter values
S = 1
w = fill(0.1, S)  # niche breadths
U = fill(1.0, S)  # total niche use
c = fill(0.002,S) # strengths of competition
Ω = sum(U)        # niche use scaling
η = fill(1e-3, S) # environmental variances
μ = fill(1e-5, S) # mutation rates
V = fill(2.0, S)  # magnitudes of drift
R = fill(0.1, S)  # innate rate of growth
a = fill(1e-2,S)  # strengths of abiotic selection
θ = fill(0.0, S)  # phenotypic optima
n = fill(10.0, S) # scaling parameter

# equilibrium abundance an the absence of interspecific interactions
# we use this as the initial abundance

C = c.*(U.^2)./.√(4 .*π.*w)
N₀ = Int64.(floor.( (n.^2).*( (R./n).-0.5*( η.*a .+ .√(μ.*a) ) )./C ) )

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

# set up initial population
X = community(S=S, x=x₀, g=g₀, N=N₀, n=n, x̄=mean.(x₀), σ²=var.(x₀),
	G=var.(g₀), R=R, a=a, θ=θ, c=c, w=w, U=U, η=η, μ=μ, V=V)

# always a good idea to inspect a single iteration
rescaled_lower(X)

# number of generations to halt at
T = 100

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

# container for individuals
#individualₕ = zeros(2)

# fill them in
for i in 1:S
	for j in 1:T
		Nₕ[i,j] =Xₕ[j].N[i]
		x̄ₕ[i,j] =Xₕ[j].x̄[i]
		σ²ₕ[i,j]=Xₕ[j].σ²[i]
		Gₕ[i,j] =Xₕ[j].G[i]
	end
end

# rescaled time
resc_time = (1:T)./N₀[1]

# total number of individuals across entire simulation
total_inds = Int64(sum(Nₕ[1,:]))

# traits of each individual across entire simulation
inds = zeros(2,total_inds)

ind = 0
for i in 1:T
	for j in 1:Int64(Nₕ[1,i])

		global ind += 1
		inds[1,ind] = i
		inds[2,ind] = Xₕ[i].x[1][j]

	end
end

scatter(inds[1,:], inds[2,:], legend=false, ms=1)

plot(resc_time,Nₕ[1,:]./n[1])
plot(resc_time,x̄ₕ[1,:]./√(n[1]))
plot(resc_time,σ²ₕ[1,:]./n[1])


histogram(Xₕ[500].x[1])
