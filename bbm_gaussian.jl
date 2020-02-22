using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim, KernelDensity,
	QuadGK, HypothesisTests

include("/home/bob/Research/Branching Brownian Motion/bbm_functions_structs.jl")

#
# demonstrating approximate Gaussian distribution for small V/N
#
# corresponds
#

# number of generations to halt at
T = 200

# initial trait values
x₀ = rand(Normal(0.0,1.0),n)

# parameter values
r = 2.0
a = 0.01
θ = 0.0
c = 0.0001
μ = 1.0
V = 20.0

# equilibrium abundance
# we use this as the initial abundance
N̂ = floor((r-.5*√(μ*a))/c)

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
X = population(x=x₀, x̄=mean(x₀), σ²=var(x₀), N=N̂,
	r=r, a=a, θ=θ, c=c, μ=μ, V=V)

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

# rescale time
Tₛ = fill(1/N₀,T)
for i in 2:T
	Tₛ[i] = Tₛ[i-1]+1/n
end

# plotting the rescaled process

plot(Tₛ,Nₕ)
plot(Tₛ,x̄ₕ)
plot(Tₛ,σ²ₕ)

# our approximate
phen_dist = kde(Xₕ[T].x)

# expected
x̄_exp = θ
σ²_exp = √(μ/a)

# integral of kernel density estimate
one = quadgk(x->pdf(phen_dist,x), -Inf, Inf, rtol=1e-3)[1]

# plotting the approximate versus expected
plot(range(-15,15),x->pdf(phen_dist,x)/one)
plot!(x->pdf(Normal(x̄_exp,√σ²_exp),x))

# define culmulative density functions (cdf's)
function phen_cdf(x,X)
	n = length(X)
	num = length(findall(X.<=x))
	return Float64(num) / Float64(n)
end

# plot the cdf's
plot(range(-15,15),x->phen_cdf(x,Xₕ[T].x))
plot!(x->cdf(Normal(x̄_exp,√σ²_exp),x))

# calculate kolmogorov statistic
function Kst_fct(x)
	abs( cdf(Normal(x̄_exp,√σ²_exp),x) - phen_cdf(x,Xₕ[T].x) )
end
findKst = maximize(Kst_fct,-20,20)
KS_stat = -findKst.res.minimum

cdf(KSOneSided(length(Xₕ[T].x)),KS_stat)

plot(range(0,1),x->cdf(KSOneSided(length(Xₕ[T].x)),x))

function KSONEcdf(x)
	return cdf(KSOneSided(length(Xₕ[T].x)),x)
end

plot(range(0,.,length=1000),y->KSONEcdf(y))
