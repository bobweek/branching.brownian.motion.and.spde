#
# Adapted from
#
# http://www.stochasticlifestyle.com/solving-systems-stochastic-pdes-using-gpus-julia/
#

using OrdinaryDiffEq, StochasticDiffEq, RecursiveArrayTools, LinearAlgebra, Plots, Parameters

const N = 100

const Mx = Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1])

# data type that holds model parameters
@with_kw mutable struct ModelParameters
    μ::Float64	# mutation rate
    R::Float64	# intrinsic rates of growth
    a::Float64	# strength of abiotic selection
    θ::Float64	# phenotypic optima
    c::Float64	# strength of competition
    V::Float64	# variance of reproductive output
    X::Vector{Float64}	# variance of reproductive output
end

# Define the constants for the PDE
θ = 30.0
ω = 100.0
c = 0.001
a = 0.001
μ = 0.1
V = 1
R = 1
X = Vector([i for i in 1:N])

pars = ModelParameters(μ=μ, R=R, a=a, θ=θ, c=c, V=V, X=X)

# Define the initial condition as normal arrays
u0 = 1000*exp.(-(X.-50).^2 / (2*ω))/√(2*π*ω)

# Define the discretized PDE as an ODE function
function f(u,p,t)
    @unpack μ, R, a, θ, c = p
    du = 0.5*μ*Mx*u + u.*(R .- 0.5*a*(X.-θ).^2 .- c*norm(u,1))
    return(du)
end

function g(u,p,t)
    @unpack V, X = p
    du = zeros(length(X))
    for i in 1:length(X)
        if u[i]>0
            du[i] = √(V*u[i])
        end
    end
    return(du)
end

#######################################################
### Solve the PDE
#######################################################

prob = ODEProblem(f,u0,(0.0,50.0),pars)
sol = solve(prob,ROCK2())

time = sol.t

abun = zeros(N,length(time))
for i in 1:length(time)
    for j in 1:N
        abun[j,i] = sol.u[i][j]
    end
end
contour(time,X,abun,fill = true)

#######################################################
### Solve the SPDE
#######################################################

prob2 = SDEProblem(f,g,u0,(0.0,50.0),pars)
sol = solve(prob2,SOSRI())

time = sol.t

abun = zeros(N,length(time))
for i in 1:length(time)
    for j in 1:N
        abun[j,i] = sol.u[i][j]
    end
end
contour(time,X,abun,fill = true)
