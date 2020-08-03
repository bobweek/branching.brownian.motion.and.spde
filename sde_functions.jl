###############################
#                             #
# data and parameter structs  #
#                             #
###############################


# data type that holds data on the community
@with_kw mutable struct CommData
	x̄::Vector{Float64}  # mean traits
	G::Vector{Float64}  # additive genetic variances
	N::Vector{Float64}  # abundances
end

# data type that holds model parameters
@with_kw mutable struct ModelParameters
	S::Int64            # species richness
	w::Vector{Float64}  # niche breadths, always zero for now...
	U::Vector{Float64}  # total niche use
	η::Vector{Float64}  # segregation variances
	c::Vector{Float64}	# strengths of competition
	a::Vector{Float64}	# strengths of abiotic selection
	μ::Vector{Float64}	# mutation rates
	V::Vector{Float64}	# variance of reproductive output
	R::Vector{Float64}	# intrinsic rates of growth
	θ::Vector{Float64}	# phenotypic optima
	Ω::Float64        	# total niche-use scaling
end

# data type that holds known parameters
@with_kw mutable struct FixedParameters
	S::Int64            # species richness
	w::Vector{Float64}  # niche breadths, always zero for now...
	U::Vector{Float64}  # total niche use
	x̃::Vector{Float64}  # initial condition (≈ deterministic equilibrium)
	ε::Float64          # convergence threshold
	T::Float64          # duration of time to integrate over
	dT::Float64					# duration of time to find deterministic equilibrium
	ΔT::Float64					# length of time to check convergence over
	bT::Float64					# duration of burnin time (=0 removes burnin)
	bins::Int64					# number of bins to use for histograms
end

# data type that holds parameters we need to estimate
@with_kw mutable struct LooseParameters
	η::Float64         # segregation variances
	c::Float64         # strengths of competition
	a::Float64         # strengths of abiotic selection
	μ::Float64         # mutation rates
	V::Float64         # variance of reproductive output
	R::Float64         # intrinsic rates of growth
	θ::Float64         # phenotypic optima
	Ω::Float64         # total niche-use scaling
end


###############################
#                             #
# drift & diffusion functions #
#                             #
###############################

#
# untransformed
#

# deterministic drift  (not to be confused with random genetic drift)
function f(u,p,t)

  # unpack model parameters
  @unpack S, w, U, η, c, a, μ, V, R, θ = p

  # unpack model variables
  x = u[(0*S+1):(1*S)] # mean traits
  G = u[(1*S+1):(2*S)] # additive genetic variances
  N = u[(2*S+1):(3*S)] # abundances

  # initialize containers for differentials
  dx = zeros(S)
  dG = zeros(S)
  dN = zeros(S)

  b = Array{Float64,2}(undef,S,S)
  for i = 1:S
    for j = 1:S
      b[i,j] = 1 / ( w[i] + w[j] + G[i] + G[j] + η[i] + η[j] )
    end
  end

  for i = 1:S

    B = 0.0
    for j = 1:S
      B += N[j] * b[i,j] * (x[j]-x[i]) * √( b[i,j] / (2.0*π) ) *
	  	exp( -b[i,j] * ( x[j]-x[i] )^2 / 2.0 )
    end

    dx[i] =  a[i] * G[i] * ( θ[i]-x[i] ) - c[i] * G[i] * B

    B = 0.0
    for j = 1:S
      B += N[j] * b[i,j] * (1.0 - b[i,j] * ( x[i]-x[j] )^2) *
	  	√( b[i,j] / (2.0*π) ) * exp( -b[i,j] * ( x[i]-x[j] )^2 / 2.0 )
    end

	B += N[i]*b[i,i]*√(b[i,i]/(2.0*π))

    dG[i] = μ[i] + ( c[i]*B - a[i] ) * G[i]^2 - V[i] * G[i] / N[i]

    B = 0.0
    for j = 1:S
      B += N[j] * √( b[i,j] / (2.0*π) ) *
	  	exp( -b[i,j] * ( x[j]-x[i] )^2 / 2.0 )
    end

    dN[i] = (R[i] - a[i]*(  ( θ[i] - x[i] )^2 + G[i] + η[i]  )/2.0 -
		c[i] * B)*N[i]

  end

  du = dx
  append!(du,dG)
  append!(du,dN)

  return(du)

end

# diffusion
function g(u,p,t)


  # unpack model parameters
  @unpack S, w, U, η, c, a, μ, V, R, θ = p

  # unpack model parameters
#  S::Int = p[1]
#  w      = p[(0*S+2):(1*S+1)]
#  U      = p[(1*S+2):(2*S+1)]
#  η      = p[(2*S+2):(3*S+1)]
#  c      = p[(3*S+2):(4*S+1)]
#  a      = p[(4*S+2):(5*S+1)]
#  μ      = p[(5*S+2):(6*S+1)]
#  v      = p[(6*S+2):(7*S+1)]
#  r      = p[(7*S+2):(8*S+1)]
#  θ      = p[(8*S+2):(9*S+1)]

  # unpack model variables
  x = u[(0*S+1):(1*S)] # mean traits
  G = u[(1*S+1):(2*S)] # additive genetic variances
  N = u[(2*S+1):(3*S)] # abundances

  # initialize containers for differentials
  dx = zeros(S)
  dG = zeros(S)
  dN = zeros(S)

  for i = 1:S

  	dN[i] = √( V[i] * N[i] )

    dx[i] = √( V[i] * G[i] / N[i] )

    dG[i] = G[i] * √( 2.0 * V[i] / N[i] )

  end

  du = dx
  append!(du,dG)
  append!(du,dN)

  return(du)

end

function n_ints(S)

  convert( Int64, ((S-1)*S/2) )

end

function rf_ρ(x,σ²,s,U,w)

  S = length(x)
  M = zeros(S,S)

  for i in 1:S
    for j in 1:S
	    b = 1/(w[i]+w[j]+σ²[i]+σ²[j])
      if i ≠ j
	      M[i,j] = U[i] * U[j] * s[i] * s[j] * sqrt(b/(2*π)) * exp(-b*(x[i]-x[j])^2/2)
      end
    end
  end

  return M

end

function vect_ρ(ρ,S)

  v = zeros(convert(Int64,S*(S-1)//2))
  k = 0

  for i in 1:S
    for j in 1:S
      if i>j
        k += 1
        v[k] = ρ[i,j]
      end
    end
  end

  return v

end

function pair2single(i,j,S)

    s = 0

    for k in 1:(S-1)
        for l in (k+1):S
            s += 1
            if i==k && j==l
                return s
            end
        end
    end

end

function single2pair(k,S)

    for i in 1:(S-1)
        for j in (i+1):S
            if pair2single(i,j,S)==k
                return [i,j]
            end
        end
    end

end

function Dₖₗ(ρ₁,ρ₂,S)

  D = 0

  for i in 1:n_ints(S)
    if ρ₁[i]>0 && ρ₂[i]>0
      D += ρ₁[i]*log(ρ₁[i]/ρ₂[i])
    end
  end

  return D

end

function p_norm(x,p)
	norm = (sum(abs.(x).^p))^(1/p)
	return norm
end

function in_box(x,y,b)
	return prod( x.-y .< b )
end

function findEcolParsₚ(qₚ,Δ,p)
	#
	S = convert(Int64, length(qₚ)/2)
	#
	w  = qₚ[(0+1):(1*S)]
  	U  = qₚ[(S+1):(2*S)]
	#
	x̄  = Δ[1][:,1]
	σ² = Δ[1][:,2]
	s  = Δ[1][:,3]
	#
	ρₒ = Δ[2]
	ρₚ = rf_ρ(x̄,σ²,s,U,w)
	ρₚ = LinearAlgebra.normalize( vect_ρ(ρₚ,S), 1)
	#
	return p_norm(ρₚ.-ρₒ,n_ints(S),p)
	#
end

function findEcolPars(qₚ,Δ)
  #
  # qₚ is a vector containing the proposed ecol pars
  # qₚ = cat( w, U, dims=1)
  #
  # Δ is a tuple containing observed data
  # Δ[1][:,1] = [ x[1], ...,x[S]  ]
  # Δ[1][:,2] = [ σ²[1],...,σ²[S] ]
  # Δ[1][:,3] = [ N[1], ...,N[S]  ]
  # Δ[2] = ρₒ
  #
  S = convert(Int64, length(qₚ)/2)
  #
  w  = qₚ[(0+1):(1*S)]
  U  = qₚ[(S+1):(2*S)]
  #
  x̄  = Δ[1][:,1]
  σ² = Δ[1][:,2]
  s  = Δ[1][:,3]
  #
  ρₒ = Δ[2]
  ρₚ = rf_ρ(x̄,σ²,s,U,w)
  ρₚ = LinearAlgebra.normalize( vect_ρ(ρₚ,S), 1 )
  #
  return Dₖₗ(ρₒ,ρₚ,S)
  #
end

function equilibrium(Δ,p)

	@unpack x̄, G, N = Δ
	@unpack S, w, c, η, r, a, μ, v = p
	val1 = 0.0
	val2 = 0.0

	for i in 1:S

		b = 1/(2.0*w[i]+2.0*G[i]+2.0*η[i]*N[i])
		val1 += r[i] + c[i]*sqrt(b/(2.0*π)) - a[i]*G[i]/2.0
		val2 += μ[i] - a[i]*G[i]^2 -v[i]^2*G[i]/N[i]
		subval1 = 0.0
		subval2 = 0.0

		for j in 1:S

			b = 1 / ( w[i] + G[i] + η[i] + w[j] + G[j] + η[j] )
			subval1 += N[j]*sqrt(b/(2*π))
 			subval2 += N[j]*b*sqrt(b/(2*π))
		end

		val1 -= c[i]*subval1
		val2 += c[i]*G[i]*subval2

	end

	val = abs(val1) + abs(val2)

	return val

end


function alpha(sol,pars)

	@unpack S, w, U, η, c, a, μ, V, R, θ, Ω = pars

	x̄ = sol[(0*S+1):(1*S)]
	G = sol[(1*S+1):(2*S)]
	N = sol[(2*S+1):(3*S)]

	α = zeros(S,S)

	for i in 1:S
		for j in 1:S
			b = 1 / ( w[i] + w[j] + G[i] + G[j] + η[i] + η[j] )
			α[i,j] = c[i] * √(b/(2*π)) * exp(-b*(x̄[i]-x̄[j])^2/2)
		end
	end

	return α

end

function beta(sol,pars)

	@unpack S, w, U, η, c, a, μ, V, R, θ, Ω = pars

	x̄ = sol[(0*S+1):(1*S)]
	G = sol[(1*S+1):(2*S)]
	N = sol[(2*S+1):(3*S)]

	β = zeros(S,S)

	for i in 1:S
		for j in 1:S
			b = 1 / ( w[i] + w[j] + G[i] + G[j] + η[i] + η[j] )
			β[i,j] = c[i] * N[j] * b * (x̄[i]-x̄[j]) * √(b/(2*π)) * exp(-b*(x̄[i]-x̄[j])^2/2)
		end
	end

	return β

end

function gamma(sol,pars)

	@unpack S, w, U, η, c, a, μ, V, R, θ, Ω = pars

	x̄ = sol[(0*S+1):(1*S)]
	G = sol[(1*S+1):(2*S)]
	N = sol[(2*S+1):(3*S)]

	γ = zeros(S,S)

	for i in 1:S
		for j in 1:S
			b = 1 / ( w[i] + w[j] + G[i] + G[j] + η[i] + η[j] )
			γ[i,j] = c[i] * N[j] * b * ( 1 - b*(x̄[i]-x̄[j])^2 ) * √(b/(2*π)) * exp(-b*(x̄[i]-x̄[j])^2/2)
		end
	end

	return γ

end

function coevolution(sol,pars)

	β = beta(sol,pars)
	γ = gamma(sol,pars)

	ℭ = zeros(pars.S,pars.S)

	for i in 1:S
		for j in 1:S
			ℭ[i,j] = sqrt(abs(β[i,j]*β[j,i])) + sqrt(abs(γ[i,j]*γ[j,i]))
		end
	end

	return ℭ
end
