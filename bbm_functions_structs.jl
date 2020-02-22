# data type that holds population and parameters

@with_kw mutable struct population
    x::Vector{Float64}  # trait values
    N::Int64            # poulation size
    x̄::Float64          # mean trait
    σ²::Float64         # phenotypic variance
    r::Float64          # intrinsic rate of growth
    a::Float64          # strength of abiotic selection
    θ::Float64          # abiotic optimum
    c::Float64          # strength of competition
    μ::Float64          # rate of diffusion (mutation)
    V::Float64          # variance in reproductive output
end

function update(X)

    @unpack x, N, x̄, σ², r, a, θ, c, μ, V = X

    W = fill(0,N)

    for j in 1:N

        # mean fitness of individual j
        w = exp( r - a*(θ-x[j])^2/2.0 - c*N )

        # parameterizing the NegativeBinomial
        q = w/V
        s = w^2/(V-w)

        # draw random number of offspring
        W[j] = rand( NegativeBinomial( s, q ), 1)[1]

    end

    # total number of offspring
    Nₚ = sum(W)

    # container for locations of offspring
    xₚ = fill(0.0,Nₚ)

    # keeps track of which individual is being born
    ct = 0

    # loop throug parents
    for j in 1:N

        # birth each offspring
        for k in 1:W[j]

            # consider next individual
            ct += 1

            # draw random trait for this individual
            xₚ[ct] = rand( Normal( x[j], √μ ), 1)[1]

        end

    end

    x̄ₚ = mean(xₚ)
    σₚ²= var(xₚ)

    Xₚ = population(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,r=r,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end
