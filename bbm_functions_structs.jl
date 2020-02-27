# data type that holds population and parameters

@with_kw mutable struct population
    x::Vector{Float64}  # trait values
    N::Int64            # poulation size
    x̄::Float64          # mean trait
    σ²::Float64         # phenotypic variance
    R::Float64          # innate rate of growth
    a::Float64          # strength of abiotic selection
    θ::Float64          # abiotic optimum
    c::Float64          # strength of competition
    μ::Float64          # rate of diffusion (mutation)
    V::Float64          # variance in reproductive output
end

@with_kw mutable struct community
    S::Int64            # number of species
    x::Matrix{Float64}  # trait values
    N::Vector{Int64}    # poulation sizes
    x̄::Vector{Float64}  # mean traits
    σ²::Vector{Float64} # phenotypic variances
    R::Vector{Float64}  # innate rates of growth
    a::Vector{Float64}  # strengths of abiotic selection
    θ::Vector{Float64}  # abiotic optima
    c::Vector{Float64}  # strengths of competition
    w::Vector{Float64}  # individual niche widths
    U::Vector{Float64}  # total niche uses
    μ::Vector{Float64}  # rates of diffusion (mutation)
    V::Vector{Float64}  # variances in reproductive output
end

# update for community with complex competition term
function comm_update(X)

    @unpack S, x, N, x̄, σ², R, a, θ, c, w, U, μ, V = X

    for i in 1:S

        W = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            # this follows exactly from SM §5.6
            #

            # container for aggregating effects of competition
            B = 0.0

            # collect effects of competition with other individuals
            # within the same population
            except_j = 1:N[i] # somehow remove j
            for k in except_j
                B += U[i]^2*exp( (x[i,j] - x[i,k])^2 / (4*w[i]) )
                    / √(4*π*w[i])
            end

            # collect effects of competition with other individuals
            # in other populations
            except_i = 1:S # somehow remove i
            for k in except_i
                for l in 1:N[k]
                    B += U[i]*U[k]*exp( (x[i,j] - x[k,l])^2 / (2*(w[i]+w[k])) )
                        / √(2*π*(w[i]+w[k]))
                end
            end

            w = exp( R[i] - a[i]*(θ[i]-x[i,j])^2/2.0 - c[i]*B )

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
        for j in 1:N[i]

            # birth each offspring
            for k in 1:W[j]

                # consider next individual
                ct += 1

                # draw random trait for this individual
                xₚ[ct] = rand( Normal( x[i,j], √μ ), 1)[1]

            end

        end

        x̄ₚ[i] = mean(xₚ)
        σₚ²[i]= var(xₚ)

    end

    Xₚ = population(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,R=R,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end

# update for community with simple competition term
function update_simple(X)

    @unpack x, N, x̄, σ², R, a, θ, c, μ, V = X

    W = fill(0,N)

    for j in 1:N

        # mean fitness of individual j
        w = exp( R - a*(θ-x[j])^2/2.0 - c*N )

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

    Xₚ = population(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,R=R,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end
