priors(::ModelD) = [
    truncated(Pareto(0.5),0,40),
    truncated(Normal(0.025,0.1),-2,2),
    truncated(Exponential(1),0,2)
    ]
priors(::ModelF) = [
    truncated(Pareto(0.5),0,40),
    truncated(Normal(0.025,0.1),-2,2)
    ]
# priors(::ModelO) = [truncated(Exponential(3),0,20)]

struct BayesModel
	m::AbstractModel
	ts
	t₀
	I₀
	t::AbstractVector
	I::AbstractVector
    result
end

function bayes(M::AbstractModel,t::AbstractVector{<:Real},I::AbstractVector;num_samples=1000,threads=1)
    ivp = makeivp(M)
    Ifun = intensity(constants(M))
    ϕ = priors(M)
    N = length(ϕ)
    syms = [Turing.@varname(θ[i]) for i in 1:N]
    # likelihood = (u,p,t,σ) -> MvLogNormal(MvNormal(log.(u),σ*ones(length(u))))
    likelihood = (u,p,t,σ) -> MvNormal(u,σ*ones(length(u)))
    
    Turing.@model function mf(x,::Type{T}=Float64) where {T<:Real}
        θ = Vector{T}(undef,N)
        for i in 1:N
            θ[i] ~ NamedDist(ϕ[i], syms[i])
        end
        p = nondimensionalize(M,θ)
        sol = OrdinaryDiffEq.solve(remake(ivp,p=p),AutoVern7(KenCarp4()),saveat=t,reltol=1e-8,abstol=1e-7)
        if (sol.retcode !== :Success) 
            @info sol.retcode
            Turing.acclogp!(_varinfo, -Inf)
            return
        end

        σ ~ InverseGamma(2,0.04)
        hs,cs = eachrow(Array(sol))
        Is = Ifun.(hs,cs)
        h_inc = Inf#something(findfirst(diff(hs) .> 0),Inf)
        I_inc = Inf#something(findfirst(diff(Is) .> 0),Inf)
        k_inc = (sol.retcode !== :Success) ? 1 : min(h_inc,I_inc)

        if !isinf(k_inc)
            k_inc = Int(k_inc)
            #@info "ODE failed at time $(t[k_inc])"
            Turing.acclogp!(_varinfo, log(t[k_inc+1]))
            return
        end
        
        # σ = nothing
        #try 
            x ~ likelihood(Is,θ,Inf,σ)
        #catch except 
         #   @show Float64.(hs)
          #  @show Float64.(Is)
           # @show k_inc
           # throw(except)
        #end
        return
    end false
    
    model = mf(I)
    if threads > 1
        result = sample(model,Turing.NUTS(0.65),MCMCThreads(),num_samples÷threads,threads; progress=true)
    else
        result = sample(model,Turing.NUTS(0.65),num_samples; progress=true)
    end
    return BayesModel(M,1,0,1,t,I,result)
end

function bayes(M::AbstractModel,t::AbstractVector{<:Quantity},I::AbstractVector;kwargs...)
    t₀ = t[1]
    ts = t[end]-t₀
    I₀ = I[1]
    tn = (t.-t₀)/ts
	In = I/I₀
    mb = bayes(M,tn,In;kwargs...)
    @show mb.result
    return BayesModel(M,ts,t₀,I₀,tn,In,mb.result)
 end
