# h(t): film thickness (relative to initial)
# c(t): film osmolarity (relative to initial)

# <:Quantity types have dimensional units, <:Real types are dimensionless


# Create a function that computes FL intensity as a function of h and c.
# Depends on initial FL and Naperian constant
function intensity(C::ModelConstants)
	f₀ = C.f₀
	ϕ = C.ϕ
	II(h,c) = (1 - exp(-ϕ*f₀*h*c)) / (1 + (f₀*c)^2)
	I₀ = II(1,1)
	return (h,c) -> II(h,c)/I₀
end

#######################################################################
# Parent type for all the ODE models
#######################################################################

abstract type AbstractModel end 

# Create an IVP for a given model instance.
function makeivp(M::AbstractModel)
	Pc = M.constants.Pc
	function ode!(du,u,p,t)
		h, c = u
		g = strain(M,p)
		gt = g(t)
		du[1] = -gt*h + Pc*(c-1) - p[1]
		du[2] = -gt*c - du[1]*c / h
	end
	return ODEProblem(ode!,[1.,1.],(0.,1.),M.parameters)
end

# Convenience access functions for all models
constants(M::AbstractModel) = M.constants
timescale(M::AbstractModel) = M.constants.ts
parameters(M::AbstractModel) = uconvert.(units(M),M.parameters)
nondimensionalize(M::AbstractModel) = nondimensionalize(M,parameters(M)) 
isknown(M::AbstractModel) = !isnothing(M.solution)
strain(M::AbstractModel) = strain(M,nondimensionalize(M,M.parameters))

# Solve a particular model type at given dimensional parameters
function solve(M::AbstractModel,p̂::AbstractVector{<:Quantity})
	ivp = makeivp(M)
	p = nondimensionalize(M,p̂)
	solution = OrdinaryDiffEq.solve(remake(ivp,p=p),Tsit5(),reltol=1e-10,abstol=1e-11)
	if solution.retcode !== :Success 
		@warn "Solution was terminated without success, final time = $(solution.t[end])"
	end
	return typeof(M)(M.constants,p̂,solution)
end

# Convert nondimensional parameters and solve
solve(M::AbstractModel,p::AbstractVector{<:Real}) = solve(M,dimensionalize(M,p))

# Access solution and intensity functions for either dimensional or dimensionless time
solution(M::AbstractModel,t::Real;idxs=nothing) = M.solution(t;idxs)
solution(M::AbstractModel) = (t;kwargs...) -> solution(M,t;kwargs...)
intensity(M::AbstractModel) = t -> intensity(M.constants)(solution(M,t)...)
intensity(M::AbstractModel,t) = intensity(M.constants)(solution(M,t)...)

# Compact display
function show(io::IO,M::AbstractModel) 
	if isempty(M.parameters)
		print("$(typeof(M)) with unknown parameters")
	else
		print(io,"$(typeof(M)) with parameters")
		for (n,p,u) in zip(names(M),M.parameters,units(M))
			print(io," $n = $(round(u,p,digits=4)),")
		end
		print(io,"\b")
	end
end

# Convert model to a dict of constants and model parameter values
function Dict(m::AbstractModel)
	p = uconvert.(units(m),parameters(m))
	Dict( ("constants"=>Dict(m.constants), "parameters"=>p) )
end

# Universal constructor
# Kludge: use the number of parameters to determine model type on the fly
function Model(con::ModelConstants,p::AbstractVector{<:Quantity})
	if length(p)==1
		ModelO(con,p)
	elseif length(p)==2
		ModelF(con,p)
	elseif length(p)==3
		ModelD(con,p)
	else
		ModelM(con,p)
	end
end

# Proper constructors
Model(::Val{'O'},c::ModelConstants) = ModelO(c)
Model(::Val{'F'},c::ModelConstants) = ModelF(c)
Model(::Val{'D'},c::ModelConstants) = ModelD(c)
Model(::Val{'M'},c::ModelConstants) = ModelM(c)
Model(::Val{'Q'},c::ModelConstants) = ModelQ(c)

#######################################################################
# Specific model types
#######################################################################

## Model M 
# Strain can relax from one nonzero value to another 
struct ModelM <: AbstractModel 
	constants::ModelConstants
	parameters::AbstractVector{<:Quantity}
	solution
end

# Constructors
ModelM(c::ModelConstants) = ModelM(c,Quantity[],nothing)
ModelM(c::ModelConstants,p̂::AbstractVector{<:Quantity}) = solve(ModelM(c),p̂)

# Parameter names and units
names(::ModelM) = ["v","a","b₁","b₂"]
units(::ModelM) = [u"μm/minute", u"s^-1", u"s^-1", u"s^-1"]

# Access strain function
strain(M::ModelM,p::AbstractVector) = t -> p[2] + p[3] * exp(-p[4]*t)

# Convert to dimensional parameters
function dimensionalize(M::ModelM,p)
	c = M.constants
	return [p[1]*c.h₀/c.ts, p[2]/c.ts, p[3]/c.ts, p[4]/c.ts]
end

# Convert to nondimensional parameters
function nondimensionalize(M::ModelM,p̂)
	c = M.constants
	return uconvert.( Unitful.NoUnits, [p̂[1]*c.ts/c.h₀, p̂[2]*c.ts, p̂[3]*c.ts, p̂[4]*c.ts] )
end

function nondimensionalize(M::ModelM,p̂::AbstractVector{<:Real})
	c = M.constants
	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
	b = ustrip(uconvert(unit(1/units(M)[2]),c.ts))
	return [ p̂[1]*a, p̂[2]*b, p̂[3]*b, p̂[4]*b ]
end

## Model Q
# experimental (not useful)

struct ModelQ <: AbstractModel 
	constants::ModelConstants
	parameters::AbstractVector{<:Quantity}
	solution
end

ModelQ(c::ModelConstants) = ModelQ(c,Quantity[],nothing)
ModelQ(c::ModelConstants,p̂::AbstractVector{<:Quantity}) = solve(ModelQ(c),p̂)

names(::ModelQ) = ["v","g₀","gₘ","g₁"]
units(::ModelQ) = [u"μm/minute", u"s^-1", u"s^-1", u"s^-1"]

function strain(M::ModelQ,p::AbstractVector{<:Real})
	return t -> p[2] + p[3]*tanh(p[4]*t)
end

function dimensionalize(M::ModelQ,p)
	c = M.constants
	return [p[1]*c.h₀/c.ts, p[2]/c.ts, p[3]/c.ts, p[4]/c.ts]
end

function nondimensionalize(M::ModelQ,p̂)
	c = M.constants
	return uconvert.( Unitful.NoUnits, [p̂[1]*c.ts/c.h₀, p̂[2]*c.ts, p̂[3]*c.ts, p̂[4]*c.ts] )
end

function nondimensionalize(M::ModelQ,p̂::AbstractVector{<:Real})
	c = M.constants
	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
	b = ustrip(uconvert(unit(1/units(M)[2]),c.ts))
	return [ p̂[1]*a, p̂[2]*b, p̂[3]*b, p̂[4]*b ]
end

## Model D
# Strain decays from nonzero to zero

struct ModelD <: AbstractModel 
	constants::ModelConstants
	parameters::AbstractVector{<:Quantity}
	solution
end

# Constructors
ModelD(c::ModelConstants) = ModelD(c,Quantity[],nothing)
ModelD(c::ModelConstants,p̂::AbstractVector{<:Quantity}) = solve(ModelD(c),p̂)

# Parameter names and units
names(::ModelD) = ["v","b₁","b₂"]
units(::ModelD) = [u"μm/minute", u"s^-1", u"s^-1"]

# Access to strain function
strain(::ModelD,p::AbstractVector) = t -> p[2] * exp(-p[3]*t)

# Convert to dimensional parameters
function dimensionalize(M::ModelD,p)
	c = M.constants
	return [p[1]*c.h₀/c.ts, p[2]/c.ts, p[3]/c.ts]
end

# Convert to nondimensional parameters
function nondimensionalize(M::ModelD,p̂::AbstractVector{<:Quantity})
	c = M.constants
	return uconvert.( Unitful.NoUnits, [ p̂[1]*c.ts/c.h₀, p̂[2]*c.ts, p̂[3]*c.ts ] )
end
function nondimensionalize(M::ModelD,p̂::AbstractVector{<:Real})
	c = M.constants
	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
	b = ustrip(uconvert(unit(1/units(M)[2]),c.ts))
	return [ p̂[1]*a, p̂[2]*b, p̂[3]*b ]
end

## Model F
# Constant strain rate
struct ModelF <: AbstractModel 
	constants::ModelConstants
	parameters::AbstractVector{<:Quantity}
	solution
end

# Constructors
ModelF(c::ModelConstants) = ModelF(c,Quantity[],nothing)
ModelF(c::ModelConstants,p̂::AbstractVector{<:Quantity}) = solve(ModelF(c),p̂)

# Parameter names and units
names(::ModelF) = ["v","a"]
units(::ModelF) = [u"μm/minute", u"s^-1"]

# Access to strain function
strain(::ModelF,p::AbstractVector) = t -> p[2]

# Convert to dimensional parameters
function dimensionalize(M::ModelF,p)
	c = M.constants
	return [p[1]*c.h₀/c.ts, p[2]/c.ts ]
end

# Convert to nondimensional parameters
function nondimensionalize(M::ModelF,p̂::AbstractVector{<:Quantity})
	c = M.constants
	return uconvert.( Unitful.NoUnits, [ p̂[1]*c.ts/c.h₀, p̂[2]*c.ts ] )
end

function nondimensionalize(M::ModelF,p̂::AbstractVector{<:Real})
	c = M.constants
	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
	b = ustrip(uconvert(unit(1/units(M)[2]),c.ts))
	return [ p̂[1]*a, p̂[2]*b ]
end

## Model O
# Evaporation only (no tangential strain or flow)
struct ModelO <: AbstractModel 
	constants::ModelConstants
	parameters::AbstractVector{<:Quantity}
	solution
end

# Constructors
ModelO(c::ModelConstants) = ModelO(c,Quantity[],nothing)
ModelO(c::ModelConstants,p̂::AbstractVector{<:Quantity}) = solve(ModelO(c),p̂)

# Parameter names and units
names(::ModelO) = ["v"]
units(::ModelO) = [u"μm/minute"]

# Access to strain function
strain(::ModelO,p::AbstractVector) = t -> 0

# Convert to dimensional parameters
function dimensionalize(M::ModelO,p)
	c = M.constants
	return [p[1]*c.h₀/c.ts ]
end

# Convert to nondimensional parameters
function nondimensionalize(M::ModelO,p̂)
	c = M.constants
	return uconvert.( Unitful.NoUnits, [ p̂[1]*c.ts/c.h₀] )
end

function nondimensionalize(M::ModelO,p̂::AbstractVector{<:Real})
	c = M.constants
	a = ustrip(uconvert(unit(1/units(M)[1]),c.ts/c.h₀))
	return [ p̂[1]*a ]
end
