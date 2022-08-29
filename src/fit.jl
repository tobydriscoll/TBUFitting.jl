# Lower and upper bounds for model parameters
lower(::ModelQ) = [0u"μm/minute",-2u"s^-1",-2u"s^-1",0u"s^-1"]
upper(::ModelQ) = [40u"μm/minute",2u"s^-1",2u"s^-1",2u"s^-1"] 
lower(::ModelM) = [0u"μm/minute",-1u"s^-1",-2u"s^-1",0u"s^-1"]
upper(::ModelM) = [40u"μm/minute",2u"s^-1",5u"s^-1",2u"s^-1"] 
lower(m::ModelD) = lower(ModelM(m.constants))[[1,3,4]]
upper(m::ModelD) = upper(ModelM(m.constants))[[1,3,4]]
lower(m::ModelF) = lower(ModelM(m.constants))[[1,2]]
upper(m::ModelF) = upper(ModelM(m.constants))[[1,2]]
lower(m::ModelO) = lower(ModelM(m.constants))[[1]]
upper(m::ModelO) = upper(ModelM(m.constants))[[1]]

# Measure residual using the trapezoid rule to approximate 2-norm
trap(x,y) = 0.5*sum((x[i+1]-x[i])*(y[i+1]+y[i]) for i in 1:length(x)-1)
function misfit(sol,Ifun,t,I) 
	Im = [Ifun(sol(t)...) for t in t]
	return trap(t,(I-Im).^2)
end

# Create callback to halt IVP solver if dh/dt > 0 or I increases at a step
function solvercb(Ifun)
	CBFOO = [0.,0.]
	function increaseI(u,t,integrator) 
		return (u[1] > integrator.uprev[1]) || (Ifun(u[1],u[2]) > Ifun(integrator.uprev...))
	end
	function posdhdt(u,t,integrator) 
		integrator.f(CBFOO,u,integrator.p,t)
		return CBFOO[1]-1e-3
	end
	affect!(integrator) = terminate!(integrator)
	solvercb = CallbackSet(DiscreteCallback(increaseI,affect!),ContinuousCallback(posdhdt,affect!))
end

# Object for a model that has been fit to data
struct FittedModel
	m::AbstractModel
	ts
	t₀
	I₀
	t::AbstractVector
	I::AbstractVector
	residual::AbstractFloat
end

# Create a fitted model from nondimensional data.
# It's assumed that the model given has already been optimized, and we are just calculating the quality of fit.
function FittedModel(m::AbstractModel,t::AbstractVector{<:Real},I::AbstractVector)
	inten = intensity(m)
	y = [(I-inten(t))^2 for (t,I) in zip(t,I)]
	misfit = sqrt(trap(t,y))
	return FittedModel(m,1,0,1,t,I,misfit)
end

# Create a fitted model from dimensional data.
# It's assumed that the model given has already been optimized, and we are just calculating the quality of fit.
function FittedModel(m::AbstractModel,t::AbstractVector{<:Quantity},I::AbstractVector)
	t₀ = t[1]
	ts = t[end]-t₀
	I₀ = I[1]
	c = constants(m)
	c = ModelConstants(c.h₀,ts,c.f̂₀)
	m = typeof(m)(c,m.parameters,m.solution)
	# get residual
	tn = (t.-t₀)/ts
	In = I/I₀
	model_n = FittedModel(m,tn,In)
	return FittedModel(m,ts,t₀,I₀,t,I,model_n.residual)
end

# Convenience accessors
constants(M::FittedModel) = constants(model(M))
parameters(f::FittedModel) = parameters(f.m)
names(f::FittedModel) = names(f.m)
model(f::FittedModel) = f.m

# Access the solution of a fitted model 
solution(f::FittedModel,t::Real;kwargs...) = solution(model(f),t;kwargs...)
function solution(f::FittedModel,t::Unitful.Time;kwargs...) 
	solution(model(f),(t-f.t₀)/f.ts;kwargs...)
end
solution(f::FittedModel,t::AbstractVector;kwargs...) = [ solution(f,t;kwargs...) for t in t ]

# Access the intensity of a fitted model
intensity(f::FittedModel,t::Real) = f.I₀*intensity(model(f),t)
intensity(f::FittedModel,t::Unitful.Time) = f.I₀*intensity(model(f),(t-f.t₀)/f.ts)
intensity(f::FittedModel,t::AbstractVector) = [ intensity(f,t) for t in t ]

# Compact display
function show(io::IO,M::FittedModel) 
	res = round(M.residual,sigdigits=4)
	print(io,"Fitted ",M.m)
	print(";  with residual $res")
end

# Convert main parts to a dictionary
function Dict(f::FittedModel)
	Dict( ("model"=>Dict(model(f)), "t"=>f.t, "I"=>f.I, "residual"=>f.residual) )
end

#######################################################################
# Optimization for fitting to data
#######################################################################

# Fit a single model to nondimensional data, given initial guesses at the parameters
function fit(
	M::AbstractModel,
	t::AbstractVector{<:Real},
	I::AbstractVector{<:Real},
	initpar;
	lower=lower(M), upper=upper(M), method=:LN_NELDERMEAD)

	@assert all(@. 0≤t≤1) "Invalid time vector."
	ivp = makeivp(M)
	inten = intensity(M.constants)

	# Define the loss function.
	misfitfail(sol) = sol.retcode==:Success ? misfit(sol,inten,t,I) : 1/sol.t[end]
	lossfun = build_loss_objective(ivp, Tsit5(), misfitfail, maxiters = 2000, verbose=false,verbose_opt=false,verbose_steps=10,reltol=1e-10,abstol=1e-11,callback=solvercb(inten))

	# Set up optimization.
	n = length(names(M))
	opt = Opt(method,n)
	opt.lower_bounds = nondimensionalize(M,lower)
	opt.upper_bounds = nondimensionalize(M,upper)
	opt.xtol_rel = 1e-5
	opt.xtol_abs = 1e-7
	opt.maxeval = 10000
	opt.min_objective = lossfun
	
	# Optimize over all initializations.
	bestmin,bestpar = Inf,nondimensionalize(M,initpar[1])
	bestret = []
	for p̂ in initpar
		p = nondimensionalize(M,p̂)
		minval,minp,ret = NLopt.optimize(opt,p)
		#(typeof(M)==ModelM) && (@show minval,p,minp,ret)  # debugging
		if minval < bestmin
			bestmin,bestpar,bestret = minval,minp,ret
		end
	end

	if bestret == :FAILURE
		@warn "optimization failed for $(typeof(M))"
	end

	#@show typeof(M),bestmin,bestpar,dimensionalize(M,bestpar)  # debugging

	return FittedModel(solve(M,bestpar),t,I)
end

# Fit all the models to nondimensional data, using simpler ones to help initialize the more complex ones.
function fit(con::ModelConstants,t::AbstractVector{<:Real},I)
	function safe(m::AbstractModel)
		# get parameters that are safely pushed away from the bounds
		p̂ = parameters(m)
		p̂ = min.(p̂,0.98*upper(m))
		return max.(p̂,0.98*lower(m))
	end

	## Start with model O
	init = [ [v*u"μm/minute"] for v in [1;5:5:35] ]
	modO = fit(ModelO(con),t,I,init)

	## Model F
	init = [ [2u"μm/minute",0.06u"s^-1"], 
		[2u"μm/minute",-0.06u"s^-1"], 
		[10u"μm/minute",0.06u"s^-1"], 
		[10u"μm/minute",-0.06u"s^-1"],
		[0u"μm/minute",-0.06u"s^-1"],
		[0u"μm/minute",0.06u"s^-1"],
		]
	
	# Use O result to initialize model F 
	p̂ = safe(modO.m) 
	append!(init,[ 
		[p̂[1],0u"s^-1"], [p̂[1],0.06u"s^-1"], [p̂[1],-0.06u"s^-1"], 
	 ] )
	modF = fit(ModelF(con),t,I,init)

	## Model D
	init = [ 
		[1u"μm/minute",0.2u"s^-1",0.2u"s^-1"], 
		[6u"μm/minute",0.2u"s^-1",0.2u"s^-1"], 
		[0u"μm/minute",0.1u"s^-1",0.2u"s^-1"],
		[15u"μm/minute",0.1u"s^-1",0.2u"s^-1"],
		[0u"μm/minute",-0.1u"s^-1",0.2u"s^-1"],
		]
	p̂ = safe(modO.m)
	append!(init,[ 
		[p̂[1],0u"s^-1",0u"s^-1"], [p̂[1],0.1u"s^-1",0.1u"s^-1"], [p̂[1],-0.1u"s^-1",0.1u"s^-1"], 
	 ] )
	p̂ = safe(modF.m)
	append!(init,[ 
		 [p̂[1],p̂[2],0u"s^-1"], [p̂[1],p̂[2],0.5u"s^-1"], 
	  ] )
	modD = fit(ModelD(con),t,I,init)

	## Model M
	init = [ 
		[1u"μm/minute",0.1u"s^-1",0.2u"s^-1",0.5u"s^-1"], 
		[6u"μm/minute",-0.1u"s^-1",0.2u"s^-1",0.5u"s^-1"], 
		[15u"μm/minute",0.1u"s^-1",0.2u"s^-1",0.5u"s^-1"], 
		[24u"μm/minute",0.1u"s^-1",-0.2u"s^-1",0.8u"s^-1"], 
		[20u"μm/minute",0.1u"s^-1",-0.2u"s^-1",0.8u"s^-1"], 
		[0u"μm/minute",0.2u"s^-1",0.2u"s^-1",0.5u"s^-1"], 
		]
	p̂ = safe(modO.m)
	append!(init,[ 
		[p̂[1],0u"s^-1",0u"s^-1",0u"s^-1"], [p̂[1],0.2u"s^-1",-0.1u"s^-1",0.5u"s^-1"],
	] )
	p̂ = safe(modF.m)
	append!(init,[ 
		 [p̂[1],p̂[2],0u"s^-1",0u"s^-1"], [p̂[1],p̂[2],-0.1u"s^-1",0.5u"s^-1"], 
	] )
	p̂ = safe(modD.m)
	append!(init,[ 
		[p̂[1],0u"s^-1",p̂[2],p̂[3]], 
	] )
	modM = fit(ModelM(con),t,I,init)

	## Model Q
	modQ = fit(ModelQ(con),t,I,init)

	## Re-try model F with what D and M found, to avoid a poor local min.
	init = [ safe(modF.m) ]
	p̂ = safe(ModelF(con,modD.m.parameters[1:2]))
	push!(init,p̂)
	pM = modM.m.parameters[1:3]
	p̂ = safe(ModelF(con,[pM[1],pM[2]+pM[3]]))
	push!(init,p̂)
	try
		modF = fit(ModelF(con),t,I,init)
	catch
		@warn "There was a problem with final model F fit"
	end

	return (O=modO,F=modF,D=modD,M=modM,Q=modQ)
end

# Fit all the models to dimensional data.
function fit(con::ModelConstants,t::AbstractVector{<:Quantity},I)
	# Rescale time and intensity.
	ts = t[end]-t[1]
	tt = (t.-t[1])/ts
	II = I/I[1]
	cc = ModelConstants(con.h₀,ts,con.f̂₀)
	return fit(cc,tt,II)  # use the nondimensional code
end
