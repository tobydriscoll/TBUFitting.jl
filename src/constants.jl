const ρ = 1e3u"g/L"       # density of water, gm/liter
const P̂c = 12.1u"μm/s"
const Vw = 18e-6u"m^3/mol"      # molar volume of water, m^3/mol
const c₀ = 300u"mol/m^3"        # serum osmolarity, mOsM or Osm/m^3, same thing; 300 ok too
const σ₀ = 0.045u"N/m";
const μ = 1.3e-3u"Pa*s";

## solute parameters
# fcr is critical concentation for max value in intensity at fixed d 
const fcr = 0.002; # material property: critical FL concentration in mole frac
const Mw_FL = 376u"g/mol"; # molecular weight of Na_2FL, g/mol
# Naperian extinction coefficient
const eps_f = 1.75e7u"1/(m*M)";   # m^-1 M^-1
# transmittance parameter
const fcr2 = ρ*fcr/Mw_FL

struct ModelConstants
	h₀::Unitful.Length
	ts::Unitful.Time
	f̂₀::Float64
	Pc::Float64
	f₀::Float64
	ϕ::Float64
end

function ModelConstants(h₀::Unitful.Length,ts::Unitful.Time,f̂₀::Number) 
	Pc,f₀,ϕ = deriveconstants(h₀,ts,f̂₀)
	ModelConstants(h₀,ts,f̂₀,Pc,f₀,ϕ)
end

function deriveconstants(h₀::Unitful.Length,ts::Unitful.Time,f̂₀::Number)
	# dimensional parameters computed for trials
	ℓ = uconvert(u"μm",(ts*σ₀*h₀^3/μ)^0.25)
	U = ℓ/ts;
	ϵ = h₀/ℓ;
	V = ϵ*U; 
	f0p2 = f̂₀*ρ/Mw_FL;

	## dimensionless parameters
	# dimensionless extinction coefficient
	ϕ = uconvert(Unitful.NoUnits, eps_f*fcr2*h₀)
	# corneal permeability (non dimnl) 
	Pc = P̂c*Vw*c₀/V;
	f₀ = f̂₀/fcr
	return Pc,f₀,ϕ
end

function show(io::IO,c::ModelConstants)
	println(io,"ModelConstants:")
	println(io,"    h₀ = $(c.h₀), ts = $(c.ts), f̂₀ = $(c.f̂₀)")
	println(io,"    Pc = $(c.Pc), f₀ = $(c.f₀), ϕ = $(c.ϕ)")
end

function Dict(c::ModelConstants)
	Dict( ("h0"=>c.h₀, "ts"=>c.ts, "f0"=>c.f̂₀) )
end