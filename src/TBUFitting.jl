module TBUFitting

using Polynomials
using Unitful
using NLopt 
using OrdinaryDiffEq
using DiffEqParamEstim
using Turing
using Statistics
using StatsBase
using Distributions: Uniform,Normal,Exponential,InverseGamma,truncated,Pareto
using Infiltrator

import Base: show,names,Dict
import Polynomials: fit 
import OrdinaryDiffEq: solve

export ModelConstants
include("constants.jl")

export Model,ModelO,ModelF,ModelD,ModelM,ModelQ,dimensionalize,nondimensionalize,isknown,strain,solve,names,units,parameters,solution,constants,timescale
include("models.jl")

export intensity,FittedModel,fit,model,adjustfit
include("fit.jl")

export initialtrim
include("trim.jl")

export BayesModel,bayes 
include("bayes.jl")

end # module
