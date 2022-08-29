module TBUFitting

using Polynomials
using Unitful
using NLopt 
using OrdinaryDiffEq
using DiffEqParamEstim
using Turing
using Statistics
using StatsBase
using Distributions: Uniform,Normal,Exponential,InverseGamma,truncated

import Base: show,names,Dict
import Polynomials: fit 
import OrdinaryDiffEq: solve

# Deal with dimensional and nondimensional constants in the models
export ModelConstants
include("constants.jl")

# Define the ODE models
export Model,ModelO,ModelF,ModelD,ModelM,ModelQ,dimensionalize,nondimensionalize,isknown,strain,solve,names,units,parameters,solution,constants,timescale
include("models.jl")

# Fit data to the models
export intensity,FittedModel,fit,model,adjustfit
include("fit.jl")

# Trim beginning and end of data to get decreasing period only
export initialtrim
include("trim.jl")

# Bayesian parameter fits (experimental)
export BayesModel,bayes 
include("bayes.jl")

end # module
