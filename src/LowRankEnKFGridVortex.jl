module LowRankEnKFGridVortex

using GridPotentialFlow
using RigidBodyTools
using LinearAlgebra
using Interpolations
using ProgressMeter
using Statistics
using NamedColors
using Distributions
using Combinatorics
using UnPack


include("ensemble.jl")
include("forecast.jl")
include("observation.jl")

include("vortex/initial.jl")
include("vortex/vortex.jl")
include("vortex/vortex_forecast.jl")

include("pressure/vortex_pressure_observation.jl")

include("DA/types.jl")
include("DA/generate_twin_experiment.jl")
include("DA/enkf.jl")
include("DA/state_utilities.jl")
include("DA/classification.jl")
include("DA/MCMC.jl")

include("experiments/routines.jl")

end
