using LowRankEnKFGridVortex
using Test

const GROUP = get(ENV, "GROUP", "All")
if GROUP == "All" || GROUP == "ForwardDiff"
    include("forwarddiff.jl")
end
