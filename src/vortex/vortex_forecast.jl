

export forecast, VortexForecast, time_advancement!, createsheddedvortices

#### FORECAST OPERATORS ####


mutable struct VortexForecast{withfreestream,Nb,Ne} <: AbstractForecastOperator

    "vortex model from GridPotentialFlow.jl"
    vvm :: Vector{VortexModel}

    intermediate_vm :: VortexModel

    "Number of vortices"
    Nv :: Int64

end

"""
    VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)

Allocate the structure for forecasting of vortex dynamics
"""
function VortexForecast(vvm::Vector{VortexModel{Nb,Ne,TS,TU,TE,TF,TX,ILS}}) where {Nb,Ne,TS,TU,TE,TF,TX,ILS}

    intermediate_bodies = deepcopy(vvm[1].bodies)
    for j=1:Nb
        intermediate_bodies[j].edges = Int64[]
    end
    intermediate_vm = VortexModel(vvm[1].g,vortices=[vvm[1].vortices...],bodies=intermediate_bodies,U∞=vvm[1].U∞)
    withfreestream = vvm[1].U∞ == 0.0 ? false : true
    Nv = length(vvm[1].vortices)
    Nx = 3*Nv
    VortexForecast{withfreestream,Nb,Ne}(vvm,intermediate_vm,Nv)
end



function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{withfreestream,Nb,Ne},i::Int64) where {withfreestream,Nb,Ne}
    @unpack vvm, intermediate_vm = fdata
    vm = vvm[i] #i-th ensemble member
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb
    
    states_to_vortices!(vm,x)
    time_advancement!(vm,intermediate_vm,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm.vortices[end-1:end])
    pushvortices!(vm,vLEnew,vTEnew)
    sol = solve(vm)

    xnew = similar(x[1:end])
    # New vortices released from the two edges augment the state vector by 3*Ne
    append!(xnew,zeros(3*Ne))
    vortices_to_states!(xnew,vm,sol,Δt)
    fdata.Nv = length(vm.vortices)
    return xnew
end

function time_advancement!(vm::VortexModel{Nb,Ne},intermediate_vm::VortexModel{Nb,0},Δt) where {Nb,Ne}
    construct_intermediate_model!(intermediate_vm,vm)
    sol = ConstrainedIBPoissonSolution(intermediate_vm._ψ, intermediate_vm._f, zeros(Float64,Nb), zeros(Float64,Ne))
    advect_vortices!(intermediate_vm,sol,Δt)
    vm.vortices = deepcopy(intermediate_vm.vortices)
    for i=1:Nb
        vm.bodies[i].Γ = deepcopy(intermediate_vm.bodies[i].Γ)
    end
end

"""construct_intermediate_model!(intermediate_vm::VortexModel{Nb,0},vm::VortexModel{Nb,Ne}) where {Nb,Ne} --> VortexModel{Nb,0}
an intermediate vortex model with all fields the same as vm except that intermediate_vm has no regularized edge. This allows
the solution of the system for the existing vortices which solves solve!(sol::ConstrainedIBPoissonSolution, 
vm::VortexModel{Nb,0,ConstrainedIBPoisson{Nb,TU,TF}}) in the vortexmodel.jl file. 
"""
function construct_intermediate_model!(intermediate_vm::VortexModel{Nb,0},vm::VortexModel{Nb,Ne}) where {Nb,Ne}
    intermediate_vm.bodies = deepcopy(vm.bodies)
    intermediate_vm.vortices = deepcopy(vm.vortices)
    intermediate_vm.U∞ = deepcopy(vm.U∞)
    for i=1:Nb
        intermediate_vm.bodies[i].edges = Int64[]
    end
end

"""Advances the motion of vortices in one time step for the existing vortices in the domain and a body with no regularized edge.
Used in the foreward model."""
function advect_vortices!(vm::VortexModel,sol::ConstrainedIBPoissonSolution,Δt)
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm,sol)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

"""Advances the motion of vortices in one time step for the existing vortices in the domain and a body with Ne regularized edge.
Used in the observation model."""
function advect_vortices!(vm::VortexModel,Δt)
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

# vortices released at one-third of the way from the edge to the last released vortex from that edge
function createsheddedvortices(plate::Plate,oldvortices)

    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE, vTE
end