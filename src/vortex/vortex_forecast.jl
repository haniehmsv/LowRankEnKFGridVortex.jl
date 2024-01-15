

export forecast, VortexForecast, time_advancement!, createsheddedvortices

#### FORECAST OPERATORS ####


mutable struct VortexForecast{withfreestream,Nb,Ne} <: AbstractForecastOperator

    "vortex model from GridPotentialFlow.jl"
    vvm :: Vector{VortexModel}

    intermediate_vvm :: Vector{VortexModel}

    "Number of vortices"
    Nv :: Int64

end

"""
    VortexForecast(vm::VortexModel,pfb::PotentialFlowBody)

Allocate the structure for forecasting of vortex dynamics
"""
function VortexForecast(vvm::Vector{VortexModel{Nb,Ne,TS,TU,TE,TF,TX,ILS}}) where {Nb,Ne,TS,TU,TE,TF,TX,ILS}

    intermediate_vvm = VortexModel[]
    for i in 1:length(vvm)
        intermediate_bodies = deepcopy(vvm[i].bodies)
        for j=1:Nb
            intermediate_bodies[j].edges = Int64[]
        end
        push!(intermediate_vvm, VortexModel(vvm[i].g,vortices=[vvm[i].vortices...],bodies=intermediate_bodies,U∞=vvm[i].U∞))
    end

    withfreestream = vvm[1].U∞ == 0.0 ? false : true
    Nv = length(vvm[1].vortices)
    Nx = 3*Nv
    VortexForecast{withfreestream,Nb,Ne}(vvm,intermediate_vvm,Nv)
end



function forecast(x::AbstractVector,t,Δt,fdata::VortexForecast{Nb,Ne},i::Int64) where {Nb,Ne}
    @unpack vvm, intermediate_vvm = fdata
    vm = vvm[i] #i-th ensemble member
    int_vm = intermediate_vvm[i]
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb
    
    states_to_vortices!(vm,x,Δt)

    #solve the system for the existing vortices (vortexmode.jl)
    #which solves solve!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,0,ConstrainedIBPoisson{Nb,TU,TF}})
    int_vm.bodies = deepcopy(vm.bodies)
    int_vm.vortices = deepcopy(vm.vortices)
    int_vm.U∞ = deepcopy(vm.U∞)
    for i=1:Nb
        int_vm.bodies[i].edges = Int64[]
    end

    sol = ConstrainedIBPoissonSolution(int_vm._ψ, int_vm._f, zeros(Float64,Nb), zeros(Float64,Ne))
    time_advancement!(int_vm,sol,Δt)
    vm.vortices = deepcopy(int_vm.vortices)
    for i=1:Nb
        vm.bodies[i].Γ = deepcopy(int_vm.bodies[i].Γ)
    end

    vLEnew, vTEnew = createsheddedvortices(points,vm.vortices[end-1:end])
    pushvortices!(vm,vLEnew,vTEnew)
    newsol = solve(vm)

    xnew = similar(x[1:end])
    # New vortices released from the two edges augment the state vector by 3*Ne
    append!(xnew,zeros(6))
    vortices_to_states!(xnew,vm,newsol,Δt)
    fdata.Nv = length(vm.vortices)
    return xnew
end

function time_advancement!(vm::VortexModel,sol::ConstrainedIBPoissonSolution,Δt)
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm,sol)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

function time_advancement!(vm::VortexModel,Δt)
    X = getvortexpositions(vm)
    Ẋ = vortexvelocities!(vm,sol)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)
end

# vortices released at one-third of the way from the edge to the last released vortex from that edge
function createsheddedvortices(plate::Plate,oldvortices)

    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE, vTE
end