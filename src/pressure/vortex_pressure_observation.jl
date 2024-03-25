
export observations, VortexPressure, Sensor, calculate_impulse

#### OBSERVATION OPERATORS ####

struct Sensor{T}
    x :: Vector{T}
    y :: Vector{T}
    Nsens :: Int64
end

mutable struct VortexPressure{Ny,withfreestream,Nb,Ne,TS<:Union{AbstractPotentialFlowSystem,Laplacian},DST,TVG<:Edges,TX<:VectorData,TVF<:VectorData,TF<:ScalarData,TP<:Nodes} <: AbstractObservationOperator{Ny,true}
    sens::Sensor
    config::VortexForecast{withfreestream,Nb,Ne,TS}
    # intermediate_vm :: VortexModel
    Δs::DST

    """for pressure calculations"""
    v̄::TVG
    Xs::TX
    v̄s::TVF
    dp::TF
    p̄::TP
end

function VortexPressure(sens::Sensor,config::VortexForecast{withfreestream,Nb,Ne,TS}) where {withfreestream,Nb,Ne,TS}
    @unpack vm = config
    Nv = config.Nv
    Nx = 3*Nv
    Ny = sens.Nsens
    Δs = dlengthmid(vm.bodies[1].points)
    v̄ = Edges(Primal,vm._ψ)
    Xs = VectorData(collect(vm.bodies))
    v̄s = VectorData(Xs)
    dp = ScalarData(Xs)
    p̄ = Nodes(Primal,vm._ψ)
    return VortexPressure{Ny,withfreestream,Nb,Ne,TS,typeof(Δs),typeof(v̄),typeof(Xs),typeof(v̄s),typeof(dp),typeof(p̄)}(sens,config,Δs,v̄,Xs,v̄s,dp,p̄)
end




function observations(x::AbstractVector,t,Δt,obs::VortexPressure{Ny,true,Nb,Ne,<:UnsteadyRegularizedIBPoisson{Nb,Ne}}) where {Ny,Nb,Ne}
    @unpack sens, config, Δs = obs
    @unpack vm = config
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    states_to_vortices!(vm,x)
    vm.bodies[1].Γ = -sum(vm.vortices.Γ[1:end-2])
    vmn = deepcopy(vm)
    soln = solve(vmn)
    γn = soln.f./Δs
    setvortexstrengths!(vmn, soln.δΓ_vec, length(vmn.vortices)-1:length(vmn.vortices))
    subtractcirculation!(vmn.bodies, soln.δΓ_vec)

    #solution at the next time step n+1
    vm1 = deepcopy(vmn)
    solve!(soln, vmn)
    # advect_vortices!(vm1,soln,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm1.vortices)
    pushvortices!(vm1,vLEnew,vTEnew)
    solnp1 = solve(vm1)
    γnp1 = solnp1.f./Δs
    # setvortexstrengths!(vm1, solnp1.δΓ_vec, length(vm1.vortices)-1:length(vm1.vortices))
    # subtractcirculation!(vm1.bodies, solnp1.δΓ_vec)

    velocity!(obs.v̄,soln.ψ,vmn.ilsys)
    GridPotentialFlow.surface_velocity!(obs.v̄s,obs.v̄,vmn.ilsys)
    dp2 = deepcopy(obs.dp)
    obs.dp, dp2 = pressurejump!(obs.dp,γn,γnp1,obs.v̄s,Δt,vmn.ilsys)
    # pressurejump!(obs.dp,soln.f,vm1,solnp1.f,obs.v̄s,Δt,Δs,vmn.ilsys)  #another approach for computing dp
    # pressure!(obs.p̄,obs.v̄,obs.dp,vmn.ilsys)
    # p⁺, p⁻ = sided_pressures(obs.p̄,obs.dp,vmn.ilsys)

    # dp_sens = surface_interpolation(obs.dp,pfb,sens)

    return obs.dp, dp2
end

function observations(x::AbstractVector,t,Δt,obs::VortexPressure{Ny,true,Nb,Ne,<:ConstrainedIBPoisson{Nb}}) where {Ny,Nb,Ne}
    @unpack sens, config, Δs = obs
    @unpack vm = config
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    states_to_vortices!(vm,x) #i-th ensemble member
    vm.bodies[1].Γ = -sum(vm.vortices.Γ)
    vmn = deepcopy(vm)
    soln = solve(vmn)
    γn = soln.f./Δs

    #solution at the next time step n+1
    vm1 = deepcopy(vmn)
    solve!(soln, vmn)
    advect_vortices!(vm1,soln,Δt)
    solnp1 = solve(vm1)
    γnp1 = solnp1.f./Δs

    velocity!(obs.v̄,soln.ψ,vmn.ilsys)
    GridPotentialFlow.surface_velocity!(obs.v̄s,obs.v̄,vmn.ilsys)
    # dp2 = deepcopy(obs.dp)
    pressurejump!(obs.dp,γn,γnp1,obs.v̄s,Δt,vmn.ilsys)
    # pressurejump!(obs.dp,soln.f,vm1,solnp1.f,obs.v̄s,Δt,Δs,vmn.ilsys)  #another approach for computing dp
    pressure!(obs.p̄,obs.v̄,obs.dp,vmn.ilsys)
    # p⁺, p⁻ = sided_pressures(obs.p̄,obs.dp,vmn.ilsys)

    # dp_sens = surface_interpolation(obs.dp,pfb,sens)

    return dp_sens, obs.p̄
end

"""impulse"""
function calculate_impulse(config::VortexForecast,Δt)
    @unpack vm = config
    @unpack bodies = vm
    @unpack points = bodies[1]
    plate = bodies[1]
    X = getvortexpositions(vm)
    Ẋ = deepcopy(X)
    sol = solve(vm)
    setvortexstrengths!(vm, sol.δΓ_vec, length(X.u)-1:length(X.u))
    subtractcirculation!(vm.bodies, sol.δΓ_vec)
    Px, Py = impulse(vm)
    vortexvelocities!(Ẋ, vm, sol.ψ)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(vm, X)

    vLE, vTE = createsheddedvortices(points,vm.vortices)
    pushvortices!(vm,vLE,vTE)
    return Px, Py
end

