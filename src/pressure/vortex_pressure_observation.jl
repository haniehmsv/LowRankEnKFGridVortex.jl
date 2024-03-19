
export observations, VortexPressure, Sensor

#### OBSERVATION OPERATORS ####

struct Sensor{T}
  x :: Vector{T}
  y :: Vector{T}
  Nsens :: Int64
end

mutable struct VortexPressure{Ny,withfreestream,Nb,Ne,TS<:Union{AbstractPotentialFlowSystem,Laplacian},DST,TVG<:Edges,TX<:VectorData,TVF<:VectorData,TF<:ScalarData,TP<:Nodes} <: AbstractObservationOperator{Ny,true}
    sens::Sensor
    config::VortexForecast{withfreestream,Nb,Ne,TS}
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
    Ny = length(sens.x)
    Δs = dlengthmid(vm.bodies[1].points)
    v̄ = Edges(Primal,vm._ψ,dtype=Real)
    Xs = VectorData(collect(vm.bodies))
    v̄s = VectorData(Xs,dtype=Real)
    dp = ScalarData(Xs,dtype=Real)
    p̄ = Nodes(Primal,vm._ψ,dtype=Real)
    return VortexPressure{Ny,withfreestream,Nb,Ne,TS,typeof(Δs),typeof(v̄),typeof(Xs),typeof(v̄s),typeof(dp),typeof(p̄)}(sens,config,Δs,v̄,Xs,v̄s,dp,p̄)
end


function observations(x::AbstractVector,t,Δt,obs::VortexPressure{Ny,true,Nb,Ne,<:UnsteadyRegularizedIBPoisson{Nb,Ne}}) where {Ny,Nb,Ne}
    @unpack sens, config, Δs = obs
    @unpack vm = config
    @unpack bodies = vm
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    #Solution at the current time step
    states_to_vortices!(vm,x) #i-th ensemble member
    vm.bodies[1].Γ = -sum(vm.vortices.Γ[1:end-2])
    vmn = deepcopy(vm)
    if !isempty(FD.partials(vmn.vortices.Γ[end-1]))
        T = get_tag(vmn.vortices.Γ[end-1])
        vmn.vortices.Γ[end-1] = FD.Dual{T}(FD.value(x[end])*Δt,FD.partials.(vmn.vortices.Γ[end-1])...)
    else
        vmn.vortices.Γ[end-1] = FD.value(x[end])*Δt
    end
    subtractcirculation!(vmn.bodies, [vmn.vortices.Γ[end-1]])
    soln = solve(vmn)
    vmn.vortices.Γ[end] = soln.δΓ_vec[1]
    subtractcirculation!(vmn.bodies, soln.δΓ_vec)
    γn = similar(soln.f)
    γn .= ScalarData(soln.f.data./Δs)

    #solution at the next time step n+1
    vm1 = deepcopy(vmn)
    solve!(soln, vmn)
    advect_vortices!(vm1,soln,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm1.vortices,DT=Real)
    pushvortices!(vm1,vLEnew,vTEnew)
    vm1.vortices.Γ[end-1] = FD.value(x[end])*Δt
    subtractcirculation!(vmn.bodies, [vmn.vortices.Γ[end-1]])
    solnp1 = solve(vm1)
    vm1.vortices.Γ[end] = solnp1.δΓ_vec[1]
    subtractcirculation!(vm1.bodies, solnp1.δΓ_vec)
    γnp1 = similar(solnp1.f)
    γnp1 .= ScalarData(solnp1.f.data./Δs)

    velocity!(obs.v̄,soln.ψ,vmn.ilsys)
    GridPotentialFlow.surface_velocity!(obs.v̄s,obs.v̄,vmn.ilsys)

    pressurejump!(obs.dp,γn,γnp1,obs.v̄s,Δt,vmn.ilsys)
    # pressure!(obs.p̄,obs.v̄,obs.dp,vmn.ilsys)
    # p⁺, p⁻ = sided_pressures(obs.p̄,obs.dp,vmn.ilsys)

    dp_sens = surface_interpolation(obs.dp,pfb,sens)

    return dp_sens
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
    γn = similar(soln.f)
    γn .= ScalarData(soln.f.data./Δs)

    #solution at the next time step n+1
    vm1 = deepcopy(vmn)
    solve!(soln, vmn)
    advect_vortices!(vm1,soln,Δt)
    solnp1 = solve(vm1)
    γnp1 = similar(solnp1.f)
    γnp1 .= ScalarData(solnp1.f.data./Δs)

    velocity!(obs.v̄,soln.ψ,vmn.ilsys)
    GridPotentialFlow.surface_velocity!(obs.v̄s,obs.v̄,vmn.ilsys)

    pressurejump!(obs.dp,γn,γnp1,obs.v̄s,Δt,vmn.ilsys)
    # pressurejump!(obs.dp,soln.f,vm1,solnp1.f,obs.v̄s,Δt,Δs,vmn.ilsys)  #another approach for computing dp
    # pressure!(obs.p̄,obs.v̄,obs.dp,vmn.ilsys)
    # p⁺, p⁻ = sided_pressures(obs.p̄,obs.dp,vmn.ilsys)

    dp_sens = surface_interpolation(obs.dp,pfb,sens)

    return dp_sens
end

