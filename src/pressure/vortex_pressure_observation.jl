
export observations, VortexPressure, Sensor

#### OBSERVATION OPERATORS ####

struct Sensor
  x :: Vector{Float64}
  y :: Vector{Float64}
  Nsens :: Int64
end

mutable struct VortexPressure{Ny,withfreestream,DST,TVG,TX,TVF,TF,TP} <: AbstractObservationOperator{Ny,true}
    sens::Sensor
    config::VortexForecast
    # intermediate_vm :: VortexModel
    Δs::DST

    """for pressure calculations"""
    v̄::TVG
    Xs::TX
    v̄s::TVF
    dp::TF
    p̄::TP
end

function VortexPressure(sens::Sensor,config::VortexForecast)
    @unpack vvm= config
    vm = vvm[1]
    withfreestream = vm.U∞ == 0.0 ? false : true
    Nv = config.Nv
    Nx = 3*Nv
    Ny = length(sens.x)
    Δs = dlengthmid(vm.bodies[1].points)
    v̄ = Edges(Primal,vm._ψ)
    Xs = VectorData(collect(vm.bodies))
    v̄s = VectorData(Xs)
    dp = ScalarData(Xs)
    p̄ = Nodes(Primal,vm._ψ)
    return VortexPressure{Ny,withfreestream,typeof(Δs),typeof(v̄),typeof(Xs),typeof(v̄s),typeof(dp),typeof(p̄)}(sens,config,Δs,v̄,Xs,v̄s,dp,p̄)
end


function observations(x::AbstractVector,t,Δt,obs::VortexPressure,i::Int64)
    @unpack sens, config, Δs = obs
    @unpack vvm = config
    @unpack bodies = vvm[i] #i-th ensemble member
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    # states_to_vortices!(vvm[i],x) #i-th ensemble member
    # vvm[i].bodies[1].Γ = -sum(vvm[i].vortices.Γ)
    # #solution at the current time step n
    # vmn = deepcopy(vvm[i])
    # construct_intermediate_model!(intermediate_vm,vmn)
    # soln = solve(intermediate_vm)
    # γn = soln.f./Δs;

    states_to_vortices!(vvm[i],x) #i-th ensemble member
    vvm[i].bodies[1].Γ = -sum(vvm[i].vortices.Γ[1:end-2])
    vmn = deepcopy(vvm[i])
    vmn.vortices.Γ[end-1] = x[end]*Δt
    subtractcirculation!(vmn.bodies, [vmn.vortices.Γ[end-1]])
    soln = solve(vmn)
    vmn.vortices.Γ[end] = soln.δΓ_vec[1]
    subtractcirculation!(vmn.bodies, soln.δΓ_vec)
    soln = solve(vmn)
    γn = soln.f./Δs

    #solution at the next time step n+1
    vm1 = deepcopy(vmn)
    # states_to_vortices!(vm1,x,Δt)
    # subtractcirculation!(vm1.bodies, [vm1.vortices.Γ[end-1]])
    # # solving for v_TE
    # solnp1 = solve(vm1)
    # vm1.vortices.Γ[end] = solnp1.δΓ_vec[1]
    # subtractcirculation!(vm1.bodies, solnp1.δΓ_vec)
    # advecting all vortices
    advect_vortices!(vm1,soln,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm1.vortices[end-1:end])
    pushvortices!(vm1,vLEnew,vTEnew)
    vm1.vortices.Γ[end-1] = x[end]*Δt
    subtractcirculation!(vmn.bodies, [vmn.vortices.Γ[end-1]])
    solnp1 = solve(vm1)
    vm1.vortices.Γ[end] = solnp1.δΓ_vec[1]
    subtractcirculation!(vm1.bodies, solnp1.δΓ_vec)
    solnp1 = solve(vm1)
    γnp1 = solnp1.f./Δs

    velocity!(obs.v̄,soln.ψ,vmn.ilsys)
    surface_velocity!(obs.v̄s,obs.v̄,vmn.ilsys)

    pressurejump!(obs.dp,γn,γnp1,obs.v̄s,Δt,vmn.ilsys)
    pressure!(obs.p̄,obs.v̄,obs.dp,vmn.ilsys)
    p⁺, p⁻ = sided_pressures(obs.p̄,obs.dp,vmn.ilsys)

    dp_sens = surface_interpolation(obs.dp,pfb,sens)

    return dp_sens, obs.p̄, p⁺, p⁻
end

