
export observations, VortexPressure, Sensor

#### OBSERVATION OPERATORS ####

struct Sensor
  x :: Vector{Float64}
  y :: Vector{Float64}
  Nsens :: Int64
end

mutable struct VortexPressure{Ny,withfreestream,DST} <: AbstractObservationOperator{Ny,true}
    sens::Sensor
    config::VortexForecast
    Δs::DST
end

function VortexPressure(sens::Sensor,config::VortexForecast)
    @unpack vvm, pfb = config
    vm = vvm[1]
    withfreestream = vm.U∞ == 0.0 ? false : true
    Nv = config.Nv
    Nx = 3*Nv
    Ny = length(sens.x)
    Δs = dlengthmid(pfb.points)
    return VortexPressure{Ny,withfreestream,typeof(Δs)}(sens,config,Δs)
end


function observations(x::AbstractVector,t,Δt,obs::VortexPressure,i::Int64)
    @unpack sens, config, Δs = obs
    @unpack vvm = config
    @unpack bodies = vvm[i] #i-th ensemble member
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    states_to_vortices!(vvm[i],x) #i-th ensemble member

    #solution at the current time step n
    vmn = deepcopy(vvm[i])
    vm1 = deepcopy(vvm[i])
    soln = solve(vmn)
    γn = soln.f./Δs;

    #advance the solution to the next time step by advection and release of new vortices
    time_advancement!(vm1,Δt)
    vLEnew, vTEnew = createsheddedvortices(points,vm1.vortices[end-1:end])
    pushvortices!(vm1,vLEnew,vTEnew)
    solnp1 = solve(vm1)
    γnp1 = solnp1.f./Δs

    v̄ = Edges(Primal,soln.ψ)
    Xs = VectorData(collect(vmn.bodies))
    v̄s = VectorData(Xs)
    velocity!(v̄,soln.ψ,vmn.ilsys)
    surface_velocity!(v̄s,v̄,vmn.ilsys)

    dp = ScalarData(Xs)
    pressurejump!(dp,γn,γnp1,v̄s,Δt,vmn.ilsys)
    # p̄ = Nodes(Primal,soln.ψ)
    # pressure!(p̄,v̄,dp,vmn.ilsys)
    # p⁺, p⁻ = sided_pressures(p̄,dp,vmn.ilsys)

    dp_sens = surface_interpolation(dp,pfb,sens)

    return dp_sens
end