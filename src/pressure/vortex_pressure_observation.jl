
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
    intermediate_vm :: VortexModel
    Δs::DST
end

function VortexPressure(sens::Sensor,config::VortexForecast)
    @unpack vvm= config
    vm = vvm[1]
    intermediate_bodies = deepcopy(vm.bodies)
    Nb = length(vm.bodies)
    for j=1:Nb
        intermediate_bodies[j].edges = Int64[]
    end
    intermediate_vm = VortexModel(vm.g,vortices=[vm.vortices...],bodies=intermediate_bodies,U∞=vm.U∞)
    withfreestream = vm.U∞ == 0.0 ? false : true
    Nv = config.Nv
    Nx = 3*Nv
    Ny = length(sens.x)
    Δs = dlengthmid(vm.bodies[1].points)
    return VortexPressure{Ny,withfreestream,typeof(Δs)}(sens,config,intermediate_vm,Δs)
end


function observations(x::AbstractVector,t,Δt,obs::VortexPressure,i::Int64)
    @unpack sens, config, intermediate_vm, Δs = obs
    @unpack vvm = config
    @unpack bodies = vvm[i] #i-th ensemble member
    #for 1 body for now
    pfb = bodies[1]
    @unpack points = pfb

    states_to_vortices!(vvm[i],x) #i-th ensemble member
    vvm[i].bodies[1].Γ = -sum(vvm[i].vortices.Γ)
    #solution at the current time step n
    vmn = deepcopy(vvm[i])
    construct_intermediate_model!(intermediate_vm,vmn)
    soln = solve(intermediate_vm)
    γn = soln.f./Δs;

    #solution at the next time step n+1
    vm1 = deepcopy(vvm[i])
    # setting the value for v_LE
    states_to_vortices!(vm1,x,Δt)
    subtractcirculation!(vm1.bodies, vm1.vortices.Γ[end-1])
    # solving for v_TE
    solnp1 = solve(vm1)
    vm1.vortices.Γ[end] = solnp1.δΓ_vec[1]
    subtractcirculation!(vm.bodies, sol.δΓ_vec)
    # advecting all vortices
    Xv = getvortexpositions(vm1)
    Ẋv = deepcopy(Xv)
    vortexvelocities!(Ẋv, vm1, solnp1.ψ)
    Xv .= Xv .+ Ẋv*Δt
    setvortexpositions!(vm1, Xv)
    solnp1 = solve(vm1)
    γnp1 = solnp1.f./Δs

    v̄ = Edges(Primal,soln.ψ)
    Xs = VectorData(collect(vmn.bodies))
    v̄s = VectorData(Xs)
    velocity!(v̄,soln.ψ,vmn.ilsys)
    surface_velocity!(v̄s,v̄,vmn.ilsys)

    dp = ScalarData(Xs)
    pressurejump!(dp,γn,γnp1,v̄s,Δt,vmn.ilsys)
    p̄ = Nodes(Primal,soln.ψ)
    pressure!(p̄,v̄,dp,vmn.ilsys)
    p⁺, p⁻ = sided_pressures(p̄,dp,vmn.ilsys)

    dp_sens = surface_interpolation(dp,pfb,sens)

    return dp_sens, p̄, p⁺, p⁻
end

