using ForwardDiff
FD = ForwardDiff
using CartesianGrids
using GridPotentialFlow

const TOL = 0.01

Δx = 0.01
Lx = 2.0
Δt = 2e-2
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx)
c = Lx/2    #chord length
U = 1   #freestream velocity
α = -π/3    #angle of attack
plate = Plate(c,2*cellsize(g))
Tr = RigidTransform((0.0,0.0),α)
update_body!(plate,Tr)
Δs = dlengthmid(plate);
pfb = PotentialFlowBody(plate,edges=[length(plate)])

#initial point vortices
vLE = Vortex(plate.x[1]+3Δt*U*cos(plate.α+π/2),plate.y[1]+3Δt*U*sin(plate.α+π/2),0.0,DT=Real)
vTE = Vortex(plate.x[end]+3Δt*U*cos(α+π/2),plate.y[end]+3Δt*U*sin(α+π/2),0.0,DT=Real)

t = 0
state = zeros(7,1)
Γ̇1= -0.751121/Δt    #vorticity flux at the LE
state[:,1] = [vLE.x, vLE.y, vLE.Γ, vTE.x, vTE.y, vTE.Γ, Γ̇1]
X = BasicEnsembleMatrix(state)
vvm = [VortexModel(g,bodies=deepcopy([pfb]),vortices=[vLE,vTE],U∞=(1.0,0.0)) for i in 1:size(state,2)]
vforec = VortexForecast(vvm)
Nsens = length(plate)
sens = setup_sensors(pfb,Nsens)
obs = VortexPressure(sens,vforec);

@testset "Checking Jacobian of pressure sensors wrt states" begin
    # creating state vector xd containing FD.Dual numbers
    x = X(1)
    cfg = FD.JacobianConfig(forecast, x)
    xd = cfg.duals
    FD.seed!(xd, x)
    seeds = cfg.seeds
    FD.seed!(xd, x, 1, seeds)

    # calculate dp, the output is of type FD.Dual
    dp_FD = observations(xd,t,Δt,obs,1)

    # Finite Difference
    dp_Findiff = zeros(Real,Nsens,length(X(1)))
    epsil = 10^(-6)
    for i in 1:length(X(1))
        vvm = [VortexModel(g,bodies=deepcopy([pfb]),vortices=[vLE,vTE],U∞=(1.0,0.0)) for i in 1:size(state,2)]
        vforec = VortexForecast(vvm)
        obs = VortexPressure(sens,vforec)
        x = deepcopy(X(1))
        dp_Findiff[:,i] .= observations(x,t,Δt,obs,1)   # dp(x)
        x[i] += epsil
        dp_Findiff[:,i] .-= observations(x,t,Δt,obs,1)  # dp(x+ε)
        dp_Findiff[:,i] ./= -epsil                      # [dp(x+ε)-dp(x)]/ε
        if i in [1,2,4,5]
            @test isapprox(dp_Findiff[:,i],FD.partials.(dp_FD,i),atol=TOL)
        end
    end
    vvm = [VortexModel(g,bodies=deepcopy([pfb]),vortices=[vLE,vTE],U∞=(1.0,0.0)) for i in 1:size(state,2)]
    vforec = VortexForecast(vvm)
    obs = VortexPressure(sens,vforec)
    dp_val = observations(X(1),t,Δt,obs,1)

    @test FD.value.(dp_FD) ≈ dp_val

end
