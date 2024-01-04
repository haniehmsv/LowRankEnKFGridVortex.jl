
export observation, VortexPressure, setup_sensors

#### OBSERVATION OPERATORS ####


mutable struct VortexPressure{Nx,Ny,withfreestream,ST} <: AbstractObservationOperator{Nx,Ny,true}
    sens::ST
    config::VortexForecast
end

function VortexPressure(sens::AbstractVector,config::VortexForecast)
    @unpack vm = config
    withfreestream = vm.U∞ == 0.0 ? false : true
    Nv = config.Nv
    Nx = 3*Nv
    Ny = length(sens)

    return VortexPressure{Nx,Ny,withfreestream,typeof(sens)}(sens,config)
end


function observations!(Y::EnsembleMatrix{Ny,Ne},X::EnsembleMatrix{Nx,Ne},t,obs::VortexPressure) where {Nx,Ny,Ne}
    for j in 1:Ne
        Y(j) .= observations(X(j),t,obs)
    end
    return Y
end

function observations(x::AbstractVector,t,obs::VortexPressure) end



"""Setting up sensors"""
function setup_sensors(Nsens;layout=(:line,1.0))

    layout_type, len = layout
  
    if layout_type == :circle
      rsens = len
      θsens = range(0,2π,length=Nsens+1)
      sens = rsens*exp.(im*θsens[1:end-1])
    elseif layout_type == :line
      ϵsens = 0.0
      lowerrow = range(-len,len,length=Nsens) .+ (-0.5ϵsens .+ ϵsens*rand(Nsens))*im
      #upperrow = range(-2.0,2.0,length=Nsens) .+ 1.0*im
      #leftside = im*range(-1.0,3.0,length=Nsens) .- 1.0
      #rightside = im*range(-1.0,3.0,length=Nsens) .+ 1.0
      sens = vcat(lowerrow,)  #upperrow);
    elseif layout_type == :dline
      ϵsens = 0.02
      lowerrow1 = range(-len,len,length=Nsens÷2) .- 0.5*ϵsens
      lowerrow2 = range(-len,len,length=Nsens÷2) .+ 0.5*ϵsens
      sens = sort(vcat(lowerrow1,lowerrow2)) .+ 0.0im
    end
    return sens
  end