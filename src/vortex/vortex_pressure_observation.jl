export observation



mutable struct VortexPressure{ST} <: AbstractObservationOperator{Nx,Ny,true}
    sens::ST
    config::VortexForecast
end

function VortexPressure(sens::AbstractVector,config::VortexForecast)
    
end