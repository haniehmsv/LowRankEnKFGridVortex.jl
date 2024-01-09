export sens, setup_sensors

struct sens
    x :: Vector{Float64}
    y :: Vector{Float64}
    Nsens :: Int64
end
  
  """Setting up sensors only for body::Plate for now"""
  function setup_sensors(pfb::PotentialFlowBody,Nsens)
    @unpack points, U, Ω = pfb
  
    xsens = range(points.x[1],points.x[end],length=Nsens)
    ysens = range(points.y[1],points.y[end],length=Nsens)
  
    return sens(vcat(xsens),vcat(ysens),Nsens)
  end