export setup_sensors
  
  """Setting up sensors only for body::Plate for now"""
  function setup_sensors(pfb::PotentialFlowBody,Nsens)
    @unpack points, U, Ω = pfb
  
    xsens = range(points.x[1],points.x[end],length=Nsens)
    ysens = range(points.y[1],points.y[end],length=Nsens)
  
    return Sensor(vcat(xsens),vcat(ysens),Nsens)
  end