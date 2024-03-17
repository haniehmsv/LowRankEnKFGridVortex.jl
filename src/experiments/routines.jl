# TODO: setup_sensors(body::Ellipse,Nsens) now works only for cicles, modify it for ellipses

import BasicInterpolators: RBFInterpolator

export setup_sensors, surface_interpolation

function setup_sensors(pfb::PotentialFlowBody,Nsens)
  return _setup_sensors(pfb.points,Nsens)
end

  """Setting up sensors for body::Plate"""
function _setup_sensors(body::Polygon,Nsens)
    
    xsens = range(body.x[1],body.x[end],length=Nsens)
    ysens = range(body.y[1],body.y[end],length=Nsens)
  
    return Sensor(collect(xsens),collect(ysens),Nsens)
end

"""Setting up sensors for body::Circle"""
function _setup_sensors(body::Ellipse,Nsens)
  ds = collect(range(0,2π,Nsens+1))
  pop!(ds)
  xsens = body.cent[1] .+ body.a .* cos.(ds)
  ysens = body.cent[2] .+ body.b .* sin.(ds)

  return Sensor(xsens,ysens,Nsens)
end

"""Interpolate scalar data p defined on the body points to the sensor locations"""
function surface_interpolation(p::ScalarData,pfb::PotentialFlowBody,sens::Sensor)
    return _surface_interpolation(p,pfb.points,sen)
end

function _surface_interpolation(p::ScalarData,body::Polygon,sens::Sensor)
  xsens = deepcopy(sens.x)
  ysens = deepcopy(sens.y)
  xs, ys = deepcopy(RigidBodyTools.collect(body))
  x_cent, y_cent = body.cent

  @. xs -= x_cent
  @. ys -= y_cent
  rs = (xs .> 0.0) .* sqrt.(xs.^2 .+ ys.^2) .+ (xs .<= 0.0) .* -sqrt.(xs.^2 .+ ys.^2)
  itp = linear_interpolation(rs, p.data)

  @. xsens -= x_cent
  @. ysens -= y_cent
  rsens = (xsens .> 0.0) .* sqrt.(xsens.^2 .+ ysens.^2) .+ (xsens .<= 0.0) .* -sqrt.(xsens.^2 .+ ysens.^2)
  psens = itp(rsens)
  return psens
end

function _surface_interpolation(p::ScalarData,body::Ellipse,sens::Sensor;ε=0.1)
  Xp = hcat(body.x, body.y)
  itp = RBFInterpolator(Xp, p.data, ε)
  psens = itp.(sens.x,sens.y)
  return psens
end