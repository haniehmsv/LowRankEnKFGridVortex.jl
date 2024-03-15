# TODO: setup_sensors(body::Ellipse,Nsens) now works only for cicles, modify it for ellipses
export setup_sensors, surface_interpolation
  
  """Setting up sensors for body::Plate"""
function setup_sensors(body::Polygon,Nsens)
    
    xsens = range(body.x[1],body.x[end],length=Nsens)
    ysens = range(body.y[1],body.y[end],length=Nsens)
  
    return Sensor(collect(xsens),collect(ysens),Nsens)
end

"""Setting up sensors for body::Circle"""
function setup_sensors(body::Ellipse,Nsens)
  ds = collect(range(0,2π,Nsens+1))
  pop!(ds)
  xsens = body.cent[1] .+ body.a .* cos.(ds)
  ysens = body.cent[2] .+ body.b .* sin.(ds)

  return Sensor(xsens,ysens,Nsens)
end

function surface_interpolation(p::ScalarData,pfb::PotentialFlowBody,sens::Sensor)
    xsens = deepcopy(sens.x)
    ysens = deepcopy(sens.y)
    xs, ys = deepcopy(RigidBodyTools.collect(pfb.points))
    x_cent, y_cent = pfb.points.cent

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