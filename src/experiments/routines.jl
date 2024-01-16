export setup_sensors, surface_interpolation
  
  """Setting up sensors only for body::Plate for now"""
  function setup_sensors(pfb::PotentialFlowBody,Nsens)
    @unpack points, U, Ω = pfb
  
    xsens = range(points.x[1],points.x[end],length=Nsens)
    ysens = range(points.y[1],points.y[end],length=Nsens)
  
    return Sensor(vcat(xsens),vcat(ysens),Nsens)
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