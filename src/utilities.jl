import CartesianGrids: Regularize

export Regularize

function CartesianGrids.Regularize(x::AbstractVector{T},y::AbstractVector{T},dx::U;
    ddftype::DataType=Yang3,graddir::Int=0,
    I0::Tuple{Int,Int}=(1,1),
    weights::Union{U,Vector{U}}=1.0,
    filter::Bool = false,
    issymmetric::Bool = false) where {T<:Real,U<:Real}

    _issymmetric = (filter ? false : issymmetric)

    n = length(x)
    @assert length(y)==n
    if !_issymmetric
    if typeof(weights) == T
    wtvec = similar(x)
    fill!(wtvec,weights/(dx*dx))
    else
    @assert length(weights)==n
    wtvec = deepcopy(weights)./(dx*dx)
    end
    else
    # if the regularization and interpolation are symmetric, then the
    # weights are automatically set to be the cell area in order to cancel it
    # in the denominator of the regularization operator.
    wtvec = similar(x)
    fill!(wtvec,1.0)
    end

    baseddf = DDF(ddftype=ddftype,dx=1.0)
    if graddir == 0
    ddf = baseddf
    else
    ddf = GradDDF(graddir,ddftype=ddftype,dx=1.0)
    end

    Regularize{length(x),filter}(x./dx.+I0[1],y./dx.+I0[2],1.0/(dx*dx),
      wtvec,ddf,_get_regularization_radius(baseddf),_issymmetric)
end