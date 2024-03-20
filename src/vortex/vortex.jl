export vortices_to_states!, states_to_vortices!, create_covariance


function vortices_to_states!(x::AbstractVector,vm::VortexModel{Nb,Ne},xold::AbstractVector) where {Nb,Ne}
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        x[3i-2:3i] .= (vortex.x, vortex.y, vortex.Γ)
    end
    x[end] = xold[end]
end

function vortices_to_states!(x::AbstractVector,vm::VortexModel{Nb,Ne}) where {Nb,Ne}
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        x[3i-2:3i] .= (vortex.x, vortex.y, vortex.Γ)
    end
end

function states_to_vortices!(vm::VortexModel{Nb,Ne},x::AbstractVector,Δt) where {Nb,Ne}
    vm.vortices.x .= x[1:3:end-1]
    vm.vortices.y .= x[2:3:end-1]
    vm.vortices.Γ .= x[3:3:end-1]
    vm.vortices.Γ[end-1] = x[end]*Δt
end

function states_to_vortices!(vm::VortexModel{Nb,Ne},x::AbstractVector) where {Nb,Ne}
    vm.vortices.x .= x[1:3:end-Ne]
    vm.vortices.y .= x[2:3:end-Ne]
    vm.vortices.Γ .= x[3:3:end-Ne]
end

"""
create_covariance(var::Real, N::Int) -> Matrix with dimensions (N*N)
creates a covariance matrix with variance var and dimention N
"""
function create_covariance(var::Real, N::Int)
  return Diagonal(var^2*ones(N))
end

"""
    state_covariance(varx, vary, varΓ, config::VortexForecast)

Create a state covariance matrix with variances `varx`, `vary` and `varΓ`
for the x, y, and strength entries for every vortex.
"""
function state_covariance(varx, vary, varΓ, config::VortexForecast)
  @unpack Nv = config

  Σx_diag = zeros(Float64,state_length(config))
  for j = 1:Nv
    Σx_diag[3i-2] = varx
    Σx_diag[3i-1] = vary
    Σx_diag[3i] = varΓ
  end

  return Diagonal(Σx_diag)
end