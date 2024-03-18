export vortices_to_states!, states_to_vortices!


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