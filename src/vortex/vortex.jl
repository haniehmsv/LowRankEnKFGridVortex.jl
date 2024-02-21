export vortices_to_states!, states_to_vortices!


function vortices_to_states!(x::AbstractVector,vm::VortexModel{Nb,Ne},xold::AbstractVector,Δt) where {Nb,Ne}
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        x[3i-2:3i] .= (vortex.x, vortex.y, vortex.Γ)
    end
    x[end] = xold[end]
end

function states_to_vortices!(vm::VortexModel{Nb,Ne},x::AbstractVector,Δt) where {Nb,Ne}
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        [vortex.x, vortex.y, vortex.Γ] .= x[3i-2:3i]
    end
    vm.vortices.Γ[end-1] = x[end]*Δt
end

function states_to_vortices!(vm::VortexModel{Nb,Ne},x::AbstractVector) where {Nb,Ne}
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        [vortex.x, vortex.y, vortex.Γ] .= x[3i-2:3i]
    end
end