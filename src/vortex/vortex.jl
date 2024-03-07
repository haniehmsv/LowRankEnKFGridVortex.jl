export vortices_to_states!, states_to_vortices!


function vortices_to_states!(x::AbstractVector,vm::VortexModel{Nb,Ne}) where {Nb,Ne}
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        x[3i-2:3i] .= (vortex.x, vortex.y, vortex.Γ)
    end
end

function states_to_vortices!(vm::VortexModel{Nb,Ne},x::AbstractVector,Δt) where {Nb,Ne}
    vm.vortices.x .= x[1:3:end]
    vm.vortices.y .= x[2:3:end]
    vm.vortices.Γ .= x[3:3:end]
    @inbounds for k=1:Ne
        vm.vortices.Γ[end-Ne+k] = x[end-Ne+k]*Δt
    end
end

function states_to_vortices!(vm::VortexModel{Nb,Ne},x::AbstractVector) where {Nb,Ne}
    vm.vortices.x .= x[1:3:end]
    vm.vortices.y .= x[2:3:end]
    vm.vortices.Γ .= x[3:3:end]
end