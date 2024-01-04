export vortices_to_states!, states_to_vortices!


function vortices_to_states!(x::AbstractVector,vm::VortexModel)
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        x[3i-2:3i] .= (vortex.x, vortex.y, vortex.Γ)
    end
end

function states_to_vortices!(vm::VortexModel,x::AbstractVector)
    @inbounds for (i, vortex) in enumerate(vm.vortices)
        [vortex.x, vortex.y, vortex.Γ] .= x[3i-2:3i]
    end
end