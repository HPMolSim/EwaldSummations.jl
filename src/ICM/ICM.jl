function IcmEwald2DInteraction(n_atoms::Int, s::T, α::T, γ::Tuple{T, T}, L::NTuple{3,T}, N_image::Int; ϵ::T = one(T)) where{T}
    r_c = s / α
    k_c = 2 * α * s
    k_set = k_set_2D(k_c, L)
    return IcmEwald2DInteraction(n_atoms, ϵ, α, r_c, k_c, γ, L, N_image, k_set)
end

function ICM_Ewald2D_energy(interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}

    ref_pos, ref_q = ICM_reflect(interaction.γ, interaction.L, interaction.N_image, position, q)

    Es = ICM_Ewald_short_energy(interaction, ref_pos, ref_q)
    El = ICM_Ewald2D_long_energy(interaction, ref_pos, ref_q)

    @debug "ICM Ewald2D, El = $El, Es = $Es"

    return Es + El
end

function ICM_Ewald2D_force(interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}

    ref_pos, ref_q = ICM_reflect(interaction.γ, interaction.L, interaction.N_image, position, q)

    force_short = ICM_Ewald_short_force(interaction, ref_pos, ref_q)
    force_long = ICM_Ewald2D_long_force(interaction, ref_pos, ref_q)

    return force_short + force_long
end

function IcmEwald3DInteraction(n_atoms::Int, s::T, α::T, γ::Tuple{T, T}, L::NTuple{3,T}, N_image::Int, N_pad::TP; ϵ::T = one(T)) where{T, TP}
    r_c = s / α
    k_c = 2 * α * s
    k_set = k_set_3D(k_c, (L[1], L[2], (2 * N_pad + 1) * L[3]))
    return IcmEwald3DInteraction(n_atoms, ϵ, α, r_c, k_c, γ, L, N_image, N_pad, k_set)
end

function ICM_Ewald3D_energy(interaction::IcmEwald3DInteraction{T, TP}, position::Vector{Point{3, T}}, q::Vector{T}) where{T, TP}

    ref_pos, ref_q = ICM_reflect(interaction.γ, interaction.L, interaction.N_image, position, q)

    Es = ICM_Ewald_short_energy(interaction, ref_pos, ref_q)
    El = ICM_Ewald3D_long_energy(interaction, ref_pos, ref_q)
    E_slab = ICM_Ewald3D_long_energy_slab(interaction, ref_pos, ref_q)

    @debug "ICM Ewald3D ELC, El = $El, Es = $Es, E_slab = $E_slab"

    return Es + El + E_slab
end

function ICM_Ewald3D_force(interaction::IcmEwald3DInteraction{T, TP}, position::Vector{Point{3, T}}, q::Vector{T}) where{T, TP}

    ref_pos, ref_q = ICM_reflect(interaction.γ, interaction.L, interaction.N_image, position, q)

    force_short = ICM_Ewald_short_force(interaction, ref_pos, ref_q)
    force_long = ICM_Ewald3D_long_force(interaction, ref_pos, ref_q)
    force_slab = ICM_Ewald3D_long_force_slab(interaction, ref_pos, ref_q)

    return force_short + force_long + force_slab
end