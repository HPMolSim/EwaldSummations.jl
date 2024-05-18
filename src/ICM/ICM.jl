function IcmEwald2DInteraction(n_atoms::Int, s::T, α::T, γ::Tuple{T, T}, L::NTuple{3,T}, N_image::Int; ϵ::T = one(T)) where{T}
    r_c = s / α
    k_c = 2 * α * s
    k_set = k_set_2D(k_c, L)
    return IcmEwald2DInteraction(n_atoms, ϵ, α, r_c, k_c, γ, L, N_image, k_set)
end

function ICMEwald2D_energy(interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}

    ref_pos, ref_q = ICM_reflect(interaction.γ, interaction.L, interaction.N_image, position, q)

    Es = ICMEwald2D_short_energy(interaction, ref_pos, ref_q)
    El = ICMEwald2D_long_energy(interaction, ref_pos, ref_q)

    return Es + El
end