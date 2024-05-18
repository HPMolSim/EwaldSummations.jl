function ICMEwald2D_long_energy(interaction::IcmEwald2DInteraction{T}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T}
    energy = Atomic{T}(zero(T))
    @threads for i in 1:interaction.n_atoms
        t = zero(T)
        t += ICMEwald2D_long_energy_k0(i, interaction, ref_position, ref_charge)
        for K in interaction.k_set
            t += ICMEwald2D_long_energy_k(i, K, interaction, ref_position, ref_charge)
        end
        atomic_add!(energy, t)
    end
    return energy[] / interaction.ϵ
end

function ICMEwald2D_long_energy_k(i::Int, K::Tuple{T, T, T}, interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    k_x, k_y, k = K
    α = interaction.α
    L = interaction.L
    t = zero(T)
    for j in 1:length(position)
        x_ij, y_ij, z_ij = position[i] - position[j]
        if !(erfc(k / (2α) - α * z_ij) ≈ zero(T))
            t += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * z_ij) * erfc(k / (2α) - α * z_ij) / (8 * L[1] * L[2] * k)
        end

        if !(erfc(k / (2α) + α * z_ij) ≈ zero(T))
            t += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp(k * z_ij) * erfc(k / (2α) + α * z_ij) / (8 * L[1] * L[2] * k)
        end
    end
    return t
end

function ICMEwald2D_long_energy_k0(i::Int, interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    α = interaction.α
    L = interaction.L
    t = zero(T)
    for j in 1:length(position)
        x_ij, y_ij, z_ij = position[i] - position[j]
        t += - q[i] * q[j] * (1 / (α * sqrt(π)) * exp(-(α * z_ij)^2) + z_ij * erf(α * z_ij)) / (4 *  L[1] * L[2])
    end
    return t
end
