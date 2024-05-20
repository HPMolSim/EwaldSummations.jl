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


function ICMEwald3D_long_energy(interaction::IcmEwald3DInteraction{T}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T}

    energy = Atomic{T}(zero(T))

    @threads for K in interaction.k_set
        t = ICMEwald3D_long_energy_k(K, interaction, ref_position, ref_charge)
        atomic_add!(energy, t)
    end

    return energy[] / (4π * interaction.ϵ)
end

function ICMEwald3D_long_energy_k(K::Tuple{T, T, T, T}, interaction::IcmEwald3DInteraction{T}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T}
    k_x, k_y, k_z, k = K
    L_x, L_y, L_z = interaction.L
    α = interaction.α
    N_pad = interaction.N_pad

    ρ_k1 = zero(Complex{T})
    for j in 1:interaction.n_atoms
        x_j, y_j, z_j = ref_position[j]
        t = ref_charge[j] * exp(1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
        ρ_k1 += t
    end
    
    ρ_k2 = zero(Complex{T})
    for j in 1:length(ref_position)
        x_j, y_j, z_j = ref_position[j]
        t = ref_charge[j] * exp( - 1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
        ρ_k2 += t
    end

    sum_k = 2π / (L_x * L_y * (2 * N_pad + 1) * L_z) * T(real(ρ_k2[] * ρ_k1[])) * exp(- k^2 / (4 * α^2)) / k^2

    return sum_k
end

function ICMEwald3D_long_energy_slab(interaction::IcmEwald3DInteraction{T}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T}

    L_x, L_y, L_z = interaction.L
    energy = zero(T)
    N_pad = interaction.N_pad

    for i in 1:interaction.n_atoms
        x_i, y_i, z_i = ref_position[i]
        q_i = ref_charge[i]
        t = zero(T)
        for j in 1:length(ref_position)
            x_j, y_j, z_j = ref_position[j]
            q_j = ref_charge[j]
            t += q_j * (z_i - z_j)^2
        end
        energy += - π / (L_x * L_y * (2 * N_pad + 1) * L_z) * q_i * t
    end

    return energy / (4π * interaction.ϵ)
end