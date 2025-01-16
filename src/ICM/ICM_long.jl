function ICM_Ewald2D_long_energy(interaction::IcmEwald2DInteraction{T}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T}
    energy = Atomic{T}(zero(T))
    @threads for i in 1:interaction.n_atoms
        t = zero(T)
        t += ICM_Ewald2D_long_energy_k0(i, interaction, ref_position, ref_charge)
        for K in interaction.k_set
            t += ICM_Ewald2D_long_energy_k(i, K, interaction, ref_position, ref_charge)
        end
        atomic_add!(energy, t)
    end
    return energy[] / interaction.ϵ
end

function ICM_Ewald2D_long_energy_k(i::Int, K::Tuple{T, T, T}, interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
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

function ICM_Ewald2D_long_energy_k0(i::Int, interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    α = interaction.α
    L = interaction.L
    t = zero(T)
    for j in 1:length(position)
        x_ij, y_ij, z_ij = position[i] - position[j]
        t += - q[i] * q[j] * (1 / (α * sqrt(π)) * exp(-(α * z_ij)^2) + z_ij * erf(α * z_ij)) / (4 *  L[1] * L[2])
    end
    return t
end


function ICM_Ewald2D_long_force(interaction::IcmEwald2DInteraction{T}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T}
    
    n_atoms = interaction.n_atoms
    force_long = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]

    ICM_Ewald2D_long_force_k0!(interaction, ref_position, ref_charge, force_long)
    for K in interaction.k_set
        ICM_Ewald2D_long_force_k!(K, interaction, ref_position, ref_charge, force_long)
    end

    return force_long
end

function ICM_Ewald2D_long_force_k!(K::Tuple{T, T, T}, interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force::Vector{Point{3, T}}) where{T}
    n_atoms = interaction.n_atoms
    α = interaction.α
    k_x, k_y, k = K

    for i in 1:n_atoms
        fxi = zero(T)
        fyi = zero(T)
        fzi = zero(T)
        for j in 1:length(charge)
            x_ij, y_ij, z_ij = position[i] - position[j]

            t1 = abs(k * z_ij) > 650.0 ? zero(T) : exp(k * z_ij)
            t2 = abs(k * z_ij) > 650.0 ? zero(T) : exp( - k * z_ij)

            sum_xy = - charge[i] * charge[j] * sin(k_x * x_ij + k_y * y_ij) * (t1 * erfc(k / (2α) + α * z_ij) + t2 * erfc(k / (2α) - α * z_ij))
            sx = k_x * sum_xy / k
            sy = k_y * sum_xy / k

            sz = charge[i] * charge[j] * cos(k_x * x_ij + k_y * y_ij) * (
                k * t1 * erfc(k / (2α) + α * z_ij) - 
                k * t2 * erfc(k / (2α) - α * z_ij) -
                2α / sqrt(π) * t1 * exp(-(k / (2α) + α * z_ij)^2) +
                2α / sqrt(π) * t2 * exp(-(k / (2α) - α * z_ij)^2) ) / k
            
            fxi += sx
            fyi += sy
            fzi += sz
        end
        s = Point(fxi, fyi, fzi)
        force[i] += - s / (4 * interaction.L[1] * interaction.L[2] * interaction.ϵ)
    end
    return force
end

function ICM_Ewald2D_long_force_k0!(interaction::IcmEwald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force::Vector{Point{3, T}}) where{T}
    α = interaction.α
    for i in 1:interaction.n_atoms
        f0i = Atomic{T}(zero(T))
        @threads for j in 1:length(charge)
            z_ij = position[i][3] - position[j][3]
            f0z = charge[i] * charge[j] * (erf(α * z_ij))  / (2 * interaction.L[1] * interaction.L[2] * interaction.ϵ)
            atomic_add!(f0i, f0z)
        end
        force[i] += Point(zero(T), zero(T), f0i[])
    end
    return force
end


function ICM_Ewald3D_long_energy(interaction::IcmEwald3DInteraction{T, TP}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T, TP}

    energy = Atomic{T}(zero(T))

    @threads for K in interaction.k_set
        t = ICM_Ewald3D_long_energy_k(K, interaction, ref_position, ref_charge)
        atomic_add!(energy, t)
    end

    return energy[] / (4π * interaction.ϵ)
end

function ICM_Ewald3D_long_energy_k(K::Tuple{T, T, T, T}, interaction::IcmEwald3DInteraction{T, TP}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T, TP}
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

function ICM_Ewald3D_long_force(interaction::IcmEwald3DInteraction{T, TP}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T, TP}
    
    n_atoms = interaction.n_atoms
    force_long = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]

    for K in interaction.k_set
        ICM_Ewald3D_long_force_k!(K, interaction, position, charge, force_long)
    end

    return force_long
end

function ICM_Ewald3D_long_force_k!(K::Tuple{T, T, T, T}, interaction::IcmEwald3DInteraction{T, TP}, position::Vector{Point{3, T}}, charge::Vector{T}, force_long::Vector{Point{3, T}}) where{T, TP}
    k_x, k_y, k_z, k = K
    L_x, L_y, L_z = interaction.L
    α = interaction.α
    N_pad = interaction.N_pad

    ρ_k = zero(Complex{T})
    for j in 1:length(charge)
        x_j, y_j, z_j = position[j]
        ρ_k += charge[j] * exp(1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
    end

    for i in 1:interaction.n_atoms
        x_i, y_i, z_i = position[i]
        force_long[i] += - charge[i] * T(imag(ρ_k * exp( - 1.0im * (k_x * x_i + k_y * y_i + k_z * z_i)))) * exp(- k^2 / (4 * α^2)) / k^2 /interaction.ϵ * Point(k_x, k_y, k_z) / (L_x * L_y * L_z * (2 * N_pad + 1))
    end

    return nothing
end

function ICM_Ewald3D_long_energy_slab(interaction::IcmEwald3DInteraction{T, TP}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T, TP}

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

function ICM_Ewald3D_long_force_slab(interaction::IcmEwald3DInteraction{T, TP}, ref_position::Vector{Point{3, T}}, ref_charge::Vector{T}) where{T, TP}

    L_x, L_y, L_z = interaction.L
    N_pad = interaction.N_pad

    force_slab = [Point(zero(T), zero(T), zero(T)) for _=1:interaction.n_atoms]

    for i in 1:interaction.n_atoms
        x_i, y_i, z_i = ref_position[i]
        q_i = ref_charge[i]
        t = zero(T)
        for j in 1:length(ref_position)
            x_j, y_j, z_j = ref_position[j]
            q_j = ref_charge[j]
            t += 4 * q_j * (z_i - z_j)
        end
        force_slab[i] += Point(zero(T), zero(T), π / (L_x * L_y * (2 * N_pad + 1) * L_z) * q_i * t)
    end

    return force_slab ./ (4π * interaction.ϵ)
end