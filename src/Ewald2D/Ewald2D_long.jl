function Ewald2D_long_energy(interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}
    energy = zero(T)
    for i in 1:interaction.n_atoms
        energy += Ewald2D_long_energy_k0(i, interaction, position, charge)
        for K in interaction.k_set
            energy += Ewald2D_long_energy_k(i, K, interaction, position, charge)
        end
    end
    return energy / interaction.ϵ
end

function Ewald2D_long_energy_k(i::Int, K::Tuple{T, T, T}, interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    k_x, k_y, k = K
    α = interaction.α
    sum_ki = zero(T)
    for j in 1:interaction.n_atoms
        x_ij, y_ij, z_ij = position[i] - position[j]
        sum_ki += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij)) / (8 * interaction.L[1] * interaction.L[2] * k)
    end
    return sum_ki
end

function Ewald2D_long_energy_k0(i::Int, interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    α = interaction.α
    sum_k0i = zero(T)
    for j in 1:interaction.n_atoms
        x_ij, y_ij, z_ij = position[i] - position[j]
        sum_k0i += - q[i] * q[j] * (1 / (α * sqrt(π)) * exp(-(α * z_ij)^2) + z_ij * erf(α * z_ij)) / (4 *  interaction.L[1] * interaction.L[2])
    end
    return sum_k0i
end

function Ewald2D_long_force(interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T<:Number}

    n_atoms = interaction.n_atoms
    
    force_long = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]

    Ewald2D_long_force_k0!(interaction, position, charge, force_long)

    for K in interaction.k_set
        Ewald2D_long_force_k!(K, interaction, position, charge, force_long)
    end

    return force_long
end

function Ewald2D_long_force_k0!(interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force_long::Vector{Point{3, T}}) where {T<:Number}
    α = interaction.α
    for i in 1:interaction.n_atoms
        for j in 1:interaction.n_atoms
            z_ij = position[i][3] - position[j][3]
            force_long[i] += Point(zero(T), zero(T), charge[i] * charge[j] * (erf(α * z_ij))  / (2 * interaction.L[1] * interaction.L[2] * interaction.ϵ))
        end
    end
    return nothing
end

function Ewald2D_long_force_k!(K::Tuple{T, T, T}, interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force_long::Vector{Point{3, T}}) where {T<:Number}
    n_atoms = interaction.n_atoms
    α = interaction.α
    k_x, k_y, k = K

    for i in 1:n_atoms
        sum_x = zero(T)
        sum_y = zero(T)
        sum_z = zero(T)
        for j in 1:n_atoms
            x_ij, y_ij, z_ij = position[i] - position[j]

            sum_xy = - charge[i] * charge[j] * sin(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij))
            sum_x += k_x * sum_xy / k
            sum_y += k_y * sum_xy / k

            sum_z += charge[i] * charge[j] * cos(k_x * x_ij + k_y * y_ij) * (
                k * exp(k * z_ij) * erfc(k / (2α) + α * z_ij) - 
                k * exp( - k * z_ij) * erfc(k / (2α) - α * z_ij) -
                2α / sqrt(π) * exp(k * z_ij) * exp(-(k / (2α) + α * z_ij)^2) +
                2α / sqrt(π) * exp(- k * z_ij) * exp(-(k / (2α) - α * z_ij)^2) ) / k
        end
        force_long[i] += - Point(sum_x, sum_y, sum_z) / (4 * interaction.L[1] * interaction.L[2] * interaction.ϵ)
    end
    return nothing
end