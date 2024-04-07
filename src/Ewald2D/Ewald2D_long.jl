function Ewald2D_long_energy(interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}
    energy = Atomic{T}(zero(T))
    @threads for i in 1:interaction.n_atoms
        t = zero(T)
        t += Ewald2D_long_energy_k0(i, interaction, position, charge)
        for K in interaction.k_set
            t += Ewald2D_long_energy_k(i, K, interaction, position, charge)
        end
        atomic_add!(energy, t)
    end
    return energy[] / interaction.ϵ
end

function Ewald2D_long_energy_k(i::Int, K::Tuple{T, T, T}, interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    k_x, k_y, k = K
    α = interaction.α
    L = interaction.L
    t = zero(T)
    for j in 1:interaction.n_atoms
        x_ij, y_ij, z_ij = position[i] - position[j]
        t += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij)) / (8 * L[1] * L[2] * k)
    end
    return t
end

function Ewald2D_long_energy_k0(i::Int, interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    α = interaction.α
    L = interaction.L
    t = zero(T)
    for j in 1:interaction.n_atoms
        x_ij, y_ij, z_ij = position[i] - position[j]
        t += - q[i] * q[j] * (1 / (α * sqrt(π)) * exp(-(α * z_ij)^2) + z_ij * erf(α * z_ij)) / (4 *  L[1] * L[2])
    end
    return t
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
        force_i = @distributed (+) for j in 1:interaction.n_atoms
            z_ij = position[i][3] - position[j][3]
            Point(zero(T), zero(T), charge[i] * charge[j] * (erf(α * z_ij))  / (2 * interaction.L[1] * interaction.L[2] * interaction.ϵ))
        end
        force_long[i] += force_i
    end
    return nothing
end

function Ewald2D_long_force_k!(K::Tuple{T, T, T}, interaction::Ewald2DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force_long::Vector{Point{3, T}}) where {T<:Number}
    n_atoms = interaction.n_atoms
    α = interaction.α
    k_x, k_y, k = K

    for i in 1:n_atoms
        s = @distributed (+) for j in 1:n_atoms
            x_ij, y_ij, z_ij = position[i] - position[j]

            sum_xy = - charge[i] * charge[j] * sin(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij))
            sx = k_x * sum_xy / k
            sy = k_y * sum_xy / k

            sz = charge[i] * charge[j] * cos(k_x * x_ij + k_y * y_ij) * (
                k * exp(k * z_ij) * erfc(k / (2α) + α * z_ij) - 
                k * exp( - k * z_ij) * erfc(k / (2α) - α * z_ij) -
                2α / sqrt(π) * exp(k * z_ij) * exp(-(k / (2α) + α * z_ij)^2) +
                2α / sqrt(π) * exp(- k * z_ij) * exp(-(k / (2α) - α * z_ij)^2) ) / k
            
            Point(sx, sy, sz)
        end
        force_long[i] += - s / (4 * interaction.L[1] * interaction.L[2] * interaction.ϵ)
    end
    return nothing
end