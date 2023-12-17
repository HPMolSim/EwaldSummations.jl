function Ewald3D_long_energy(interaction::Ewald3DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}
    energy = zero(T)

    energy += Ewald3D_long_energy_k0(interaction, position, charge)

    for K in interaction.k_set
        energy += Ewald3D_long_energy_k(K, interaction, position, charge)
    end
    return energy / (4π * interaction.ϵ)
end

function Ewald3D_long_energy_k(K::Tuple{T, T, T, T}, interaction::Ewald3DInteraction{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    k_x, k_y, k_z, k = K
    L_x, L_y, L_z = interaction.L
    α = interaction.α

    ρ_k = zero(ComplexF64)
    for j in 1:interaction.n_atoms
        x_j, y_j, z_j = position[j]
        ρ_k += q[j] * exp(1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
    end

    sum_k = 2π / (L_x * L_y * L_z) * T(ρ_k' * ρ_k) * exp(- k^2 / (4 * α^2)) / k^2

    return sum_k
end

function Ewald3D_long_energy_k0(interaction::Ewald3DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}

    L_x, L_y, L_z = interaction.L

    P = Point(zero(T), zero(T), zero(T))
    for i in 1:interaction.n_atoms
        P += charge[i] * position[i]
    end

    energy_k0 = 2π * dist2(P, Point(zero(T), zero(T), zero(T))) / (2 * interaction.ϵ_inf + one(T)) / (L_x * L_y * L_z)

    return energy_k0
end

function Ewald3D_long_force(interaction::Ewald3DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}
    n_atoms = interaction.n_atoms
    
    force_long = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]

    Ewald3D_long_force_k0!(interaction, position, charge, force_long)

    for K in interaction.k_set
        Ewald3D_long_force_k!(K, interaction, position, charge, force_long)
    end

    return force_long
end

function Ewald3D_long_force_k!(K::Tuple{T, T, T, T}, interaction::Ewald3DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force_long::Vector{Point{3, T}}) where{T}
    k_x, k_y, k_z, k = K
    L_x, L_y, L_z = interaction.L
    α = interaction.α

    ρ_k = zero(ComplexF64)
    for j in 1:interaction.n_atoms
        x_j, y_j, z_j = position[j]
        ρ_k += charge[j] * exp(1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
    end

    for i in 1:interaction.n_atoms
        x_i, y_i, z_i = position[i]
        force_long[i] += - charge[i] * T(imag(ρ_k * exp( - 1.0im * (k_x * x_i + k_y * y_i + k_z * z_i)))) * exp(- k^2 / (4 * α^2)) / k^2 /interaction.ϵ * Point(k_x, k_y, k_z) / (L_x * L_y * L_z)
    end

    return nothing
end

function Ewald3D_long_force_k0!(interaction::Ewald3DInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}, force_long::Vector{Point{3, T}}) where{T}

    L_x, L_y, L_z = interaction.L

    P = Point(zero(T), zero(T), zero(T))
    for i in 1:interaction.n_atoms
        P += charge[i] * position[i]
    end

    para = - one(T) / (2 * interaction.ϵ_inf + one(T)) / (2 * interaction.ϵ * L_x * L_y * L_z)

    for i in 1:interaction.n_atoms
        force_long[i] += para * T(2) * P * charge[i]
    end

    return nothing
end