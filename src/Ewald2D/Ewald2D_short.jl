function Ewald2D_short_energy(interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    neighbor_list = neighbor.neighbor_list

    boundary = Q2dBoundary(interaction.L[1], interaction.L[2], interaction.L[3])

    energy_short = Atomic{T}(zero(T))

    @threads for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, interaction.r_c)
        if iszero(r_sq)
            E = zero(T)
        else
            q_1 = q[i]
            q_2 = q[j]
            E = Ewald2D_Es_pair(q_1, q_2, interaction.α, r_sq)
        end
        atomic_add!(energy_short, E)
    end

    @threads for i in 1:interaction.n_atoms
        t = Ewald2D_Es_self(q[i], interaction.α)
        atomic_add!(energy_short, t)
    end

    return energy_short[] / (4π * interaction.ϵ)
end

"""
    Ewald2D_short_energy_i(i::Int, interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where {T}

Calculate the short-range part of the Ewald summation energy for a given particle `i` in a 2D periodic system.

# Arguments
- `i::Int`: Index of the particle for which the short-range energy is calculated.
- `interaction::Ewald2DInteraction{T}`: Struct containing parameters for the Ewald summation, including the system dimensions and the Ewald parameter `α`.
- `neighbor::CellList3D{T}`: Neighbor list containing pairs of interacting particles and their distances.
- `position::Vector{Point{3, T}}`: Vector of particle positions in 3D space.
- `q::Vector{T}`: Vector of particle charges.

# Returns
- `energy_short::T`: The short-range part of the Ewald summation energy for particle `i`.

# Details
The function iterates over the neighbor list to compute the pairwise short-range interaction energy between particle `i` and its neighbors. It also includes the self-interaction energy for particle `i`. The total short-range energy is normalized by the dielectric constant `ϵ` and a factor of `4π`.

"""
function Ewald2D_short_energy_i(i::Int, interaction::Ewald2DInteraction{T}, neighbor_dict::Dict{Int, Vector{Int}}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    boundary = Q2dBoundary(interaction.L[1], interaction.L[2], interaction.L[3])

    energy_short = zero(T)

    for j in neighbor_dict[i]
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, interaction.r_c)
        if !iszero(r_sq)
            q_1 = q[i]
            q_2 = q[j]
            E = Ewald2D_Es_pair(q_1, q_2, interaction.α, r_sq)
            energy_short += E / 2
        end
    end

    energy_short += Ewald2D_Es_self(q[i], interaction.α)

    return energy_short / (4π * interaction.ϵ)
end

"""
    Ewald2D_short_energy_N(N::Int, interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where {T}

Calculate the short-range part of the Ewald summation energy of the first N particles for a 2D system.

# Arguments
- `N::Int`: The number of particles.
- `interaction::Ewald2DInteraction{T}`: The Ewald interaction parameters.
- `neighbor::CellList3D{T}`: The cell list for neighbor searching.
- `position::Vector{Point{3, T}}`: The positions of the particles.
- `q::Vector{T}`: The charges of the particles.

# Returns
- `energy_short_N::Vector{T}`: A vector containing the short-range energy for each particle.
"""
function Ewald2D_short_energy_N(N::Int, interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    energy_short_N = zeros(T, N)
    neighbor_dict = Dict{Int, Vector{Int}}()
    for i in 1:N
        neighbor_dict[i] = Int[]
    end
    for (i, j, _) in neighbor.neighbor_list
        (i ≤ N) && push!(neighbor_dict[i], j)
        (j ≤ N) && push!(neighbor_dict[j], i)
    end
    @threads for i in 1:N
        energy_short_N[i] = Ewald2D_short_energy_i(i, interaction, neighbor_dict, position, q)
    end
    return energy_short_N
end

function Ewald2D_Es_pair(q_1::T, q_2::T, α::T, r_sq::T) where{T}
    return q_1 * q_2 * erfc(α * sqrt(r_sq)) / sqrt(r_sq)
end

function Ewald2D_Es_self(q::T, α::T) where{T}
    return - q^2 * α / sqrt(π)
end

function Ewald2D_short_force(interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    neighbor_list = neighbor.neighbor_list
 
    force_short = [[Point(zero(T), zero(T), zero(T)) for _=1:interaction.n_atoms] for _ in 1:nthreads()]
    boundary = Q2dBoundary(interaction.L[1], interaction.L[2], interaction.L[3])

    @threads for (i, j, ρ) in neighbor_list
        n = threadid()
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = q[i]
            q_2 = q[j]
            F_ij = Ewald2D_Fs_pair(q_1, q_2, interaction.α, coord_1, coord_2)
            force_short[n][i] += F_ij / (4π * interaction.ϵ)
            force_short[n][j] -= F_ij / (4π * interaction.ϵ)
        end
    end

    return sum(force_short)
end

function Ewald2D_Fs_pair(q_1::T, q_2::T, α::T, coord_1::Point{3, T}, coord_2::Point{3, T}) where{T}
    energy = r -> q_1 * q_2 * erfc(α * r) / r
    Δr = sqrt(dist2(coord_1, coord_2))
    if iszero(Δr)
        return Point(zero(T), zero(T), zero(T))
    else
        force = - ForwardDiff.derivative(energy, Δr)
        return force * (coord_1 - coord_2) / Δr
    end
end