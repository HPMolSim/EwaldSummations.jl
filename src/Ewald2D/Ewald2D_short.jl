function Ewald2D_short_energy(interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    neighbor_list = neighbor.neighbor_list

    boundary = Q2dBoundary(interaction.L[1], interaction.L[2], interaction.L[3])

    energy_short = @distributed (+) for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, interaction.r_c)
        if iszero(r_sq)
            E = zero(T)
        else
            q_1 = q[i]
            q_2 = q[j]
            E = Ewald2D_Es_pair(q_1, q_2, interaction.α, r_sq)
        end
        E
    end

    energy_short += @distributed (+) for i in 1:interaction.n_atoms
        Ewald2D_Es_self(q[i], interaction.α)
    end

    return energy_short / (4π * interaction.ϵ)
end

function Ewald2D_Es_pair(q_1::T, q_2::T, α::T, r_sq::T) where{T}
    return q_1 * q_2 * erfc(α * sqrt(r_sq)) / sqrt(r_sq)
end

function Ewald2D_Es_self(q::T, α::T) where{T}
    return - q^2 * α / sqrt(π)
end

function Ewald2D_short_force(interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    neighbor_list = neighbor.neighbor_list
 
    force_short = [Point(zero(T), zero(T), zero(T)) for _=1:interaction.n_atoms]
    boundary = Q2dBoundary(interaction.L[1], interaction.L[2], interaction.L[3])

    for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = q[i]
            q_2 = q[j]
            F_ij = Ewald2D_Fs_pair(q_1, q_2, interaction.α, coord_1, coord_2)
            force_short[i] += F_ij / (4π * interaction.ϵ)
            force_short[j] -= F_ij / (4π * interaction.ϵ)
        end
    end

    return force_short
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