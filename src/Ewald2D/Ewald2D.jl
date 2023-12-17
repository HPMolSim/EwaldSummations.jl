function Ewald2DInteraction(n_atoms::Int, s::T, α::T, L::NTuple{3,T}; ϵ::T = one(T)) where{T}
    r_c = s / α
    k_c = 2 * α * s
    k_set = k_set_2D(k_c, L)
    return Ewald2DInteraction(n_atoms, ϵ, α, r_c, k_c, L, k_set)
end

function k_set_2D(k_c::T, L::NTuple{3,T}) where{T}
    L_x, L_y, L_z = L

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1

    k_set = Vector{NTuple{3, T}}()

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            k_x = m_x * 2π / L_x
            k_y = m_y * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if 0 < k <= k_c
                push!(k_set, (k_x, k_y, k))
            end
        end
    end

    return k_set
end

function energy(interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    charge = [atoms[info.particle_info[i].id].charge for i in 1:interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:interaction.n_atoms]

    energy_long = Ewald2D_long_energy(interaction, position, charge)
    energy_short = Ewald2D_short_energy(interaction, neighbor, position, charge)

    return energy_long + energy_short
end

function force(interaction::Ewald2DInteraction{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    charge = [atoms[info.particle_info[i].id].charge for i in 1:interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:interaction.n_atoms]

    force_long = Ewald2D_long_force(interaction, position, charge)
    force_short = Ewald2D_short_force(interaction, neighbor, position, charge)

    return force_long + force_short
end