function Ewald3DInteraction(n_atoms::Int, s::T, α::T, L::NTuple{3,T}; ϵ::T = one(T), ϵ_inf::T = Inf) where{T}
    r_c = s / α
    k_c = 2 * α * s
    k_set = k_set_3D(k_c, L)
    return Ewald3DInteraction(n_atoms, ϵ, ϵ_inf, α, r_c, k_c, L, k_set)
end

function k_set_3D(k_c::T, L::NTuple{3,T}) where{T}
    L_x, L_y, L_z = L

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1
    mz_max = ceil(Int, k_c * L_z / 2π) + 1

    k_set = Vector{NTuple{4, T}}()

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            for m_z in - mz_max : mz_max
                k_x = m_x * 2π / L_x
                k_y = m_y * 2π / L_y
                k_z = m_z * 2π / L_z
                k = sqrt(k_x^2 + k_y^2 + k_z^2)
                if 0 < k <= k_c
                    push!(k_set, (k_x, k_y, k_z, k))
                end
            end
        end
    end

    return k_set
end

function ExTinyMD.energy(interaction::Ewald3DInteraction{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    charge = [atoms[info.particle_info[i].id].charge for i in 1:interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:interaction.n_atoms]

    energy_long = Ewald3D_long_energy(interaction, position, charge)
    energy_short = Ewald3D_short_energy(interaction, neighbor, position, charge)

    @debug "Ewald3D El = $energy_long, Es = $energy_short"

    return energy_long + energy_short
end

function force(interaction::Ewald3DInteraction{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    charge = [atoms[info.particle_info[i].id].charge for i in 1:interaction.n_atoms]
    position = [info.particle_info[i].position for i in 1:interaction.n_atoms]

    force_long = Ewald3D_long_force(interaction, position, charge)
    force_short = Ewald3D_short_force(interaction, neighbor, position, charge)

    return force_long + force_short
end