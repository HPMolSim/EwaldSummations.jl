function CellListICM(reflect_position::Vector{Point{3, T}}, L::NTuple{3, T}, r_c::T, N_image::Int, n_atoms::Int) where T
    pos = [[p...] for p in reflect_position]
    nlist = neighborlist(pos, r_c; unitcell=[L[1], L[2], (2 * N_image + 1) * L[3] + 2 * r_c])

    rr_list = Vector{Tuple{Int64, Int64, T}}()
    ri_list = Vector{Tuple{Int64, Int64, T}}()
    for (i, j, r) in nlist
        if i <= n_atoms && j <= n_atoms
            push!(rr_list, (i, j, r))
        elseif i <= n_atoms && j > n_atoms
            push!(ri_list, (i, j, r))
        elseif i > n_atoms && j <= n_atoms
            push!(ri_list, (i, j, r))
        end
    end
    return rr_list, ri_list
end

function ICMEwald2D_short_energy(interaction::IcmEwald2DInteraction{T}, ref_pos::Vector{Point{3, T}}, ref_q::Vector{T}) where{T}

    rr_list, ri_list = CellListICM(ref_pos, interaction.L, interaction.r_c, interaction.N_image, interaction.n_atoms)

    energy_short = Atomic{T}(zero(T))

    @threads for (i, j, ρ) in rr_list
        q_1 = ref_q[i]
        q_2 = ref_q[j]
        E = Ewald2D_Es_pair(q_1, q_2, interaction.α, ρ^2)
        atomic_add!(energy_short, E)
    end

    @threads for (i, j, ρ) in ri_list
        q_1 = ref_q[i]
        q_2 = ref_q[j]
        E = Ewald2D_Es_pair(q_1, q_2, interaction.α, ρ^2) / T(2)
        atomic_add!(energy_short, E)
    end

    @threads for i in 1:interaction.n_atoms
        t = Ewald2D_Es_self(ref_q[i], interaction.α)
        atomic_add!(energy_short, t)
    end

    return energy_short[] / (4π * interaction.ϵ)
end