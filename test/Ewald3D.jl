using ExTinyMD, EwaldSummations
using Test

@testset "Compare Ewald3D energy with direct sum" begin
    @info "Test for Ewald3D energy"
    n_atoms = 100
    L = 100.0
    boundary = Boundary((L, L, L), (1, 1, 1))

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    Ewald3D_interaction = Ewald3DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 0.3, ϵ_inf = 1.0)

    neighbor = CellList3D(info, Ewald3D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald3D_interaction, neighbor, info, atoms)

    sys_3d = Sys3D((L, L, L), (30, 30, 30), ϵ = 0.3)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    energy_direct = Energy_3D(sys_3d, coords, charge)

    @test isapprox(energy_ewald, energy_direct, atol = 1e-2)
end

@testset "Compare Ewald3D force with direct sum" begin
    @info "Test for Ewald3D force"
    n_atoms = 100
    L = 100.0
    boundary = Boundary((L, L, L), (1, 1, 1))

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    Ewald3D_interaction = Ewald3DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 0.3, ϵ_inf = 1.0)

    neighbor = CellList3D(info, Ewald3D_interaction.r_c, boundary, 1)
    force_ewald = force(Ewald3D_interaction, neighbor, info, atoms)

    sys_3d = Sys3D((L, L, L), (30, 30, 30), ϵ = 0.3)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    force_direct = Force_3D(sys_3d, coords, charge)

    for i in 1:n_atoms
        @test isapprox(force_direct[i][1], force_ewald[i][1], atol = 1e-4)
        @test isapprox(force_direct[i][2], force_ewald[i][2], atol = 1e-4)
        @test isapprox(force_direct[i][3], force_ewald[i][3], atol = 1e-4)
    end
end