using ExTinyMD, EwaldSummations
using Test
using Random
Random.seed!(1234)

@testset "ICM_Ewald3D_ELC vs Ewald2D" begin
    @info "Test for ICM_Ewald3D energy and force"
    n_atoms = 32
    L = 10.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]

    γ = (0.0, 0.0)

    ICMEwald3D_interaction = IcmEwald3DInteraction(n_atoms, 6.0, 2.0, γ, (L, L, L), 0, 20)
    energy_icmewald = ICM_Ewald3D_energy(ICMEwald3D_interaction, coords, charge)
    force_icmewald = ICM_Ewald3D_force(ICMEwald3D_interaction, coords, charge)

    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 6.0, 2.0, (L, L, L), ϵ = 1.0)
    neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, neighbor, info, atoms)
    force_ewald = force(Ewald2D_interaction, neighbor, info, atoms)
    
    @test energy_icmewald ≈ energy_ewald
    for i in 1:n_atoms
        @test force_icmewald[i][1] ≈ force_ewald[i][1]
        @test force_icmewald[i][2] ≈ force_ewald[i][2]
        @test force_icmewald[i][3] ≈ force_ewald[i][3]
    end
end

@testset "ICM_Ewald3D_ELC vs ICM_Ewald2D" begin
    n_atoms = 32
    L = 10.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]

    for γ in [(0.5, 0.5), (0.5, -0.5), (-0.5, -0.5)]

        ICMEwald3D_interaction = IcmEwald3DInteraction(n_atoms, 6.0, 2.0, γ, (L, L, L), 5, 10)
        energy_icmewald3d = ICM_Ewald3D_energy(ICMEwald3D_interaction, coords, charge)
        force_icmewald3d = ICM_Ewald3D_force(ICMEwald3D_interaction, coords, charge)

        ICMEwald2D_interaction = IcmEwald2DInteraction(n_atoms, 6.0, 2.0, γ, (L, L, L), 5)
        energy_icmewald2d = ICM_Ewald2D_energy(ICMEwald2D_interaction, coords, charge)
        force_icmewald2d = ICM_Ewald2D_force(ICMEwald2D_interaction, coords, charge)
        
        @test energy_icmewald3d ≈ energy_icmewald2d
        for i in 1:n_atoms
            @test force_icmewald3d[i][1] ≈ force_icmewald2d[i][1]
            @test force_icmewald3d[i][2] ≈ force_icmewald2d[i][2]
            @test force_icmewald3d[i][3] ≈ force_icmewald2d[i][3]
        end
    end
end