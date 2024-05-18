using ExTinyMD, EwaldSummations
using Test

@testset "ICM_Ewald2D vs Ewald2D" begin
    @info "Test for ICM_Ewald2D energy"
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

    ICMEwald2D_interaction = IcmEwald2DInteraction(n_atoms, 6.0, 2.0, γ, (L, L, L), 1)
    energy_icmewald = ICMEwald2D_energy(ICMEwald2D_interaction, coords, charge)

    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 6.0, 2.0, (L, L, L), ϵ = 1.0)
    neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, neighbor, info, atoms)
    
    @test energy_icmewald ≈ energy_ewald
end

@testset "ICM_Ewald2D vs ICM_direct_sum" begin
    @info "Test for ICM_Ewald2D energy"
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

    N_real = 100
    N_img = 20
    for γ in [(0.5, 0.5), (-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5)]

        ICMEwald2D_interaction = IcmEwald2DInteraction(n_atoms, 3.0, 1.0, γ, (L, L, L), N_img)
        energy_ewald = ICMEwald2D_energy(ICMEwald2D_interaction, coords, charge)

        sys_q2d = SysQ2D(γ, (L, L, L), N_real, N_img)
        ref_pos, ref_charge = SysQ2DInit(sys_q2d, coords, charge)
        energy_direct_sum = Energy_Q2D(sys_q2d, coords, charge, ref_pos, ref_charge)

        @test isapprox(energy_direct_sum, energy_ewald, atol = 1e-2)
    end
end