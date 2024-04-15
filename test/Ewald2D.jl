@testset "compare Ewald2D energy with ICM" begin
    @info "Test for Ewald2D energy"
    n_atoms = 100
    L = 100.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 0.3)

    neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, neighbor, info, atoms)

    N_real = 200
    N_img = 0
    sys_q2d = SysQ2D((0.0, 0.0), (L, L, L), N_real, N_img, ϵ = 0.3)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    ref_pos, ref_charge = SysQ2DInit(sys_q2d, coords, charge)
    energy_icm = Energy_Q2D(sys_q2d, coords, charge, ref_pos, ref_charge)

    @test isapprox(energy_icm, energy_ewald, atol = 1e-2)
end

@testset "compare Ewald2D force with ICM" begin
    @info "Test for Ewald2D force"
    n_atoms = 100
    L = 100.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 0.3)

    neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    force_ewald = force(Ewald2D_interaction, neighbor, info, atoms)

    N_real = 300
    N_img = 0
    sys_q2d = SysQ2D((0.0, 0.0), (L, L, L), N_real, N_img, ϵ = 0.3)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    ref_pos, ref_charge = SysQ2DInit(sys_q2d, coords, charge)
    force_icm = Force_Q2D(sys_q2d, coords, charge, ref_pos, ref_charge)

    for i in 1:n_atoms
        @test isapprox(force_icm[i][1], force_ewald[i][1], atol = 1e-3)
        @test isapprox(force_icm[i][2], force_ewald[i][2], atol = 1e-3)
        @test isapprox(force_icm[i][3], force_ewald[i][3], atol = 1e-3)
    end
end