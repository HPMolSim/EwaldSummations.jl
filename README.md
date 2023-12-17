# EwaldSummations

[![Build Status](https://github.com/HPMolSim/EwaldSummations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HPMolSim/EwaldSummations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/HPMolSim/EwaldSummations.jl.svg?branch=main)](https://travis-ci.com/HPMolSim/EwaldSummations.jl)
[![Coverage](https://codecov.io/gh/HPMolSim/EwaldSummations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HPMolSim/EwaldSummations.jl)


`EwaldSummations.jl` is an implementation of the basic Ewald summation algorithm, including Ewald3D and Ewald2D. This package also include the naive approach: sum directly over periodic directions.

This package is designed to provide standard/correct results of electrostatic interaction in charged particle systems, and can be used as a tool box when developing new algorithms.

## Getting Started

Add this package and its deps `ExTinyMD.jl` in julia by typing `]` in Julia REPL and then
```julia
pkg> add ExTinyMD, EwaldSummations
```
to install the package.

## Examples of Usage

Here are some examples of usage, for more details please see the test files.

### Triply periodic systems

Here is an example for calculating energy of 3D systems, via Ewald3D or direction summation.

```julia
# Init the system
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

# randomly generate position of particles
info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

# take alpha = 0.2 and s = 3.0
Ewald3D_interaction = Ewald3DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 1.0, ϵ_inf = 1.0)

# setup cell list for short range interaction
neighbor = CellList3D(info, Ewald3D_interaction.r_c, boundary, 1)

# calculation interaction energy
energy_ewald = energy(Ewald3D_interaction, neighbor, info, atoms)

# direction summation, take 60 * 60 * 60 periodic images into account
sys_3d = Sys3D((L, L, L), (30, 30, 30), ϵ = 1.0)
coords = [p_info.position for p_info in info.particle_info]
charge = [atoms[p_info.id].charge for p_info in info.particle_info]
energy_direct = Energy_3D(sys_3d, coords, charge)
```

### Doubly periodic systems

Here is an example for calculating energy of Q2D systems, via Ewald3D or direction summation.

```julia
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

Ewald2D_interaction = Ewald2DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 1.0)

neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
energy_ewald = energy(Ewald2D_interaction, neighbor, info, atoms)

# consider 400 * 400 periodic images, with no boundary mismatch
sys_q2d = SysQ2D((0.0, 0.0), (L, L, L), 200, 0, ϵ = 1.0)
coords = [p_info.position for p_info in info.particle_info]
charge = [atoms[p_info.id].charge for p_info in info.particle_info]
ref_pos, ref_charge = SysQ2DInit(sys_q2d, coords, charge)
energy_icm = Energy_Q2D(sys_q2d, coords, charge, ref_pos, ref_charge)
```

### Doubly periodic systems with dielectric mismatch

Here is an example for calculating energy of Q2D systems confined by dielectric mismatch, via direction summation together with image charge method.

```julia
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

# consider 100 * 100 periodic images, and 20 layers of image charges
sys_q2d = SysQ2D((0.9, 0.9), (L, L, L), 50, 10, ϵ = 1.0)
coords = [p_info.position for p_info in info.particle_info]
charge = [atoms[p_info.id].charge for p_info in info.particle_info]
ref_pos, ref_charge = SysQ2DInit(sys_q2d, coords, charge)
energy_icm = Energy_Q2D(sys_q2d, coords, charge, ref_pos, ref_charge)
```

## How to Contribute

If you find any bug or have any suggestion, please open an issue.