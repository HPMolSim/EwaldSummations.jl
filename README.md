# EwaldSummations

[![Build Status](https://github.com/HPMolSim/EwaldSummations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HPMolSim/EwaldSummations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HPMolSim/EwaldSummations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HPMolSim/EwaldSummations.jl)


`EwaldSummations.jl` is an implementation of the basic Ewald summation algorithm, including Ewald3D and Ewald2D. This package also include the naive approach: sum directly over periodic directions.
Specially, for dielectric confined quasi-2D systems, this package provide a set of tools based on the image charge method, combining with the Ewald2D or Ewald3D-ELC method.

This package is designed to provide standard/correct results of electrostatic interaction in charged particle systems, and can be used as a tool box when developing new algorithms.
For some of the implemented methods, the package support parallel computation by `Base.Threads`, start that by `julia -t 8 your_script.jl`.

## Getting Started

Add this package and its deps `ExTinyMD.jl` in julia by typing `]` in Julia REPL and then
```julia
pkg> add ExTinyMD, EwaldSummations
```
to install the package.

## Numerical Methods

In this package, we provide various methods to calculate the electrostatic interaction energy of charged particle systems. For system without dielectric mismatch, the following methods are available:
* Direction summation: sum directly over periodic directions (3D and quasi-2D)
* Ewald3D: the standard Ewald summation for 3D systems
* Ewald2D: the standard Ewald summation for 2D systems

Then for quasi-2D systems confined by dielectric mismatch, we combine the image charge method and the methods above, given by
* Image charge reflection: generating the image charge series according to the dielectric mismatch
* ICM + direct summation: combining the image charge method and the direct summation
* ICM + Ewald2D: combining the image charge method and the Ewald2D method
* ICM + Ewald3D + ELC: extending the quasi-2D system as triply periodic, and calculate the energy via Ewald3D and electrostatic layer correction (ELC)

For detailed numerical formula, please refer to the reference[^yuan].

## Examples of Usage

Here are some examples of usage, for more details please see the test files.

### Triply periodic systems

Here is an example for calculating energy of 3D systems, via Ewald3D or direction summation.

```julia
# Init the system
using ExTinyMD, EwaldSummations

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

### Doubly periodic systems without dielectric mismatch

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

s = 3.0
alpha = 0.2
γ = (0.9, 0.9)
N_img = 10
N_pad = 20
N_real = 50

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

# direction summations, consider 100 * 100 periodic images, and 20 layers of image charges
sys_q2d = SysQ2D(γ, (L, L, L), N_real, N_img, ϵ = 1.0)
coords = [p_info.position for p_info in info.particle_info]
charge = [atoms[p_info.id].charge for p_info in info.particle_info]
ref_pos, ref_charge = SysQ2DInit(sys_q2d, coords, charge)
energy_icm = Energy_Q2D(sys_q2d, coords, charge, ref_pos, ref_charge)

# ICM-Ewald2D, consider 20 layers of image charges
ICMEwald2D_interaction = IcmEwald2DInteraction(n_atoms, s, alpha, γ, (L, L, L), N_img)
energy_icmewald2d = ICM_Ewald2D_energy(ICMEwald2D_interaction, coords, charge)

# ICM-Ewald3D, consider 20 layers of image charges and 40 layers of padding
ICMEwald3D_interaction = IcmEwald3DInteraction(n_atoms, s, alpha, γ, (L, L, L), N_img, N_pad)
energy_icmewald3d = ICM_Ewald3D_energy(ICMEwald3D_interaction, coords, charge)
```

## Citation

If you use this package in your research, please cite at least one of the following paper:

```


```

## How to Contribute

If you find any bug or have any suggestion, please open an issue.


[^yuan]: J. Yuan, H. S. Antila, and E. Luijten, Particle–particle particle–mesh algorithm for electrolytes between charged dielectric interfaces, J. Chem. Phys., 154 (2021), p. 094115.
