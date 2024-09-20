module EwaldSummations

using SpecialFunctions, ForwardDiff, CellListMap
using Base.Threads, Distributed
using OhMyThreads: tmapreduce

using Reexport
@reexport using ExTinyMD

export CoulumbEnergy, CoulumbForce, force
export ICM_reflect, SysQ2D, SysQ2DInit, Sys3D, Energy_Q2D, Force_Q2D, Force_self_Q2D, Energy_3D, Force_3D
export Ewald2DInteraction, Ewald2D_long_energy, Ewald2D_short_energy, Ewald2D_short_force, Ewald2D_long_force, Ewald2D_long_energy_k, Ewald2D_long_energy_k0
export Ewald2D_short_energy_i, Ewald2D_short_energy_N, Ewald2D_long_energy_N
export Ewald3DInteraction, Ewald3D_long_energy, Ewald3D_long_energy_k, Ewald3D_long_energy_k0, Ewald3D_long_force, Ewald3D_long_force_k!, Ewald3D_long_force_k0!

export IcmEwald2DInteraction, CellListICM, ICM_Ewald2D_long_energy, ICM_Ewald_short_energy, ICM_Ewald2D_energy
export IcmEwald3DInteraction, ICM_Ewald3D_long_energy, ICM_Ewald3D_slab_energy, ICM_Ewald3D_energy
export ICM_Ewald2D_force, ICM_Ewald3D_force

include("types.jl")

include("direct_sum/Coulumb.jl")
include("direct_sum/direct_sum_Q2D.jl")
include("direct_sum/direct_sum_3D.jl")

include("Ewald2D/Ewald2D.jl")
include("Ewald2D/Ewald2D_long.jl")
include("Ewald2D/Ewald2D_short.jl")

include("Ewald3D/Ewald3D.jl")
include("Ewald3D/Ewald3D_long.jl")
include("Ewald3D/Ewald3D_short.jl")

include("ICM/ICM.jl")
include("ICM/ICM_short.jl")
include("ICM/ICM_long.jl")

end
