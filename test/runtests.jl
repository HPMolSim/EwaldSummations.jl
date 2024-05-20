using EwaldSummations
using Test
using ExTinyMD

@testset "EwaldSummations.jl" begin
    include("Ewald2D.jl")
    include("Ewald3D.jl")
    include("ICM_Ewald2D.jl")
    include("ICM_Ewald3D.jl")
end
