using EwaldSummations
using Test
using ExTinyMD

@testset "EwaldSummations.jl" begin
    include("Ewald2D.jl")
    include("Ewald3D.jl")
end
