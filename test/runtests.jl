using BVHFiles
using Test


@testset "Graph" begin
    include("graph.jl")
end

@testset "File" begin
    include("file.jl")
end

@testset "Interpolation" begin
    include("interpolation.jl")
end