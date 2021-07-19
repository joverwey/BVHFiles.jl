using BVHFiles
using Test


@testset "graph" begin
    println("Starting test of 'graph'")

    g = BVHGraph(3)
    @test add_edge!(g, 1, 2) == true
    @inferred add_edge!(g, 2, 3, offset = [1.0, 3.0, 5.0])
    @inferred add_vertex!(g, name = "JOINT L3", sequence = :ZXY)
    
    @inferred has_vertex(g, 2)
    @test has_vertex(g, 2) == true
    @test has_vertex(g, 5) == false

    @inferred has_edge(g, 2, 3)
    @test has_edge(g, 2, 3) == true
    @test has_edge(g, 3, 2) == false

    @inferred props(g, 3)
    
    @inferred AbstractString name(g)
    @inferred AbstractString name(g, 3)
    @inferred AbstractString name!(g, "Test.bvh")
    @inferred AbstractString name!(g, 1, "JOINT hip")
    @inferred AbstractString name!(g, 2, "JOINT L5")

    @inferred sequence(g)
    @inferred sequence(g, 3)
    sequence!(g, :ZXY)
    sequence!(g, 1, :ZXY)

    @inferred nframes(g)
    @test nframes(g) == 1
    @inferred frames(g)

    @inferred frametime(g)
    @test frametime(g) == 0.01

    @inferred Matrix{Float64} rotations(g, 1)
    @inferred Matrix{Float64} rotations!(g, 1, rand(Float64, 10, 3))
    @inferred Matrix{Float64} rotations(g, 1, 5)

    @inferred Matrix{Float64} positions(g)
    @inferred Matrix{Float64} positions(g, 1)
    @inferred Matrix{Float64} positions!(g, rand(Float64, 10, 3))
    @inferred Matrix{Float64} positions!(g, 1, rand(Float64, 10, 3))
    @inferred Matrix{Float64} positions(g, 1, 5)
    
    @inferred Vector{Float64} offset(g)
    @inferred Vector{Float64} offset(g, 1, 2)
    @inferred Vector{Float64} offset!(g, rand(Float64, 3))
    @inferred Vector{Float64} offset!(g, 1, 2, rand(Float64, 3))

    @inferred find(g, "JOINT L5")
    @test find(g, "JOINT hip") == 1
    @test find(g, "JOINT L5") == 2
    
    @inferred eltype(g)
    @inferred edgetype(g)
    @test is_directed(g) == true
    @inferred print(g)

    @inferred nv(g)
    @test nv(g) == 4
    @inferred ne(g)
    @test ne(g) == 2

    @inferred vertices(g)
    @inferred edges(g)

    @inferred inneighbors(g, 3)
    @test inneighbors(g, 3) == [2]
    @inferred outneighbors(g, 1)
    @test outneighbors(g, 1) == [2]

    @inferred rem_edge!(g, 2, 3)
    @inferred rem_vertex!(g, 3)
    @inferred fadj(g)
    @inferred badj(g)

    @inferred BVHGraph(3)
end