using BVHFiles
using Test


@testset "interpolation" begin
    println("Starting test of 'interpolation'")

    g = load("Example.bvh") |>
        interpolate!(1)

    @test frames(g) == 481


    g = load("Example.bvh") |>
    interpolate!(2)

    @test frames(g) == 721

    
    g = load("Example.bvh") |>
    interpolate!(3)

    @test frames(g) == 961
end