using BVHFiles
using Test


@testset "interpolation" begin
    println("Starting test of 'interpolation'")

    g = load("Example.bvh") |>
        interpolate!(1)

    @test nframes(g) == 481


    g = load("Example.bvh") |>
    interpolate!(2)

    @test nframes(g) == 721

    
    g = load("Example.bvh") |>
    interpolate!(3)

    @test nframes(g) == 961
end