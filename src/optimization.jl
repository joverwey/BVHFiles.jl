export squared_errors, total_squared_errors, optimize_offsets!



"""
    squared_errors(g::BVHGraph, f::Integer)

Calculate the sum of the squared differences between the current and stored positions 
of the vertices for a given frame `f`.

Only those vertices are taken into account that have got positions for every frame.

See also: [`total_squared_errors`](@ref)
"""
squared_errors(g::BVHGraph, f::Integer) = sum(norm(positions(g, v, f) - global_position(g, v, f))^2 for v in filter(v -> size(positions(g, v), 1) == nframes(g), vertices(g)))

squared_errors(f::Integer) = g -> squared_errors(g, f)


"""
    total_squared_errors(g::BVHGraph)

Calculate the sum of the squared differences between the current and stored positions 
of the vertices for all frames.

Only those vertices are taken into account that have got positions for every frame.

See also: [`squared_errors`](@ref)
"""
total_squared_errors(g::BVHGraph) = sum(squared_errors(g, f) for f in frames(g))

total_squared_errors() = g -> total_squared_errors(g)


"""
    optimize_offsets!(g::BVHGraph)

For each edge calculate the average norm of the global positions between its source 
and destination, then adjust the offset vector accordingly.

After removing vertices this function can reduce the deviations of the vertices from 
their original global positions since the new offsets are usually biased.

# Examples
```julia
julia> g = load("Test.bvh") |>
           global_positions! |>
           remove_joint!(7) |>
           remove_joint!(13) |>
           remove_joints!("J_L_Bale", "J_R_Bale", "J_L4", "J_L3", "J_L1", "J_T12", "J_T10", "J_T9", "J_T8", "J_T6", "J_T5", "J_T4", "J_T3", "J_T2") |>
           optimize_offsets!
BVHGraph
Name: Test.bvh
[...]
```

See also: [`total_squared_errors`](@ref)
"""
function optimize_offsets!(g::BVHGraph)
    for v in vertices(g)
        v == 1 && continue
        v₋₁ = inneighbors(g, v)[1]
        off = offset(g, v₋₁, v)
        avg = 0.0

        for f in frames(g)
            d = positions(g, v, f) - positions(g, v₋₁, f)
            avg += norm(d) / nframes(g)
        end
        
        scale = avg / norm(off)
        offset!(g, v₋₁, v, scale * off)
    end

    return g
end

optimize_offsets!() = g -> optimize_offsets!(g)