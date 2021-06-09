export squared_errors, total_squared_errors, optimize_offsets!



squared_errors(g::BVHGraph, f::Integer) = sum(norm(positions(g, v, f) - global_position(g, v, f))^2 for v in filter(v -> size(positions(g, v), 1) == nframes(g), vertices(g)))

squared_errors(f::Integer) = g -> squared_errors(g, f)


total_squared_errors(g::BVHGraph) = sum(squared_errors(g, f) for f in frames(g))

total_squared_errors() = g -> total_squared_errors(g)


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