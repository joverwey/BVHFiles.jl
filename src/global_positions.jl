export global_positions!, global_position



"""
    global_positions!(g::BVHGraph)

Calculate the global positions for every vertex and frame in `g` and store them.

See also: [`global_position`](@ref)
"""
function global_positions!(g::BVHGraph)
    for v in vertices(g)
        positions!(g, v, zeros(Float64, nframes(g), 3))
    end

    function calculate_positions(v::Integer, f::Integer, P::Matrix{Float64})
        off = offset(g, v)

        if outneighbors(g, v) != []
            R = rotation(g, v, f)
            M = P * [   R[1, 1] R[1, 2] R[1, 3] off[1]; 
                        R[2, 1] R[2, 2] R[2, 3] off[2]; 
                        R[3, 1] R[3, 2] R[3, 3] off[3]; 
                        0.0 0.0 0.0 1.0]
    
            positions(g, v)[f, :] = M[1:3, 4]
    
            for n in outneighbors(g, v)
                calculate_positions(n, f, M)
            end
        else
            vec = P * [off[1:3]..., 1.0]
            positions(g, v)[f, :] = vec[1:3]
        end
    end

    for f in frames(g)
        off = offset(g)
        pos = positions(g)[f, :]
        R = rotation(g, 1, f)
        M = [   R[1, 1] R[1, 2] R[1, 3] pos[1] + off[1]; 
                R[2, 1] R[2, 2] R[2, 3] pos[2] + off[2]; 
                R[3, 1] R[3, 2] R[3, 3] pos[3] + off[3]; 
                0.0 0.0 0.0 1.0]

        positions(g, 1)[f, :] = M[1:3, 4]

        for n in outneighbors(g, 1)
            calculate_positions(n, f, M)
        end
    end

    return g
end

global_positions!() = g -> global_positions!(g)


"""
    global_position(g::BVHGraph, v::Integer, f::Integer, N::Matrix{Float64} = Matrix(1.0I, 4, 4))

Return the global position of vertex `v` for frame `f`.

This function does not store the result in `g`. `N` should not be changed.

See also: [`global_positions!`](@ref)
"""
function global_position(g::BVHGraph, v::Integer, f::Integer, N::Matrix{Float64} = Matrix(1.0I, 4, 4))
    if outneighbors(g, v) != []
        R = rotation(g, v, f)
    else
        R = one(RotMatrix{3, Float64})
    end
    
    if v != 1
        v₋₁ = inneighbors(g, v)[1]
        off = offset(g, v₋₁, v)
        A = [   R[1, 1] R[1, 2] R[1, 3] off[1]; 
                R[2, 1] R[2, 2] R[2, 3] off[2]; 
                R[3, 1] R[3, 2] R[3, 3] off[3]; 
                0.0 0.0 0.0 1.0] * N
        return global_position(g, v₋₁, f, A)
    else
        off = offset(g)
        pos = positions(g)[f, :]
        A = [   R[1, 1] R[1, 2] R[1, 3] pos[1] + off[1]; 
                R[2, 1] R[2, 2] R[2, 3] pos[2] + off[2]; 
                R[3, 1] R[3, 2] R[3, 3] pos[3] + off[3]; 
                0.0 0.0 0.0 1.0] * N
        return A[1:3, 4]
    end
end