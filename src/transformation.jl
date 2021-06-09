export add_frames!, add_joint!, change_sequence!, change_sequences!, project!, 
remove_joint!, remove_joints!, rename!, replace_offset!, replace_offsets!, scale!, zero!



function add_frames!(g::BVHGraph, frames::Integer)
    positions!(g, [positions(g); zeros(Float64, frames, 3)])

    for v in vertices(g)
        if outneighbors(g, v) != []
            rotations!(g, v, [rotations(g, v); zeros(Float64, frames, 3)])
        end
    end

    nframes!(g, nframes(g) + frames)
    return g
end

add_frames!(frames::Integer) = g -> add_frames!(g, frames)


function add_joint!(g::BVHGraph, v₋₁::Integer, v₊₁::Integer, name::AbstractString; fraction::Float64 = 0.5)
    v = add_vertex!(g, name = "JOINT $(name)", sequence = sequence(g, v₋₁), rotations = zeros(Float64, nframes(g), 3))
    off = offset(g, v₋₁, v₊₁)
    add_edge!(g, v₋₁, v, offset = fraction * off)
    add_edge!(g, v, v₊₁, offset = (1 - fraction) * off)
    rem_edge!(g, v₋₁, v₊₁)
    return g
end

add_joint!(v₋₁::Integer, v₊₁::Integer, name::AbstractString; fraction::Float64 = 0.5) = g -> add_joint!(g, v₋₁, v₊₁, name, fraction = fraction)
add_joint!(v₋₁::AbstractString, v₊₁::AbstractString, name::AbstractString; fraction::Float64 = 0.5) = g -> add_joint!(g, find(g, v₋₁), find(g, v₊₁), name, fraction = fraction)

function add_joint!(g::BVHGraph, v₋₁::Integer, name::AbstractString, off::Vector{Float64}, nb::Vector = outneighbors(g, v₋₁))
    v = add_vertex!(g, name = "JOINT $(name)", sequence = sequence(g, v₋₁), rotations = zeros(Float64, nframes(g), 3))

    for n in nb
        if n in outneighbors(g, v₋₁)
            add_edge!(g, v, n, offset = offset(g, v₋₁, n) - off)
            rem_edge!(g, v₋₁, n)
        end
    end

    add_edge!(g, v₋₁, v, offset = off)
    return g
end

add_joint!(v₋₁::Integer, name::AbstractString, off::Vector{Float64}, nb::Vector = outneighbors(g, v₋₁)) = g -> add_joint!(g, v₋₁, name, off, nb)
add_joint!(v₋₁::AbstractString, name::AbstractString, off::Vector{Float64}) = g -> add_joint!(g, find(g, v₋₁), name, off)
add_joint!(g::BVHGraph, v₋₁::AbstractString, name::AbstractString, off::Vector{Float64}, nb::Vector = outneighbors(g, find(g, v₋₁))) = add_joint!(g, find(g, v₋₁), name, off, nb)


"""
    change_sequence!(g::BVHGraph, v::Integer, sym::Symbol)

Change the rotation order of vertex `v` to `sym`. The euler angles are adjusted accordingly.

Valid Symbols are `:XYZ`, `:XYX`, `:XZX`, `:XZY`, `:YXZ`, `:YZX`, `:YXY`, `:YZY`, 
`:ZXY`, `:ZYX`, `:ZXZ`, `:ZYZ`.

See also: [`change_sequences!`](@ref)
"""
function change_sequence!(g::BVHGraph, v::Integer, sym::Symbol)
    r = constructor(sym)
    
    for f in frames(g)
        R = rotation(g, v, f) |> r
        rotations(g, v)[f, :] = degrees(R)
    end

    sequence!(g, v, sym)
    return g
end

change_sequence!(v::Integer, sym::Symbol) = g -> change_sequence!(g, v, sym)


"""
    change_sequences!(g::BVHGraph, sym::Symbol)

Change the rotation order of all vertices to `sym`. The euler angles are adjusted accordingly.

Valid Symbols are `:XYZ`, `:XYX`, `:XZX`, `:XZY`, `:YXZ`, `:YZX`, `:YXY`, `:YZY`, 
`:ZXY`, `:ZYX`, `:ZXZ`, `:ZYZ`.

See also: [`change_sequence!`](@ref)
"""
function change_sequences!(g::BVHGraph, sym::Symbol)
    for v in vertices(g)
        outneighbors(g, v) != [] && change_sequence!(g, v, sym)
    end

    return g
end

change_sequences!(sym::Symbol) = g -> change_sequences!(g, sym)


function project!(g::BVHGraph, h::BVHGraph, T::Matrix = Matrix(1.0I, 3, 3))
    for hv in vertices(h)

        if outneighbors(h, hv) != []
            gv = find(g, name(h, hv))

            for f in frames(h)
                R = rotation(h, hv, f)
                rotation!(g, gv, f, T * R * inv(T))
            end
        end
    end

    for f in frames(h)
        position!(g, f, T * positions(h)[f, :])
    end

    frametime!(g, frametime(h))
    return g
end

project!(x...) = g -> project!(g, x...)


function remove_joint!(g::BVHGraph, v::Integer, v₊₁::Integer = outneighbors(g, v) != [] ? outneighbors(g, v)[1] : 0)
    v₋₁ = inneighbors(g, v)[1]

    if outneighbors(g, v) != []
        oᵥ = offset(g, v₋₁, v)
        oᵥ₊₁ = offset(g, v, v₊₁)
        o = oᵥ + oᵥ₊₁

        for f in frames(g)
            Rᵥ = rotation(g, v, f)
            Rᵥ₋₁ = rotation(g, v₋₁, f)
            B = rotation_between(o, inv(Rᵥ) * oᵥ + oᵥ₊₁)
            rotation!(g, v₋₁, f, Rᵥ₋₁ .* Rᵥ .* B)

            for n in outneighbors(g, v)

                if outneighbors(g, n) != []
                    rotation!(g, n, f, inv(B) * rotation(g, n, f))
                end
            end

            for n in outneighbors(g, v₋₁)

                if outneighbors(g, n) != []
                    rotation!(g, n, f, inv(B) .* inv(Rᵥ) .* rotation(g, n, f))
                end
            end
        end

        for n in outneighbors(g, v)
            add_edge!(g, v₋₁, n)
            offset!(g, v₋₁, n, oᵥ + offset(g, v, n))
        end

    elseif outneighbors(g, v₋₁) |> length == 1
        name!(g, v₋₁, "End Site")
        rotations!(g, v₋₁, zeros(Float64, 1, 3))
    end

    rem_vertex!(g, v)
    return g
end

remove_joint!(v::Integer, x...) = g -> remove_joint!(g, v, x...)
remove_joint!(nameᵥ::AbstractString, nameᵥ₊₁::AbstractString) = g -> remove_joint!(g, find(g, nameᵥ), find(g, nameᵥ₊₁))


function remove_joints!(g::BVHGraph, names::AbstractString...)
    for nameᵥ in names
        remove_joint!(g, find(g, nameᵥ))
    end

    return g
end

remove_joints!(names::AbstractString...) = g -> remove_joints!(g, names...)


function rename!(g::BVHGraph, dict::Dict{String,String})
    for k in keys(dict)
        name!(g, find(g, k), "JOINT $(dict[k])")
    end

    return g
end

rename!(dict::Dict{String,String}) = g -> rename!(g, dict)


function replace_offset!(g::BVHGraph, h::BVHGraph, gv₊₁::Integer, T::Matrix{Float64} = Matrix(1.0I, 3, 3); change_rotation::Bool = true)
    gv = inneighbors(g, gv₊₁)[1]
    hv = find(h, name(g, gv))
    hv₊₁ = find_outneighbor(h, hv, name(g, gv₊₁))
    goffᵥ₊₁ = offset(g, gv, gv₊₁)
    hoffᵥ₊₁ = T * offset(h, hv, hv₊₁)
    scale = norm(goffᵥ₊₁) / norm(hoffᵥ₊₁)
    offset!(g, gv, gv₊₁, scale * hoffᵥ₊₁)

    if change_rotation
        B = rotation_between(hoffᵥ₊₁, goffᵥ₊₁)

        for f in frames(g)
            Rᵥ = rotation(g, gv, symᵥ, f)
            rotation!(g, gv, f, B * Rᵥ)
        
            for n in outneighbors(g, gv)
        
                if outneighbors(g, n) != []
                    Rᵥ₊₁ = rotation(g, n, symᵥ₊₁, f)
                    rotation!(g, n, f, inv(Rᵥ) * inv(B) * Rᵥ * Rᵥ₊₁)
                end
            end
        end
    end

    return g
end


function replace_offsets!(g::BVHGraph, h::BVHGraph, exclude::Vector, T::Matrix{Float64} = Matrix(1.0I, 3, 3))
    push!(exclude, 1)

    for n in outneighbors(g, 1)
        push!(exclude, n)
    end

    for v in vertices(g)
        v in exclude || replace_offset!(g, h, v, T)
    end

    return g
end

replace_offsets!(x...) = g -> replace_offsets!(g, x...)


function scale!(g::BVHGraph, scale::Float64)
    for e in edges(g)
        offset!(g, src(e), dst(e), offset(g, src(e), dst(e)) * scale)
    end

    offset!(g, offset(g) * scale)
    positions!(g, positions(g) * scale)
    return g
end

scale!(scale::Float64) = g -> scale!(g, scale)


function zero!(g::BVHGraph)
    frames = nframes(g)
    positions!(g, zeros(Float64, frames, 3))

    for v in vertices(g)
        neighbors(g, v) != [] && rotations!(g, v, zeros(Float64, frames, 3))
    end

    return g
end

zero!() = g -> zero!(g)