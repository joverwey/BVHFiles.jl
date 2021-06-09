export load, save



"""
    load(filename::AbstractString)

Read a BVH file with name `filename` and return a `BVHGraph`.

# Examples
```julia
julia> g = load("Test.bvh")
BVHGraph
Name: Test.bvh
[...]
```
"""
function load(filename::AbstractString)
    list = read(filename, String) |> split((' ', '\t', '\n', '\r'))
    g = BVHGraph(1, name = filename, 
                    offset = [parse(Float64, e) for e in list[6:8]],
                    sequence = short(list[11:13]))

    name!(g, 1, "ROOT" * ' ' * list[3])
    sequence!(g, 1, short(list[14:16]))
    i = 17

    function add_joint(v₋₁::Integer)
        while true
            v = add_vertex!(g, name = list[i] * ' ' * list[i + 1])
            add_edge!(g, v₋₁, v, offset = [parse(Float64, e) for e in list[i + 4:i + 6]])
        
            if list[i] == "JOINT"
                sequence!(g, v, short(list[i + 9:i + 11]))
                i += 12
                add_joint(v)
                i += 1
            else
                i += 8
            end
        
            list[i] != "JOINT" && list[i] != "End" && break
        end
    end

    add_joint(1)
    frames = parse(Int64, list[i + 3])
    nframes!(g, frames)
    frametime!(g, parse(Float64, list[i + 6]))
    positions!(g, zeros(Float64, frames, 3))

    for v in vertices(g)
        outneighbors(g, v) != [] && rotations!(g, v, zeros(Float64, frames, 3))
    end

    i += 7

    function add_frame(v::Integer, f::Integer)
        rotations(g, v)[f, :] = [parse(Float64, e) for e in list[i:i + 2]]
        i += 3
        
        for n in outneighbors(g, v)
            outneighbors(g, n) != [] && add_frame(n, f)
        end
    end

    for f in 1:frames
        positions(g)[f, :] = [parse(Float64, e) for e in list[i:i + 2]]
        rotations(g, 1)[f, :] = [parse(Float64, e) for e in list[i + 3:i + 5]]
        i += 6

        for n in outneighbors(g, 1)
            outneighbors(g, n) != [] && add_frame(n, f)
        end
    end

    return g
end


"""
    save(g::BVHGraph, filename::AbstractString)

Save a BVHGraph `g` to `filename` as a BVH file.
"""
function save(g::BVHGraph, filename::AbstractString)

    function add_hierarchy(io::IO, v₋₁::Integer, v::Integer, tabs::Integer)
        println(io, '\t'^tabs, name(g, v))
        println(io, '\t'^tabs, "{")
        println(io, '\t'^(tabs + 1), "OFFSET", ['\t' * string(e) for e in offset(g, v₋₁, v)]...)

        if outneighbors(g, v) != []
            println(io, '\t'^(tabs + 1), "CHANNELS 3 ", sequence(g, v) |> long)

            for n in outneighbors(g, v)
                add_hierarchy(io, v, n, tabs + 1)
            end
        end

        println(io, '\t'^tabs, "}")
    end

    function add_numbers(io::IO, f::Integer, v::Integer)
        print(io, ['\t' * string(e) for e in rotations(g, v, f)]...)

        for n in outneighbors(g, v)
            outneighbors(g, n) != [] && add_numbers(io, f, n)
        end
    end

    open(filename, "w") do io
        println(io, "HIERARCHY")
        println(io, name(g, 1))
        println(io, "{")
        println(io, "\tOFFSET", ['\t' * string(e) for e in offset(g)]...)
        println(io, "\tCHANNELS 6 ", sequence(g) |> long_position, " ", sequence(g, 1) |> long)

        for n in outneighbors(g, 1)
            add_hierarchy(io, 1, n, 1)
        end

        println(io, "}")
        println(io, "MOTION")
        println(io, "Frames: ", nframes(g))
        println(io, "Frame Time: ", frametime(g))

        for f in frames(g)
            pos = positions(g)[f, :]
            print(io, pos[1], "\t", pos[2], "\t", pos[3])
            print(io, ['\t' * string(e) for e in rotations(g, 1, f)]...)

            for n in outneighbors(g, 1)
                add_numbers(io, f, n)
            end

            print(io, "\n")
        end
    end
end