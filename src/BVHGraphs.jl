import Base: eltype, ==, issubset, zero, show

import LightGraphs: AbstractGraph, vertices, inneighbors, outneighbors, 
add_vertex!, add_edge!, rem_vertex!, rem_edge!, edges, is_directed, 
eltype, edgetype, nv, ne, has_vertex, has_edge, neighbors, src, dst

import LightGraphs.SimpleGraphs: SimpleDiGraph, SimpleEdge, Edge, fadj, badj

export BVHGraph, vertices, inneighbors, outneighbors, add_vertex!, add_edge!, 
rem_vertex!, rem_edge!, props, edges, is_directed, eltype, edgetype, nv, ne, 
has_vertex, has_edge, ==, zero, fadj, badj, issubset, neighbors, 
name, name!, frames, frames!, frametime, frametime!, offset, offset!, 
sequence, sequence!, rotations, rotations!, positions, positions!, 
find, find_outneighbor, show, src, dst



Base.@kwdef mutable struct GProps
    name::AbstractString = "BVHGraph.bvh"
    offset::Vector{Float64} = zeros(Float64, 3)
    sequence::Symbol = :XYZ
    positions::Matrix{Float64} = zeros(Float64, 1, 3)
    frames::Int64 = 1
    frametime::Float64 = 0.01
end


Base.@kwdef mutable struct VProps
    name::AbstractString = "End Site"
    sequence::Symbol = :XYZ
    rotations::Matrix{Float64} = zeros(Float64, 1, 3)
    positions::Matrix{Float64} = zeros(Float64, 1, 3)
end


Base.@kwdef mutable struct EProps
    offset::Vector{Float64} = zeros(Float64, 3)
end


"""
    struct BVHGraph{T <: Integer} <: AbstractGraph{T}

Used to store motion capture data. Usually constructed by using `load` on a BVH file. 
"""
struct BVHGraph{T <: Integer} <: AbstractGraph{T}
    graph::SimpleDiGraph{T}
    gprops::GProps
    vprops::Dict{T,VProps}
    eprops::Dict{SimpleEdge{T},EProps}
end

function BVHGraph(n::Integer; 
    name::AbstractString = "BVHGraph.bvh", 
    offset::Vector{Float64} = zeros(Float64, 3),
    sequence::Symbol = :XYZ,
    positions::Matrix{Float64} = zeros(Float64, 1, 3),
    frames::Int64 = 1,
    frametime::Float64 = 0.01)

    T = eltype(n)
    g = SimpleDiGraph(n)
    gprops = GProps(name = name, 
                    offset = offset, 
                    sequence = sequence, 
                    positions = positions, 
                    frames = frames, 
                    frametime = frametime)
    vprops = Dict{T,VProps}()
    eprops = Dict{SimpleEdge{T},EProps}()

    for v in vertices(g)
        vprops[v] = VProps()
    end

    return BVHGraph(g, gprops, vprops, eprops)
end


SimpleDiGraph(g::BVHGraph) = g.graph


vertices(g::BVHGraph) = vertices(g.graph)

inneighbors(g::BVHGraph, v::Integer) = inneighbors(g.graph, v)
outneighbors(g::BVHGraph, v::Integer) = fadj(g.graph, v)


function add_vertex!(g::BVHGraph; 
    name::AbstractString = "End Site", 
    sequence::Symbol = :XYZ, 
    rotations::Matrix{Float64} = zeros(Float64, 1, 3), 
    positions::Matrix{Float64} = zeros(Float64, 1, 3))

    add_vertex!(g.graph) || return 0
    v = nv(g)
    g.vprops[v] = VProps(name, sequence, rotations, positions)
    return v
end

function add_edge!(g::BVHGraph, s::Integer, d::Integer; offset::Vector{Float64} = zeros(Float64, 3))
    add_edge!(g.graph, s, d) || return false
    e = Edge(s, d)
    g.eprops[e] = EProps(offset = offset)
    return true
end

function rem_vertex!(g::BVHGraph, v::Integer)
    v in vertices(g) || return false

    lastv = nv(g)
    lastvprops = props(g, lastv)
    lasteoutprops = Dict(n => props(g, lastv, n) for n in outneighbors(g, lastv))
    lasteinprops = Dict(n => props(g, n, lastv) for n in inneighbors(g, lastv))

    for n in outneighbors(g, lastv)
        delete!(g.eprops, Edge(lastv, n))
    end

    for n in inneighbors(g, lastv)
        delete!(g.eprops, Edge(n, lastv))
    end

    if v != lastv

        for n in outneighbors(g, v)
            delete!(g.eprops, Edge(v, n))
        end

        for n in inneighbors(g, v)
            delete!(g.eprops, Edge(n, v))
        end

    end

    delete!(g.vprops, lastv)
    retval = rem_vertex!(g.graph, v)
    retval || return false

    if v != lastv
        g.vprops[v] = lastvprops

        for n in outneighbors(g, v)
            g.eprops[Edge(v, n)] = lasteoutprops[n]
        end

        for n in inneighbors(g, v)
            g.eprops[Edge(n, v)] = lasteinprops[n]
        end
    end

    return true
end

function rem_edge!(g::BVHGraph, s::Integer, d::Integer)
    delete!(g.eprops, Edge(s, d))
    rem_edge!(g.graph, s, d)
end


props(g::BVHGraph) = getproperty(g, :gprops)
props(g::BVHGraph, v::Integer) = get(g.vprops, v, VProps())
props(g::BVHGraph, s::Integer, d::Integer) = get(g.eprops, Edge(s, d), EProps())



edges(g::BVHGraph) = edges(g.graph)

is_directed(::Type{BVHGraph}) = true
is_directed(::Type{BVHGraph{T}}) where T = true
is_directed(g::BVHGraph) = true

eltype(g::BVHGraph) = eltype(g.graph)
edgetype(g::BVHGraph) = edgetype(g.graph)

nv(g::BVHGraph) = nv(g.graph)
ne(g::BVHGraph) = ne(g.graph)

has_vertex(g::BVHGraph, x...) = has_vertex(g.graph, x...)
@inline has_edge(g::BVHGraph, x...) = has_edge(g.graph, x...)

==(g::BVHGraph, h::BVHGraph) = g.graph == h.graph

zero(::BVHGraph) = BVHGraph(0)
zero(g::BVHGraph) = BVHGraph(0, name = g.name)



@inline fadj(g::BVHGraph, x...) = fadj(g.graph, x...)
@inline badj(g::BVHGraph, x...) = badj(g.graph, x...)


issubset(g::T, h::T) where T <: BVHGraph = issubset(g.graph, h.graph)


neighbors(g::BVHGraph, v::Integer) = outneighbors(g, v)



name(g::BVHGraph) = getproperty(props(g), :name)
name(g::BVHGraph, v::Integer) = getproperty(props(g, v), :name)
name(v::Integer) = g -> name(g, v)

name!(g::BVHGraph, name::AbstractString) = setproperty!(props(g), :name, name)
name!(name::AbstractString) = g -> name!(g, name)
name!(g::BVHGraph, v::Integer, name::AbstractString) = setproperty!(props(g, v), :name, name)
name!(v::Integer, name::AbstractString) = g -> name!(g, v, name)


frames(g::BVHGraph) = getproperty(props(g), :frames)

frames!(g::BVHGraph, f::Int64) = setproperty!(props(g), :frames, f)
frames!(f::Int64) = frames!(g, f)
# frames(g::BVHGraph) = 1:getproperty(props(g), :frames)


frametime(g::BVHGraph) = getproperty(props(g), :frametime)

frametime!(g::BVHGraph, t::Float64) = setproperty!(props(g), :frametime, t)
frametime!(t::Float64) = g -> frametime!(g, t)


offset(g::BVHGraph) = getproperty(props(g), :offset)
offset(g::BVHGraph, s::Integer, d::Integer) = getproperty(props(g, s, d), :offset)
offset(g::BVHGraph, d::Integer) = offset(g, inneighbors(g, d)[1], d)


offset!(g::BVHGraph, o::Vector{Float64}) = setproperty!(props(g), :offset, o)
offset!(g::BVHGraph, s::Integer, d::Integer, o::Vector{Float64}) = setproperty!(props(g, s, d), :offset, o)
offset!(g::BVHGraph, d::Integer, o::Vector{Float64}) = offset!(g, inneighbors(g, d)[1], d, o)


sequence(g::BVHGraph) = getproperty(props(g), :sequence)
sequence(g::BVHGraph, v::Integer) = getproperty(props(g, v), :sequence)
sequence(v::Integer) = g -> sequence(g, v)


sequence!(g::BVHGraph, sym::Symbol) = setproperty!(props(g), :sequence, sym)
sequence!(sym::Symbol) = g -> sequence!(g, sym)
sequence!(g::BVHGraph, v::Integer, sym::Symbol) = setproperty!(props(g, v), :sequence, sym)
sequence!(v::Integer, sym::Symbol) = g -> sequence!(g, v, sym)


rotations(g::BVHGraph, v::Integer) = getproperty(props(g, v), :rotations)
rotations(g::BVHGraph, v::Integer, f::Integer) = rotations(g, v)[f, :]


rotations!(g::BVHGraph, v::Integer, matrix::Matrix{Float64}) = setproperty!(props(g, v), :rotations, matrix)


positions(g::BVHGraph) = getproperty(props(g), :positions)
positions(g::BVHGraph, v::Integer) = getproperty(props(g, v), :positions)
positions(g::BVHGraph, v::Integer, f::Integer) = positions(g, v)[f, :]


positions!(g::BVHGraph, matrix::Matrix{Float64}) = setproperty!(props(g), :positions, matrix)
positions!(g::BVHGraph, v::Integer, matrix::Matrix{Float64}) = setproperty!(props(g, v), :positions, matrix)



function find(g::BVHGraph, vertexname::AbstractString)
    for v in vertices(g)
        name(g, v) |> endswith(vertexname) && return v
    end

    return 0
end


function find_outneighbor(g::BVHGraph, v, vertexname::AbstractString)
    endswith(vertexname, "End Site") && return outneighbors(g, v)[1]

    for n in outneighbors(g, v)
        name(g, n) |> endswith(vertexname) && return n
    end

    return 0
end



function show(io::IO, g::BVHGraph)
    println(io, "BVHGraph")
    println(io, "Name: $(name(g))")
    println(io, "Frames: $(frames(g))")
    println(io, "Frame Time: $(frametime(g))")
    println(io, "Vertices: $(nv(g))\n")
    println(io, get_names_formatted(g))
    return nothing
end

function get_names_formatted(g::BVHGraph, v::Integer = 1, s::AbstractString = "", spaces::Integer = 0)
    s *= "$(' '^spaces)$(v) $(name(g, v))\n"

    for n in outneighbors(g, v)
        s = get_names_formatted(g, n, s, spaces + 1)
    end

    return s
end