



Base.split(dlm) = s -> split(s, dlm, keepempty = false)

shorten(a::Vector{SubString{String}}) = Symbol(a[1][1], a[2][1], a[3][1])

extend_rotation(sym::Symbol) = join([c * "rotation" for c ∈ string(sym)[1:3]], ' ')
extend_position(sym::Symbol) = join([c * "position" for c ∈ string(sym)[1:3]], ' ')

constructor(sym::Symbol) = getfield(BVHFiles, Symbol("Rot", sym))

degrees(R::Rotation) = [getfield(R, θ) |> rad2deg for θ ∈ [Symbol("theta", n) for n ∈ 1:3]]

radians(v::Vector{Float64}, sym::Symbol, constr::Type = constructor(sym)) = constr([deg2rad(θ) for θ ∈ v]...)


rotation(g::BVHGraph, v::Integer, sym::Symbol, f::Integer) = radians(rotations(g, v, f), sym)
rotation(g::BVHGraph, v::Integer, f::Integer) = radians(rotations(g, v, f), sequence(g, v))

function rotation!(g::BVHGraph, v::Integer, f::Integer, rot)
    sym = sequence(g, v)
    constr = constructor(sym)
    rotations(g, v)[f, :] = rot |> constr |> degrees
    return nothing
end

Rx(ψ) = [1.0 0.0 0.0; 0.0 cos(ψ) -sin(ψ); 0.0 sin(ψ) cos(ψ)]
Ry(θ) = [cos(θ) 0.0 sin(θ); 0.0 1.0 0.0; -sin(θ) 0.0 cos(θ)]
Rz(φ) = [cos(φ) -sin(φ) 0.0; sin(φ) cos(φ) 0.0; 0.0 0.0 1.0]

for name in ("Rxyz", "Rxyx", "Rxzy", "Rxzx", "Ryxz", "Ryzx", "Ryxy", "Ryzy", "Rzxy", "Rzyx", "Rzxz", "Rzyz")
    s = name * "(ψ, θ, φ) = " * "R" * name[2] * "(ψ) * " * "R" * name[3] * "(θ) * " * "R" * name[4] * "(φ)"
    eval(Meta.parse(s))
end

function rot(g::BVHGraph, v::Integer, ψ, θ, φ)
    sym = sequence(g, v)
    f = getfield(BVHFiles, Symbol("R", string(sym) |> lowercase))
    return f(ψ, θ, φ)
end

rot(g::BVHGraph, v::Integer, vec) = rot(g, v, vec...)

function rot(sym::Symbol, ψ, θ, φ)
    f = getfield(BVHFiles, Symbol("R", string(sym) |> lowercase))
    return f(ψ, θ, φ)
end

rot(sym::Symbol, vec) = rot(sym, vec...)