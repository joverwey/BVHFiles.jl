



Base.split(dlm) = s -> split(s, dlm, keepempty = false)

short(a::Vector{SubString{String}}) = Symbol(a[1][1], a[2][1], a[3][1])

long(sym::Symbol) = join([c * "rotation" for c ∈ string(sym)[1:3]], ' ')
long_position(sym::Symbol) = join([c * "position" for c ∈ string(sym)[1:3]], ' ')

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

function position!(g::BVHGraph, f::Integer, pos)
    positions(g)[f, :] = pos
    return nothing
end

Rx(ψ) = [1.0 0.0 0.0; 0.0 cos(ψ) -sin(ψ); 0.0 sin(ψ) cos(ψ)]
Ry(θ) = [cos(θ) 0.0 sin(θ); 0.0 1.0 0.0; -sin(θ) 0.0 cos(θ)]
Rz(φ) = [cos(φ) -sin(φ) 0.0; sin(φ) cos(φ) 0.0; 0.0 0.0 1.0]

Rxyz(ψ, θ, φ) = Rx(ψ) * Ry(θ) * Rz(φ)
Rxyz(v) = Rxyz(v...)