module VoronoiController

using LinearAlgebra
using Polyhedra
using JuMP

const Polyhedron = DefaultPolyhedron{
    Float64,
    Polyhedra.Intersection{Float64,Vector{Float64},Int64},
    Polyhedra.Hull{Float64,Vector{Float64},Int64}
}

_hot(α, i, N) = [j == i ? α : zero(α) for j = 1:N]

struct Dynamic
    A::Matrix{Float64}
    b::Vector{Float64}
end

struct Piece
    domain::Polyhedron
    dynamic::Dynamic
end

struct PWASystem
    pieces::Vector{Piece}
    initials::Vector{Polyhedron}
    unsafes::Vector{Polyhedron}
end

struct Edge
    source::Int
    target::Int
end

struct Graph
    edges::Vector{Edge}
end

include("polyhedra.jl")
include("graph.jl")
include("pwalyapunov.jl")
include("pwclyapunov.jl")

end # module VoronoiController
