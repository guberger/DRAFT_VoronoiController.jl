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

struct Edge
    source::Int
    target::Int
end

struct Dynamic
    A::Matrix{Float64}
    b::Vector{Float64}
end

struct Piece
    domain::Polyhedron
    dynamic::Dynamic
end

struct Graph
    pieces::Vector{Piece}
    unsafes::Vector{Polyhedron}
    initials::Vector{Polyhedron}
    piece_edges::Vector{Edge}
    unsafe_edges::Vector{Edge}
    initial_edges::Vector{Edge}
end

include("partition.jl")
include("graph.jl")
include("lyapunov.jl")

end # module VoronoiController
