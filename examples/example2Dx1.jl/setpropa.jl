module Example2Dx1

using LinearAlgebra
using Printf
using Polyhedra

using JuMP
using Gurobi
const GUROBI_ENV = Gurobi.Env()
const optimizer = optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
)
const lib = DefaultLibrary{Float64}(optimizer)

using YAML
using CSV
using DataFrames

include("../../src/VoronoiController.jl")
VC = VoronoiController

env = YAML.load_file(string(@__DIR__, "/env.yaml"))
controller = CSV.read(string(@__DIR__, "/counterExamples.csv"), DataFrame)

N = env["numStateDim"]

using Plots

safe = VC.rectangle(env["workspace"]["lb"], env["workspace"]["ub"], lib)

samples = [
    (collect(x), reshape(collect(a), 2, 2), reshape(collect(b), 2))
    for (x, a, b) in zip(
        zip(controller.x1, controller.x2),
        zip(controller.A11, controller.A21, controller.A12, controller.A22),
        zip(controller.b1, controller.b2)
    )
]
# for sample in samples
#     sample[2] .= [0.95 -0.2; 0.2 0.75]
#     sample[3] .= zeros(N)
# end
points = getindex.(samples, 1)

cells = VC.voronoi_partition(points, lib)
for cell in cells
    intersect!(cell, hrep(safe))
    removehredundancy!(cell)
    @assert !isempty(cell)
end

pieces = VC.Piece[]
for (point, A, b) in samples
    iok = 0
    for (i, cell) in enumerate(cells)
        point ∉ cell && continue
        iok = i; break
    end
    @assert iok > 0
    push!(pieces, VC.Piece(cells[iok], VC.Dynamic(A, b)))
end

terminal = VC.rectangle(
    env["termSet"]["lb"] .- 0.01, env["termSet"]["ub"] .+ 0.01, lib
)

# transients
acc = HalfSpace{Float64,Vector{Float64}}[]
transients = Polyhedron[]
for h in halfspaces(terminal)
    H = hrep([HalfSpace(-h.a, -h.β)])
    if !isempty(acc)
        H = intersect(H, hrep(acc))
    end
    push!(acc, h)
    transient = polyhedron(H, lib)
    removehredundancy!(transient)
    isempty(transient) && continue
    push!(transients, transient)
end

sets = [VC.rectangle(env["initSet"]["lb"], env["initSet"]["ub"], lib)]

nstep = 100
for istep = 1:nstep
    isempty(sets) && break
    plot(xlabel="x1", ylabel="x2")
    for cell in cells
        plot!(cell, fc=nothing, lc=:red)
    end
    plot!(safe, fc=nothing)
    plot!(terminal, fc=nothing, ls=:dot)
    old_sets = sets
    new_sets = VC.Polyhedron[]
    for old_set in old_sets
        hit = intersect(old_set, hrep(terminal))
        if isempty(hit)
            push!(new_sets, old_set)
            continue
        end
        for transient in transients
            set = intersect(old_set, hrep(transient))
            removehredundancy!(set)
            isempty(set) && continue
            push!(new_sets, set)
        end
    end
    for set in new_sets
        plot!(set, fa=0.1, fc=:green, lc=nothing)
    end
    old_sets = new_sets
    new_sets = VC.Polyhedron[]
    for old_set in old_sets
        for piece in pieces
            dynamic = piece.dynamic
            domain = intersect(old_set, hrep(piece.domain))
            isempty(domain) && continue
            set = Polyhedra.translate(dynamic.A*domain, dynamic.b)
            removehredundancy!(set)
            push!(new_sets, set)
        end
    end
    for set in new_sets
        plot!(set, fc=nothing, lc=:blue)
    end
    fname = string("setpropa/", @sprintf("%03i", istep), ".png")
    savefig(string(@__DIR__, "/data/", fname))
    global sets = new_sets
end

end # module