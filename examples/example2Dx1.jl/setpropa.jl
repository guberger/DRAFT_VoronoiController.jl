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

safe = VC.rectangle(
    float(env["workspace"]["lb"]), float(env["workspace"]["ub"]), lib
)

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
    float(env["termSet"]["lb"]) .- 0.01, float(env["termSet"]["ub"]) .+ 0.01, lib
)
transients = VC.complement(terminal, lib)

sets = [VC.rectangle(
    float(env["initSet"]["lb"]), float(env["initSet"]["ub"]), lib
)]

unsafes = VC.complement(safe, lib)

nstep = 100
for istep = 1:nstep
    plt = plot(xlabel="x1", ylabel="x2")
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
    savefig(string(@__DIR__, "/set_", @sprintf("%03i", istep), ".png"))
    global sets = new_sets
end

# old_pieces = pieces
# pieces = VC.Piece[]
# for old_piece in old_pieces
#     domain = intersect(old_piece.domain, hrep(terminal))
#     if isempty(domain)
#         push!(pieces, old_piece)
#         continue
#     end
#     for transient in transients
#         domain = intersect(old_piece.domain, hrep(transient))
#         removehredundancy!(domain)
#         isempty(domain) && continue
#         push!(pieces, VC.Piece(domain, old_piece.dynamic))
#     end
# end


end # module