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
solver() = Model(optimizer)
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
for sample in samples
    sample[2] .= [0.95 -0.2; 0.2 0.75]
    sample[3] .= zeros(N)
end
points = getindex.(samples, 1)

cells = VC.voronoi_partition(points, lib)
for cell in cells
    intersect!(cell, hrep(safe))
    removehredundancy!(cell)
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

terminal = VC.rectangle(env["termSet"]["lb"], env["termSet"]["ub"], lib)

transients = VC.complement(terminal, lib)
old_pieces = pieces
pieces = VC.Piece[]
for old_piece in old_pieces
    domain = intersect(old_piece.domain, hrep(terminal))
    for transient in transients
        domain = intersect(old_piece.domain, hrep(transient))
        removehredundancy!(domain)
        isempty(domain) && continue
        push!(pieces, VC.Piece(domain, old_piece.dynamic))
    end
end

initials = [VC.rectangle(env["initSet"]["lb"], env["initSet"]["ub"], lib)]

unsafes = VC.complement(safe, lib)
frame = VC.rectangle([-2, -2], [2, 2], lib)
for unsafe in unsafes
    for unsafe in unsafes
        intersect!(unsafe, hrep(frame))
        removehredundancy!(unsafe)
    end
end

nstep = 10
for istep = 1:nstep
    plot(xlabel="x1", ylabel="x2")
    plot!(xlims=(-2.05, 2.05), ylims=(-2.05, 2.05))
    for initial in initials
        plot!(initial, fa=0.1, lc=nothing)
    end
    for piece in pieces
        plot!(piece.domain, fc=nothing, lc=:red)
    end
    plot!(safe, fc=nothing)
    plot!(terminal, fc=nothing, ls=:dot)
    for unsafe in unsafes
        plot!(unsafe, fa=0.2, fc=:red, ls=:dash)
    end
    fname = string("pwclyap_", @sprintf("%03i", istep), ".png")
    savefig(string(@__DIR__, "/data/", fname))
end

end # module