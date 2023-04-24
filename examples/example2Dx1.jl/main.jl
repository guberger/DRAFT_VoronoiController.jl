module Example2Dx1

using LinearAlgebra
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

plt = plot(xlabel="x1", ylabel="x2")

safe = VC.rectangle(
    float(env["workspace"]["lb"]), float(env["workspace"]["ub"]), lib
)
plot!(safe, fc=nothing)

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

terminal = VC.rectangle(
    float(env["termSet"]["lb"]), float(env["termSet"]["ub"]), lib
)
plot!(terminal, fa=0.1)

transients = VC.complement(terminal, lib)
old_pieces = pieces
pieces = VC.Piece[]
for old_piece in old_pieces
    domain = intersect(old_piece.domain, hrep(terminal))
    if isempty(domain)
        push!(pieces, old_piece)
        continue
    end
    for transient in transients
        domain = intersect(old_piece.domain, hrep(transient))
        removehredundancy!(domain)
        isempty(domain) && continue
        push!(pieces, VC.Piece(domain, old_piece.dynamic))
    end
end

np1 = 4
np2 = 4
x1min = minimum(x -> getindex(x, 1), Polyhedra.points(vrep(safe)))
x1max = maximum(x -> getindex(x, 1), Polyhedra.points(vrep(safe)))
x1s = range(x1min, x1max, length=np1)
x2min = minimum(x -> getindex(x, 2), Polyhedra.points(vrep(safe)))
x2max = maximum(x -> getindex(x, 2), Polyhedra.points(vrep(safe)))
x2s = range(x2min, x2max, length=np2)
points = collect(map(collect, Iterators.product(x1s, x2s)))[:]
cells = VC.voronoi_partition(points, lib)

old_pieces = pieces
pieces = VC.Piece[]
for old_piece in old_pieces
    for cell in cells
        domain = intersect(old_piece.domain, hrep(cell))
        removehredundancy!(domain)
        isempty(domain) && continue
        push!(pieces, VC.Piece(domain, old_piece.dynamic))
    end
end

colors = palette(:default).colors
for piece in pieces
    plot!(piece.domain, fc=nothing, lc=:red)
    x = center_of_mass(piece.domain)
    c = colors[mod(hash(piece.dynamic) - 1, length(colors)) + 1]
    scatter!([x[1]], [x[2]], mc=c)
end

initials = [VC.rectangle(
    float(env["initSet"]["lb"]), float(env["initSet"]["ub"]), lib
)]
for initial in initials
    plot!(initial, fa=0.1, lc=nothing)
end

unsafes = VC.complement(safe, lib)
frame = VC.rectangle(float([-2, -2]), float([2, 2]), lib)
for unsafe in unsafes
    for unsafe in unsafes
        intersect!(unsafe, hrep(frame))
        removehredundancy!(unsafe)
    end
end
for unsafe in unsafes
    plot!(unsafe, fa=0.2, fc=:red, ls=:dash)
end

plot!(xlims=(-2.05, 2.05), ylims=(-2.05, 2.05))

graph = VC.transition_graph(pieces, unsafes, initials)

iplot = rand(1:length(graph.pieces))*0
for edge in graph.piece_edges
    i1, i2 = edge.source, edge.target
    i1 != iplot && continue
    cell1, cell2 = graph.pieces[i1].domain, graph.pieces[i2].domain
    dynamic1 = graph.pieces[i1].dynamic
    postcell1 = Polyhedra.translate(dynamic1.A*cell1, dynamic1.b)
    plot!(postcell1, fc=nothing)
    x1 = center_of_mass(cell1)
    continue
    if i1 == i2
        scatter!([x1[1]], [x1[2]], shape=:xcross, ms=5, msw=3)
    else
        x2 = center_of_mass(cell2)
        plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:closed))
    end
end

iplot = rand(1:length(graph.pieces))*0
for edge in graph.unsafe_edges
    i1, i2 = edge.source, edge.target
    i1 != iplot && continue
    cell1, cell2 = graph.pieces[i1].domain, graph.unsafes[i2]
    dynamic1 = graph.pieces[i1].dynamic
    postcell1 = Polyhedra.translate(dynamic1.A*cell1, dynamic1.b)
    plot!(postcell1, fc=nothing)
    x1 = center_of_mass(cell1)
    x2 = center_of_mass(cell2)
    plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:open))
end

iplot = rand(1:length(graph.pieces))*0
for edge in graph.initial_edges
    i1, i2 = edge.source, edge.target
    i1 != iplot && continue
    cell1, cell2 = graph.pieces[i1].domain, graph.initials[i2]
    x1 = center_of_mass(cell1)
    x2 = center_of_mass(cell2)
    plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:open))
end

model, as_var, βs_var, r_var = VC.pwa_lyapunov_model(graph, 1e3, N, solver)

optimize!(model)

display(solution_summary(model))
display(objective_value(model))

as = map(a -> value.(a), as_var)
βs = map(β -> value(β), βs_var)

ymin = +Inf
ymax = -Inf
for (i, piece) in enumerate(graph.pieces)
    points = Polyhedra.points(vrep(piece.domain))
    ylocmax = maximum(x -> dot(as[i], x) + βs[i], points)
    ylocmin = minimum(x -> dot(as[i], x) + βs[i], points)
    global ymin, ymax = min(ymin, ylocmin), max(ymax, ylocmax)
end
display((ymin, ymax))

nlev = 50
for (i, piece) in enumerate(graph.pieces)
    # continue
    points = Polyhedra.points(vrep(piece.domain))
    ylocmax = maximum(x -> dot(as[i], x) + βs[i], points)
    ylocmin = minimum(x -> dot(as[i], x) + βs[i], points)
    x1min, x2min = minimum(x -> x[1], points), minimum(x -> x[2], points)
    x1max, x2max = maximum(x -> x[1], points), maximum(x -> x[2], points)
    rad = max(x1max - x1min, x2max - x2min)
    # nlevloc = ylocmin < ylocmax - (ymax - ymin)/1e3 ? ceil(Int, rad/0.1) : 1
    # levs = range(ylocmin, ylocmax, length=nlevloc + 1)
    nlevloc = nlev
    levs = range(ymin, ymax, length=nlevloc + 1)
    for ilev in 1:nlevloc
        H = hrep([
            HalfSpace(+as[i], +levs[ilev + 1] - βs[i]),
            HalfSpace(-as[i], -levs[ilev + 0] + βs[i])
        ])
        p = intersect(piece.domain, H)
        isempty(p) && continue
        midlev = 0.5*(levs[ilev] + levs[ilev + 1])
        plot!(p, fill_z=midlev, fc=cgrad(scale=(ymin, ymax)), line_z=midlev)
    end
end

for (i, piece) in enumerate(graph.pieces)
    # positive: hatched
    H = hrep([HalfSpace(-as[i], 0 + βs[i])])
    p = intersect(piece.domain, H)
    isempty(p) && continue
    plot!(p, lw=0, fc=:green, fillstyle=:/)
end

for (i, piece) in enumerate(graph.pieces)
    # negative: image
    H = hrep([HalfSpace(+as[i], 0 - βs[i])])
    p = intersect(piece.domain, H)
    isempty(p) && continue
    dynamic = graph.pieces[i].dynamic
    postp = Polyhedra.translate(dynamic.A*p, dynamic.b)
    plot!(postp, fc=nothing)
end

for edge in graph.unsafe_edges
    display((as[edge.source], βs[edge.source]))
end

savefig(string(@__DIR__, "/img.png"))
display(plt)

end # module