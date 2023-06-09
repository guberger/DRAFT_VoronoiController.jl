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

safe = VC.rectangle(env["workspace"]["lb"], env["workspace"]["ub"], lib)
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

terminal = VC.rectangle(env["termSet"]["lb"], env["termSet"]["ub"], lib)
plot!(terminal, fa=0.1)

transients = VC.complement(terminal)
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

np1 = 3
np2 = 3
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

initial = VC.rectangle(env["initSet"]["lb"], env["initSet"]["ub"], lib)
plot!(initial, fa=0.1, lc=nothing)

unsafes = VC.complement(safe)
display(length(unsafes))
frame = VC.rectangle([-2, -2], [2, 2], lib)
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

sys = VC.PWASystem(pieces, initial, unsafes)
prob = VC.pwa_lyapunov_problem(sys)

iplot = rand(1:length(sys.pieces))*0
for edge in prob.graph_piece.edges
    i1, i2 = edge.source, edge.target
    i1 != iplot && continue
    domain1, domain2 = prob.domains[i1], prob.domains[i2]
    dynamic = prob.dynamics[i1]
    post = Polyhedra.translate(dynamic.A*domain1, dynamic.b)
    plot!(post, fc=nothing)
    x1 = center_of_mass(domain1)
    # continue
    if i1 == i2
        scatter!([x1[1]], [x1[2]], shape=:xcross, ms=5, msw=3)
    else
        x2 = center_of_mass(domain2)
        plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:closed))
    end
end

iplot = rand(1:length(sys.pieces))*0
for edge in prob.graph_unsafe.edges
    i1, i2 = edge.source, edge.target
    i1 != iplot && continue
    domain1, domain2 = prob.domains[i1], prob.unsafes[i2]
    dynamic = prob.dynamics[i1]
    post = Polyhedra.translate(dynamic.A*domain1, dynamic.b)
    plot!(post, fc=nothing)
    x1 = center_of_mass(domain1)
    x2 = center_of_mass(domain2)
    plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:open))
end

iplot = rand(1:length(sys.pieces))*0
for i in prob.support_initial
    i != iplot && continue
    domain1, domain2 = prob.domains[i], initial
    x1 = center_of_mass(domain1)
    x2 = center_of_mass(domain2)
    plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:open))
end

model, as_var, βs_var, r_var = VC.pwa_lyapunov_model(prob, 1e3, N, solver)

optimize!(model)

display(solution_summary(model))
display(objective_value(model))

as = map(a -> value.(a), as_var)
βs = map(β -> value(β), βs_var)

ymin = +Inf
ymax = -Inf
for (i, domain) in enumerate(prob.domains)
    points = Polyhedra.points(vrep(domain))
    ylocmax = maximum(x -> dot(as[i], x) + βs[i], points)
    ylocmin = minimum(x -> dot(as[i], x) + βs[i], points)
    global ymin, ymax = min(ymin, ylocmin), max(ymax, ylocmax)
end
display((ymin, ymax))

nlev = 50
for (i, domain) in enumerate(prob.domains)
    # break
    points = Polyhedra.points(vrep(domain))
    ylocmax = maximum(x -> dot(as[i], x) + βs[i], points)
    ylocmin = minimum(x -> dot(as[i], x) + βs[i], points)
    nlevloc = nlev
    levs = range(ymin, ymax, length=nlevloc + 1)
    for ilev in 1:nlevloc
        H = hrep([
            HalfSpace(+as[i], +levs[ilev + 1] - βs[i]),
            HalfSpace(-as[i], -levs[ilev + 0] + βs[i])
        ])
        p = intersect(domain, H)
        isempty(p) && continue
        midlev = 0.5*(levs[ilev] + levs[ilev + 1])
        plot!(p, fill_z=midlev, fc=cgrad(scale=(ymin, ymax)), line_z=midlev)
    end
end

for (i, domain) in enumerate(prob.domains)
    # positive: hatched
    # break
    H = hrep([HalfSpace(-as[i], 0 + βs[i])])
    p = intersect(domain, H)
    isempty(p) && continue
    plot!(p, lw=0, fc=:green, fillstyle=:/)
end

for (i, domain) in enumerate(prob.domains)
    # negative: image
    # break
    H = hrep([HalfSpace(+as[i], 0 - βs[i])])
    p = intersect(domain, H)
    isempty(p) && continue
    dynamic = prob.dynamics[i]
    post = Polyhedra.translate(dynamic.A*p, dynamic.b)
    plot!(post, fc=nothing)
end

for edge in prob.graph_unsafe.edges
    display((as[edge.source], βs[edge.source]))
end

savefig(string(@__DIR__, "/data/pwalyap/001.png"))
display(plt)

end # module