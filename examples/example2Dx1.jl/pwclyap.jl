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

terminal = VC.rectangle(
    env["termSet"]["lb"] .- 0.1, env["termSet"]["ub"] .+ 0.1, lib
)

transients = VC.complement(terminal)
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

initial = VC.rectangle(env["initSet"]["lb"], env["initSet"]["ub"], lib)

unsafes = VC.complement(safe)
frame = VC.rectangle([-2, -2], [2, 2], lib)
for unsafe in unsafes
    for unsafe in unsafes
        intersect!(unsafe, hrep(frame))
        removehredundancy!(unsafe)
    end
end

sys = VC.PWASystem(pieces, initial, unsafes)
prob = VC.pwc_lyapunov_problem(sys)

nstep = 500
for istep = 1:nstep
    print(@sprintf("%3i", istep), " - ")
    # plt = plot(xlabel="x1", ylabel="x2")
    # plot!(xlims=(-2.05, 2.05), ylims=(-2.05, 2.05))

    # plot!(safe, fc=nothing)
    # plot!(initial, fa=0.2, fc=:yellow, lc=nothing)

    # for unsafe in unsafes
    #     plot!(unsafe, fa=0.2, fc=:red, ls=:dash)
    # end

    # for (i, domain) in enumerate(prob.domains)
    #     plot!(domain, fc=nothing, lc=:blue)
    #     x = prob.points[i]
    #     scatter!([x[1]], [x[2]], c=:black)
    #     dynamic = prob.dynamics[i]
    #     post = Polyhedra.translate(dynamic.A*domain, dynamic.b)
    #     plot!(post, fc=nothing, lc=:red)
    # end

    dataprob = VC.data_pwc_lyapunov_problem(prob)

    # for edge in dataprob.graph.edges
    #     i1, i2 = edge.source, edge.target
    #     x1, x2 = prob.points[i1], prob.points[i2]
    #     plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:closed))
    # end

    # for i in dataprob.support_initial
    #     x = prob.points[i]
    #     scatter!([x[1]], [x[2]], shape=:xcross, c=:orange, msw=2)
    # end

    # for i in dataprob.support_unsafe
    #     x = prob.points[i]
    #     scatter!([x[1]], [x[2]], shape=:xcross, c=:red)
    # end

    # fname = string("pwclyap/", @sprintf("%03i", istep), "_graph.png")
    # savefig(string(@__DIR__, "/data/", fname))

    model, cs_var, r_var = VC.pwc_lyapunov_model(dataprob, 1e3, solver)

    optimize!(model)

    @assert termination_status(model) == MOI.OPTIMAL
    @assert primal_status(model) == MOI.FEASIBLE_POINT
    @assert dual_status(model) == MOI.FEASIBLE_POINT
    r = objective_value(model)
    println("r generator: ", r)
    @assert r > 0

    cs = map(c -> value(c), cs_var)

    ymin = minimum(cs)
    ymax = maximum(cs)

    # plt = plot(xlabel="x1", ylabel="x2")
    # plot!(xlims=(-2.05, 2.05), ylims=(-2.05, 2.05))

    # for (i, domain) in enumerate(prob.domains)
    #     plot!(domain, fill_z=cs[i], fc=cgrad(scale=(ymin, ymax)), lc=:blue)
    #     if cs[i] > 0
    #         x = prob.points[i]
    #         scatter!([x[1]], [x[2]], shape=:xcross, c=:red)
    #     end
    # end

    print("      verify safe: ")
    r, i1, i2 = VC.verify_pwc_lyapunov_safe(prob.graph_unsafe, cs)
    println((r, i1, i2))

    if r ≤ 0
        domain1, domain2 = prob.domains[i1], prob.unsafes[i2]
        dynamic = prob.dynamics[i1]
        pre = VC._prehs(halfspaces(domain2), dynamic)
        hit = intersect(domain1, hrep(pre))
        x1 = center_of_mass(hit)
        x2 = dynamic.A*x1 + dynamic.b
        # plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:closed))
        inew = VC.split_piece!(prob, i1, x1, 1e-9)
        @assert inew > 0
    else
        print("      verify decrease: ")
        r, i1, i2 = VC.verify_pwc_lyapunov_decrease(prob.overgraph_piece, cs)
        println((r, i1, i2))
        if r ≤ 0
            domain1, domain2 = prob.domains[i1], prob.domains[i2]
            dynamic = prob.dynamics[i1]
            pre = VC._prehs(halfspaces(domain2), dynamic)
            hit = intersect(domain1, hrep(pre))
            try
                x1 = center_of_mass(hit)
                x2 = dynamic.A*x1 + dynamic.b
                # plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:closed))
                inew1 = VC.split_piece!(prob, i1, x1, 1e-9)
                @assert inew1 > 0
                inew2 = VC.split_piece!(prob, i2, x2, 1e-9)
                @assert inew2 > 0
                push!(prob.undergraph_piece.edges, VC.Edge(inew1, inew2))
            catch err
                @warn(err)
                # plot!(domain1, fc=nothing, lc=:green)
                # plot!(domain2, fc=nothing, lc=:white, ls=:dot)
                # post = Polyhedra.translate(dynamic.A*domain1, dynamic.b)
                # plot!(post, fc=nothing, lc=:red)
                r = 1
            end
        end
    end

    # for edge in prob.overgraph_piece.edges
    #     i1, i2 = edge.source, edge.target
    #     x1, x2 = prob.points[i1], prob.points[i2]
    #     plot!([x1[1], x2[1]], [x1[2], x2[2]], arrow=(style=:closed))
    # end

    # fname = string("/pwclyap/", @sprintf("%03i", istep), "_lyap.png")
    # savefig(string(@__DIR__, "/data/", fname))

    if r > 0
        display(r); break
    end
end

end # module