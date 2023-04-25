struct PWALyapunovProblem
    npiece::Int
    domains::Vector{Polyhedron}
    dynamics::Vector{Dynamic}
    initial::Polyhedron
    unsafes::Vector{Polyhedron}
    graph_piece::Graph
    support_initial::BitSet
    graph_unsafe::Graph
end

function pwa_lyapunov_problem(sys::PWASystem)
    domains = map(piece -> piece.domain, sys.pieces)
    dynamics = map(piece -> piece.dynamic, sys.pieces)
    graph_piece = transition_graph(sys.pieces, domains)
    support_initial = intersection_support(domains, sys.initial)
    graph_unsafe = transition_graph(sys.pieces, sys.unsafes)
    return PWALyapunovProblem(
        length(sys.pieces), domains, dynamics, sys.initial, sys.unsafes,
        graph_piece, support_initial, graph_unsafe
    )
end

_avar(model, N) = @variable(model, [1:N], lower_bound=-1, upper_bound=1)

function add_to_expressions!(expr, λ, v, N)
    for i = 1:N
        add_to_expression!(expr[i], λ, v[i])
    end
end

function add_lagrange_negative!(model, da, dβ, a, β, N)
    λ = @variable(model, lower_bound=0)
    add_to_expressions!(da, λ, -a, N)
    add_to_expression!(dβ, λ, β)
end

function add_lagrange_set!(model, da, dβ, hs, N)
    for h in hs
        add_lagrange_negative!(model, da, dβ, h.a, h.β, N)
    end
end

function add_lagrange_set_image!(model, da, dβ, hs, A, b, N)
    for h in hs
        add_lagrange_negative!(model, da, dβ, A'*h.a, h.β - dot(h.a, b), N)
    end
end

function add_pwa_lyapunov_piece_edge(
        model, domains, dynamics, as, βs, r, edge, N
    )
    i1, i2 = edge.source, edge.target
    domain1, domain2 = domains[i1], domains[i2]
    A, b = dynamics[i1].A, dynamics[i1].b
    a1, a2 = as[i1], as[i2]
    β1, β2 = βs[i1], βs[i2]
    da::Vector{AffExpr} = A'*a2 - a1
    dβ::AffExpr = β2 - β1 + r
    add_lagrange_set!(model, da, dβ, halfspaces(domain1), N)
    add_lagrange_set_image!(model, da, dβ, halfspaces(domain2), A, b, N)
    @constraint(model, da .== 0)
    @constraint(model, dβ == 0)
end

function add_pwa_lyapunov_unsafe_edge(
        model, domains, dynamics, unsafes, as, βs, r, edge, N
    )
    i1, i2 = edge.source, edge.target
    domain1, domain2 = domains[i1], unsafes[i2]
    A, b = dynamics[i1].A, dynamics[i1].b
    a1 = as[i1]
    β1 = βs[i1]
    da::Vector{AffExpr} = -a1
    dβ::AffExpr = r - β1
    add_lagrange_set!(model, da, dβ, halfspaces(domain1), N)
    add_lagrange_set_image!(model, da, dβ, halfspaces(domain2), A, b, N)
    @constraint(model, da .== 0)
    @constraint(model, dβ == 0)
end

function add_pwa_lyapunov_initial_support(
        model, domains, initial, as, βs, i, N
    )
    domain1, domain2 = domains[i], initial
    a1 = as[i]
    β1 = βs[i]
    da::Vector{AffExpr} = a1
    dβ::AffExpr = β1
    add_lagrange_set!(model, da, dβ, halfspaces(domain1), N)
    add_lagrange_set!(model, da, dβ, halfspaces(domain2), N)
    @constraint(model, da .== 0)
    @constraint(model, dβ == 0)
end

function pwa_lyapunov_model(prob::PWALyapunovProblem, rmax, N, solver)
    domains = prob.domains
    dynamics = prob.dynamics
    initial = prob.initial
    unsafes = prob.unsafes
    model = solver()
    as = [_avar(model, N) for i = 1:prob.npiece]
    βs = [@variable(model) for i = 1:prob.npiece]
    r = @variable(model, upper_bound=rmax)

    for edge in prob.graph_piece.edges
        add_pwa_lyapunov_piece_edge(
            model, domains, dynamics, as, βs, r, edge, N)
    end
    for edge in prob.graph_unsafe.edges
        add_pwa_lyapunov_unsafe_edge(
            model, domains, dynamics, unsafes, as, βs, r, edge, N)
    end
    for i in prob.support_initial
        add_pwa_lyapunov_initial_support(model, domains, initial, as, βs, i, N)
    end

    @objective(model, Max, r)

    return model, as, βs, r
end