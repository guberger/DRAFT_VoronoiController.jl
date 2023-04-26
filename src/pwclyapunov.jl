struct DataPWCLyapunovProblem
    nnode::Int
    points::Vector{Vector{Float64}}
    graph::Graph
    support_initial::BitSet
    support_unsafe::BitSet
end

mutable struct PWCLyapunovProblem
    npiece::Int
    domains::Vector{Polyhedron}
    dynamics::Vector{Dynamic}
    points::Vector{Vector{Float64}}
    unsafes::Vector{Polyhedron}
    overgraph_piece::Graph
    undergraph_piece::Graph
    support_initial::BitSet
    graph_unsafe::Graph
    support_default::BitSet
end

function pwc_lyapunov_problem(sys::PWASystem)
    noninits = complement(sys.initial)
    pieces = Piece[]
    support_initial = BitSet()
    for piece in sys.pieces
        domain = intersect(piece.domain, hrep(sys.initial))
        removehredundancy!(domain)
        if isempty(domain)
            push!(pieces, piece)
            continue
        end
        push!(pieces, Piece(domain, piece.dynamic))
        push!(support_initial, length(pieces))
        for noninit in noninits
            domain = intersect(piece.domain, hrep(noninit))
            removehredundancy!(domain)
            isempty(domain) && continue
            push!(pieces, Piece(domain, piece.dynamic))
        end
    end
    domains = map(piece -> piece.domain, pieces)
    dynamics = map(piece -> piece.dynamic, pieces)
    points = map(piece -> center_of_mass(piece.domain), pieces)
    overgraph_piece = transition_graph(pieces, domains)
    graph_unsafe = transition_graph(pieces, sys.unsafes)
    npiece = length(pieces)
    return PWCLyapunovProblem(
        npiece, domains, dynamics, points, sys.unsafes,
        overgraph_piece, Graph(Set{Edge}()),
        support_initial, graph_unsafe, BitSet(1:npiece)
    )
end

function split_piece!(prob::PWCLyapunovProblem, i, xnew, δ)::Int
    if i ∈ prob.support_default
        prob.points[i] = xnew
        delete!(prob.support_default, i)
        return i
    end
    xprev = prob.points[i]
    # (xnew - xprev)'*(x - xprev) ≤ || xnew - xprev ||^2 / 2
    v = xnew - xprev
    nv = norm(v)
    nv < δ && return i
    a, β = v/nv, nv/2 + dot(v, xprev)/nv
    domain1 = intersect(prob.domains[i], hrep([HalfSpace(+a, +β)]))
    removehredundancy!(domain1)
    isempty(domain1) && return 0
    domain2 = intersect(prob.domains[i], hrep([HalfSpace(-a, -β)]))
    removehredundancy!(domain2)
    isempty(domain2) && return 0
    prob.domains[i] = domain1
    push!(prob.domains, domain2)
    push!(prob.dynamics, prob.dynamics[i])
    push!(prob.points, xnew)
    prob.npiece += 1
    inew = prob.npiece
    if i ∈ prob.support_initial
        push!(prob.support_initial, inew)
    end
    edges_to_remove = Set{Edge}()
    edges_to_add = Set{Edge}()
    for edge in prob.overgraph_piece.edges
        i1, i2 = edge.source, edge.target
        if i1 == i
            pre = _prehs(halfspaces(prob.domains[i2]), prob.dynamics[i])
            hit = intersect(prob.domains[i], hrep(pre))
            if isempty(hit)
                push!(edges_to_remove, Edge(i, i2))
            end
            hit = intersect(prob.domains[inew], hrep(pre))
            if !isempty(hit)
                push!(edges_to_add, Edge(inew, i2))
            end
        end
        if i2 == i
            pre = _prehs(halfspaces(prob.domains[i]), prob.dynamics[i1])
            hit = intersect(prob.domains[i1], hrep(pre))
            if isempty(hit)
                push!(edges_to_remove, Edge(i1, i))
            end
            pre = _prehs(halfspaces(prob.domains[inew]), prob.dynamics[i1])
            hit = intersect(prob.domains[i1], hrep(pre))
            if !isempty(hit)
                push!(edges_to_add, Edge(i1, inew))
            end
        end
    end
    union!(prob.overgraph_piece.edges, edges_to_add)
    setdiff!(prob.overgraph_piece.edges, edges_to_remove)
    empty!(edges_to_add)
    empty!(edges_to_remove)
    for edge in prob.graph_unsafe.edges
        i1, i2 = edge.source, edge.target
        if i1 == i
            pre = _prehs(halfspaces(prob.unsafes[i2]), prob.dynamics[i])
            hit = intersect(prob.domains[i], hrep(pre))
            if isempty(hit)
                push!(edges_to_remove, Edge(i, i2))
            end
            hit = intersect(prob.domains[inew], hrep(pre))
            if !isempty(hit)
                push!(edges_to_add, Edge(inew, i2))
            end
        end
    end
    union!(prob.graph_unsafe.edges, edges_to_add)
    setdiff!(prob.graph_unsafe.edges, edges_to_remove)
    return inew
end

function data_pwc_lyapunov_problem(prob::PWCLyapunovProblem)
    support_unsafe = BitSet()
    for edge in prob.graph_unsafe.edges
        edge.source ∈ prob.support_default && continue
        push!(support_unsafe, edge.source)
    end
    return DataPWCLyapunovProblem(
        prob.npiece, prob.points, prob.undergraph_piece,
        prob.support_initial, support_unsafe
    )    
end

function pwc_lyapunov_model(prob::DataPWCLyapunovProblem, rmax, solver)
    model = solver()
    cs = [@variable(model) for i = 1:prob.nnode]
    r = @variable(model, upper_bound=rmax)

    for edge in prob.graph.edges
        i1, i2 = edge.source, edge.target
        @constraint(model, cs[i2] ≤ cs[i1] - r)
    end
    for i in prob.support_initial
        @constraint(model, cs[i] ≤ 0)
    end
    for i in prob.support_unsafe
        @constraint(model, cs[i] ≥ r)
    end
    for i1 = 1:prob.nnode
        for i2 = 1:(i1 - 1)
            dx = norm(prob.points[i1] - prob.points[i2])
            @constraint(model, cs[i1] - cs[i2] ≤ dx)
            @constraint(model, cs[i2] - cs[i1] ≤ dx)
        end
    end

    @objective(model, Max, r)

    return model, cs, r    
end

function verify_pwc_lyapunov_safe(graph_unsafe, cs)
    rmin::Float64 = Inf
    i1opt::Int, i2opt::Int = 0, 0
    for edge in graph_unsafe.edges
        i1, i2 = edge.source, edge.target
        r = cs[i1]
        r ≥ rmin && continue
        rmin, i1opt, i2opt = r, i1, i2
    end
    return rmin, i1opt, i2opt
end

function verify_pwc_lyapunov_decrease(graph_piece, cs)
    rmin::Float64 = Inf
    i1opt::Int, i2opt::Int = 0, 0
    for edge in graph_piece.edges
        i1, i2 = edge.source, edge.target
        r = cs[i1] - cs[i2]
        r ≥ rmin && continue
        rmin, i1opt, i2opt = r, i1, i2
    end
    return rmin, i1opt, i2opt
end