_avar(model, N) = @variable(model, [1:N], lower_bound=-1, upper_bound=1)

function pwa_lyapunov_model(graph, rmax, N, solver)
    M = length(graph.pieces)
    model = solver()
    as = [_avar(model, N) for i = 1:M]
    βs = [@variable(model) for i = 1:M]
    r = @variable(model, upper_bound=rmax)

    for edge in graph.piece_edges
        i1, i2 = edge.source, edge.target
        domain1, domain2 = graph.pieces[i1].domain, graph.pieces[i2].domain
        dynamic = graph.pieces[i1].dynamic
        At, b = transpose(dynamic.A), dynamic.b
        a1, a2 = as[i1], as[i2]
        β1, β2 = βs[i1], βs[i2]
        da = At*a2 - a1
        dβ = β2 - β1 + r
        for h in halfspaces(domain1)
            λ = @variable(model, lower_bound=0)
            da -= λ*h.a
            dβ += λ*h.β
        end
        for h in halfspaces(domain2)
            λ = @variable(model, lower_bound=0)
            da -= λ*(At*h.a)
            dβ += λ*(h.β - dot(h.a, b))
        end
        @constraint(model, da .== 0)
        @constraint(model, dβ == 0)
    end
    for edge in graph.unsafe_edges
        i1, i2 = edge.source, edge.target
        domain1, domain2 = graph.pieces[i1].domain, graph.unsafes[i2]
        dynamic = graph.pieces[i1].dynamic
        At, b = transpose(dynamic.A), dynamic.b
        a1 = as[i1]
        β1 = βs[i1]
        da = -a1
        dβ = r - β1
        for h in halfspaces(domain1)
            λ = @variable(model, lower_bound=0)
            da -= λ*h.a
            dβ += λ*h.β
        end
        for h in halfspaces(domain2)
            λ = @variable(model, lower_bound=0)
            da -= λ*(At*h.a)
            dβ += λ*(h.β - dot(h.a, b))
        end
        @constraint(model, da .== 0)
        @constraint(model, dβ == 0)
    end
    for edge in graph.initial_edges
        i1, i2 = edge.source, edge.target
        domain1, domain2 = graph.pieces[i1].domain, graph.initials[i2]
        a1 = as[i1]
        β1 = βs[i1]
        da = a1
        dβ = β1
        for h in halfspaces(domain1)
            λ = @variable(model, lower_bound=0)
            da -= λ*h.a
            dβ += λ*h.β
        end
        for h in halfspaces(domain2)
            λ = @variable(model, lower_bound=0)
            da -= λ*h.a
            dβ += λ*h.β
        end
        @constraint(model, da .== 0)
        @constraint(model, dβ == 0)
    end

    @objective(model, Max, r)

    return model, as, βs, r
end