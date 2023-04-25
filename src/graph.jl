function _prehs(hs, dynamic)
    At, b = transpose(dynamic.A), dynamic.b
    return [HalfSpace(At*h.a, h.β - dot(h.a, b)) for h in hs]
end

function transition_graph(pieces, sets)
    edges = Edge[]
    for (i1, piece) in enumerate(pieces)
        for (i2, set) in enumerate(sets)
            pre = _prehs(halfspaces(set), piece.dynamic)
            hit = intersect(piece.domain, hrep(pre))
            isempty(hit) && continue
            push!(edges, Edge(i1, i2))
        end
    end
    return Graph(edges)
end

function intersection_graph(sets1, sets2)
    edges = Edge[]
    for (i1, set1) in enumerate(sets1)
        for (i2, set2) in enumerate(sets2)
            hit = intersect(set1, hrep(set2))
            isempty(hit) && continue
            push!(edges, Edge(i1, i2))
        end
    end
    return Graph(edges)
end