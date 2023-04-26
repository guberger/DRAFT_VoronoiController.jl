function _prehs(hs, dynamic)
    At, b = transpose(dynamic.A), dynamic.b
    return [HalfSpace(At*h.a, h.β - dot(h.a, b)) for h in hs]
end

function transition_graph(pieces, sets)
    edges = Set{Edge}()
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
    edges = Set{Edge}()
    for (i2, set2) in enumerate(sets2)
        supp = intersection_support(sets1, set2)
        for i1 in supp
            push!(edges, Edge(i1, i2))
        end
    end
    return Graph(edges)
end

function intersection_support(sets1, set2)
    supp = BitSet()
    H = hrep(set2)
    for (i1, set1) in enumerate(sets1)
        hit = intersect(set1, H)
        isempty(hit) && continue
        push!(supp, i1)
    end
    return supp
end