function _pre_image(p, dynamic)
    At, b = transpose(dynamic.A), dynamic.b
    hs = [
        HalfSpace(At*hs.a, hs.β - dot(hs.a, b)) for hs in halfspaces(p)
    ]
    return hrep(hs)
end

function transition_graph(pieces, unsafes, initials)
    piece_edges = Edge[]
    unsafe_edges = Edge[]
    initial_edges = Edge[]
    for (i1, piece1) in enumerate(pieces)
        for (i2, piece2) in enumerate(pieces)
            pre = _pre_image(piece2.domain, piece1.dynamic)
            hit = intersect(piece1.domain, pre)
            isempty(hit) && continue
            push!(piece_edges, Edge(i1, i2))
        end
        for (i2, unsafe) in enumerate(unsafes)
            pre = _pre_image(unsafe, piece1.dynamic)
            hit = intersect(piece1.domain, pre)
            isempty(hit) && continue
            push!(unsafe_edges, Edge(i1, i2))
        end
        for (i2, initial) in enumerate(initials)
            hit = intersect(piece1.domain, initial)
            isempty(hit) && continue
            push!(initial_edges, Edge(i1, i2))
        end
    end
    return Graph(
        pieces, unsafes, initials, piece_edges, unsafe_edges, initial_edges
    )
end