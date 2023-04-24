function rectangle(lb, ub, lib)
    bb = zip(lb, ub)
    N = length(bb)
    hs = HalfSpace{Float64,Vector{Float64}}[]
    for (i, (xl, xu)) in enumerate(bb)
        push!(hs, HalfSpace(_hot(+1.0, i, N), +xu))
        push!(hs, HalfSpace(_hot(-1.0, i, N), -xl))
    end
    return polyhedron(hrep(hs), lib)::Polyhedron
end

function voronoi_partition(points, lib)
    cells = Polyhedron[]
    for (i1, point1) in enumerate(points)
        hs = HalfSpace{Float64,Vector{Float64}}[]
        for (i2, point2) in enumerate(points)
            i1 == i2 && continue
            v = point2 - point1
            nv = norm(v)
            push!(hs, HalfSpace(v/nv, dot(v, point1)/nv + nv/2))
        end
        push!(cells, polyhedron(hrep(hs), lib))
    end
    return cells
end

function complement(p, lib)
    acc = HalfSpace{Float64,Vector{Float64}}[]
    ps = Polyhedron[]
    for h in halfspaces(p)
        H = hrep([HalfSpace(-h.a, -h.β)])
        if !isempty(acc)
            H = intersect(H, hrep(acc))
        end
        push!(ps, polyhedron(H, lib))
        acc = push!(acc, h)
    end
    return ps
end