function rectangle(lb, ub, lib)::Polyhedron
    bb = zip(lb, ub)
    N = length(bb)
    hs = HalfSpace{Float64,Vector{Float64}}[]
    for (i, (xl, xu)) in enumerate(bb)
        push!(hs, HalfSpace(_hot(+1.0, i, N), +xu))
        push!(hs, HalfSpace(_hot(-1.0, i, N), -xl))
    end
    return polyhedron(hrep(hs), lib)
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

function complement(p)
    ps = Polyhedron[]
    hs = halfspaces(p)
    for id in Iterators.product([(0, 1) for i in eachindex(hs)]...)
        all(isone, id) && continue
        H = HalfSpace{Float64,Vector{Float64}}[]
        for (i, h) in enumerate(hs)
            if id[i] == 0
                push!(H, HalfSpace(-h.a, -h.β))
            else
                push!(H, HalfSpace(+h.a, +h.β))
            end
        end
        p = polyhedron(hrep(H), library(p))
        removehredundancy!(p)
        isempty(p) && continue
        push!(ps, p)
    end
    return ps
end