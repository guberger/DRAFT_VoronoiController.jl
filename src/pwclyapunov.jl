struct DataPWCLyapunovProblem
    nnode::Int
    points::Vector{Vector{Float64}}
    graph::Graph
    set_initial::BitSet
    set_unsafe::BitSet
end