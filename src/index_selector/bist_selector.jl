
""" BiST strategy """
struct BiST <: IndexSelector
    name::String
    L::Int
    K::Int
    M::Int
    patience::Int
    data::SelectorData
    function BiST(L::Int, K::Int; patience = K * L)
        name = "BiST_$(L)_$(K)"
        new(name, L, K, 1, patience, SelectorData(L, K))
    end
end

""" Get next index to optimize over """
function pre(f::BiST, X::Matrix{Int})
    if length(f.data.prev_visited) > 0
        i_prev, j_prev = f.data.prev_visited[end]
    else
        i_prev, j_prev = 0, 1
    end
    i = i_prev + 1
    j = j_prev
    if i == f.L + 1
        i = 1
        j += 1
    end
    if j == f.K + 1
        j = 1
    end
    return Vector{Tuple{Int,Int}}([(i, j)])
end

""" Generate log data """
function generate_log(f::BiST)
    log = Dict(
        "index_selector_name" => f.name,
        "L" => f.L,
        "K" => f.K,
        "patience" => f.patience,
    )
    return merge(log, generate_log(f.data))
end
