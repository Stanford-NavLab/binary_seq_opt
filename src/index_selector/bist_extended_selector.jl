
""" BiST strategy with extra indices """
struct BiSTExtended <: IndexSelector
    name::String
    L::Int
    K::Int
    M::Int
    max_columns::Int
    columnwise_limit::Int
    patience::Int
    data::SelectorData
    prev_visited::Vector{Tuple{Int,Int}}
    randomize_M::Bool
    function BiSTExtended(
        L::Int,
        K::Int,
        M::Int;
        patience = K * L,
        max_columns::Int = 2,
        columnwise_limit::Int = L,
        randomize_M::Bool = true,
    )
        name = "BiSTExtended_$(L)_$(K)_$(M)_$(max_columns)_$(columnwise_limit)"
        new(
            name,
            L,
            K,
            M,
            max_columns,
            columnwise_limit,
            patience,
            SelectorData(L, K),
            Vector{Tuple{Int,Int}}([(0, 1)]),
            randomize_M,
        )
    end
end

""" Get next index to optimize over """
function pre(f::BiSTExtended, X::Matrix{Int})
    # add extra indices
    column_set = Set{Int}()
    inds = Dict{Int,Vector{Int}}(i => collect(1:f.L) for i = 1:f.K)
    index_list = Vector{Tuple{Int,Int}}()
    M = f.randomize_M ? rand(1:f.M) : f.M
    for m = 1:M
        if m == 1
            # BiST sequential indices
            i_prev, j_prev = f.prev_visited[end]
            i = i_prev + 1
            j = j_prev
            if i == f.L + 1
                i = 1
                j += 1
            end
            if j == f.K + 1
                j = 1
            end
            f.prev_visited[end] = (i, j)
        else
            # extra indices
            cols = [k for (k, s) in inds if length(s) > f.L - f.columnwise_limit]
            j = rand(cols)
            i = rand(1:length(inds[j]))
        end
        # add index
        push!(column_set, j)
        push!(index_list, (inds[j][i], j))
        deleteat!(inds[j], i)
        if length(column_set) >= f.max_columns
            for col = 1:f.K
                if !(col in column_set)
                    delete!(inds, col)
                end
            end
        end
    end
    return index_list
end

""" Generate log data """
function generate_log(f::BiSTExtended)
    log = Dict(
        "index_selector_name" => f.name,
        "L" => f.L,
        "K" => f.K,
        "M" => f.M,
        "patience" => f.patience,
        "columnwise_limit" => f.columnwise_limit,
        "max_columns" => f.max_columns,
    )
    return merge(log, generate_log(f.data))
end
