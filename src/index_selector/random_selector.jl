
""" Select M random indices """
struct RandomSampler <: IndexSelector
    L::Int
    K::Int
    M::Int
    columnwise_limit::Int
    max_columns::Int
    patience::Int
    data::SelectorData
    function RandomSampler(L::Int, K::Int, M::Int; 
        patience=typemax(Int), 
        columnwise_limit::Int = typemax(Int),
        max_columns::Int = K,
    )
        new(L, K, M, columnwise_limit, max_columns,
            patience, SelectorData(L, K))
    end
end

""" Generate indices to optimize over """
function pre(f::RandomSampler)
    column_set = Set{Int}()
    inds = Dict{Int,Vector{Int}}(i => collect(1:f.L) for i=1:f.K)
    index_list = Vector{Tuple{Int,Int}}()
    for _=1:f.M
        j = rand([k for (k, s) in inds if length(s) > f.L - f.columnwise_limit])
        i = rand(1:length(inds[j]))
        push!(column_set, j)
        push!(index_list, (inds[j][i], j))
        deleteat!(inds[j], i)
        if length(column_set) >= f.max_columns
            for col=1:f.K
                if !(col in column_set)
                    delete!(inds, col)
                end
            end
        end
    end
    return index_list
end

""" update with results, return true if terminated """
function post(f::RandomSampler, 
    new_objective::Float64, 
    index_list::Vector{Tuple{Int,Int}}
)::Bool
    update!(f.data, new_objective, index_list)
    stop = num_iters_same(f.data) >= f.patience
    return stop
end

""" Generate log data """
function generate_log(f::RandomSampler)
    log = Dict(
        "index_selector_name" => "RandomSampler",
        "L" => f.L,
        "K" => f.K,
        "M" => f.M,
        "patience" => f.patience,
        "columnwise_limit" => columnwise_limit,
    )
    return merge(log, generate_log(f.data))
end
