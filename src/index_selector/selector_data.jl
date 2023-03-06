
""" Track data used by index selectors """
struct SelectorData
    not_improved_global::Vector{Int}
    not_improved_columnwise::Vector{Int}
    best_objective::Vector{Float64}
    total_visits::Matrix{Int}
    prev_visited::Vector{Tuple{Int,Int}}
    function SelectorData(L::Int, K::Int)
        new([0], zeros(Int, K), [Inf], zeros(Int, L, K), Vector{Tuple{Int,Int}}())
    end
end

""" Update data with index list and new objective value """
function update!(
    f::SelectorData,
    new_objective::Float64,
    index_list::Vector{Tuple{Int,Int}},
)
    columns = sort(collect(Set([j for (_, j) in index_list])))
    if new_objective < f.best_objective[1]
        f.best_objective[1] = new_objective
        f.not_improved_global[1] = 0
        f.not_improved_columnwise[columns] .= 0
    else
        f.not_improved_global[1] += 1
        f.not_improved_columnwise[columns] .+= 1
    end

    for (i, j) in index_list
        f.total_visits[i, j] += 1
    end

    empty!(f.prev_visited)
    append!(f.prev_visited, index_list)
end

""" Get number of iterations best objective not improved """
function num_iters_same(f::SelectorData)::Int
    return f.not_improved_global[1]
end

""" Get number of iterations best objective not improved for a given column """
function num_col_iters_same(f::SelectorData, column::Int)::Int
    return f.not_improved_columnwise[column]
end

""" Get indices visited in the previous iteration """
function get_previous_indices(f::SelectorData)::Vector{Tuple{Int,Int}}
    return f.prev_visited
end

""" Get total visit counts for all indices """
function get_visit_counts(f::SelectorData)::Matrix{Int}
    return f.total_visits
end

""" Generate log dictionary """
function generate_log(f::SelectorData)
    return Dict(
        "not_improved_global" => f.not_improved_global,
        "not_improved_columnwise" => f.not_improved_columnwise,
        "best_objective" => f.best_objective,
        "total_visits" => f.total_visits,
        "prev_visited" => f.prev_visited,
    )
end
