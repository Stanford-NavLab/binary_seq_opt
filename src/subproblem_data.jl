
""" return useful data structures """
function parse_index_list(L::Int, K::Int, index_list::Vector{Tuple{Int,Int}})
    # variable column => variable indices
    variable_dict = Dict{Int,Vector{Int}}()
    for (i, j) in index_list
        if j in keys(variable_dict)
            push!(variable_dict[j], i)
        else
            variable_dict[j] = [i]
        end
    end

    # variable and fixed columns
    variable_columns = Set(keys(variable_dict))
    parameter_columns = Set([i for i = 1:K if !(i in variable_columns)])

    # ensure lists are sorted
    for col in variable_columns
        sort!(variable_dict[col])
    end

    # variable column => fixed indices
    parameter_dict = Dict{Int,Vector{Int}}(k => [i for i = 1:L] for k in variable_columns)
    for j in variable_columns
        deleteat!(parameter_dict[j], variable_dict[j])
    end

    return variable_dict, parameter_dict, variable_columns, parameter_columns
end

""" quadratic index set: () """
function form_quadratic_indices(
    variable_columns::Set{Int},
    variable_dict::Dict{Int,Vector{Int}},
)
    # quadratic auto and cross correlation terms
    quad_indices = Vector{Tuple{Int,Int,Int,Int}}()
    for j in variable_columns
        rows_j = variable_dict[j]
        append!(
            quad_indices,
            [(r1, j, r2, j) for r1 in rows_j for r2 in rows_j if r2 >= r1],
        )

        for k in [_k for _k in variable_columns if _k > j]
            rows_k = variable_dict[k]
            append!(quad_indices, [(r1, j, r2, k) for r1 in rows_j for r2 in rows_k])
        end
    end
    return Set(quad_indices)
end

""" calculate objective and fixed correlation indices """
function form_correlation_indices(L::Int, variable_columns::Set{Int})
    correlation_indices = Vector{Tuple{Int,Int,Int}}()
    for i in variable_columns
        append!(correlation_indices, [(i, i, k) for k = 1:Int(floor(L / 2))])
        for j in variable_columns
            if j > i
                append!(correlation_indices, [(i, j, k) for k = 0:L-1])
            end
        end
    end
    return Set(correlation_indices)
end

""" Convert index list into more useful data structure 
    # query index set
    # query variable columns
    # query purely fixed columns
    # query column => indices of variables
    # query column => indices of fixed parameters
    # query quadratic z indices
    # index set (i,j,k) for all auto and cross correlations
    # index set (i,j,k) for all variable auto and cross correlations
    # index set (i,j,k) for all fixed auto and cross correlations
"""
struct SubproblemData
    L::Int
    K::Int
    index_set::Set{Tuple{Int,Int}}
    variable_cols::Set{Int}
    fixed_cols::Set{Int}
    variable_rows::Dict{Int,Vector{Int}}
    fixed_rows::Dict{Int,Vector{Int}}
    quad_index_set::Set{Tuple{Int,Int,Int,Int}}
    correlation_set::Set{Tuple{Int,Int,Int}}
    variable_correlation_set::Set{Tuple{Int,Int,Int}}
    fixed_correlation_set::Set{Tuple{Int,Int,Int}}
    function SubproblemData(X::Matrix{Int}, index_list::Vector{Tuple{Int,Int}})
        L, K = size(X)

        v_dict, p_dict, v_cols, p_cols = parse_index_list(L, K, index_list)
        quad_index_set = form_quadratic_indices(v_cols, v_dict)
        correlation_indices = form_correlation_indices(L, v_cols)
        fixed_correlation_set =
            Set([(i, j, k) for i in v_cols for j in p_cols for k = 0:L-1])
        correlation_set = union(correlation_indices, fixed_correlation_set)

        new(
            L,
            K,
            Set(index_list),
            v_cols,
            p_cols,
            v_dict,
            p_dict,
            quad_index_set,
            correlation_set,
            correlation_indices,
            fixed_correlation_set,
        )
    end
end
