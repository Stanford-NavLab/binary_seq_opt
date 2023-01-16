
""" Select M random indices """
struct RandomSampler <: IndexSelector
    name::String
    L::Int
    K::Int
    M::Int
    columnwise_limit::Int
    max_columns::Int
    patience::Int
    data::SelectorData
    boost_col_probs::Bool
    function RandomSampler(L::Int, K::Int, M::Int; 
        patience=typemax(Int), 
        columnwise_limit::Int = typemax(Int),
        max_columns::Int = K,
        boost_col_probs::Bool = false,
    )
        name = "RandomSampler_$(L)_$(K)_$(M)"
        new(name, L, K, M, columnwise_limit, max_columns,
            patience, SelectorData(L, K), boost_col_probs)
    end
end

""" Generate indices to optimize over """
function pre(f::RandomSampler, X::Matrix{Int})
    column_set = Set{Int}()
    inds = Dict{Int,Vector{Int}}(i => collect(1:f.L) for i=1:f.K)
    index_list = Vector{Tuple{Int,Int}}()
    for _=1:f.M
        cols = [k for (k, s) in inds if length(s) > f.L - f.columnwise_limit]

        # boost probability of columns with peak correlation values
        if f.boost_col_probs
            psl = PSL(X)
            FX = [fft(X[:, k]) for k = 1:f.K]

            # find columns achieving peak correlation level
            boost_list = []
            for i in cols
                corrs = hcat([
                    real(ifft(FX[i] .* conj.(FX[j]))) for j = i:f.K]...)
                if any(abs.(corrs) .> psl - 0.5)
                    push!(boost_list, i)
                end
            end

            # boost probability of choosing those columns
            append!(cols, vcat([boost_list 
                for _=1:Int(floor(length(cols) / max(1, length(boost_list))))
                ]...)
            )
        end

        # sample index
        j = rand(cols)
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
function post(f::IndexSelector, 
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
        "index_selector_name" => f.name,
        "L" => f.L,
        "K" => f.K,
        "M" => f.M,
        "patience" => f.patience,
        "columnwise_limit" => f.columnwise_limit,
    )
    return merge(log, generate_log(f.data))
end
