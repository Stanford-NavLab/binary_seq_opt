
""" Create and solve JuMP model for solving BCD subproblem
    X: binary code matrix, size (sequence length, number of sequences)
    index_list: list of tuples of indices of X which are variable
    solver: JuMP-compatible MIP solver
        - Example instantiation for solver:
            solver = optimizer_with_attributes(
                Gurobi.Optimizer,
                "OutputFlag" => 1,
                "MIPGap" => 1e-10,
            )
        - to solve using brute-force, pass solver = nothing
"""
function solve_bcd_subproblem(
    t::Int,
    X::Matrix{Int},
    index_list::Vector{Tuple{Int,Int}},
    obj::Function,
    solver,
    stop_if_improved::Bool,
    disallow_shifts::Bool,
    balanced::Bool,
)
    L, K = size(X)

    # preprocess index list
    prob_data = SubproblemData(X, index_list)

    # form JuMP model and binary variables
    model = Model(solver; add_bridges = false)
    @variable(model, _x[prob_data.index_set], Bin)
    @expression(model, x, 2_x .- 1)
    @variable(model, z[prob_data.quad_index_set])

    # disallow shifted versions of sequences
    if disallow_shifts
        for j in prob_data.variable_cols
            vrows = prob_data.variable_rows[j]
            frows = prob_data.fixed_rows[j]
            Xjk = X[:, j]
            for l = j:j
                for k in (j == l ? (1:Int(floor(L / 2))) : (0:L-1))
                    Xlk = circshift(X[:, l], -k)
                    @constraint(
                        model,
                        sum([Xlk[i] * x[(i, j)] for i in vrows]) + sum([Xlk[i] * Xjk[i] for i in frows]) <=
                        L - 1
                    )
                end
            end
        end
    end

    # handle balance constraints
    if balanced
        for j in prob_data.variable_cols
            if L % 2 == 0
                @constraint(model, sum([i in prob_data.variable_rows[j] ? x[(i, j)] : X[i,j] for i=1:L]) == 0)
            else
                @constraint(model, sum([i in prob_data.variable_rows[j] ? x[(i, j)] : X[i,j] for i=1:L]) <= 1)
                @constraint(model, sum([i in prob_data.variable_rows[j] ? x[(i, j)] : X[i,j] for i=1:L]) >= -1)
            end
        end
    end

    # generate linking constraints for auxiliary variables
    for (ij, j, ik, k) in prob_data.quad_index_set
        @constraint(model, z[(ij, j, ik, k)] <= x[(ij, j)] - x[(ik, k)] + 1)
        @constraint(model, z[(ij, j, ik, k)] <= x[(ik, k)] - x[(ij, j)] + 1)
        @constraint(model, z[(ij, j, ik, k)] >= -x[(ij, j)] - x[(ik, k)] - 1)
        @constraint(model, z[(ij, j, ik, k)] >= x[(ij, j)] + x[(ik, k)] - 1)
    end

    # correlations between fixed column j and variable column i at fixed indices
    function constant_corr_vector(i::Int, j::Int)
        YSc = hcat([circshift(reverse(X[:, j]), k) for k in prob_data.fixed_rows[i]]...)
        xsc = X[prob_data.fixed_rows[i], i]
        return YSc * xsc
    end

    Y = Dict(
        (i, j) => constant_corr_vector(i, j) for i in prob_data.variable_cols for
        j in prob_data.fixed_cols
    )

    # build internal correlation expressions, populate with current values
    @expression(
        model,
        corr[(i, j, k) in prob_data.correlation_set],
        (i, j, k) in prob_data.variable_correlation_set ?
        AffExpr(dot(X[:, i], circshift(X[:, j], k))) : AffExpr(Y[(i, j)][k+1])
    )

    # correlation between columns i and j at shift k
    for (i, j, k) in prob_data.variable_correlation_set
        for _l = 1:L
            _m = mod1(_l - k, L)
            # if autocorrelation, sort (_l, _m) to ensure no duplicate
            l, m = i == j ? (min(_l, _m), max(_l, _m)) : (_l, _m)

            # replace constant terms with affine expressions
            if (l, i, m, j) in prob_data.quad_index_set
                base = X[l, i] * X[m, j]
                add_to_expression!(corr[(i, j, k)], z[(l, i, m, j)] - base)
            elseif (l, i) in prob_data.index_set
                base = X[l, i] * X[m, j]
                add_to_expression!(corr[(i, j, k)], x[(l, i)] * X[m, j] - base)
            elseif (m, j) in prob_data.index_set
                base = X[l, i] * X[m, j]
                add_to_expression!(corr[(i, j, k)], X[l, i] * x[(m, j)] - base)
            end
        end
    end

    # build external correlation expressions (affine expressions), if any
    if length(prob_data.fixed_cols) > 0
        # column i contains variables, column j is fixed
        for i in prob_data.variable_cols
            for j in prob_data.fixed_cols
                # create reduced Y matrix
                Y = hcat(
                    [circshift(reverse(X[:, j]), k) for k in prob_data.variable_rows[i]]...,
                )
                for k = 0:L-1
                    add_to_expression!(
                        corr[(i, j, k)],
                        dot(Y[k+1, :], [x[(l, i)] for l in prob_data.variable_rows[i]]),
                    )
                end
            end
        end
    end

    # form objective and solve
    obj(model, prob_data, t, X, stop_if_improved)
    solver_time = @elapsed optimize!(model)

    # return optimized matrix
    Xnew = copy(X)
    if solution_summary(model).result_count > 0
        for (i, j) in index_list
            Xnew[i, j] = Int(round(value(x[(i, j)])))
        end
    end

    # if !(solution_summary(model).objective_value ≈ obj(Xnew))
    #     display((solution_summary(model).objective_value, obj(Xnew)))
    # end
    # @assert solution_summary(model).objective_value ≈ obj(Xnew)

    return Xnew, solver_time
end

""" solve subproblem via brute force """
function solve_bcd_subproblem(
    t::Int,
    X::Matrix{Int},
    index_list::Vector{Tuple{Int,Int}},
    obj::Function,
    solver::Nothing,
    stop_if_improved::Bool,
    disallow_shifts::Bool,
    balanced::Bool,
)
    start = time()
    L, K = size(X)
    N = length(index_list)
    index_vec = [i + L * (j - 1) for (i, j) in index_list]
    best_obj = Inf
    best_X = copy(X)

    var_cols = Set([j for (_, j) in index_list])
    FX = [fft(X[:, k]) for k = 1:K]
    J_curr = obj(X, K = K, FX = FX)
    for i in ProgressBar(0:2^N-1)
        s = bitstring(i)[end-N+1:end]
        x_temp = vec(copy(X))
        x_temp[index_vec] = [2 * Int(b == '1') - 1 for b in s]
        X_temp = reshape(x_temp, (L, K))
        for col in var_cols
            FX[col] = fft(X_temp[:, col])
        end
        J = obj(X_temp, K = K, FX = FX)
        if J < best_obj
            best_X .= X_temp
            best_obj = J
        end
        if stop_if_improved && J < J_curr
            return best_X
        end
        clear!(:s)
        clear!(:x_temp)
        clear!(:X_temp)
    end
    clear!(:var_cols)
    clear!(:FX)
    return best_X, time() - start
end
