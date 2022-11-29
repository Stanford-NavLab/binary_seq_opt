
""" Create and solve JuMP model """
function solve_bcd_subproblem(
    X::Matrix{Int}, index_list::Vector{Tuple{Int,Int}}, obj::Function, solver
)
    L, _ = size(X)

    # preprocess index list
    prob_data = SubproblemData(X, index_list)

    # form JuMP model and binary variables
    model = Model(solver; add_bridges = false)
    @variable(model, _x[prob_data.index_set], Bin)
    @expression(model, x, 2_x .- 1)
    @variable(model, z[prob_data.quad_index_set])

    # generate linking constraints for auxiliary variables
    for (ij, j, ik, k) in prob_data.quad_index_set
        @constraint(model, z[(ij, j, ik, k)] <= x[(ij, j)] - x[(ik, k)] + 1)
        @constraint(model, z[(ij, j, ik, k)] <= x[(ik, k)] - x[(ij, j)] + 1)
        @constraint(model, z[(ij, j, ik, k)] >= -x[(ij, j)] - x[(ik, k)] - 1)
        @constraint(model, z[(ij, j, ik, k)] >= x[(ij, j)] + x[(ik, k)] - 1)
    end

    # build internal correlation expressions
    @expression(model, corr[(i, j, k) in prob_data.variable_correlation_set], 
                AffExpr(dot(X[:, i], circshift(X[:, j], k))))
    for (i, j, k) in prob_data.variable_correlation_set
        for l = 1:L
            m = mod1(l - k, L)
            if (l, i, m, j) in prob_data.quad_index_set
                base = X[l, i] * X[m, j]
                add_to_expression!(corr[(i, j, k)], -base)
                add_to_expression!(corr[(i, j, k)], z[(l, i, m, j)])
            elseif (l, i) in prob_data.index_set
                base = X[l, i] * X[m, j]
                add_to_expression!(corr[(i, j, k)], -base)
                add_to_expression!(corr[(i, j, k)], x[(l, i)] * X[m, j])
            elseif (m, j) in prob_data.index_set
                base = X[l, i] * X[m, j]
                add_to_expression!(corr[(i, j, k)], -base)
                add_to_expression!(corr[(i, j, k)], X[l, i] * x[(m, j)])
            end
        end
    end

    # build external correlation expressions, if any
    if length(prob_data.fixed_cols) > 0
        function constant_corr_vector(i::Int, j::Int)
            YSc = hcat([circshift(reverse(X[:, j]), k) 
                        for k in prob_data.fixed_rows[i]]...)
            xsc = X[prob_data.fixed_rows[i], i]
            return YSc * xsc
        end
        Y = Dict((i, j) => constant_corr_vector(i, j)
            for i in prob_data.variable_cols for j in prob_data.fixed_cols
        )

        @expression(model, 
                    ecorr[(i, j, k) in prob_data.fixed_correlation_set], 
                    AffExpr(Y[(i, j)][k + 1]))

        for i in prob_data.variable_cols
            for j in prob_data.fixed_cols
                # create reduced Y matrix
                Y = hcat([circshift(reverse(X[:, j]), k) 
                          for k in prob_data.variable_rows[i]]...)

                for k=0:L-1
                    add_to_expression!(
                        ecorr[(i, j, k)], 
                        dot(Y[k + 1, :],
                        [x[(l, i)] for l in prob_data.variable_rows[i]]))
                end

            end
        end
    else
        @expression(model, ecorr, 0)
    end

    # form objective and solve
    obj(model, prob_data)
    optimize!(model)

    Xnew = copy(X)
    for (i, j) in index_list
        Xnew[i, j] = Int(round(value(x[(i, j)])))
    end

    return Xnew
end

""" solve subproblem via brute force """
function solve_bcd_subproblem(
    X::Matrix{Int}, index_list::Vector{Tuple{Int,Int}}, 
    obj::Function, solver::Nothing
)
    L, K = size(X)
    N = length(index_list)
    index_vec = [i + L * (j - 1) for (i, j) in index_list]
    best_obj = Inf
    best_X = copy(X)
    for i in ProgressBar(0:2^N-1)
        s = bitstring(i)[end-N+1:end]
        x_temp = vec(copy(X))
        x_temp[index_vec] = [2 * Int(b == '1') - 1 for b in s]
        X_temp = reshape(x_temp, (L, K))
        J = obj(X_temp)
        if J < best_obj
            best_X .= X_temp
            best_obj = J
        end
    end
    return best_X
end
