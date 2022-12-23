
# calculate objective for a given matrix
function PSL(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([
            real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return norm(vcat(vec(auto), vec(cross)), Inf)
    else
        return norm(vec(auto), Inf)
    end
end

# form JuMP expression for objective, given correlations
function PSL(model::Model, prob_data::SubproblemData, iter::Int, X::Matrix{Int})
    # solve PSL minimization problem
    @variable(model, t)
    @constraint(model, model[:corr] .<= t)
    @constraint(model, -model[:corr] .<= t)
    @objective(model, Min, t)

    # if iter % 5 == 0
    #     return
    # end

    # optimize!(model)
    # p = value(t)
    # xvals = Dict((i, j) => value.(model[:_x][(i, j)]) for (i, j) in prob_data.index_set)

    # # # maximize L1 subject to PSL constraint
    # @variable(model, _v[prob_data.correlation_set], Bin)
    # @expression(model, v, 2_v .- 1)

    # @variable(model, vub[prob_data.correlation_set])
    # @constraint(model, model[:corr] .- v * p .<= vub)
    # @constraint(model, -model[:corr] .+ v * p .<= vub)

    # @objective(model, Min, mean([
    #     (i == j ? prob_data.K : 1.0) * vub[(i, j, k)]^2 for (i, j, k) in prob_data.correlation_set
    # ]))

    # fix(t, p)
    # for (i, j) in prob_data.index_set
    #     fix(model[:_x][(i, j)], xvals[(i, j)])
    # end

    # # one index at a time?
    # for (i, j) in sample(collect(prob_data.index_set), min(8, length(prob_data.index_set)); replace=false)
    #     unfix(model[:_x][(i, j)])
    #     optimize!(model)
    #     xval_new = value(model[:_x][(i, j)])
    #     fix(model[:_x][(i, j)], xval_new)
    # end
    # @objective(model, Min, 0.0)
end
