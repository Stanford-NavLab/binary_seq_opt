
# calculate objective for a given matrix
function LPSL(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    return PSL(X)
end

# form JuMP expression for objective, given correlations
function LPSL(model::Model, prob_data::SubproblemData, iter::Int, X::Matrix{Int})
    @variable(model, t)
    @constraint(model, model[:corr] .<= t)
    @constraint(model, -model[:corr] .<= t)
    # @objective(model, Min, t)

    # if iter % 5 == 0
    #     return
    # end

    p = PSL(X)
    # @constraint(model, t <= p)

    @variable(model, _v[prob_data.correlation_set], Bin)
    @expression(model, v, 2_v .- 1)

    @variable(model, vub[prob_data.correlation_set])
    @constraint(model, model[:corr] .- v * p .<= 0.1vub)
    @constraint(model, -model[:corr] .+ v * p .<= vub)
    @objective(model, Min, t + mean(vub))

    # @objective(model, Min, mean([
    #     (i == j ? 2prob_data.K : 1.0) * vub[(i, j, k)]^2 for (i, j, k) in prob_data.correlation_set
    # ]))

end
