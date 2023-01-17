
# calculate objective for a given matrix
function MAE(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([
            real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return sum(abs.(vcat(vec(auto), vec(cross))))
    else
        return sum(abs.(auto))
    end
end

# form JuMP expression for objective, given correlations
function MAE(model::Model, prob_data::SubproblemData, t::Int, X::Matrix{Int}, stop_if_improved::Bool)
    @variable(model, abs_corr[prob_data.correlation_set])
    @constraint(model, model[:corr] .<= abs_corr)
    @constraint(model, -model[:corr] .<= abs_corr)

    if prob_data.L % 2 == 0
        # even length: double autocorrelations except at shift L / 2
        @expression(model, J, sum([
            ((i == j && k == Int(prob_data.L / 2)) ? 2 : 1) * abs_corr[(i, j, k)]
            for (i, j, k) in prob_data.correlation_set
        ]))
    else
        # odd length: double autocorrelations
        @expression(model, J, sum([
            (i == j ? 2 : 1) * abs_corr[(i, j, k)]
            for (i, j, k) in prob_data.correlation_set
        ]))
    end
    @objective(model, Min, J)
end
