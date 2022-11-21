
# calculate objective for a given matrix
function ISL(X::Matrix{Int})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    cross = hcat([
        real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

    return mean(
        vcat(
            2vec(auto) .^ 2,
            2 * vec(cross[1, :]) .^ 2,
            2 * 2vec(cross[2:end, :]) .^ 2,
        ),
    )
end

# form JuMP expression for objective, given correlations
function ISL(model::Model, prob_data::SubproblemData)
    @expression(model, corr2, sum([
        (k == 0 ? 1 : 2) * model[:corr][(i, j, k)] ^ 2 
        for (i, j, k) in prob_data.variable_correlation_set]))

    if length(prob_data.fixed_cols) > 0
        @expression(model, ecorr2, sum([
            (k == 0 ? 1 : 2) * model[:ecorr][(i, j, k)] ^ 2 
            for (i, j, k) in prob_data.fixed_correlation_set]))
    else
        @expression(model, ecorr2, 0)
    end

    @objective(model, Min, mean(vcat(corr2, ecorr2)))
end
