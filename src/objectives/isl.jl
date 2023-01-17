
# calculate objective for a given matrix
function ISL(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([
            real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return sum(
            vcat(
                2vec(auto) .^ 2,
                2vec(cross[1, :]) .^ 2,
                4vec(cross[2:end, :]) .^ 2,
            ),
        )
    else
        return sum(2vec(auto) .^ 2)
    end
end

# form JuMP expression for objective, given correlations
function ISL(model::Model, prob_data::SubproblemData, t::Int, X::Matrix{Int}, stop_if_improved::Bool)
    @expression(model, J, 2sum([
        (k == 0 ? 1 : 2) * model[:corr][(i, j, k)] ^ 2 
        for (i, j, k) in prob_data.correlation_set
    ]))

    @objective(model, Min, J)

    # if stop_if_improved
    #     @constraint(model, t/ PSL(X) <= 1 - 1e-10)
    # else
    #     @objective(model, Min, t)
    # end    
end
