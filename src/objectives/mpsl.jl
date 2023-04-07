
# calculate objective for a given matrix
function MPSL(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)

    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return sum(
            vcat(
                [norm(auto[:, k], Inf) for k = 1:K],
                [2norm(cross[:, k], Inf) for k = 1:size(cross)[2]],
            ),
        ) / (K + 2size(cross)[2])
    else
        return mean([norm(auto[:, k], Inf) for k = 1:K])
    end
end

# form JuMP expression for objective, given correlations
function MPSL(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    # average of peak correlations
    inds = collect(Set([(i, j) for (i, j, _) in prob_data.correlation_set]))
    @variable(model, t[inds])
    for (i, j, k) in prob_data.correlation_set
        @constraint(model, model[:corr][(i, j, k)] <= t[(i, j)])
        @constraint(model, -model[:corr][(i, j, k)] <= t[(i, j)])
    end

    if prob_data.K > 1
        @expression(
            model,
            J,
            sum([(i == j ? 1.0 : 2.0) * t[(i, j)] for (i, j) in inds]) / (prob_data.K + (prob_data.K - 1) * prob_data.K)
        )
    else
        @expression(model, J, mean([t[(i, i)] for i = 1:prob_data.K]))
    end

    if stop_if_improved
        @constraint(model, J / MPSL(X) <= 1 - 1e-10)
    else
        @objective(model, Min, J)
    end
end
