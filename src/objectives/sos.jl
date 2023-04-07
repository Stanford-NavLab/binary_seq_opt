
# calculate objective for a given matrix
function SOS(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)

    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return sum(vcat(vec(auto) .^ 2, vec(cross) .^ 2))
    else
        return sum(vec(auto) .^ 2)
    end
end

# form JuMP expression for objective, given correlations
function SOS(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    if prob_data.L % 2 == 0
        # even length: double autocorrelations except at shift L / 2
        @expression(
            model,
            J,
            sum([
                ((i == j && k == Int(prob_data.L / 2)) ? 2 : 1) * model[:corr][(i, j, k)]^2
                for (i, j, k) in prob_data.correlation_set
            ])
        )
    else
        # odd length: double autocorrelations
        @expression(
            model,
            J,
            sum([
                (i == j ? 2 : 1) * model[:corr][(i, j, k)]^2 for
                (i, j, k) in prob_data.correlation_set
            ])
        )
    end
    @objective(model, Min, J)
end
