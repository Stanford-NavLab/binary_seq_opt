
# calculate objective for a given matrix with autocorrelation zero property
function ACZ(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)
    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2] for k = 1:K]...)
    return sum(vec(auto .- (size(X)[1] % 2 == 0 ? 0 : -1)) .^ 2)
end

# form JuMP expression for objective, given correlations
function ACZ(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    @expression(
        model,
        J,
        sum([
            2 * model[:corr][(i, i, k)]^2
            for (i, k) in prob_data.autocorrelation_set if abs(mod(k, prob_data.L)) == 1
        ])
    )
    @objective(model, Min, J)
end
