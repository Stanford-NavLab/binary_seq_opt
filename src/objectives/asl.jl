
# calculate objective for a given matrix
function ASL(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)

    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return sum(
            vcat(2abs.(vec(auto)), 2abs.(vec(cross[1, :])), 4abs.(vec(cross[2:end, :]))),
        )
    else
        return sum(2abs.(vec(auto)))
    end
end

# form JuMP expression for objective, given correlations
function ASL(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    @variable(model, abs_corr[prob_data.correlation_set])
    @constraint(model, model[:corr] .<= abs_corr)
    @constraint(model, -model[:corr] .<= abs_corr)
    @objective(
        model,
        Min,
        2sum([
            (k == 0 ? 1 : 2) * abs_corr[(i, j, k)] for
            (i, j, k) in prob_data.correlation_set
        ])
    )
end
