
# calculate objective for a given matrix
function MSOS(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)

    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return max(mean(vec(auto) .^ 2), mean(vec(cross) .^ 2))
    else
        return mean(vec(auto) .^ 2)
    end
end

# form JuMP expression for objective, given correlations
function MSOS(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    # populate missing correlation terms
    missing_ac = []
    missing_cc = []
    missing_vals = Dict()
    fixed_cols = sort([c for c in prob_data.fixed_cols])
    for i in fixed_cols
        append!(missing_ac, [(i, i, k) for k = 1:Int(floor(prob_data.L / 2))])
        auto = acor(prob_data.X[:, i])
        for k = 1:Int(floor(prob_data.L / 2))
            missing_vals[(i, i, k)] = auto[k+1]
        end
        for j in fixed_cols
            if j > i
                append!(missing_cc, [(i, j, k) for k = 0:prob_data.L-1])
                cross = xcor(prob_data.X[:, i], prob_data.X[:, j])
                for k = 0:prob_data.L-1
                    missing_vals[(i, j, k)] = cross[k+1]
                end
            end
        end
    end
    @expression(
        model,
        corr_ac_other[(i, i, k) in missing_ac],
        AffExpr(missing_vals[(i, i, k)])
    )
    @expression(
        model,
        corr_cc_other[(i, j, k) in missing_cc],
        AffExpr(missing_vals[(i, j, k)])
    )

    # define objective
    @variable(model, t)
    if prob_data.L % 2 == 0
        # even length: double autocorrelations except at shift L / 2
        @expression(
            model,
            autocorrelation,
            sum([
                ((k == Int(prob_data.L / 2)) ? 2 : 1) * model[:corr][(i, i, k)]^2 for
                (i, k) in prob_data.autocorrelation_set
            ]) + sum(
                vcat(
                    0,
                    [
                        ((k == Int(prob_data.L / 2)) ? 2 : 1) *
                        model[:corr_ac_other][(i, i, k)]^2 for (i, i, k) in missing_ac
                    ],
                ),
            )
        )
        auto_terms = sum([
            ((k == Int(prob_data.L / 2)) ? 2 : 1) for
            (i, k) in prob_data.autocorrelation_set
        ])
        +sum([((k == Int(prob_data.L / 2)) ? 2 : 1) for (i, i, k) in missing_ac])
    else
        # odd length: double autocorrelations
        @expression(
            model,
            autocorrelation,
            sum([
                2 * model[:corr][(i, i, k)]^2 for (i, k) in prob_data.autocorrelation_set
            ]) + sum(
                vcat(
                    0,
                    [2 * model[:corr_ac_other][(i, i, k)]^2 for (i, j, k) in missing_ac],
                ),
            )
        )
        auto_terms = 2 * length(prob_data.autocorrelation_set) + 2 * length(missing_ac)
    end
    @expression(
        model,
        crosscorrelation,
        sum([model[:corr][(i, j, k)]^2 for (i, j, k) in prob_data.crosscorrelation_set]) + sum(vcat(0, [model[:corr_cc_other][(i, j, k)]^2 for (i, j, k) in missing_cc]))
    )
    cross_terms = length(prob_data.crosscorrelation_set) + length(missing_cc)

    @constraint(model, autocorrelation / auto_terms <= t)
    @constraint(model, crosscorrelation / cross_terms <= t)
    @objective(model, Min, t)
end
