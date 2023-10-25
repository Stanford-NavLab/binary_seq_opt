
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
                ((i == j && k != Int(prob_data.L / 2)) ? 2 : 1) * model[:corr][(i, j, k)]^2
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




# calculate objective for a given matrix with autocorrelation zero property
function ACZSOS(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)
    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    # weight = size(X)[1]
    # if size(X)[1] % 2 == 0
    #     auto[1, :] = sqrt(weight) * auto[1, :]
    # else
    #     auto[1, :] = sqrt(weight) * (auto[1, :] .+ 1)
    # end
    if !all(isapprox.(auto[1, :], (size(X)[1] % 2 == 0 ? 0 : -1)))
        return Inf
    end
    if K > 1
        cross = hcat([real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return sum(vcat(vec(auto) .^ 2, vec(cross) .^ 2))
    else
        return sum(vec(auto) .^ 2)
    end
end

# form JuMP expression for objective, given correlations
function ACZSOS(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    # weight = prob_data.L
    for (i, k) in prob_data.autocorrelation_set
        if abs(mod(k, prob_data.L)) == 1
            @constraint(model, model[:corr][(i, i, k)] == (prob_data.L % 2 == 0 ? 0 : -1))
        end
    end

    if prob_data.L % 2 == 0
        # even length: double autocorrelations except at shift L / 2
        # @expression(
        #     model,
        #     J,
        #     sum([
        #         (abs(mod(k, prob_data.L)) == 1 && i == j ? weight : 1) * ((i == j && k != Int(prob_data.L / 2)) ? 2 : 1) * model[:corr][(i, j, k)]^2
        #         for (i, j, k) in prob_data.correlation_set
        #     ])
        # )
        @expression(
            model,
            J,
            sum([
                ((i == j && k != Int(prob_data.L / 2)) ? 2 : 1) * model[:corr][(i, j, k)]^2
                for (i, j, k) in prob_data.correlation_set
            ])
        )
    else
        # odd length: double autocorrelations
        # @expression(
        #     model,
        #     J,
        #     sum([
        #         (i == j ? 2 : 1) * (abs(mod(k, prob_data.L)) == 1 && i == j ? weight * (model[:corr][(i, j, k)] + 1)^2 : model[:corr][(i, j, k)]^2)
        #         for (i, j, k) in prob_data.correlation_set
        #     ])
        # )
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





# create SOS objective for low correlation zone case
function SOS_LCZ_5(
    X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}};
    K = size(X)[2],
    FX = [fft(X[:, k]) for k = 1:K],
)
    lcz_width = Int(ceil(size(X)[1] * 5 / 100))
    auto = hcat([real(ifft(FX[k] .* conj.(FX[k])))[2:lcz_width] for k = 1:K]...)
    if K > 1
        cross = hcat([real(ifft(FX[i] .* conj.(FX[j])))[1:lcz_width] for i = 1:K for j = i+1:K]...)

        return sum(vcat(vec(auto) .^ 2, vec(cross) .^ 2))
    else
        return sum(vec(auto) .^ 2)
    end
end


# form JuMP expression for objective, given correlations
function SOS_LCZ_5(
    model::Model,
    prob_data::SubproblemData,
    t::Int,
    X::Matrix{Int},
    stop_if_improved::Bool,
)
    lcz_width = Int(ceil(size(X)[1] * 5 / 100))
    selected_corr_set = [(i, j, k) for (i, j, k) in prob_data.correlation_set if abs(mod(k, prob_data.L)) < lcz_width]

    if prob_data.L % 2 == 0
        # even length: double autocorrelations except at shift L / 2
        @expression(
            model,
            J,
            sum([
                model[:corr][(i, j, k)]^2
                for (i, j, k) in selected_corr_set
            ])
        )
    else
        # odd length: double autocorrelations
        @expression(
            model,
            J,
            sum([
                model[:corr][(i, j, k)]^2 for
                (i, j, k) in selected_corr_set
            ])
        )
    end
    @objective(model, Min, J)
end
