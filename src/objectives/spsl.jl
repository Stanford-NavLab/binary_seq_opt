
# calculate objective for a given matrix
function SPSL(X::Union{Matrix{Int},Adjoint{Int,Matrix{Int}}})
    return PSL(X)
end

# form JuMP expression for objective, given correlations
function SPSL(model::Model, prob_data::SubproblemData, iter::Int, X::Matrix{Int})
    # identify where peak is achieved
    K = size(X)[2]
    L_ac = Int(floor(prob_data.L / 2))
    FX = [fft(X[:, k]) for k = 1:K]

    top2psl = zeros(2)
    for i = 1:K
        ac = abs.(real(ifft(FX[i] .* conj.(FX[i])))[2:L_ac+1])
        top2 = ac[partialsortperm(ac, 1:2, rev = true)]
        top2psl .= max(top2psl, top2)
        for j = i+1:K
            cc = abs.(real(ifft(FX[i] .* conj.(FX[j])))[1:end])
            top2 = cc[partialsortperm(cc, 1:2, rev = true)]
            top2psl .= max(top2psl, top2)
        end
    end

    correlation_subset = []
    for i = 1:K
        ac = abs.(real(ifft(FX[i] .* conj.(FX[i])))[2:L_ac+1])
        shifts = collect(1:length(ac))[ac.>top2psl[1]-0.5]
        append!(correlation_subset, [(i, i, k) for k in shifts])
        for j = i+1:K
            cc = abs.(real(ifft(FX[i] .* conj.(FX[j])))[1:end])
            shifts = collect(0:prob_data.L-1)[cc.>top2psl[1]-0.5]
            append!(correlation_subset, [(i, j, k) for k in shifts])
        end
    end
    correlation_subset =
        intersect(Set(correlation_subset), prob_data.variable_correlation_set)

    display(correlation_subset)
    if length(correlation_subset) == 0
        @variable(model, t)
        @constraint(model, t >= 1)
        @constraint(model, t <= -1)
        @objective(model, Min, 0.0)
        return
    end

    # solve PSL minimization problem over correlation subset
    @variable(model, t)
    @constraint(model, model[:corr][collect(correlation_subset)] .<= t)
    @constraint(model, -model[:corr][collect(correlation_subset)] .<= t)

    @constraint(
        model,
        model[:corr][collect(
            setdiff(prob_data.variable_correlation_set, correlation_subset),
        )] .<= top2psl[1] - 2
    )
    @constraint(
        model,
        -model[:corr][collect(
            setdiff(prob_data.variable_correlation_set, correlation_subset),
        )] .<= top2psl[1] - 2
    )

    @objective(model, Min, t)
end
