
# calculate objective for a given matrix
function SOS(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([
            real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return mean(vcat(vec(auto) .^ 2, vec(cross) .^ 2))
    else
        return mean(vec(auto) .^ 2)
    end
end

# form JuMP expression for objective, given correlations
function SOS(model::Model, prob_data::SubproblemData, t::Int, X::Matrix{Int}, stop_if_improved::Bool)
    @expression(model, J, sum([
        model[:corr][(i, j, k)] ^ 2 
        for (i, j, k) in prob_data.correlation_set
    ]))

    @objective(model, Min, J)
end
