
# calculate objective for a given matrix
function PSL(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    cross = hcat([
        real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

    return norm(vcat(vec(auto), vec(cross)), Inf)
end

# form JuMP expression for objective, given correlations
function PSL(model::Model, prob_data::SubproblemData)
    @variable(model, t)
    @constraint(model, model[:corr] .<= t)
    @constraint(model, -model[:corr] .<= t)
    @constraint(model, model[:ecorr] .<= t)
    @constraint(model, -model[:ecorr] .<= t)
    @objective(model, Min, t)
end
