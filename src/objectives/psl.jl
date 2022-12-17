
# calculate objective for a given matrix
function PSL(X::Union{Matrix{Int},Adjoint{Int, Matrix{Int}}})
    K = size(X)[2]
    FX = [fft(X[:, k]) for k = 1:K]

    auto = hcat([
        real(ifft(FX[k] .* conj.(FX[k])))[2:end] for k = 1:K]...)
    if K > 1
        cross = hcat([
            real(ifft(FX[i] .* conj.(FX[j]))) for i = 1:K for j = i+1:K]...)

        return norm(vcat(vec(auto), vec(cross)), Inf)
    else
        return norm(vec(auto), Inf)
    end
end

# form JuMP expression for objective, given correlations
function PSL(model::Model, prob_data::SubproblemData)
    # solve PSL minimization problem
    @variable(model, t)
    @constraint(model, model[:corr] .<= t)
    @constraint(model, -model[:corr] .<= t)
    @objective(model, Min, t)

    # subject to current PSL constraint, push correlations towards a boundary
    optimize!(model)
    p = value(t)

    # square_signed = (a) -> sign(a) * a^2
    # form_coef = (a) -> square_signed(rand((-1,1)) * p - a)
    # form_coef = (a) -> square_signed(sign(a + 1e-9rand((-1,1))) * p - a)
    form_coef = (a) -> rand((-1, 1))
    # form_coef = (a) -> sign(a + 0.5rand((-1, 1)))

    c = Dict([
        # (i, j, k) => sign(value(model[:corr][(i,j,k)]) + 1e-9(rand()-0.5))
        (i, j, k) => form_coef(value(model[:corr][(i,j,k)]))
        for (i, j, k) in prob_data.correlation_set
    ])

    fix(t, p)
    @objective(model, Max, sum([
        c[(i,j,k)] * model[:corr][(i,j,k)] 
        for (i, j, k) in prob_data.correlation_set
    ]) / length(prob_data.correlation_set))
end


"""

    choose target: rand((-1, 1)) * p

    rand((-1, 1)) * p - a

    +p
    |
    | a1
    |
    0
    |
    | a2
    |
    -p

    -p - a1 = -|p| - |a1|
    -p - a2 = -|-p - a2|

    +p - a1 = |p-a1|
    +p - a2 = |p| + |a2|
""";