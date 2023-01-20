
""" Block coordinate descent for binary sequence set optimization """
struct BCD{I,O,S,D}
    X0::Matrix{Int}
    X::Matrix{Int}
    X_best::Matrix{Int}
    obj_best::Vector{Float64}
    iteration_times::Vector{Float64}
    index_selector::I
    objective::O
    solver::S
    stop_if_improved::Bool
    disallow_shifts::Bool
    obj_values::Vector{Float64}
    log_path::String
    log_name::String
    log_freq::Int
    function BCD(
        index_selector::IndexSelector, 
        objective::Function,
        solver;
        X0::Matrix{Int}=randb(index_selector.L, index_selector.K),
        stop_if_improved::Bool = false,
        disallow_shifts::Bool = false,
        log_path::String="", 
        log_name::String="BCD-"* string(objective) * "-" * index_selector.name * "-" * Dates.format(now(), "HH_MM_SS_MS") * ".jls",
        log_freq::Int=1,
    )
        new{typeof(index_selector),
            typeof(objective),
            typeof(solver),
            typeof(log)
        }(
            X0,
            copy(X0),
            copy(X0),
            [objective(X0)],
            Vector{Float64}([]),
            index_selector,
            objective,
            solver,
            stop_if_improved,
            disallow_shifts,
            zeros(0),
            log_path,
            log_name,
            log_freq,
        )
    end
end

""" run BCD """
function (f::BCD)(T::Int; verbose::Bool = true)
    for t=1:T
        # perform BCD step
        stop, obj_val, elapsed_time = step(f, t)
        push!(f.obj_values, obj_val)
        push!(f.iteration_times, elapsed_time)
        log!(f, t)

        # update best, if not using a descent method
        if obj_val < f.obj_best[1]
            f.obj_best .= obj_val
            f.X_best .= copy(f.X)
        end

        if verbose
            @printf "Iteration %d Objective: %f\n" t obj_val
        end

        if stop
            break
        end
    end
end

""" Select indices, perform optimization, update data """
function step(f::BCD, t::Int)
    index_list = pre(f.index_selector, f.X)

    start = time()
    Xnew = solve_bcd_subproblem(t, f.X, index_list, f.objective, f.solver, f.stop_if_improved, f.disallow_shifts)
    elapsed = time() - start
    new_obj = f.objective(Xnew)
    f.X .= Xnew

    stop = post(f.index_selector, new_obj, index_list)
    return stop, new_obj, elapsed
end

""" write log to file """
function log!(f::BCD, t::Int)
    if f.log_path != "" && t % f.log_freq == 0
        selector_log = generate_log(f.index_selector)
        bcd_log = Dict(
            "objective" => string(f.objective),
            "solver" => string(f.solver),
            "X" => f.X,
            "X0" => f.X0,
            "obj_values" => f.obj_values,
            "iteration_times" => f.iteration_times,
            "log_path" => f.log_path,
            "log_name" => f.log_name,
            "log_freq" => f.log_freq,
        )
        log = merge(bcd_log, selector_log)
        serialize(joinpath(f.log_path, f.log_name), log)
    end
end
