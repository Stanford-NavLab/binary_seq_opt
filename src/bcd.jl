
""" Block coordinate descent for binary sequence set optimization """
struct BCD{I,O,S,D}
    X::Matrix{Int}
    index_selector::I
    objective::O
    solver::S
    obj_values::Vector{Float64}
    log_path::String
    log_name::String
    log_freq::Int
    function BCD(
        index_selector::IndexSelector, 
        objective::Function,
        solver;
        X0::Matrix{Int}=randb(index_selector.L, index_selector.K),
        log_path::String="", 
        log_name::String="BCD-"* index_selector.name * "-" * Dates.format(now(), "HH_MM_SS_MS") * ".jls",
        log_freq::Int=1,
    )
        new{typeof(index_selector),
            typeof(objective),
            typeof(solver),
            typeof(log)
        }(
            X0,
            index_selector,
            objective,
            solver,
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
        stop, obj_val = step(f)
        push!(f.obj_values, obj_val)
        log!(f, t)

        if verbose
            @printf "Iteration %d Objective: %f\n" t obj_val
        end

        if stop
            break
        end
    end
end

""" Select indices, perform optimization, update data """
function step(f::BCD)
    index_list = pre(f.index_selector)

    Xnew = solve_bcd_subproblem(f.X, index_list, f.objective, f.solver)
    new_obj = f.objective(Xnew)
    f.X .= Xnew

    stop = post(f.index_selector, new_obj, index_list)
    return stop, new_obj
end

""" write log to file """
function log!(f::BCD, t::Int)
    if f.log_path != "" && t % f.log_freq == 0
        selector_log = generate_log(f.index_selector)
        bcd_log = Dict(
            "objective" => string(f.objective),
            "solver" => string(f.solver),
            "X" => f.X,
            "obj_values" => f.obj_values,
            "log_path" => f.log_path,
            "log_name" => f.log_name,
            "log_freq" => f.log_freq,
        )
        log = merge(bcd_log, selector_log)
        serialize(joinpath(f.log_path, f.log_name), log)
    end
end
