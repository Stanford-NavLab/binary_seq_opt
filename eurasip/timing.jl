
using Serialization, Dates

include("../binary_seq_opt.jl")

nargs = length(ARGS)
argnames = [
    ("seed", Int, 0),
    ("L", Int, 257),
    ("K", Int, 130),
    ("M", Int, 24),
    ("objective", String, "SOS"),
    ("num_runs", Int, 10),
    ("num_cols_min", Int, 1),
    ("num_cols_step", Int, 1),
    ("num_cols_max", Int, 50),
    ("brute_force", Bool, false),
    ("solver_procs", Int, 2),
    ("solver_time_limit", Float64, 30),
]

args = Dict()
for i in eachindex(argnames)
    if nargs >= i
        args[argnames[i][1]] =
            argnames[i][2] == String ? ARGS[i] : parse(argnames[i][2], ARGS[i])
    else
        args[argnames[i][1]] = argnames[i][3]
    end
end

# create experiment path
results_path = "../results"
if ~ispath(results_path)
    mkpath(results_path)
end

# define solver
if args["brute_force"] || args["M"] <= 4
    solver = nothing
else
    using Gurobi
    solver = optimizer_with_attributes(
        Gurobi.Optimizer,
        "Threads" => args["solver_procs"],
        "TimeLimit" => args["solver_time_limit"],
        "OutputFlag" => 1,
        "MIPGap" => 1e-12,
    )
end

# parse objective function
obj_sym = Symbol(args["objective"])
objective = @eval $obj_sym


# run experiment
results = Dict()
for cols=args["num_cols_min"]:args["num_cols_step"]:args["num_cols_max"]
    t = 0
    for r=1:args["num_runs"]
        Random.seed!(args["seed"]+r)
        X0 = randb(args["L"], args["K"])
    
        index_selector = RandomSampler(
            args["L"],
            args["K"],
            args["M"];
            max_columns = cols,
            randomize_M = false,
        )
        bcd = BCD(index_selector, objective, solver; X0 = X0)
        t += @elapsed bcd(1)
    end
    results[cols] = t / args["num_runs"]
end
display(results)

serialize(joinpath(results_path, "timing_results_$(args["L"])_$(args["K"])_$(args["M"]).jls"), results)

# col = collect(args["num_cols_min"]:args["num_cols_step"]:args["num_cols_max"])
# times = [results[c] for c in col]

# using Plots
# plot(col, times, show = true, marker=:circle)
# readline()
