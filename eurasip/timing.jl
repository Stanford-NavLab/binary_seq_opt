
using Serialization, Dates

include("../binary_seq_opt.jl")

nargs = length(ARGS)
argnames = [
    ("seed", Int, 0),
    ("L", Int, 127),
    ("K", Int, 66),
    ("M", Int, 25),
    ("objective", String, "SOS"),
    ("num_runs", Int, 1),
    ("num_cols_min", Int, 2),
    ("num_cols_step", Int, 10),
    ("num_cols_max", Int, 66),
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
using Gurobi
solver = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => args["solver_procs"],
    "TimeLimit" => args["solver_time_limit"],
    "OutputFlag" => 1,
    "MIPGap" => 1e-12,
)

# parse objective function
obj_sym = Symbol(args["objective"])
objective = @eval $obj_sym

# run experiment
elapsed_times = Dict()
solver_times = Dict()
for cols=args["num_cols_min"]:args["num_cols_step"]:args["num_cols_max"]
    t = 0
    s = 0
    for r=1:args["num_runs"]
        Random.seed!(args["seed"]+r)
        X0 = randb(args["L"], args["K"])
    
        index_selector = RandomSampler(
            args["L"],
            args["K"],
            args["M"];
            columnwise_limit = Int(ceil(args["M"] / cols)),
            max_columns = cols,
            randomize_M = false,
        )
        bcd = BCD(index_selector, objective, solver; X0 = X0)
        bcd(1)
        t += bcd.iteration_times[end]
        s += bcd.solver_times[end]
    end
    elapsed_times[cols] = t / args["num_runs"]
    solver_times[cols] = s / args["num_runs"]
end

path = joinpath(results_path, "timing_results_$(args["L"])_$(args["K"])_$(args["M"]).jls")

cols = collect(args["num_cols_min"]:args["num_cols_step"]:args["num_cols_max"])
results = Dict(
    "columns" => cols,
    "elapsed_times" => [elapsed_times[k] for k in cols],
    "solver_times" => [solver_times[k] for k in cols],
)
serialize(path, results)

# col = collect(args["num_cols_min"]:args["num_cols_step"]:args["num_cols_max"])
# times = [results[c] for c in col]

# using Plots
# plot(col, times, show = true, marker=:circle)
# readline()
