#=
# only need to run this once to install necessary packages
import Pkg

Pkg.add("Gurobi")
Pkg.build("Gurobi")

Pkg.add("ProgressBars")
Pkg.add("SparseArrays")
Pkg.add("LinearAlgebra")
Pkg.add("FFTW")
Pkg.add("JuMP")
Pkg.add("LinearOperators")
Pkg.add("StatsBase")
=#

using Serialization, Dates

include("../binary_seq_opt.jl")

# parse command line arguments
nargs = length(ARGS)
argnames = [
    ("seed", Int, 0),
    ("path", String, ""),
    ("L", Int, 127),
    ("K", Int, 66),
    ("M", Int, 20),
    ("objective", String, "ACZSOS"),
    ("randomize_M", Bool, false),
    ("max_iter", Int, 10^6),
    ("patience", Int, 10^6),
    ("brute_force", Bool, false),
    ("max_columns", Int, 66),
    ("log_freq", Int, 1),
    ("solver_procs", Int, 2),
    ("solver_time_limit", Float64, Inf),
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

if args["path"] != ""
    # load initial code from previous result
    prev_data = deserialize(args["path"])
    X0 = prev_data["X"]
else
    # generate initial code with random seed
    Random.seed!(args["seed"])
    X0 = randb(args["L"], args["K"])
end

# parse objective function
obj_sym = Symbol(args["objective"])
objective = @eval $obj_sym

# define index selector
if args["M"] == 1
    index_selector = BiST(args["L"], args["K"])
else
    index_selector = BiSTExtended(
        args["L"],
        args["K"],
        args["M"];
        max_columns = args["max_columns"],
        patience = args["patience"],
        randomize_M = args["randomize_M"],
    )
end

# set up and run BCD solver
bcd = BCD(
    index_selector,
    objective,
    solver;
    X0 = X0,
    log_path = results_path,
    log_freq = args["log_freq"],
)

bcd(args["max_iter"])
