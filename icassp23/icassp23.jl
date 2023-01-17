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
    ("L", Int, 63),
    ("K", Int, 4),
    ("M", Int, 1),
    ("objective", String, "ISL"),
    ("max_iter", Int, 63 * 4),
    ("log_freq", Int, 1),
    ("solver_procs", Int, 2),
    ("brute_force", Bool, false),
    ("max_columns", Int, 2),
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
if args["brute_force"]
    solver = nothing
else
    using Gurobi
    solver = optimizer_with_attributes(
        Gurobi.Optimizer,
        "Threads" => args["solver_procs"],
        "TimeLimit" => args["solver_time_limit"],
        "OutputFlag" => 1,
        "MIPGap" => 1e-10,
    )
end

# generate initial code
Random.seed!(args["seed"])
X0 = randb(args["L"], args["K"])

# parse objective function
obj_sym = Symbol(args["objective"])
objective = @eval $obj_sym

# define index selector
if args["M"] == 1
    index_selector = BiST(args["L"], args["K"])
else
    index_selector = BiSTExtended(args["L"], args["K"], args["M"]; max_columns =  args["max_columns"])
end

# set up and run BCD solver
bcd = BCD(index_selector, objective, solver; X0=X0,
          log_path=results_path, log_freq=args["log_freq"])

bcd(args["max_iter"])
