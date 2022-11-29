using LinearAlgebra, SparseArrays, FFTW, JuMP
using ProgressBars, Printf, Random, Statistics, StatsBase
using Distributed, Dates, Serialization

include("src/utils.jl")
include("src/gold_codes.jl")

# block coordinate descent core
include("src/subproblem_data.jl")
include("src/bcd_subproblem.jl")

# load objective functions
include("src/objectives/isl.jl")
include("src/objectives/psl.jl")

# index selection strategies
abstract type IndexSelector end
include("src/index_selector/selector_data.jl")
include("src/index_selector/random_selector.jl")

# block coordinate descent
include("src/bcd.jl")


using Gurobi


# # test

# solver = optimizer_with_attributes(
#     Gurobi.Optimizer,
#     "OutputFlag" => 1,
#     "MIPGap" => 1e-10,
# )
# # solver = nothing

# index_selector = RandomSampler(31, 2, 12, 100; columnwise_limit=6)
# bcd = BCD(index_selector, PSL, Gurobi.Optimizer)

# bcd(1000)

# plot(bcd.obj_values)
