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
include("src/objectives/sos.jl")
include("src/objectives/asl.jl")
include("src/objectives/mae.jl")
include("src/objectives/psl.jl")
include("src/objectives/spsl.jl")
include("src/objectives/mpsl.jl")
include("src/objectives/lpsl.jl")

# index selection strategies
abstract type IndexSelector end
include("src/index_selector/selector_data.jl")
include("src/index_selector/random_selector.jl")

# block coordinate descent
include("src/bcd.jl")


# using Gurobi
# using Plots

# solver = optimizer_with_attributes(
#     Gurobi.Optimizer,
#     "OutputFlag" => 1,
#     "MIPGap" => 1e-10,
# )
# # solver = nothing


# """ check that solver & brute force agree """
# function test_bcd(L::Int, K::Int, M::Int, obj::Function)
#     solver = optimizer_with_attributes(
#         Gurobi.Optimizer,
#         "OutputFlag" => 1,
#         "MIPGap" => 1e-10,
#     )
#     X0 = randb(L, K)
#     index_selector = RandomSampler(L, K, M)
#     index_list = pre(index_selector)
#     Xsolve = solve_bcd_subproblem(X0, index_list, obj, solver)
#     Xbrute = solve_bcd_subproblem(X0, index_list, obj, nothing)
#     return obj(Xsolve) ≈ obj(Xbrute)
# end
# [@assert test_bcd(31, 13, 7, ISL) for _=1:1000]


# Random.seed!(0)
# X0 = randb(31, 2)

# index_selector = RandomSampler(31, 11, 7; columnwise_limit=8)
# bcd = BCD(index_selector, ISL, solver)

# bcd(100)

# using Plots
# plot(bcd.obj_values)
