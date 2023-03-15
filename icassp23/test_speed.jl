using Serialization, Dates

include("../binary_seq_opt.jl")

nargs = length(ARGS)
argnames = [("K", Int, 4)]

args = Dict()
for i in eachindex(argnames)
    if nargs >= i
        args[argnames[i][1]] =
            argnames[i][2] == String ? ARGS[i] : parse(argnames[i][2], ARGS[i])
    else
        args[argnames[i][1]] = argnames[i][3]
    end
end

using Gurobi
solver = optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1, "MIPGap" => 1e-12)

results = Dict()
runs = 10

for L in [127, 511, 1023]
    for M in [1, 5, 10, 15, 20, 25, 30]
        if M == 1
            index_selector = BiST(L, args["K"])
        else
            index_selector =
                BiSTExtended(L, args["K"], M; max_columns = 2, randomize_M = false)
        end
        bcd = BCD(index_selector, SOS, solver;)
        results[(L, M)] = []
        for i = 1:runs+1
            push!(results[(L, M)], @elapsed bcd(1))
        end
    end
end

serialize("timing_results.jls", results)
