using ClusterManagers
using Random

np = parse(Int, ARGS[1])
Sre = parse(Int, ARGS[2])
temp = parse(Int, ARGS[3])

include("../load.jl")
include("../parameters_static.jl")



sim = many_trajectory_solver(p,saveat=10.0, seed=abs(rand(Int)), maxiters=Int(1.0e9))
save_datal("sim_threaded_data_static_pump(S=$(Sre), temp=$(temp), t=800, $np).jld2", sim)