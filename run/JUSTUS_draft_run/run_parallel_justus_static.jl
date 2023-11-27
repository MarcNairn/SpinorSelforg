using ClusterManagers
using Random

np = parse(Int,ARGS[1])

include("../load.jl")
include("../parameters_static.jl")

Sre = Real(p.S₁)

sim = many_trajectory_solver(p,saveat=10.0, seed=abs(rand(Int)), maxiters=Int(1.0e9))
save_datal("~/SpinorSelforg/cluster_data/static_pump$(Sre)_$(p.tspan)_nt=$(np).jld2", sim)
