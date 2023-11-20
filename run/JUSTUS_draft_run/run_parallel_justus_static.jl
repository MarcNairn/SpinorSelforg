using Pkg
Pkg.activate("..")

using Distributed
using ClusterManagers

np = parse(Int,ARGS[1])

addprocs(SlurmManager(np))

@everywhere push!(LOAD_PATH,"../src")
@everywhere include("../load.jl")

using Random

# Random.seed!(232344)

@everywhere include("../parameters_static.jl")

sim = many_trajectory_solver(p,saveat=10.0, seed=abs(rand(Int), maxiters=Int(1.0e9)))
save_datal("data_vs_pump_S4_part" * string(i) * ".jld2",sim)

for i in workers()
	rmprocs(i)
end
