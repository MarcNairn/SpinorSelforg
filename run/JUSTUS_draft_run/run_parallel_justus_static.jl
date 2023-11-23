using Pkg
Pkg.activate("..")

using Distributed
@everywhere using ClusterManagers

np = parse(Int,ARGS[1])


@everywhere push!(LOAD_PATH,"../src")
@everywhere include("../load.jl")

using Random

@everywhere include("../parameters_static.jl")

sim = many_trajectory_solver(p,saveat=10.0, seed=abs(rand(Int)), maxiters=Int(1.0e9))
save_datal("static_pump$(p.S₁)_$(p.S₂)_$(p.tspan).jld2", sim)


for i in workers()
	rmprocs(i)
end
