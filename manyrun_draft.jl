### THIS SCRIPT USES THE FUNCTIONALITIES DEFINED IN selforg_core.jl TO RUN THE NUMERICS OVER A RANGE OF PARAMETERS ###

using JLD2
include("src/selforg_core.jl")


p_array = vcat([fill(System_p(0.0, 0.0, i , i , 10.0, -100.0, 100.0, 10.0, 100, (0.0, 600.0), 100)) for i in 25:50]...)


sim_array = vcat([fill(many_trajectory_solver(p_array[i], saveat=10.0, seed=abs(rand(Int)), maxiters=Int(1e9))) for i in 1:length(p_array)]...)

# this object may be used for PLOTTING
sim_extract_array = vcat([extract_solution(sim_array[i]) for i in 1:length(sim_array)]...)


# SAVE THE sim_array
save_datal("sim_data/pump_range/g(25-50)_range_data(t_range=600,N_MC=100).jld", sim_extract_array)
