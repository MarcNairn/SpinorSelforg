using DifferentialEquations
using JLD2

 function simple_ode!(du, u, p, t)
    du[1] = -0.1 * u[1]
end

index = parse(Int, ARGS[1])

u0 = [1.0]

tspan = (0.0, 10.0)
prob = ODEProblem(simple_ode!, u0, tspan)

 function prob_func(prob, i, repeat)
    remake(prob, u0 = rand() * prob.u0)
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)

println("Start of computation")
@time sim = solve(ensemble_prob, Tsit5(), trajectories = 1)
println("End of computation")


# Save the results to a data file
output_file = "distributed/cluster_data/parallel_sim_example_$(index).jld2"
save(output_file, "sim", sim)

println("REACHED THE END")