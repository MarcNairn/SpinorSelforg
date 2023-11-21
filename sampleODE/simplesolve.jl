@everywhere using Distributed
@everywhere using DifferentialEquations
@everywhere using ClusterManagers
@everywhere using JLD2

# Number of trajectories
num_trajectories = 40

# Set up the SlurmManager to add worker processes
@everywhere println("Number of workers before addprocs: ", nworkers()) # for diagnostics
@everywhere addprocs(SlurmManager(num_trajectories))
@everywhere println("Number of workers after addprocs: ", nworkers())

# @everywhere macro ensures that simple_ode! is defined on all workers
@everywhere function simple_ode!(du, u, p, t)
    du[1] = -0.1 * u[1]
end

u0 = [1.0]

tspan = (0.0, 10.0)
prob = ODEProblem(simple_ode!, u0, tspan)

# @everywhere macro ensures that prob_func is defined on all workers
@everywhere function prob_func(prob, i, repeat)
    remake(prob, u0 = rand() * prob.u0)
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)

# Run the ensemble simulation in parallel
@time sim = solve(ensemble_prob, Tsit5(), EnsembleDistributed(), trajectories = num_trajectories)

# Save the results to a data file
output_file = "parallel_sim_example.jld2"
save(output_file, "sim", sim)

# Optionally, you can also save additional information or parameters if needed

# Remove worker processes
@everywhere for i in workers()
    rmprocs(i)
end
