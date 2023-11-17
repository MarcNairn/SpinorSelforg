using DifferentialEquations

# Linear ODE start at 0.5 and solves from t=0.0 to t=1.0
prob = ODEProblem((u, p, t) -> 1.01u, 0.5, (0.0, 1.0))

function prob_func(prob, i, repeat)
    remake(prob, u0 = rand() * prob.u0)
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
elt = @elapsed sim = solve(ensemble_prob, Tsit5(), EnsembleDistributed(), trajectories = 1000)
println("done in $elt seconds")

for i in workers()
	rmprocs(i)
end
