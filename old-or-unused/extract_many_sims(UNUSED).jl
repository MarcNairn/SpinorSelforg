#### FUNCTION TO RECOVER LAST STEPS 
function get_last_steps(sol::Sol)
    Sol(sol.u[:,end-1:end],sol.p,sol.t[end-1:end],sol.alg)
end

function get_last_steps(sim::Array{Sol,1})
    [get_last_steps(sol) for sol in sim]
end

function count_indices(data) # to know number of independent indeces in "data"
    if isa(data, AbstractArray)
        return 1 + count_indices(data[1])
            else
        return 0
    end
end


##### MANY_SIM STRUCTURE
##### stored as many_sims[g_index][N_MC_index][t_index, u_index]

flatten_sim = vcat(sim_array...)

sim_extract_array = vcat([extract_solution(flatten_sim[i]) for i in 1:length(flatten_sim)]...)