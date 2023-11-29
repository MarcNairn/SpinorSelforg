#THIS SCRIPT EXTRACTS ALL SIMULATIONS STORED AS .JLD FILES INTO A SIMULATION ARRAY FROM WHICH WE MAY EXTRACT THE OBSERVABLE sim_data
#ALSO USED TO CLEAN UP ALL THE USED JLD2 FILES AND STORE THEM IN N_TRAJECTORY SIZED BATCHES
using JLD2

include("selforg_core.jl")



#TO EXTRACT ALL SIMULATION DATA INTO A COMMON ARRAY
function load_datalll(directory::AbstractString)
    sim = []

    for file in readdir(directory)
        if endswith(file, ".jld2")
            filetemp = joinpath(directory, file)
            println("loading $filetemp...")
            push!(sim, load_datal(filetemp))
        end
    end

    merge_sim(sim...)
end

sim_array = load_datalll("../cluster_files/") #of type Vector{Sol} (or Array{Sol,1})

sorted_sims = split_sim_from_par(sim_array)

#RUN THROUGH ALL SIM DATA AND CATEGORISE IT IN BATCHES OF THE SAME PARAMETERS

for batch in sorted_sims 
    g = Int(batch[1].p.S₁)
    trajs = length(batch)
    save_datal("../sorted_sims/S=$(g), Ntraj=$(trajs)-sim_data.jld2", batch)
end