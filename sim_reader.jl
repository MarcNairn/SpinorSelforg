#THIS SCRIPT EXTRACTS ALL SIMULATIONS STORED AS .JLD FILES INTO A SIMULATION ARRAY FROM WHICH WE MAY EXTRACT THE OBSERVABLE sim_data
using JLD2

#include("src/custom_functions.jl")
include("src/selforg_core.jl")



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

sim_array = load_datalll("cluster_files/") #of type Vector{Sol} (or Array{Sol,1})

