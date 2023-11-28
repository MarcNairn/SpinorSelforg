#THIS SCRIPT EXTRACTS ALL SIMULATIONS STORED AS .JLD FILES INTO A SIMULATION ARRAY FROM WHICH WE MAY EXTRACT THE OBSERVABLE sim_data
using JLD2

include("src/custom_functions.jl")
include("src/selforg_core.jl")

#re-import in case struct System_p undefined when loading data
struct System_p
    U₁::Float64 
    U₂::Float64
    S₁::Complex{Float64}
    S₂::Complex{Float64}
    Δₑ::Float64
    Δc::Float64
    κ::Float64
    temp::Float64
    N::Int
    tspan::Tuple{Float64,Float64}
    N_MC::Int
end

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

sim_array = load_datalll("cluster_files/")

