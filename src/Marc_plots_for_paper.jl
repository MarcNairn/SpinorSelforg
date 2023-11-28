#HERE WE LIST ADDITIONAL PLOTS FOR THE SPIN-MOTION ORGANISATION PAPER IN PREPARATION 
#THAT ARE NOT IMMEDIATELY LISTED IN THE LUIGIPLOTTING MODULE 
 
include("selforg_core.jl")



function ai_ar_scatter(sim::Sol)
    println("found $(size(sim)[end]) trajectories")
    println("-------------------------------------")

    categories, par_list = split_sim_from_par(sim, true)

    println(par_list)

    S = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = (Float64[], Float64[])

end