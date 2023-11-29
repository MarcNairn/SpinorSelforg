#HERE WE LIST ADDITIONAL PLOTS FOR THE SPIN-MOTION ORGANISATION PAPER IN PREPARATION 
#THAT ARE NOT IMMEDIATELY LISTED IN THE LUIGIPLOTTING MODULE 
 
include("selforg_core.jl")


#=
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
=#


function plot_adaga_vs_S(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim, true)

    println(par_list)

    S = Float64[]
    y = (Float64[], Float64[])
    for x in categories
        push!(S, abs(x[1].p.S₁))

        m, s, q = expect(adaga, x)
        push!(y[1], m[end])
        push!(y[2], s[end])
    end

    matplotlib[:rc]("axes", labelpad=1)

    fig, ax = subplots(figsize=[6.2, 4.])

    ax[:set_ylabel](L"\langle a^\dag a\rangle")
    #ax[:errorbar](S, y[1], yerr=y[2], color="C0", fmt="o")
    ax[:plot](S, y[1], color = "C0")
    ax[:fill_between](S, y[1].+y[2], y[1].-y[2], color="C0", alpha=0.5)
    ax[:set_xlabel](L"{S_1}={S_2}")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_adaga_vs_S(solorsim,filename::String)
    fig, ax = plot_adaga_vs_S(solorsim)
    fig[:savefig](filename)
end


function plot_ai_ar_scatter(sim::Array{Sol,1}, selected_S::Int)
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    #categories, par_list = split_sim_from_par(sim, true)

    selected_set, par_list = split_sim_from_par(sim, false)[selected_S]

    #println(par_list)

    S = Float64[]
    y1 = Float64[]#, Float64[])
    y2 = Float64[]#, Float64[])
    
    for i in 1:length(selected_set)
        m_ai[i] = expect(absai, selected_set[i])
        push!(y1, m_ai[i])
        #push!(y1[2], s_ai[end])

        m_ar[i]= expect(absar, selected_set[i])
        push!(y2, m_ar[i])
        #push!(y2[2], s_ar[end])
    end

   # for x in categories
    #    current_S = abs(x[1].p.S₁)
        
        # Check if the current trajectory has the desired S value
     #   if current_S == selected_S
      #      push!(S, current_S)


       # end
    #end

    matplotlib[:rc]("axes", labelpad=1)
    fig, ax = subplots(figsize=[6.2, 4.])

    # Scatter plot
    ax[:scatter](y1[1], y2[1], label="m_ai vs m_ar", color="C0")

    # Add labels and title
    ax[:set_xlabel]("m_ai")
    ax[:set_ylabel]("m_ar")
    ax[:legend]()

    return fig, ax
end

function plot_ai_ar_scatter(solorsim,selected_S::Int,filename::String)
    fig, ax = plot_ai_ar_scatter(solorsim, selected_S::Int)
    fig[:savefig](filename)
end