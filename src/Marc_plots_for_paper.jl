#HERE WE LIST ADDITIONAL PLOTS FOR THE SPIN-MOTION ORGANISATION PAPER IN PREPARATION 
#THAT ARE NOT IMMEDIATELY LISTED IN THE LUIGIPLOTTING MODULE 

using Interpolations
using JLD2
using Plots

include("sim_reader.jl")



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
    ax[:fill_between](S, y[1].+y[2], y[1].-y[2], color="C0", alpha=0.2)
    ax[:set_xlabel](L"{S_1}={S_2}")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_adaga_vs_S(solorsim,filename::String)
    fig, ax = plot_adaga_vs_S(solorsim)
    fig[:savefig](filename)
end


function plot_ai_ar_scatter(sim::Array{Sol,1},bins=50)

    fig, ax = subplots()
    sorted_sim = split_sim_from_par(sim) 
    
    N::Int = sorted_sim[1][1].p.N #extract atom number of simulations

    # #end to extract SR phase 
    # #1 to extract normal phase
    for traj in 1:length(sorted_sim[end])
        a_values = sorted_sim[end][traj].u[5N+1, end]
        b_values = sorted_sim[end][traj].u[5N+2, end]
        c_values = sorted_sim[1][traj].u[5N+1, end]
        d_values = sorted_sim[1][traj].u[5N+2, end]
        scatter(a_values, b_values,color= "C0", alpha=0.5, s=25) #a_values along x, b_values along y
        scatter(c_values, d_values,color = "C2", alpha=0.5, s=25)
    end

    xlabel(L"a_r")
    ylabel(L"a_i")

    grid(true)
    xlim(-15, 15)
    ylim(-15, 15)

    # Create 2D histogram
    # plt[:xlim](-25, 25)
    # plt[:ylim](-25, 25)
    # plt[:xlabel](L"$a_r$")
    # plt[:ylabel](L"$a_i$")
    # plt[:hist2d]([sorted_sim[end][traj].u[5N+1, end] for traj in 1:length(sorted_sim[end])],
    #        [sorted_sim[end][traj].u[5N+2, end] for traj in 1:length(sorted_sim[end])],
    #        bins=bins, cmap="Blues", alpha=0.5)
    # plt[:hist2d]([sorted_sim[1][traj].u[5N+1, end] for traj in 1:length(sorted_sim[end])],
    # [sorted_sim[1][traj].u[5N+2, end] for traj in 1:length(sorted_sim[end])],
    # bins=bins, cmap="Reds", alpha=0.5)

    # colorbar()
    legend(["SR", "Normal"])
end


# function plot_ai_ar_scatter_histogram(sorted_sim, a, b, S_i, bins=20)

#     a_values = [sorted_sim[S_i][traj].u[a, end] for traj in 1:length(sorted_sim[S_i])]
#     b_values = [sorted_sim[S_i][traj].u[b, end] for traj in 1:length(sorted_sim[S_i])]

#     scatter(a_values, b_values, alpha=0.5, label="Scatter Plot")

#     xlabel(L"a_r")
#     ylabel(L"a_i")

#     # Add grid
#     grid(true)

#     # Set axis limits to center the plot
#     xlims!(-15, 15)
#     ylims!(-15, 15)

#     # Create 2D histogram
#     histogram2d(a_values, b_values, bins=bins, c=:blues, alpha=0.5)

#     #colorbar()
# end



# function plot_ai_ar_heatmap(sorted_sim, a, b, S_i, bins=20)
#     a_values = [sorted_sim[S_i][traj].u[a, end] for traj in 1:length(sorted_sim[S_i])]
#     b_values = [sorted_sim[S_i][traj].u[b, end] for traj in 1:length(sorted_sim[S_i])]

#     xlabel(L"a_r")
#     ylabel(L"a_i")

#     # Add grid
#     grid(true)

#     # Create 2D histogram heatmap with PyPlot
#     plt[:hist2d](a_values, b_values, bins=bins, cmap="coolwarm", alpha=0.7, extent=[-15, 15, -15, 15])

#     colorbar()

# end


function plot_vs_temp(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    S_list = unique([abs(p.S₁) for p in par_list])
    temp_list = unique([p.temp for p in par_list])
    if length(temp_list)*length(S_list) > length(categories)
        @info("Noplt all data have been plotted.")
        S_list = S_list[1:end-1]
        categories = reshape(categories[1:length(temp_list)*length(S_list)],length(temp_list),length(S_list))
    else
        categories = reshape(categories,length(temp_list),length(S_list))
    end

    S_int = intervallize_array(S_list)
    temp_int = intervallize_array(temp_list)

    y1 = map(x->expect(adaga,x)[1][end],categories)

    y2 = (Float64[],Float64[])
    S = Float64[]
    temp = Float64[]
    for x in categories
        push!(S,abs(x[1].p.S₁))
        push!(temp,abs(x[1].p.temp))

        m,s,q = expect(adaga,x)
        push!(y2[1],m[end])
        push!(y2[2],s[end])
    end

    matplotlib[:rc]("axes", labelpad=1)
    matplotlib[:rc]("image", cmap="cividis")

    fig, ax = subplots(2,2,figsize=[6.2, 4.])

    ax[1,1][:set_ylabel](L"temp")
    ax[1,1][:set_xlabel](L"{S_1}={S_2}")
    cs11 = ax[1,1][:pcolormesh](S_int,temp_int,y1)
    cbar11 = fig[:colorbar](cs11,ax=ax[1,1])
    cbar11[:set_label](L"\langle a^\dag a\rangle")

    ax[1,2][:set_ylabel](L"temp")
    ax[1,2][:set_xlabel](L"{S_1}={S_2}")
    cs12 = ax[1,2][:contourf](S_list,temp_list,y1,101)
    cbar12 = fig[:colorbar](cs12,ax=ax[1,2])
    cbar12[:set_label](L"\langle a^\dag a\rangle")

    ax[2,1][:set_ylabel](L"temp")
    ax[2,1][:set_xlabel](L"{S_1}={S_2}")
    cs21 = ax[2,1][:pcolormesh](S_int,temp_int,y1,norm=matplotlib[:colors][:LogNorm]())
    cbar21 = fig[:colorbar](cs21,ax=ax[2,1])
    cbar21[:set_label](L"\langle a^\dag a\rangle")

    ax[2,2][:set_ylabel](L"temp")
    ax[2,2][:set_xlabel](L"{S_1}={S_2}")
    cs22 = ax[2,2][:contourf](S_list,temp_list,y1,norm=matplotlib[:colors][:LogNorm]())
    cbar22 = fig[:colorbar](cs22,ax=ax[2,2])
    cbar22[:set_label](L"\langle a^\dag a\rangle")

    # ax[3,1][:set_ylabel](L"temp")
    # ax[3,1][:set_xlabel](L"{S_1}={S_2}")
    # cs31 = ax[3,1][:tricontourf](S,De,y2[1],101,norm=matplotlib[:colors][:LogNorm]())
    # cbar31 = fig[:colorbar](cs31,ax=ax[3,1])
    # cbar31[:set_label](L"\langle a^\dag a\rangle")

    # ax[3,2][:set_ylabel](L"temp")
    # ax[3,2][:set_xlabel](L"{S_1}={S_2}")
    # cs32 = ax[3,2][:tricontourf](S,temp,y2[1],101)
    # cbar32 = fig[:colorbar](cs32,ax=ax[3,2])
    # cbar32[:set_label](L"\langle a^\dag a\rangle")


    # ax[3,1][:set_ylabel](L"phase of cavity field $a$")
    # ax[3,1][:errorbar](S,y3[1],yerr=y3[2],color="C2",fmt="o")
    # ax[3,1][:set_yticks](pi/8*collect(0:4))
    # ax[3,1][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[end, 1][:set_xlabel](L"{S_1}={S_2}")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_vs_temp(solorsim,filename::String)
    fig, ax = plot_vs_temp(solorsim)
    fig[:savefig](filename)
end


function plot_interp_threshold(sim::Array{Sol, 1})

    sorted_sims = split_sim_from_par(sim)
    
    S = Float64[]
    y = (Float64[], Float64[])
    
    for sorted_sim in sorted_sims
        push!(S, abs(sorted_sim[1].p.S₁)) #extract coupling strengths from sorted_sims
    
        m, s, q = expect(adaga, sorted_sim) #extract mean, standard dev and 90 quantile
        push!(y[1], m[end])
        push!(y[2], s[end])
    end 
    
    
    #Need to make S into an equally spaced Tuple
    x = S[1]:1.0:S[end]
    
    itp_linear = linear_interpolation(x, y[1])
    itp_cubic = cubic_spline_interpolation(x, y[1])
    # Interpolation functions
    f_linear(x) = itp_linear(x)
    f_quad(x) = itp_quad(x)
    f_cubic(x) = itp_cubic(x)
    # Plots
    width, height = 1500, 800 
    
    x_new = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline
    

    plot(f_linear, x_new, w=5,label="Linear interpolation", color=:"#fc8d62")
    plot!(f_cubic, x_new, linestyle=:dash, w=5, label="Cubic Spline interpolation", color=:"#8da0cb")
    scatter!(S, y[1], markersize=10,label="Data points", color=:white)


    xlabel!(L"S", fontsize=20)
    ylabel!(L"⟨a^\dagger a\rangle", fontsize=20)
    
    plot!(size=(width, height), legend=:topleft, legendfontsize=15, grid=false, xtickfontsize=12, ytickfontsize=12)

    # Calculate second derivative using Central Finite Difference
    delta_x = x_new[2] - x_new[1]
    f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2
    # Plot second derivative
    #plot!(f_cubic_second_derivative, x_new[2:end-1], linestyle=:dot, w=3, label=L"Second Derivative: CFD")


    # Inset plot
    inset_x_range = 35:0.01:37.5
    inset_y_values = f_cubic_second_derivative.(inset_x_range)
    plot!(inset_x_range, inset_y_values, color="#8da0cb", inset=(1, bbox(0.1, 0.5, 0.3, 0.3)), subplot=2, linestyle=:dash, linewidth=5,xlabel=L"S", ylabel=L"\dfrac{d^2}{dS^2}⟨a^\dagger a\rangle", label=nothing, grid=false, xtickfontsize=10, ytickfontsize=10, xlabelfontsize=15, ylabelfontsize=15)
    # Draw a horizontal line at y=0 in the inset plot
    hline!([0], color=:red, alpha=0.5, linestyle=:solid, w=3, subplot=2, label=nothing)

    closest_to_zero_index = argmin(abs.(inset_y_values))
    closest_to_zero_x = inset_x_range[closest_to_zero_index]
    closest_to_zero_y = inset_y_values[closest_to_zero_index]
    scatter!([closest_to_zero_x], [closest_to_zero_y], markersize=7, marker=:circle, markercolor=:white, markerstrokecolor=:black, markerstrokealpha=0.6, subplot=2, label="Crossover")

end

function plot_interp_threshold(solorsim,filename::String)
    fig, ax = plot_interp_threshold(solorsim)
    fig[:savefig](filename)
end