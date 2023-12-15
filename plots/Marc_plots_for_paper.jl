#HERE WE LIST ADDITIONAL PLOTS FOR THE SPIN-MOTION ORGANISATION PAPER IN PREPARATION 
#THAT ARE NOT IMMEDIATELY LISTED IN THE LUIGIPLOTTING MODULE 

using Interpolations
using JLD2
using DiffEqBase
using LaTeXStrings
using LinearAlgebra
using LsqFit
using Measurements
using Random: MersenneTwister
using PyPlot
using StatsBase

include("../src/sim_reader.jl")



function plot_adaga_vs_S(sim::Array{Sol, 1})
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
    
    # Plot the line
    ax[:plot](S, y[1], color="C0")

    ax[:fill_between](S, y[1] + y[2], y[1] - y[2], color="C0", alpha=0.2)

    ax[:set_xlabel](L"{|S_1|}={|S_2|}")

    fig[:tight_layout](h_pad=0., w_pad=-0.)
    return fig, ax
end

function plot_adaga_vs_S(solorsim,filename::String)
    fig, ax = plot_adaga_vs_S(solorsim)
    fig[:savefig](filename)
end

# function plot_ai_ar_scatter(sim::Array{Sol, 1}, bins=50)
#     sorted_sim = split_sim_from_par(sim)

#     N::Int = sorted_sim[1][1].p.N  # Extract atom number of simulations

#     # Extracting values for scatter plot
#     scatter_values_sr = hcat(
#         [sorted_sim[end][traj].u[5N + 1, end] for traj in 1:length(sorted_sim[end])],
#         [sorted_sim[end][traj].u[5N + 2, end] for traj in 1:length(sorted_sim[end])]
#     )

#     scatter_values_normal = hcat(
#         [sorted_sim[1][traj].u[5N + 1, end] for traj in 1:length(sorted_sim[end])],
#         [sorted_sim[1][traj].u[5N + 2, end] for traj in 1:length(sorted_sim[end])]
#     )

#     mean_x_normal = mean(scatter_values_normal[:, 1])
#     mean_y_normal = mean(scatter_values_normal[:, 2])

#     mean_x_sr = mean(scatter_values_sr[:, 1])
#     mean_y_sr = mean(scatter_values_sr[:, 2])

#     Plots.plot(scatter_values_sr[:, 1], scatter_values_sr[:, 2], seriestype=:scatter, color=:blue, alpha=0.5, label="SR", markersize=5)
#     Plots.plot!(scatter_values_normal[:, 1], scatter_values_normal[:, 2], seriestype=:scatter, color=:red, alpha=0.5, label="Normal", markersize=5)

#     xlabel!(L"a_r", fontsize=15)
#     ylabel!(L"a_i", fontsize=15)

#     plot!(legend=true)
#     plot!(grid=true)

#     xlims!((mean_x_sr - 2*std(scatter_values_sr[:, 1]), mean_x_sr + 2*std(scatter_values_sr[:, 1])))
#     ylims!((mean_y_sr - 2*std(scatter_values_sr[:, 2]), mean_y_sr + 2*std(scatter_values_sr[:, 2])))

# end


function plot_ordering_vs_S(sims::Array{Sol,1}...)

    colorlist = ["C1","C2","C3","C4","C5","C6"]
    linelist = ["-","--",":","-.","-."]

    # fig, ax = subplots(3,1,figsize=[3.4, 5.3],sharex=true)
    fig, ax = subplots(4,1,figsize=[3.4, 7.3],sharex=true)

    ax[1].set_ylabel(L"cavity population $\langle a^\dag a \rangle$")
    ax[2].set_ylabel(L"order parameter $\vert\Phi\vert$")
    ax[3].set_ylabel(L"bunching parameter $\mathcal{B}$")
    ax[4].set_ylabel(L"final kinetic energy $E_\mathrm{kin}/\hbar\kappa$")
    ax[end].set_xlabel(L"pump strength $\sqrt{N}S/\kappa$")

    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim,true)

        S = Float64[]
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        y3 = Float64[]
        y4 = Array{Float64,1}[]
        y5 = Float64[]
        y6 = Array{Float64,1}[]
        y7 = Float64[]
        y8 = Array{Float64,1}[]
        for x in categories
            push!(S,abs(x[1].p.S₁)*sqrt(x[1].p.N)/x[1].p.κ)

            m,s,q = expect(adaga,x)
            push!(y1,m[end])
            push!(y2,q[end])

            m,s,q = expect(absX,x)
            push!(y3,m[end])
            push!(y4,q[end])

            m,s,q = expect(Cos2,x)
            push!(y5,m[end])
            push!(y6,q[end])

            m,s,q = expect(Ekin,x)./x[1].p.κ
            push!(y7,m[end])
            push!(y8,q[end])
        end

        A = sortslices(hcat(S,y1,vcat(y2'...),y3,vcat(y4'...),y5,vcat(y6'...),y7,vcat(y8'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]
        y3 = A[:,5]
        y4 = A[:,6:7]
        y5 = A[:,8]
        y6 = A[:,9:10]
        y7 = A[:,11]
        y8 = A[:,12:13]

        # matplotlib[:rc]("axes", labelpad=1)

        if par_list[1].temp == 0.0
            label = "\$temp=0\$"
        elseif par_list[1].temp < par_list[1].κ
            label = "\$temp=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].temp))*"\$"
        else
            label = "\$temp="*string(trunc(Int,par_list[1].temp/par_list[1].κ))*"\\kappa\$"
        end

        ax[1].plot(S,y1,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[1].fill_between(S,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[2].plot(S,y3,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[2].fill_between(S,y4[:,1],y4[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[3].plot(S,y5,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[3].fill_between(S,y6[:,1],y6[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[4].plot(S,y7,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[4].fill_between(S,y8[:,1],y8[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    ax[1].legend(handlelength=2.5,loc="upper left",bbox_to_anchor=(0.2, 1.08),framealpha=1)
    fig.tight_layout(h_pad=0.)

    for i in 1:4
        letter = Char(Int('a')+i-1)
        ax[i].text(0.05,0.87,"("*letter*")",transform=ax[i].transAxes)
    end


    return fig, ax
end

# function plot_ordering_vs_S(filename::String,sims::Array{Sol,1}...)
#     fig, ax = plot_ordering_vs_S(sims...)
#     fig.savefig(filename,dpi=1200)
# end


# function plot_vs_temp(sim::Array{Sol,1})
#     println("found $(size(sim)[end]) trajectories")
#     println("-----------------------------------------")

#     categories, par_list = split_sim_from_par(sim,true)

#     S_list = unique([abs(p.S₁) for p in par_list])
#     temp_list = unique([p.temp for p in par_list])
#     if length(temp_list)*length(S_list) > length(categories)
#         @info("Noplt all data have been plotted.")
#         S_list = S_list[1:end-1]
#         categories = reshape(categories[1:length(temp_list)*length(S_list)],length(temp_list),length(S_list))
#     else
#         categories = reshape(categories,length(temp_list),length(S_list))
#     end

#     S_int = intervallize_array(S_list)
#     temp_int = intervallize_array(temp_list)

#     y1 = map(x->expect(adaga,x)[1][end],categories)

#     y2 = (Float64[],Float64[])
#     S = Float64[]
#     temp = Float64[]
#     for x in categories
#         push!(S,abs(x[1].p.S₁))
#         push!(temp,abs(x[1].p.temp))

#         m,s,q = expect(adaga,x)
#         push!(y2[1],m[end])
#         push!(y2[2],s[end])
#     end

#     matplotlib[:rc]("axes", labelpad=1)
#     matplotlib[:rc]("image", cmap="cividis")

#     fig, ax = subplots(2,2,figsize=[6.2, 4.])

#     ax[1,1][:set_ylabel](L"temp")
#     ax[1,1][:set_xlabel](L"{S_1}={S_2}")
#     cs11 = ax[1,1][:pcolormesh](S_int,temp_int,y1)
#     cbar11 = fig[:colorbar](cs11,ax=ax[1,1])
#     cbar11[:set_label](L"\langle a^\dag a\rangle")

#     ax[1,2][:set_ylabel](L"temp")
#     ax[1,2][:set_xlabel](L"{S_1}={S_2}")
#     cs12 = ax[1,2][:contourf](S_list,temp_list,y1,101)
#     cbar12 = fig[:colorbar](cs12,ax=ax[1,2])
#     cbar12[:set_label](L"\langle a^\dag a\rangle")

#     ax[2,1][:set_ylabel](L"temp")
#     ax[2,1][:set_xlabel](L"{S_1}={S_2}")
#     cs21 = ax[2,1][:pcolormesh](S_int,temp_int,y1,norm=matplotlib[:colors][:LogNorm]())
#     cbar21 = fig[:colorbar](cs21,ax=ax[2,1])
#     cbar21[:set_label](L"\langle a^\dag a\rangle")

#     ax[2,2][:set_ylabel](L"temp")
#     ax[2,2][:set_xlabel](L"{S_1}={S_2}")
#     cs22 = ax[2,2][:contourf](S_list,temp_list,y1,norm=matplotlib[:colors][:LogNorm]())
#     cbar22 = fig[:colorbar](cs22,ax=ax[2,2])
#     cbar22[:set_label](L"\langle a^\dag a\rangle")

#     # ax[3,1][:set_ylabel](L"temp")
#     # ax[3,1][:set_xlabel](L"{S_1}={S_2}")
#     # cs31 = ax[3,1][:tricontourf](S,De,y2[1],101,norm=matplotlib[:colors][:LogNorm]())
#     # cbar31 = fig[:colorbar](cs31,ax=ax[3,1])
#     # cbar31[:set_label](L"\langle a^\dag a\rangle")

#     # ax[3,2][:set_ylabel](L"temp")
#     # ax[3,2][:set_xlabel](L"{S_1}={S_2}")
#     # cs32 = ax[3,2][:tricontourf](S,temp,y2[1],101)
#     # cbar32 = fig[:colorbar](cs32,ax=ax[3,2])
#     # cbar32[:set_label](L"\langle a^\dag a\rangle")


#     # ax[3,1][:set_ylabel](L"phase of cavity field $a$")
#     # ax[3,1][:errorbar](S,y3[1],yerr=y3[2],color="C2",fmt="o")
#     # ax[3,1][:set_yticks](pi/8*collect(0:4))
#     # ax[3,1][:set_yticklabels]([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

#     ax[end, 1][:set_xlabel](L"{S_1}={S_2}")

#     fig[:tight_layout](h_pad=0., w_pad=-0.)
#     return fig, ax
# end

# function plot_vs_temp(solorsim,filename::String)
#     fig, ax = plot_vs_temp(solorsim)
#     fig[:savefig](filename)
# end

# function plot_interp_threshold(sims::Vector{Vector{Sol}})
# #initiate empty plot
# plot()
#     for sim in sims
#         sorted_sims = split_sim_from_par(sim)

#         S = Float64[]
#         y = (Float64[], Float64[])
        
#         for sorted_sim in sorted_sims
#             push!(S, abs(sorted_sim[1].p.S₁)) #extract coupling strengths from sorted_sims
        
#             m, s, q = expect(adaga, sorted_sim) #extract mean, standard dev and 90 quantile
#             push!(y[1], m[end])
#             push!(y[2], s[end])
#         end 
           
#         #Need to make S into an equally spaced Tuple
#         x = S[1]:1.0:S[end]
        
#         itp_linear = linear_interpolation(x, y[1])
#         itp_cubic = cubic_spline_interpolation(x, y[1])
#         # Interpolation functions
#         f_linear(x) = itp_linear(x)
#         f_cubic(x) = itp_cubic(x)
#         # Plots
#         width, height = 1500, 800 
        
#         x_new = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline
        

#         plot(f_linear, x_new, w=5,label="Linear interpolation", color=:"#fc8d62")
#         plot!(f_cubic, x_new, linestyle=:dash, w=5, label="Cubic Spline interpolation", color=:"#8da0cb")
#         scatter!(S, y[1], markersize=10,label="Data points", color=:white)


#         xlabel!(L"S", fontsize=20)
#         ylabel!(L"⟨a^\dagger a\rangle", fontsize=20)
        
#         plot!(size=(width, height), legend=:topleft, legendfontsize=15, grid=false, xtickfontsize=12, ytickfontsize=12)

#         # Calculate second derivative using Central Finite Difference
#         delta_x = x_new[2] - x_new[1]
#         f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2



#         # Inset plot
#         inset_x_range = 35:0.01:37.5
#         inset_y_values = f_cubic_second_derivative.(inset_x_range)
#         plot!(inset_x_range, inset_y_values, color="#8da0cb", inset=(1, bbox(0.1, 0.5, 0.3, 0.3)), subplot=2, linestyle=:dash, linewidth=5,xlabel=L"S", ylabel=L"\dfrac{d^2}{dS^2}⟨a^\dagger a\rangle", label=nothing, grid=false, xtickfontsize=10, ytickfontsize=10, xlabelfontsize=15, ylabelfontsize=15)
#         # Draw a horizontal line at y=0 in the inset plot
#         hline!([0], color=:red, alpha=0.5, linestyle=:solid, w=3, subplot=2, label=nothing)

#         closest_to_zero_index = argmin(abs.(inset_y_values))
#         closest_to_zero_x = inset_x_range[closest_to_zero_index]
#         closest_to_zero_y = inset_y_values[closest_to_zero_index]
#         scatter!([closest_to_zero_x], [closest_to_zero_y], markersize=7, marker=:circle, markercolor=:white, markerstrokecolor=:black, markerstrokealpha=0.6, subplot=2, label="Crossover")
#     end
# end

# CHANGE TO AVOID USING PLOTS; CHANGE TO USE PYPLOT INSTEAD

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
    f_cubic(x) = itp_cubic(x)
    
    x_smooth = S[1]:0.1:S[end] # smoother interval, necessary for cubic spline

    fig, ax1 = plt.subplots()

    ax1.plot(x, y[1], "o", label="data", color="black", alpha=0.8)
    ax1.plot(x_smooth, f_linear.(x_smooth), linestyle=":", label="linear", color="black", alpha=0.8, markersize=0.5)
    ax1.plot(x_smooth, f_cubic.(x_smooth), linestyle="--", label="cubic", color="gray", alpha=0.5, markersize=0.5)

    # Customize plot
    ax1.set_xlabel(L"|S₁|=|S₂|")
    ax1.set_ylabel(L"\langle a^\dagger a\rangle")
    ax1.legend(loc="bottom right")
    #Inset 
    # Calculate second derivative using Central Finite Difference
    delta_x = x_smooth[2] - x_smooth[1]
    f_cubic_second_derivative(x) = (f_cubic(x - delta_x) - 2 * f_cubic(x) + f_cubic(x + delta_x)) / delta_x^2

    inset_x_range = 35:0.01:37.5
    inset_y_values = f_cubic_second_derivative.(inset_x_range)
    left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(inset_x_range, inset_y_values, color="black", "--", alpha=0.5)#, inset=(1, bbox(0.1, 0.5, 0.3, 0.3)), subplot=2, linestyle=:dash, linewidth=5,, label=nothing, grid=false, xtickfontsize=10, ytickfontsize=10, xlabelfontsize=15, ylabelfontsize=15)
    ax2.set_xlabel(L"S")
    ax2.set_ylabel(L"\partial_S^2⟨a^\dagger a\rangle")
    # Draw a horizontal line at y=0 in the inset plot
    ax2.axhline(y=0, color=:black, alpha=0.5, label=nothing)

    closest_to_zero_index = argmin(abs.(inset_y_values))
    closest_to_zero_x = inset_x_range[closest_to_zero_index]
    closest_to_zero_y = inset_y_values[closest_to_zero_index]
    ax2.scatter([closest_to_zero_x], [closest_to_zero_y], color=:black, alpha=0.8)
end

function plot_interp_threshold(solorsim,filename::String)
    fig, ax = plot_interp_threshold(solorsim)
    fig[:savefig](filename)
end






####################################################################################





####################################################################################






function plot_spinspositionhisto(sim::Array{Sol,1})
    # matplotlib[:rc]("axes", labelpad=1)
    matplotlib.rc("image", cmap="inferno")

    sim_list = Array{Sol,1}[]

    for i in split_sim(sim)
        if length(i)>0
            push!(sim_list,i)
        end
    end

    nrows = length(sim_list)

    # fig, ax = subplots(nrows,2,figsize=[4.65, 2.3*nrows],sharex=true,sharey=false)
    fig, ax = subplots(nrows,2,figsize=[3.7, 1.5*nrows],sharex=true,sharey=false)

    for (i,a) in enumerate(sim_list)
        u1 = join_trajectories(a)
        N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in a) # for backwards compatibility
        x = mod2pi.(u1[1:N])/(2pi)

        nbins = trunc(Int, ^(N,1//3))
        println("number of bins: "*string(nbins))

        ax1 = ax[i,1]
        ax1.set_ylabel(L"spin $\langle\sigma_x\rangle$")
        ax1.set_yticks(collect(-1:0.5:1))
        ax1.set_yticklabels([L"-1",L"-\frac{1}{2}",L"0",L"\frac{1}{2}",L"1"])
        ax1.hist2d(x,u1[2N+1:3N],bins=nbins,density=true)#,norm=matplotlib.colors.LogNorm())

        ax2 = ax[i,2]
        ax2.set_ylabel(L"spin $\langle\sigma_y\rangle$")
        ax2.set_yticks(collect(-1:0.5:1))
        ax2.set_yticklabels([L"-1",L"-\frac{1}{2}",L"0",L"\frac{1}{2}",L"1"])
        cs = ax2.hist2d(x,u1[3N+1:4N], bins=nbins,density=true)

        ax2.set_xticks([cs[2][1],0.5,cs[2][end]])
        ax2.set_xticklabels(["0","0.5","1"])

    end

    fig.tight_layout()

    fig.subplots_adjust(top=0.85,bottom=0.13)
    cbar_ax = fig.add_axes([0.2, 0.92, 0.7, 0.02])
    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(),cax=cbar_ax,orientation="horizontal")
    # cbar_ax.xaxis.set_ticks_position("top")
    cbar_ax.xaxis.set_label_position("top")
    cbar_ax.xaxis.set_tick_params(pad=1)
    cbar.set_label("density (a.u.)",labelpad=3.)

    fig.text(0.5, 0.03, L"atom position mod $\lambda_\mathrm{c}$ (units of $\lambda_\mathrm{c}$)", ha="center")

    return fig, ax

end

function plot_spinspositionhisto(filename::String,sim::Array{Sol,1})
    fig, ax = plot_spinspositionhisto(sim)
    fig.savefig(filename)
end