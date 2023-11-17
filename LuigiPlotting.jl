module LPlots

#using System1

using DiffEqBase
using LaTeXStrings
using LinearAlgebra
using LsqFit
using Measurements
using Random: MersenneTwister
using PyPlot
using StatsBase


struct Sol
    u::Array{Float64,2}
    p
    t::Array{Float64,1}
    alg::String
end

Sol(u,p,t) = Sol(u,p,t,"")

function extract_solution(sol::RODESolution)
    [Sol(hcat(sol.u...),sol.prob.p,sol.t,repr(sol.alg))]
end

function extract_solution(sim::EnsembleSolution)
    trajectories = Array{Sol,1}()
    for sol in sim
        append!(trajectories,extract_solution(sol))
    end
    trajectories
end

function extract_solution(sim::Array{Sol,1})
    sim
end


function intervallize_array(arr::Array{Float64,1})
    y = zeros(length(arr)+1)
    y[2:end] .= arr
    y[1:end-1] .+= arr
    y = y/2
    y[1] = 2arr[1]-y[2]
    y[end] = 2arr[end]-y[end-1]
    return y
end

model_exp(t,p) = p[1]*exp.(-p[2]*(t.-p[3])) .+ p[4]

function plot_Ekin(sim::Array{Sol,1})

    y,y_std,y_q90 = expect(Ekin,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    # p0 = [900.,1.,1.,100.]
    # fit = curve_fit(model_exp, tlist, y, p0)

    matplotlib[:rc]("axes", labelpad=1.5)

    # matplotlib[:pyplot][:xkcd]()

    fig, ax = subplots(1,1,figsize=[2.325, 2.2])

    color="C3"

    ax[:set_ylabel](L"kinetic energy (units of $\hbar\omega_\mathrm{r}$)")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[:set_xscale]("log")
    ax[:plot](tlist,y,color=color,label="mean kinetic\n energy")
    ax[:fill_between](tlist,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2,linewidth=0.1)
    # ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)
    ax[:plot]((tlist[1],tlist[end]),(100.,100.),ls=":",label=L"$\hbar\kappa$")
    # ax[:plot](tlist,model_exp(tlist,fit.param),ls="--",label="fit")
    ax.legend(loc="upper right", bbox_to_anchor=(1.1, 1.),framealpha=1)
    fig[:tight_layout]()

    return fig, ax
end

function plot_Ekin(sim::Array{Sol,1},filename::String)
    fig, ax = plot_Ekin(sim)
    fig.savefig(filename)
end

function plot_Ekin_momdistro(sim::Array{Sol,1})

    y,y_std,y_q90 = expect(Ekin,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    # p0 = [900.,1.,1.,100.]
    # fit = curve_fit(model_exp, tlist, y, p0)


    u0 = join_trajectories(sim,1)
    u1 = join_trajectories(sim,size(sim[1].t)[1])

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))


    minp0 = minimum(u0[N+1:2N])
    maxp0 = maximum(u0[N+1:2N])
    rangep0 = range(minp0,length=trunc(Int,sqrt(N)),stop=maxp0)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep0)
    pdfp0 = normalize(pdfp0)

    minp1 = minimum(u1[N+1:2N])
    maxp1 = maximum(u1[N+1:2N])
    rangep1 = range(minp1,length=trunc(Int,sqrt(N)),stop=maxp1)
    pdfp1 = fit(Histogram,u1[N+1:2N],rangep0)
    pdfp1 = normalize(pdfp1)


    # matplotlib[:rc]("axes", labelpad=1.5)


    fig, ax = subplots(1,2,figsize=[4.65, 2.2])

    color="C3"

    ax[1].set_ylabel(L"kinetic energy (units of $\hbar\omega_\mathrm{r}$)")
    ax[1].set_xlabel(L"time (units of $\omega_\mathrm{r}^{-1}$)")
    ax[1].plot(tlist,y,color=color,label="mean kinetic\n energy")
    ax[1].fill_between(tlist,y_q90[1,:],y_q90[2,:],color=color,alpha=0.2,linewidth=0.1)
    # ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)
    ax[1].plot((tlist[1],tlist[end]),(100.,100.),ls=":",label=L"$\hbar\kappa$")
    # ax[:plot](tlist,model_exp(tlist,fit.param),ls="--",label="fit")
    ax[1].legend(loc="upper right", bbox_to_anchor=(1.05, 1.),framealpha=1)


    ax[2].set_ylabel("momentum distribution")
    ax[2].set_xlabel(L"atom momentum $p$ (units of $\hbar k_\mathrm{c}$)")
    ax[2].step(pdfp0.edges[1][1:end-1],pdfp0.weights,label="initial",where="post")
    ax[2].step(pdfp1.edges[1][1:end-1],pdfp1.weights,linestyle="--",label="final",where="post")

    ax[2].legend(loc="upper right", bbox_to_anchor=(1.1, 1.),framealpha=1)

    ax[1].text(0.1,0.87,"(a)",transform=ax[1].transAxes)
    ax[2].text(0.1,0.87,"(b)",transform=ax[2].transAxes)

    fig[:tight_layout]()

    return fig, ax
end

function plot_Ekin_momdistro(sim::Array{Sol,1},filename::String)
    fig, ax = plot_Ekin_momdistro(sim)
    fig.savefig(filename)
end

function plot_Ekin_momdistro_def(sim::Array{Sol,1})

    y,y_std,y_q90 = expect(Ekin,sim)
    y_q90 = hcat(y_q90...)
    tlist = sim[1].t

    # p0 = [900.,1.,1.,100.]
    # fit = curve_fit(model_exp, tlist, y, p0)


    u0 = join_trajectories(sim,1)
    u1 = join_trajectories(sim,size(sim[1].t)[1])

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))


    minp0 = minimum(u0[N+1:2N])
    maxp0 = maximum(u0[N+1:2N])
    rangep0 = range(minp0,length=trunc(Int,sqrt(N)),stop=maxp0)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep0)
    pdfp0 = normalize(pdfp0)

    minp1 = minimum(u1[N+1:2N])
    maxp1 = maximum(u1[N+1:2N])
    rangep1 = range(minp1,length=trunc(Int,sqrt(N)),stop=maxp1)
    pdfp1 = fit(Histogram,u1[N+1:2N],rangep0)
    pdfp1 = normalize(pdfp1)


    # matplotlib[:rc]("axes", labelpad=1.5)


    fig, ax = subplots(1,2,figsize=[4.4, 2.1])

    color="C3"

    ax[1].set_ylabel(L"$E_\mathrm{kin}/\hbar\kappa$")
    ax[1].set_xlabel(L"$t/\omega_\mathrm{r}^{-1}$")
    ax[1].plot(tlist,y./sim[1].p.κ,color=color)#,label=L"E_\mathrm{kin}")
    ax[1].fill_between(tlist,y_q90[1,:]./sim[1].p.κ,y_q90[2,:]./sim[1].p.κ,color=color,alpha=0.2,linewidth=0.1)
    # ax[:fill_between](tlist.+1,y.+y_std,y.-y_std,color=color,alpha=0.5)
    ax[1].plot((tlist[1],tlist[end]),(1.,1.),ls=":",label=L"$\hbar\kappa$")
    # ax[:plot](tlist,model_exp(tlist,fit.param),ls="--",label="fit")
    ax[1].legend(loc="upper right")#, bbox_to_anchor=(1.05, 1.),framealpha=1)


    ax[2].set_ylabel("momentum distribution")
    ax[2].set_xlabel(L"$p/\hbar k_\mathrm{c}$")
    ax[2].step(pdfp0.edges[1][1:end-1],pdfp0.weights,label="initial",where="post")
    ax[2].step(pdfp1.edges[1][1:end-1],pdfp1.weights,linestyle="--",label="final",where="post")

    ax[2].legend(loc="upper right", bbox_to_anchor=(1.1, 1.),framealpha=1)

    ax[1].text(0.1,0.87,"(a)",transform=ax[1].transAxes)
    ax[2].text(0.1,0.87,"(b)",transform=ax[2].transAxes)

    fig[:tight_layout](w_pad=0.13)

    return fig, ax
end

function plot_Ekin_momdistro_def(sim::Array{Sol,1},filename::String)
    fig, ax = plot_Ekin_momdistro_def(sim)
    fig.savefig(filename)
end

function plot_cavity(sim::Array{Sol,1})

    ady,ady_std,ady_q90 = expect(adiabaticadaga,sim)
    ady_q90 = hcat(ady_q90...)
    tlist = sim[1].t

    y,y_std,y_q90 = expect(adaga,sim)
    y_q90 = hcat(y_q90...)


    matplotlib[:rc]("axes", labelpad=1.5)

    # matplotlib[:pyplot][:xkcd]()

    fig, ax = subplots(1,1,figsize=[3.4, 2.3])

    color="C3"

    ax[:set_ylabel](L"cavity population (units of $\hbar\omega_\mathrm{r}$)")
    ax[:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    ax[:plot](tlist,y,color="C1",label="cavity population")
    ax[:fill_between](tlist,y_q90[1,:],y_q90[2,:],color="C1",alpha=0.2,linewidth=0.1)
    ax[:plot](tlist,ady,color="C2",ls="--",label="adiabatic cavity population")
    # ax[:fill_between](tlist,ady_q90[1,:],ady_q90[2,:],color="C2",alpha=0.2,linewidth=0.1)
    ax[:legend](handlelength=2.5)
    fig[:tight_layout]()

    return fig, ax
end

function plot_cavity(sim::Array{Sol,1},filename::String)
    fig, ax = plot_cavity(sim)
    fig[:savefig](filename)
end

function plot_momentum_distro(sim::Array{Sol,1})

    u0 = join_trajectories(sim,1)
    u1 = join_trajectories(sim,size(sim[1].t)[1])

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))


    minp0 = minimum(u0[N+1:2N])
    maxp0 = maximum(u0[N+1:2N])
    rangep0 = range(minp0,length=trunc(Int,sqrt(N)),stop=maxp0)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep0)
    pdfp0 = normalize(pdfp0)

    minp1 = minimum(u1[N+1:2N])
    maxp1 = maximum(u1[N+1:2N])
    rangep1 = range(minp1,length=trunc(Int,sqrt(N)),stop=maxp1)
    pdfp1 = fit(Histogram,u1[N+1:2N],rangep0)
    pdfp1 = normalize(pdfp1)

    # matplotlib[:rc]("axes", labelpad=2.)

    fig, ax = subplots(1,1,figsize=[2.325, 2.2])

    ax[:set_ylabel]("momentum distribution")
    ax[:set_xlabel](L"atom momentum $p$ (units of $\hbar k_\mathrm{c}$)")
    ax[:step](pdfp0.edges[1][1:end-1],pdfp0.weights,label="initial",where="post")
    ax[:step](pdfp1.edges[1][1:end-1],pdfp1.weights,linestyle="--",label="final",where="post")

    ax.legend(loc="upper right", bbox_to_anchor=(1.1, 1.),framealpha=1)
    # ax[:legend](handlelength=2.5)

    fig[:tight_layout]()

    return fig, ax
end

function plot_momentum_distro(sim::Array{Sol,1},filename::String)
    fig, ax = plot_momentum_distro(sim)
    fig[:savefig](filename)
end

function plot_CR_vs_S(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    categories, par_list = split_sim_from_par(sim,true)

    println(par_list)
    p0 = [900.,1.,1.,100.]

    S = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = (Float64[], Float64[])
    for x in categories
        push!(S,abs(x[1].p.S₁))

        tlist = x[1].t
        m,s,q = expect(Ekin,x)
        fit = curve_fit(model_exp, tlist, m, p0)

        p0 = fit.param

        push!(y1[1],p0[2])
        push!(y1[2],stderror(fit)[2])

        push!(y2[1],m[end])
    end

    # matplotlib[:rc]("axes", labelpad=1)
    print(y1)
    fig, ax = subplots(2,1,figsize=[3.4, 4.3],sharex=true)

    ax[2][:set_xlabel](L"pump strength $\abs{S_1}=\abs{S_2}$")

    ax[1][:set_ylabel](L"cooling rate (units of $boh$)")
    ax[1][:errorbar](S,y1[1],yerr=y1[2],color="C0",fmt=".")

    ax[2][:set_ylabel](L"final kinetic energy (units of $\hbar\kappa$)")
    ax[2][:plot](S,y2[1],color="C0")


    fig[:tight_layout]()
    return fig, ax
end

function plot_CR_vs_S(solorsim,filename::String)
    fig, ax = plot_CR_vs_S(solorsim)
    fig[:savefig](filename)
end

function plot_adaga_vs_S(sims::Array{Sol,1}...)
    colorlist = ["C1","C2","C3"]
    linelist = ["-","--",":"]

    fig, ax = subplots(1,1,figsize=[3.4, 2.3])

    ax[:set_ylabel](L"cavity population $\langle a^\dag a \rangle$")
    ax[:set_xlabel](L"pump strength $\abs{S}$")

    for (i,sim) in enumerate(sims)
        println("found $(size(sim)[end]) trajectories")
        println("-----------------------------------------")

        categories, par_list = split_sim_from_par(sim,true)

        S = Float64[]
        # y1 = (Float64[], Float64[], Array{Float64,1}[])
        y1 = Float64[]
        y2 = Array{Float64,1}[]
        y3 = (Float64[], Float64[])
        y4 = (Float64[], Float64[])
        for x in categories
            push!(S,abs(x[1].p.S₁))

            m,s,q = expect(adaga,x)

            push!(y1,m[end])
            push!(y2,q[end])
        end

        A = sortslices(hcat(S,y1,vcat(y2'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]

        # matplotlib[:rc]("axes", labelpad=1)

        label = "\$\\Delta_e="*string(par_list[1].Δₑ/par_list[1].κ)*"\\kappa\$"
        ax[:plot](S,y1,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[:fill_between](S,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    ax[:legend](handlelength=2.5)
    fig[:tight_layout]()

    return fig, ax
end

function plot_adaga_vs_S(filename::String,sims::Array{Sol,1}...)
    fig, ax = plot_adaga_vs_S(sims...)
    fig[:savefig](filename)
end

function plot_ordering_vs_S(sims::Array{Sol,1}...)
    # dont remember what it does, use the one without _def

    colorlist = ["C1","C2","C3","C4"]
    linelist = ["-","--",":","-."]

    fig, ax = subplots(3,1,figsize=[2.1, 3.3],sharex=true)
    # fig, ax = subplots(4,1,figsize=[3.4, 6.3],sharex=true)

    ax[1][:set_ylabel](L"$\langle a^\dag a \rangle$")
    ax[2][:set_ylabel](L"$\vert\Phi\vert$")
    ax[3][:set_ylabel](L"$\mathcal{B}$")
    # ax[4][:set_ylabel](L"final kinetic energy $E_\mathrm{kin}$ (units of $\hbar\kappa$)")
    ax[end][:set_xlabel](L"pump strength $S$")

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
            push!(S,abs(x[1].p.S₁))

            m,s,q = expect(adaga,x)
            push!(y1,m[end])
            push!(y2,q[end])

            m,s,q = expect(absX,x)
            push!(y3,m[end])
            push!(y4,q[end])

            m,s,q = expect(Cos2,x)
            push!(y5,m[end])
            push!(y6,q[end])

            m,s,q = expect(Ekin,x)
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

        if par_list[1].Δₑ != 0.0
            label = "\$\\Delta_e=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*"\$"
        else
            label = "\$\\Delta_e=0\$"
        end
        ax[1][:plot](S,y1,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[1][:fill_between](S,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[2][:plot](S,y3,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[2][:fill_between](S,y4[:,1],y4[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        ax[3][:plot](S,y5,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[3][:fill_between](S,y6[:,1],y6[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)

        # ax[4][:plot](S,y7,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        # ax[4][:fill_between](S,y8[:,1],y8[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    ax[1][:legend](handlelength=1.5,loc="upper left",bbox_to_anchor=(-0.03, 1.12),framealpha=0.)
    fig[:tight_layout](h_pad=0.)

    # for i in 1:3
    #     letter = Char(Int('a')+i-1)
    #     ax[i].text(0.05,0.87,"("*letter*")",transform=ax[i].transAxes)
    # end


    return fig, ax
end


function plot_ordering_vs_S(filename::String,sims::Array{Sol,1}...)
    fig, ax = plot_ordering_vs_S(sims...)
    fig[:savefig](filename)
end

function plot_spinspositionhisto(sim::Array{Sol,1})
    # matplotlib[:rc]("axes", labelpad=1)
    matplotlib.rc("image", cmap="viridis")

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
        u1 = join_trajectories(a,size(a[1].t)[1])
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


    # for i in 1:nrows
    #     for j in 1:2
    #         letter = Char(Int('a')+i*2+j-3)
    #         ax[i,j].text(-0.25,0.87,"("*letter*")",transform=ax[i,j].transAxes)
    #     end
    # end


    # for i in ax
    # end

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




function plot_dynamics(sim::Array{Sol,1})

    a,b,c,d = split_sim(sim)

    tlist = sim[1].t

    # matplotlib[:rc]("axes", labelpad=1)

    colorlist = ["C0","C1","C2","C3","C4"]

    linelist = ["-",":"]

    plot_grid = [[X Y], [ar ai], [Sx Sy],  [Sz Cos2]]

    fig, ax = subplots(4,2,figsize=[4.65,5.6],sharex="col")

    ax[1,1][:set_ylabel](L"$X$")#,labelpad=-4.)
    ax[1,2][:set_ylabel](L"$Y$")#,labelpad=-4.)
    ax[2,1][:set_ylabel](L"$\langle \hat{a}_{\mathrm{r}} \rangle$")
    ax[2,2][:set_ylabel](L"$\langle \hat{a}_{\mathrm{i}} \rangle$")
    ax[3,1][:set_ylabel](L"$\langle \hat J_x \rangle$")
    ax[3,2][:set_ylabel](L"$\langle \hat J_y \rangle$")
    ax[4,1][:set_ylabel](L"$\langle \hat J_z \rangle$")
    ax[4,2][:set_ylabel](L"$\mathcal{B}$")

    ax[4,1][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    ax[4,2][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)")
    # ax[1][:set_xscale]("log")
    # ax[2][:set_xscale]("log")

    for (i,part) in enumerate([b,c])
        for (j,row) in enumerate(plot_grid)
            for (k, o) in enumerate(row)
                Xs,X_std,X_q90 = expect(o,part)
                if o == ai
                    Xs = -Xs
                    X_q90 = -X_q90
                end
                X_q90 = hcat(X_q90...)
                ax[j,k][:plot](tlist.+1,Xs,color=colorlist[i],linestyle=linelist[i])
                ax[j,k][:fill_between](tlist.+1,X_q90[1,:],X_q90[2,:],color=colorlist[i],alpha=0.2)
            end
        end
    end

    ax[1,2][:get_shared_y_axes]()[:join](ax[1,1],ax[1,2])
    ax[1,2].autoscale()

    ax[2,2].get_shared_y_axes()[:join](ax[2,1],ax[2,2])
    ax[2,2].autoscale()

    ax[3,2].get_shared_y_axes()[:join](ax[3,1],ax[3,2])
    ax[3,2].autoscale()

    for i in 1:4
        for j in 1:2
            letter = Char(Int('a')+i*2+j-3)
            ax[i,j].text(-0.3,0.87,"("*letter*")",transform=ax[i,j].transAxes)
        end
    end

    fig[:tight_layout](h_pad=0.3)#(h_pad=0.1, w_pad=0.)

    return fig, ax
end

function plot_dynamics(sim::Array{Sol,1},filename::String)
    fig, ax = plot_dynamics(sim)
    fig[:savefig](filename)
end

function plot_dynamics_def(sim::Array{Sol,1})

    a,b,c,d = split_sim(sim)

    tlist = sim[1].t

    # matplotlib[:rc]("axes", labelpad=1)

    colorlist = ["C0","C1","C2","C3","C4"]

    linelist = ["-",":"]

    labellist = [L"$\Phi<0$",L"$\Phi>0$"]

    plot_grid = [[X Y], [ar ai]]

    fig, ax = subplots(2,2,figsize=[4.2,2.4],sharex="col")

    ax[1,1][:set_ylabel](L"$X$",labelpad=-4.)
    ax[1,2][:set_ylabel](L"$Y$",labelpad=-4.)
    ax[2,1][:set_ylabel](L"$\langle \hat{a}_{\mathrm{r}} \rangle$",labelpad=-4.)
    ax[2,2][:set_ylabel](L"$\langle \hat{a}_{\mathrm{i}} \rangle$",labelpad=-4.)
    # ax[3,1][:set_ylabel](L"$\langle \hat J_x \rangle$")
    # ax[3,2][:set_ylabel](L"$\langle \hat J_y \rangle$")
    # ax[4,1][:set_ylabel](L"$\langle \hat J_z \rangle$")
    # ax[4,2][:set_ylabel](L"$\mathcal{B}$")

    ax[2,1][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)",labelpad=-1.)
    ax[2,2][:set_xlabel](L"time (units of $\omega_\mathrm{r}^{-1}$)",labelpad=-1.)
    # ax[1][:set_xscale]("log")
    # ax[2][:set_xscale]("log")

    for (i,part) in enumerate([b,c])
        for (j,row) in enumerate(plot_grid)
            for (k, o) in enumerate(row)
                Xs,X_std,X_q90 = expect(o,part)
                if o == ai
                    Xs = -Xs
                    X_q90 = -X_q90
                end
                X_q90 = hcat(X_q90...)
                ax[j,k][:plot](tlist.+1,Xs,color=colorlist[i],linestyle=linelist[i],label=labellist[i])
                ax[j,k][:fill_between](tlist.+1,X_q90[1,:],X_q90[2,:],color=colorlist[i],alpha=0.2)
            end
        end
    end

    ax[1,2][:get_shared_y_axes]()[:join](ax[1,1],ax[1,2])
    ax[1,2].autoscale()

    ax[1,1].legend()

    ax[2,2].get_shared_y_axes()[:join](ax[2,1],ax[2,2])
    ax[2,2].autoscale()

    # ax[3,2].get_shared_y_axes()[:join](ax[3,1],ax[3,2])
    # ax[3,2].autoscale()

    # for i in 1:4
    #     for j in 1:2
    #         letter = Char(Int('a')+i*2+j-3)
    #         ax[i,j].text(-0.3,0.87,"("*letter*")",transform=ax[i,j].transAxes)
    #     end
    # end

    fig[:tight_layout](h_pad=0.3,w_pad=0.3)#(h_pad=0.1, w_pad=0.)

    return fig, ax
end

function plot_dynamics_def(sim::Array{Sol,1},filename::String)
    fig, ax = plot_dynamics_def(sim)
    fig[:savefig](filename)
end

function plot_vs_phi(sim::Array{Sol,1})
    println("found $(size(sim)[end]) trajectories")
    println("-----------------------------------------")

    colorlist = ["C1","C0","C2","C3","C4"]

    linelist = [":","-"]

    categories, par_list = split_sim_from_par(sim,true)

    sp = sortperm(par_list,lt=(x,y)->isless(-1//2*angle(x.S₂/x.S₁),-1//2*angle(y.S₂/y.S₁)))

    phi = Float64[]
    # y1 = (Float64[], Float64[], Array{Float64,1}[])
    y1 = (Float64[], Float64[])
    y2 = (Float64[], Float64[])
    y3 = (Float64[], Float64[])
    y4 = ((Float64[], Float64[]),(Float64[], Float64[]),(Float64[], Float64[]),(Float64[], Float64[]))
    y5 = ((Float64[], Float64[]),(Float64[], Float64[]),(Float64[], Float64[]),(Float64[], Float64[]))
    y6 = ((Float64[], Float64[]),(Float64[], Float64[]))
    y7 = ((Float64[], Float64[]),(Float64[], Float64[]))
    for x in categories[sp]
        ϕ = -1//2*angle(x[1].p.S₂/x[1].p.S₁)

        push!(phi,ϕ)

        m1,s1,q1 = expect(absX,x)
        m2,s2,q2 = expect(absY,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y2[1],y.val)
        push!(y2[2],y.err)

        m1,s1,q1 = expect(absar,x)
        m2,s2,q2 = expect(absai,x)
        a = m1[end] ± s1[end]
        b = m2[end] ± s2[end]
        y = angle(a+im*b)
        push!(y3[1],y.val)
        push!(y3[2],y.err)

        spsim = split_sim(x)

        tmp = 0
        for c in spsim
            if !isempty(c)
                tmp+=1
            end
        end
        if tmp>2
            println(ϕ)
            println(length(spsim[1]))
            println(length(spsim[2]))
            println(length(spsim[3]))
            println(length(spsim[4]))
        end

        for (i,c) in enumerate(spsim)
            if isempty(c)
                # println(string(ϕ)*"is empty for i="*string(i))
                push!(y4[i][1],NaN)
                push!(y4[i][2],NaN)
                push!(y5[i][1],NaN)
                push!(y5[i][2],NaN)
            else

                m1,s1,q1 = expect(X,c)
                m2,s2,q2 = expect(Y,c)
                a = m1[end] ± s1[end]
                b = m2[end] ± s2[end]
                y = angle(a+im*b)
                push!(y4[i][1],y.val)
                push!(y4[i][2],y.err)

                m1,s1,q1 = expect(ar,c)
                m2,s2,q2 = expect(ai,c)
                a = m1[end] ± s1[end]
                b = m2[end] ± s2[end]
                y = angle(a-im*b)
                push!(y5[i][1],y.val)
                push!(y5[i][2],y.err)
            end
        end

        spsim = split_simPhi(x)

        println(ϕ)
        println("Φ>0: "*string(length(spsim[1])))
        println("Φ<0: "*string(length(spsim[2])))

        for (i,c) in enumerate(spsim)
            if isempty(c)
                # println(string(ϕ)*"is empty for i="*string(i))
                push!(y6[i][1],NaN)
                push!(y6[i][2],NaN)
                push!(y7[i][1],NaN)
                push!(y7[i][2],NaN)
            else

                m1,s1,q1 = expect(X,c)
                m2,s2,q2 = expect(Y,c)
                a = m1[end] ± s1[end]
                b = m2[end] ± s2[end]
                y = angle(a+im*b)
                push!(y6[i][1],y.val)
                push!(y6[i][2],y.err)

                m1,s1,q1 = expect(ar,c)
                m2,s2,q2 = expect(ai,c)
                a = m1[end] ± s1[end]
                b = m2[end] ± s2[end]
                y = angle(a-im*b)
                push!(y7[i][1],y.val)
                push!(y7[i][2],y.err)
            end
        end
    end

    fig, ax = subplots(2,1,figsize=[2.7, 3.2],sharex=true)

    ax[1].set_ylabel(L"spin phase $\phi_{\mathrm{s}}$")
    # ax[1].errorbar(phi,y2[1],yerr=y2[2],color="C1",fmt=".")
    # ax[1].set_yticks(pi/8*collect(0:4))
    # ax[1].set_yticklabels([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    ax[2].set_ylabel(L"cavity field phase $\phi_\mathrm{c}$")
    # ax[2].errorbar(phi,y3[1],yerr=y3[2],color="C2",fmt=".")
    # ax[2].set_yticks(pi/8*collect(0:4))
    # ax[2].set_yticklabels([L"0",L"\pi/8",L"\pi/4",L"3\pi/8",L"\pi/2"])

    # for c in y4
    #     # println(c)
    #     ax[1].errorbar(phi,c[1],yerr=c[2])
    # end
    # ax[1].set_yticks(pi/2*collect(-2:2))
    # ax[1].set_yticklabels([L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"])


    # for c in y5
    #     ax[2].errorbar(phi,c[1],yerr=c[2])#,fmt="o")
    # end
    # ax[2].set_yticks(pi/2*collect(-2:2))
    # ax[2].set_yticklabels([L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"])


    for (i,c) in enumerate(y6)
        # println(c)
        if i==1
            ax[1].plot(phi,c[1],color=colorlist[i],linestyle=linelist[i])
            ax[1].fill_between(phi,c[1]-c[2],c[1]+c[2],color=colorlist[i],alpha=0.2)
        else
            ax[1].plot(phi[phi.<0],c[1][phi.<0],color=colorlist[i],linestyle=linelist[i])
            ax[1].fill_between(phi[phi.<0],c[1][phi.<0]-c[2][phi.<0],c[1][phi.<0]+c[2][phi.<0],
                               color=colorlist[i],alpha=0.2)
            ax[1].plot(phi[phi.>=0],c[1][phi.>=0],color=colorlist[i],linestyle=linelist[i])
            ax[1].fill_between(phi[phi.>=0],c[1][phi.>=0]-c[2][phi.>=0],c[1][phi.>=0]+c[2][phi.>=0],
                               color=colorlist[i],alpha=0.2)
        end
        # ax[1].errorbar(phi,c[1],yerr=c[2],fmt=".")
    end
    ax[1].set_yticks(pi/2*collect(-2:2))
    ax[1].set_yticklabels([L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"])


    for (i,c) in enumerate(y7)
        if i==1
            legendlabel = L"$\Phi>0$"
            ax[2].plot(phi[phi.<π/4],c[1][phi.<π/4],color=colorlist[i],linestyle=linelist[i],label=legendlabel)
            ax[2].fill_between(phi[phi.<π/4],c[1][phi.<π/4]-c[2][phi.<π/4],c[1][phi.<π/4]+c[2][phi.<π/4],
                               color=colorlist[i],alpha=0.2)
            ax[2].plot(phi[phi.>=π/4],c[1][phi.>=π/4],color=colorlist[i],linestyle=linelist[i])
            ax[2].fill_between(phi[phi.>=π/4],c[1][phi.>=π/4]-c[2][phi.>=π/4],c[1][phi.>=π/4]+c[2][phi.>=π/4],
                               color=colorlist[i],alpha=0.2)
        else
            legendlabel = L"$\Phi<0$"
            ax[2].plot(phi,c[1],color=colorlist[i],linestyle=linelist[i],label=legendlabel)
            ax[2].fill_between(phi,c[1]-c[2],c[1]+c[2],color=colorlist[i],alpha=0.2)
        end

        # ax[2].errorbar(phi,c[1],yerr=c[2],label=legendlabel,fmt=".")
    end
    ax[2].set_yticks(pi/2*collect(-2:2))
    ax[2].set_yticklabels([L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"])



    ax[2].set_xlabel(L"phase $\phi$")
    ax[2].set_xticks(pi/4*collect(-2:2))
    ax[2].set_xticklabels([L"-\pi/2",L"-\pi/4",L"0",L"\pi/4",L"\pi/2"])


    fig.tight_layout(h_pad=0.4)#(h_pad=0,w_pad=0.5)
    ax[2].legend(loc="upper left", bbox_to_anchor=(0.02, 1.3),framealpha=1)

    return fig, ax
end

function plot_vs_phi(solorsim,filename::String)
    fig, ax = plot_vs_phi(solorsim)
    fig.savefig(filename)
end

function plot_initial_conditions(sim::Array{Sol,1})

    u0 = join_trajectories(sim,1)

    N::Int = sum(try sol.p.N catch; sol.p[10] end for sol in sim) # for backwards compatibility
    nbins = trunc(Int, sqrt(N))

    rangex = range(0,length=trunc(Int,sqrt(N)),stop=2pi)
    pdfx0 = fit(Histogram,mod2pi.(u0[1:N]),rangex)
    pdfx0 = normalize(pdfx0)

    minp = minimum(u0[N+1:2N])
    maxp = maximum(u0[N+1:2N])
    rangep = range(minp,length=trunc(Int,sqrt(N)),stop=maxp)
    pdfp0 = fit(Histogram,u0[N+1:2N],rangep)
    pdfp0 = normalize(pdfp0)

    ranges = range(-1.1,length=trunc(Int,sqrt(N)),stop=1.1)
    pdfsx = fit(Histogram,u0[2N+1:3N],ranges)
    pdfsy = fit(Histogram,u0[3N+1:4N],ranges)
    pdfsz = fit(Histogram,u0[4N+1:5N],ranges)
    pdfsx = normalize(pdfsx)
    pdfsy = normalize(pdfsy)
    pdfsz = normalize(pdfsz)


    matplotlib[:rc]("axes")

    fig, ax = subplots(2,2,figsize=[4.65,3.6])

    ax[1, 1][:set_ylabel]("distribution")
    ax[1, 1][:set_xlabel](L"atom position $x$")
    ax[1, 1][:set_ylim]([0.,maximum(pdfx0.weights)*1.1])
    ax[1, 1][:set_xticks](pi/2*collect(0:4))
    ax[1, 1][:set_xticklabels]([L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"])
    ax[1, 1][:step](pdfx0.edges[1][1:end-1],pdfx0.weights,label="initial",where="post")

    ax[1, 2][:set_ylabel]("distribution")
    ax[1, 2][:set_xlabel](L"atom momentum $p$")
    ax[1, 2][:step](pdfp0.edges[1][1:end-1],pdfp0.weights,label="initial",where="post")

    ax[2, 1][:set_ylabel]("distribution")
    ax[2, 1][:set_xlabel](L"spins $\sigma^x$ and $\sigma^y$")
    ax[2, 1][:set_xlabel](L"spins $\sigma^x$ and $\sigma^y$")
    ax[2, 1][:step](pdfsx.edges[1][1:end-1],pdfsx.weights,label=L"\sigma^x",where="post")
    ax[2, 1][:step](pdfsy.edges[1][1:end-1],pdfsy.weights,label=L"\sigma^y",where="post",linestyle="--")
    ax[2, 1][:legend](handlelength=2.5)

    ax[2, 2][:set_ylabel]("distribution")
    ax[2, 2][:set_xlabel](L"spins $\sigma^z$")
    ax[2, 2][:step](pdfsz.edges[1][1:end-1],pdfsz.weights,label=L"\sigma^z",where="post")


    # fig[:tight_layout](h_pad=0., w_pad=-0.)
    fig[:tight_layout](h_pad=0.7,w_pad=0.7)

    return fig, ax
end

function plot_initial_conditions(sim::Array{Sol,1},filename::String)
    fig, ax = plot_initial_conditions(sim)
    fig[:savefig](filename)
end



function plot_single_ordering_vs_S(sims::Array{Sol,1}...)

    colorlist = ["C1","C2","C3","C4"]
    linelist = ["-","--",":","-."]

    fig, ax = plt[:subplots](1, 1, figsize=[3.5, 2.2])

    ax[:set_ylabel](L"$\langle a^\dag a \rangle$")
    ax[:set_xlabel](L"pump strength $S$")

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
            push!(S,abs(x[1].p.S₁))

            m,s,q = expect(adaga,x)
            push!(y1,m[end])
            push!(y2,q[end])
        end

        A = sortslices(hcat(S,y1,vcat(y2'...)),dims=1)
        S = A[:,1]
        y1 = A[:,2]
        y2 = A[:,3:4]

        if par_list[1].Δₑ != 0.0
            label = "\$\\Delta_e=\\kappa/"*string(trunc(Int,par_list[1].κ/par_list[1].Δₑ))*"\$"
        else
            label = "\$\\Delta_e=0\$"
        end
        ax[:plot](S,y1,ls=linelist[i],color=colorlist[i],label=latexstring(label))
        ax[:fill_between](S,y2[:,1],y2[:,2],color=colorlist[i],alpha=0.2,linewidth=0.1)
    end

    ax[:legend](handlelength=1.5,loc="upper left",bbox_to_anchor=(-0.03, 1.12),framealpha=0.)
    fig[:tight_layout](h_pad=0.)

    # for i in 1:3
    #     letter = Char(Int('a')+i-1)
    #     ax[i].text(0.05,0.87,"("*letter*")",transform=ax[i].transAxes)
    # end


    return fig, ax
end


function plot_single_ordering_vs_S(filename::String,sims::Array{Sol,1}...)
    fig, ax = plot_single_ordering_vs_S(sims...)
    fig[:savefig](filename)
end

end



