
#define parameter struct
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

function vect_f_det(du, u, p, t)
    N::Int = p.N #Atom number
    N_MC::Int = p.N_MC  # Number of Monte Carlo samples
#=
    xv = u[1:N, :] # extract random positisions
    pv = u[N+1:2N, :] # extract random momenta
    
    sxv = u[2N+1:3N, :]
    syv = u[3N+1:4N, :] # extract σˣⱼ and σʸⱼ 
    szv = u[4N+1:5N, :] # extract σᶻⱼ = -1, atoms in the ground state
    
    axv = u[5N+1, :]
    apv = u[end, :]
    
    dax = du[1:N, :]
    dap = du[N+1:2N, :]

    dsx = du[2N+1:3N, :]
    dsy = du[2N+1:3N, :]
    dsz = du[4N+1:5N, :]

    dx = du[end-1, :]
    dp = du[end, :]


    dsx=(-p.Δₑ.*syv)
    dsy=(p.Δₑ.*sxv .- 2*p.S₁*cos(xv).*szv.*axv)
    dsz=(2*p.S₁*cos(xv).*syv.*axv)           

    dx=2*pv; 
    dp=p.S₁*sxv.*sin(xv).*axv


    dax=(p.Δc*apv .- p.κ*axv)
    dap=(-p.Δc*axv .- p.κ*apv .- 2*p.S₁*sum(cos(xv).*sxv))

=#
    REP::Float64 = real(p.S₁) + real(p.S₂) # ℜ(S₁) + ℜ(S₂)
    REM::Float64 = real(p.S₁) - real(p.S₂) # ℜ(S₁) - ℜ(S₂)
    IMP::Float64 = imag(p.S₁) + imag(p.S₂) # ℑ(S₁) + ℑ(S₂)
    IMM::Float64 = imag(p.S₁) - imag(p.S₂) # ℑ(S₁) - ℑ(S₂)
    aa::Vector{Float64} = (u[5N+1, :].^2 + u[5N+2, :].^2 .- 0.5) # aᵣ² + aᵢ² - 1/2
    bb::Float64 = 2p.Δc
    sinuj = sin.(u[1:N, :])
    cosuj = cos.(u[1:N, :])
    cc= (IMM * u[2N+1:3N, :] + REM * u[3N+1:4N, :]) .* cosuj
    dd= (REP * u[2N+1:3N, :] - IMP * u[3N+1:4N, :]) .* cosuj
    # sincos.(u[j, :])

    #bb = 0((1 - u[4N+j, :]) * p.U₁ + (1 + u[4N+j, :]) * p.U₂) * cosuj^2
    #cc = 
    #dd = (REP * u[2N+1:3N, :] - IMP * u[3N+1:4N, :]) * cosuj

    du[1:N, :] .= 2u[N+1:2N, :]
    du[N+1:2N, :] .= sinuj .* (cosuj .* ((1 .- u[4N+1:5N, :]) * p.U₁ .+ (1 .+ u[4N+1:5N, :]) * p.U₂) .* aa .+
        (REP * u[2N+1:3N, :] .* u[5N+1, :] .+ REM * u[3N+1:4N, :] .* u[5N+2, :] .+
        IMM * u[2N+1:3N, :] .* u[5N+2, :] .- IMP * u[3N+1:4N, :] .* u[5N+1, :]))
    du[2N+1:3N, :] .= -u[3N+1:4N, :] * (p.Δₑ .- (p.U₁ .- p.U₂) .* cosuj^2 .* aa) .-
        2cosuj .* u[4N+1:5N, :] .* (IMP * u[5N+1, :] - REM * u[5N+2, :])
    du[3N+1:4N, :] .= u[2N+1:3N, :] .* (p.Δₑ .- (p.U₁ - p.U₂) * cosuj.^2 .* aa) .-
        2cosuj .* u[4N+1:5N, :] .* (REP * u[5N+1, :] + IMM * u[5N+2, :])
    du[4N+1:5N, :] .= 2cosuj .* (IMP * u[2N+1:3N, :] .* u[5N+1, :] +
        IMM * u[3N+1:4N, :] .* u[5N+2, :] - REM * u[2N+1:3N, :] .* u[5N+2, :] +
        REP * u[3N+1:4N, :] .* u[5N+1, :])   
    du[5N+1, :] = bb/2 * u[5N+2, :] + cc/2 #- p.κ * u[5N+1, :]
    du[5N+2, :] .= -bb/2 * u[5N+1, :] + dd/2 #- p.κ * u[5N+2, :]
end

function vect_f_noise(du, u, p, t)
    N::Int = p.N
    N_MC::Int = p.N_MC
    du[1:5N, :] .= 0.0
    du[5N+1, :] .= sqrt(2*p.κ)
    du[5N+2, :] .= sqrt(2*p.κ)
end

function vect_initial_conditions(p::System_p, seed=abs(rand(Int)))
    N::Int = p.N
    N_MC::Int = p.N_MC
    Random.seed!(seed) # random number generator
    u0 = zeros(5N + 2, N_MC)

    u0[1:N, :] .= 2pi.*rand(N, N_MC) # generate random positions
    u0[N+1:2N, :] .= p.temp .* randn(N, N_MC) # generate random momenta

    u0[2N+1:4N, :] .= 0.1 .*randn(2N, N_MC) # σˣⱼ and σʸⱼ 
    u0[4N+1:5N, :] .= -1. # σᶻⱼ = -1, atoms in the ground state

    u0[5N+1:end, :] .= 0.0 # cavity empty

    return u0
end


function vect_define_prob_from_parameters(p::System_p,seed=abs(rand(Int)))
    # initial conditions
    Random.seed!(seed) # random number generator
    u0 = vect_initial_conditions(p,abs(rand(Int)))

    prob = SDEProblem(vect_f_det,vect_f_noise,u0,p.tspan,p)

    return prob
end

function vect_traj_solver(p::System_p;saveat::Float64=1.0,seed::Int=abs(rand(Int)),maxiters::Int=Int(1e9))
    prob = vect_define_prob_from_parameters(p,seed)
    elt = @elapsed sim = solve(prob::SDEProblem, SOSRA(), saveat=saveat, maxiters=maxiters, progress=true)
    println("done in $elt seconds")
    return sim
end
