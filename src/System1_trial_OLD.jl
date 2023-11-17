module System1_trial
# this module contains the equations describing system1

export f  # this is the last chosen function - see at the bottom
export f_noise

N = 100

function f1(du,u,p,t)
    # t is in unit of w_R

    # u = [ xⱼ,     pⱼ,      σˣⱼ,       σʸⱼ,      σᶻⱼ,    aᵣ,    aᵢ]
    #     [1:N, N+1:2N, 2N+1..3N, 3N+1..4N, 4N+1..5N, 5N+1, 5N+2]
    # x_j in units of 1/k (x'_j = k x_j)
    # p_j in units of hbar*k (p'_j = p_j/(hbar*k))

    # parameters p are:
    # p[1] = U₁
    # p[2] = U₂
    # p[3] = ℜ(S₁) + ℜ(S₂)
    # p[4] = ℜ(S₁) - ℜ(S₂)
    # p[5] = ℑ(S₁) + ℑ(S₂)
    # p[6] = ℑ(S₁) - ℑ(S₂)
    # p[7] = Δₑ
    # p[8] = Δ\_c
    # p[9] = κ

    # ancilla variables
    aa = (u[5N+1]^2 + u[5N+2]^2 - 0.5)::Float64 # aᵣ² + aᵢ² - 1/2
    bb = 2p[8]::Float64
    cc = 0.0
    dd = 0.0

    for j = 1:N
        # positions x_j
        du[j] = 2u[N+j]
        # momenta p_j
        du[N+j] = cos(u[j])sin(u[j]) * ( (1-u[4N+j])p[1] + (1+u[4N+j])p[2] ) *
            aa + sin(u[j]) * ( p[3]u[2N+j]u[5N+1] + p[4]u[3N+j]u[5N+2] +
                               p[6]u[2N+j]u[5N+2] - p[5]u[3N+j]u[5N+1] )
        # sigmax_j
        du[2N+j] = - u[3N+j] * ( p[7] - (p[1]-p[2])*cos(u[j])^2*aa ) -
            2cos(u[j])*u[4N+j] * ( p[5]u[5N+1] - p[4]u[5N+2] )
        # sigmay_j
        du[3N+j] = u[2N+j] * ( p[7] - (p[1]-p[2])*cos(u[j])^2*aa ) -
            2cos(u[j])*u[4N+j] * ( p[3]u[5N+1] + p[6]u[5N+2] )
        # sigmaz_j
        du[4N+j] = 2cos(u[j]) * ( p[5]u[2N+j]u[5N+1] + p[6]u[3N+j]u[5N+2] -
                                  p[4]u[2N+j]u[5N+2] + p[3]u[3N+j]u[5N+1] )

        # ancilla variables
        bb -= ((1-u[4N+j])p[1] + (1+u[4N+j])p[2])cos(u[j])^2
        cc -= ( p[6]u[2N+j] + p[4]u[3N+j] )cos(u[j])
        dd += ( p[3]u[2N+j] - p[5]u[3N+j] )cos(u[j])
    end

    # a_r
    du[5N+1] = bb/2 * u[5N+2] + cc/2 - p[9]u[5N+1]
    # a_i
    du[5N+2] = -bb/2 * u[5N+1] + dd/2 - p[9]u[5N+2]

end

function f2(du,u,p,t)
    # t is in unit of w_R

    # u = [ xⱼ,     pⱼ,      σˣⱼ,       σʸⱼ,      σᶻⱼ,    aᵣ,    aᵢ]
    #     [1:N, N+1:2N, 2N+1..3N, 3N+1..4N, 4N+1..5N, 5N+1, 5N+2]
    # x_j in units of 1/k (x'_j = k x_j)
    # p_j in units of hbar*k (p'_j = p_j/(hbar*k))

    # parameters p are:
    # p[1] = U₁
    # p[2] = U₂
    # p[3] = ℜ(S₁) + ℜ(S₂)
    # p[4] = ℜ(S₁) - ℜ(S₂)
    # p[5] = ℑ(S₁) + ℑ(S₂)
    # p[6] = ℑ(S₁) - ℑ(S₂)
    # p[7] = Δₑ
    # p[8] = Δ\_c
    # p[9] = κ

    # ancilla variables
    aa = (u[5N+1]^2 + u[5N+2]^2 - 0.5)::Float64 # aᵣ² + aᵢ² - 1/2
    bb = 2p[8]::Float64
    cc = 0.0
    dd = 0.0

    for j = 1:N
        # ancilla variables
        sinuj, cosuj = sincos(u[j])


        # positions x_j
        du[j] = 2u[N+j]
        # momenta p_j
        du[N+j] = sinuj * ( cosuj*( (1-u[4N+j])p[1] + (1+u[4N+j])p[2] ) * aa +
                            ( p[3]u[2N+j]u[5N+1] + p[4]u[3N+j]u[5N+2] +
                              p[6]u[2N+j]u[5N+2] - p[5]u[3N+j]u[5N+1] ) )
        # sigmax_j
        du[2N+j] = - u[3N+j] * ( p[7] - (p[1]-p[2])*cosuj^2*aa ) -
            2cosuj*u[4N+j] * ( p[5]u[5N+1] - p[4]u[5N+2] )
        # sigmay_j
        du[3N+j] = u[2N+j] * ( p[7] - (p[1]-p[2])*cosuj^2*aa ) -
            2cosuj*u[4N+j] * ( p[3]u[5N+1] + p[6]u[5N+2] )
        # sigmaz_j
        du[4N+j] = 2cosuj * ( p[5]u[2N+j]u[5N+1] + p[6]u[3N+j]u[5N+2] -
                                  p[4]u[2N+j]u[5N+2] + p[3]u[3N+j]u[5N+1] )

        # ancilla variables
        bb -= ((1-u[4N+j])p[1] + (1+u[4N+j])p[2])cosuj^2
        cc -= ( p[6]u[2N+j] + p[4]u[3N+j] )cosuj
        dd += ( p[3]u[2N+j] - p[5]u[3N+j] )cosuj
    end

    # a_r
    du[5N+1] = bb/2 * u[5N+2] + cc/2 - p[9]u[5N+1]
    # a_i
    du[5N+2] = -bb/2 * u[5N+1] + dd/2 - p[9]u[5N+2]

end

function f(du,u,p,t)
    # t is in unit of w_R

    # u = [ xⱼ,     pⱼ,      σˣⱼ,       σʸⱼ,      σᶻⱼ,    aᵣ,    aᵢ]
    #     [1:N, N+1:2N, 2N+1..3N, 3N+1..4N, 4N+1..5N, 5N+1, 5N+2]
    # x_j in units of 1/k (x'_j = k x_j)
    # p_j in units of hbar*k (p'_j = p_j/(hbar*k))

    # parameters p are:
    # p[1] = U₁
    # p[2] = U₂
    # p[3] = ℜ(S₁) + ℜ(S₂)
    # p[4] = ℜ(S₁) - ℜ(S₂)
    # p[5] = ℑ(S₁) + ℑ(S₂)
    # p[6] = ℑ(S₁) - ℑ(S₂)
    # p[7] = Δₑ
    # p[8] = Δ\_c
    # p[9] = κ
    # p[10] = N

    # ancilla variables
    N::Int = p[10]
    aa::Float64 = (u[5N+1]^2 + u[5N+2]^2 - 0.5) # aᵣ² + aᵢ² - 1/2
    bb::Float64 = 2p[8]
    cc::Float64 = 0.0
    dd::Float64 = 0.0

    for j = 1:N
        # ancilla variables
        sinuj::Float64, cosuj::Float64 = sincos(u[j])
        bb -= ((1-u[4N+j])p[1] + (1+u[4N+j])p[2])cosuj^2
        cc -= ( p[6]u[2N+j] + p[4]u[3N+j] )cosuj
        dd += ( p[3]u[2N+j] - p[5]u[3N+j] )cosuj

        # positions x_j
        du[j] = 2u[N+j]
        # momenta p_j
        du[N+j] = sinuj * ( cosuj*( (1-u[4N+j])p[1] + (1+u[4N+j])p[2] ) * aa +
                            ( p[3]u[2N+j]u[5N+1] + p[4]u[3N+j]u[5N+2] +
                              p[6]u[2N+j]u[5N+2] - p[5]u[3N+j]u[5N+1] ) )
        # sigmax_j
        du[2N+j] = - u[3N+j] * ( p[7] - (p[1]-p[2])*cosuj^2*aa ) -
            2cosuj*u[4N+j] * ( p[5]u[5N+1] - p[4]u[5N+2] )
        # sigmay_j
        du[3N+j] = u[2N+j] * ( p[7] - (p[1]-p[2])*cosuj^2*aa ) -
            2cosuj*u[4N+j] * ( p[3]u[5N+1] + p[6]u[5N+2] )
        # sigmaz_j
        du[4N+j] = 2cosuj * ( p[5]u[2N+j]u[5N+1] + p[6]u[3N+j]u[5N+2] -
                                  p[4]u[2N+j]u[5N+2] + p[3]u[3N+j]u[5N+1] )
    end

    # a_r
    du[5N+1] = bb/2 * u[5N+2] + cc/2 - p[9]u[5N+1]
    # a_i
    du[5N+2] = -bb/2 * u[5N+1] + dd/2 - p[9]u[5N+2]

end

function f_noise(du,u,p,t)
    N::Int = p[10]
    du[1:5N] .= 0.0
    du[5N+1] = sqrt(2*p[9])
    du[5N+2] = sqrt(2*p[9])
end

end
