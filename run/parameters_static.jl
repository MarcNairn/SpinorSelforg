# parameters
N = 100 # number of particles
U₁ = U₂ = 0.0
κ = 100. # decay rate of the cavity field (units of ω\_R)
S₁ = 5κ/sqrt(N)
S₂ = 5κ/sqrt(N)
Δₑ = 10
Δc = -κ
num_monte = 50

tspan = (0.0, 400.0)

temp = 10
p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)