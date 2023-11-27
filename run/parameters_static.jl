#Read pump strength from parameters
g = parse(Int, ARGS[2]) #in the case kappa=!100 or N=!100 need to renormalise the coupling. 

# parameters
N = 100 # number of particles
U₁ = U₂ = 0.0
κ = 100. # decay rate of the cavity field (units of ω\_R)
S₁ = g
S₂ = g
Δₑ = 10
Δc = -κ
num_monte = 1 

tspan = (0.0, 800.0)

temp = 10
p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)
