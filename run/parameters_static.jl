#Read pump strength from parameters
g = parse(Int, ARGS[2]) #in the case kappa=!100 or N=!100 need to renormalise the coupling. 
temp = parse(Int, ARGS[3]) #same goes for here

# parameters
N = 100 # number of particles
U₁ = U₂ = 0.0
κ = 100. # decay rate of the cavity field (units of ω\_R)
S₁ = g
S₂ = g
Δₑ = 10
Δc = -κ
num_monte = 1 #increased momentarily to optimise full parameter array in slurm queue 
#for num_monte>>1 need ~6-12hrs on cluster
#recall per trajectory need around 20 minutes

tspan = (0.0, 800.0)



p = System_p(U₁,U₂,S₁,S₂,Δₑ,Δc,κ,temp,N,tspan,num_monte)
