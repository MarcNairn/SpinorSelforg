using Pkg
Pkg.activate("..")

using DifferentialEquations: solve, SOSRA

include("src/selforg_core.jl") 
