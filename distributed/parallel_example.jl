using Distributed   
@everywhere using ClusterManagers

println("Workers before: ", nworkers())
println("Processes before: ", nprocs())


np = parse(Int,ARGS[1])
array = parse(Int, ARGS[2])

addprocs(SlurmManager(np))

println("Workers after: ", nworkers())
println("Processes afer: ", nprocs())

println("Passed array size: ", array)