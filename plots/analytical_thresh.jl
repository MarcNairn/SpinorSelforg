using QuadGK
using PyPlot


# Define parameters
selected_deltaDv = [5, 10, 20, 100]
omegauv = range(0, stop=100, length=201)

gc_values = Vector{Vector{Float64}}(undef, length(selected_deltaDv))

opacities = range(1.0, stop=0.2, length=length(selected_deltaDv)) # Vary opacities

for (idx, deltaDv) in enumerate(selected_deltaDv)
    gc = zeros(Float64, length(omegauv))

    for (bb, omegau) in enumerate(omegauv)
        N = 100
        sim = 10
        kappa = 100
        Deltac = 100
        omegar = 1 / 100 * deltaDv / (2 * π)

        fun(t) = exp(-deltaDv^2 * t^2 / 2) * (sin(omegau * t) + cos(omegau * t) * omegar * t)
        
        # Numerical integration using QuadGK
        Int, _ = quadgk(fun, 0, Inf;order=9)

        gc[bb] = sqrt((kappa^2 + Deltac^2 / (2 * Deltac) / Int) / N)
    end

    gc_values[idx] = gc
end

# Create the plot with varying opacities
figure()

for (idx, deltaDv) in enumerate(selected_deltaDv)
    plot(omegauv, gc_values[idx], label="δ_D= $deltaDv", linewidth=1, color="#a6cee3" ,alpha=opacities[end-idx+1])
end

xlabel(L"Atomic frequency $\omega_\uparrow$", fontsize=20)
ylabel(L"Coupling strength $g_c$", fontsize=20)
grid(false)
legend(loc="best", fontsize=14)
ylim(9, 40)
fontsize = 18
xticks(fontsize=fontsize)
yticks(fontsize=fontsize)
gca().set_position([0.125, 0.125, 0.775, 0.775])  # Adjust the position of the axis as needed

# Set the DPI for the figure
PyPlot.gcf().set_size_inches(8, 6)  # Adjust the width and height as needed
# savefig("threshold_Doppler_inset.png", dpi=600)  # Save the plot with x DPI
show()
