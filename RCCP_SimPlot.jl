# ====================================
# RCCP_SimPlot.jl
# RTCS Simulation plotting functions.
# ====================================

using PyPlot
using DataFrames, IndexedTables

# Include RCCP Simulation utility functions
include("RCCP_SimUtil.jl")

# https://github.com/JuliaPy/PyPlot.jl
function plot_simulation_data()
    test_set = ["A_instance2_11scen_1contract.txt"]  # "A_instance2_11scen_1NDU.txt" , "A_instance2_11scen.txt"]  # "A_instance2_11scen_1contract.txt"
    test_name = "Antoine"
    for instance_name in test_set
        println("Plotting simulation data for instance $(instance_name)")
        trace_df = read_trace_dataframe(test_name, instance_name)

    end
    println("\nDone.")


    x = linspace(0,2*pi,1000)
    y = sin.(3 * x + 4 * cos.(2 * x))
    #x = range(0; stop=2*pi, length=1000); y = sin.(3 * x + 4 * cos.(2 * x));
    plot(x, y, color="red", linewidth=2.0, linestyle="--")
    title("A sinusoidally modulated sinusoid")
end

#plot_simulation_data()
