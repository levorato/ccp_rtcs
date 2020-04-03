# =================================
# RCCP_det_Tests.jl
# RCCP Deterministic Tests
# =================================

# Include RCCP Deterministic model
include("RCCP_det.jl")
include("RCCP_Robust.jl")

time_limit = 7200
base_folder = pwd() * "/instances/"
test_set = ["Example_3.1_A.txt", "Example_3.1_B.txt"]
output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "log"])
df = DataFrame(Instance = String[], TestName = String[], ModelName = String[],
    Value = Float64[], Time = Float64[], IsOptimal = Bool[], Solution = String[])
test_name = "Example"

for instance_name in test_set
    instance_as_dict = read_tabulated_data(joinpath(normpath(base_folder), instance_name))
    output_file_log = joinpath(normpath(output_path), "RCCP_Robust_Tests" * "_" * instance_name * ".log")
    output_file_inf_log = joinpath(normpath(output_path), "RCCP_Robust_Tests" * "_" * instance_name * "_Infeasible.log")
    scenario_logger = open(output_file_log, "w+")
    infeasible_logger = open(output_file_inf_log, "w+")
    println("***** Testing instance = $(instance_name) *****")
    sol_value, solution, time, status = solve_robust_model(instance_as_dict, scenario_logger, infeasible_logger, time_limit, true, 1)
    y = solution["y"]
    println("Solution: y = $(y)")
    push!(df, [instance_name, test_name, "Robust_RCCP_Example", sol_value, time, status, string(y)])
    close(scenario_logger)
    close(infeasible_logger)
end

# Create an output file to write the dataframe to disk
output_path = joinpath(normpath(pwd()), "output")
if ! isdir(output_path)
    mkdir(output_path)
end
output_file = joinpath(normpath(output_path), test_name * "_Robust_RCCP_Test_Example_3.1.csv")
println("\nSaving to results file to $(output_file)...")
CSV.write(output_file, df)
