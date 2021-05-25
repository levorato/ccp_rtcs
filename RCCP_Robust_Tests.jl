# RCCP_Robust_Tests.jl

include("RCCP_Robust.jl")

time_limit = 7200
base_folder = antoine_instances_folder
test_set = ["A_instance2_11scen_1NDU.txt"]  # "A_instance2_11scen.txt"
#test_set = ["A_instance2_11scen_noNDU.txt", "A_instance2_11scen_noNDU.txt"]
test_name = "Antoine"
output_path = create_full_dir(normpath(EXPERIMENT_OUTPUT_FOLDER), ["output", "simulation", "log"])
df = DataFrame(Instance = String[], TestName = String[], ModelName = String[],
    Value = Float64[], Time = Float64[], IsOptimal = Bool[])
for instance_name in test_set
    instance_as_dict = read_tabulated_data(joinpath(normpath(base_folder), instance_name))
    output_file_log = joinpath(normpath(output_path), "RCCP_Robust_Tests" * "_" * instance_name * ".log")
    output_file_inf_log = joinpath(normpath(output_path), "RCCP_Robust_Tests" * "_" * instance_name * "_Infeasible.log")
    scenario_logger = open(output_file_log, "w+")
    infeasible_logger = open(output_file_inf_log, "w+")
    for rel_gap_tol in [0.0005] #, 0.01, 0.005] # [0.10, 0.075, 0.05]  # 10%, 7,5%, 5%
        println("***** Testing rel_gap_tol = $(rel_gap_tol) *****")
        sol_value, y, time, status = solve_robust_model(instance_as_dict, scenario_logger, infeasible_logger, time_limit, true, 1, rel_gap_tol)
        push!(df, [instance_name, test_name, "Robust_RCCP_GapTol" * string(rel_gap_tol), sol_value, time, status])
    end
    # Create an output file to write the dataframe to disk
    output_path = joinpath(normpath(EXPERIMENT_OUTPUT_FOLDER), "output")
    mkpath(output_path)
    output_file = joinpath(normpath(output_path), test_name * "_Robust_RCCP_Test_GapTol.csv")
    println("\nSaving to results file to $(output_file)...")
    CSV.write(output_file, df)
    close(scenario_logger)
    close(infeasible_logger)
end
