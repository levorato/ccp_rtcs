include("config.jl")
include("RCCP_Robust_Simulation.jl")

EXPERIMENT_NAME = "run_sim_forecast_avg"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

function concatenation()
    scenario_logger = open(joinpath(tempdir(), "rccp_concat_log.txt"), "w+")
    #concatenate_trace_df_list("Antoine", "A_instance2_100s.txt", scenario_logger, base_path)
    #concatenate_trace_df_list("Antoine", "A_instance2_1NDU_100s.txt", scenario_logger, base_path)
    concatenate_trace_df_list("deterministic", "antoine-skew", "A_instance2_1000s_skewed-right.txt", scenario_logger)
    #concatenate_trace_df_list("Antoine", "A_instance2_1NDU_1000s.txt", scenario_logger, base_path)
    close(scenario_logger)
end

concatenation()
