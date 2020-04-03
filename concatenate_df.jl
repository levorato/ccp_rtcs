include("RCCP_Robust_Simulation.jl")

function concatenation()
    scenario_logger = open("/tmp/rccp_concat_log.txt", "w+")
    base_path = "/home/mlevorato/rccp_experiments"
    #concatenate_trace_df_list("Antoine", "A_instance2_100s.txt", scenario_logger, base_path)
    concatenate_trace_df_list("Antoine", "A_instance2_1NDU_100s.txt", scenario_logger, base_path)
    concatenate_trace_df_list("Antoine", "A_instance2_1000s.txt", scenario_logger, base_path)
    #concatenate_trace_df_list("Antoine", "A_instance2_1NDU_1000s.txt", scenario_logger, base_path)
    close(scenario_logger)
end

concatenation()
