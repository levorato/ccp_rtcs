# =============================================================================
# run_sim.jl
# Run CCP RTCS simulation for deterministic, robust-box and
# robust-budget models. Stores results in experiment output-folder.
# =============================================================================

include("cli_arguments.jl")
include("RCCP_Facade.jl")

EXPERIMENT_NAME = "run_sim_japan_forecast_avg"

parsed_args = parse_commandline()
# If possible, do not modify the lines below
parsed_args["output-folder"] = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

# Solve models
#parsed_args["what"] = "solve"
#parsed_args["solver-verbose"] = true
#parsed_args["solver-time-limit"] = 28800.0  # 8h time-limit
#parsed_args["max-cores"] = 16
#parsed_args["instances"] = "antoine-skew"

# Simulate based on existing model results
models_to_simulate = ["robust-budget", "robust-box", "deterministic"]
instance_group_list = ["japan-10"]
parsed_args["what"] = "simulate"
parsed_args["forecast-type"] = "average"  # average-based RTCS forecast
for instance_group in instance_group_list
	for model in models_to_simulate
		parsed_args["model"] = model
		parsed_args["instances"] = instance_group
		parsed_args["instance-name"] = "instance_deltamin10_winter.txt"
		parsed_args["simulation-nthreads"] = 12
		parsed_args["final-scenario"] = 100
		run_experiment(parsed_args)
	end
end
