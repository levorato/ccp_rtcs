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
###models_to_simulate = ["deterministic"]  # ["robust-budget", "robust-box", "deterministic"]
instance_group_list = ["japan-10"]
parsed_args["what"] = "simulate"
parsed_args["forecast-type"] = "average"  # average-based RTCS forecast
###parsed_args["gamma-values"] = [0, 50, 100]  # [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
for instance_group in instance_group_list
	for model in ["deterministic", "robust-budget"]
		if model == "deterministic"
			parsed_args["gamma-values"] = [0, 50, 100]
		else
			parsed_args["gamma-values"] = [0, 20, 40, 60, 80, 100]  # [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
		end
		parsed_args["model"] = model
		parsed_args["instances"] = instance_group
		parsed_args["instance-name"] = "instance_deltamin10_spring.txt"
		parsed_args["simulation-nthreads"] = 8
		parsed_args["final-scenario"] = 1000
		parsed_args["model-policy"] = ["ignore_model"]
		parsed_args["sim-strategy"] = ["naive"]
		run_experiment(parsed_args)
	end
end
