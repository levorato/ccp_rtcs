# =============================================================================
# solve_models.jl
# Solve all CCP models at time period t = 1: deterministic, robust-box and
# robust-budget. Stores model results in experiment output-folder.
# =============================================================================

include("cli_arguments.jl")
include("RCCP_Facade.jl")

EXPERIMENT_NAME = "run_sim_japan_forecast_avg_2"

parsed_args = parse_commandline()
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))
parsed_args["output-folder"] = EXPERIMENT_OUTPUT_FOLDER

parsed_args["what"] = "solve"
parsed_args["solver-verbose"] = true
parsed_args["solver-time-limit"] = 28800.0  # 8h time-limit
parsed_args["max-cores"] = 8
parsed_args["gamma-values"] = [0, 20, 40, 50, 60, 80, 100]

# Solve models
models_to_solve = ["robust-budget"]  # ["deterministic"]  # , "robust-box", "robust-budget"]
instance_group_list = ["japan-10"]

for instance_group in instance_group_list
	for model in models_to_solve
		parsed_args["model"] = model
		parsed_args["instances"] = instance_group
		parsed_args["instance-name"] = "instance_deltamin10_summer.txt"
		parsed_args["relative-gap-tol"] = 1e-5
		run_experiment(parsed_args)
	end
end

# Simulate based on model results
#parsed_args["what"] = "simulate"
#parsed_args["model"] = "robust-box"
#run_experiment(parsed_args)
#parsed_args["model"] = "deterministic"
#run_experiment(parsed_args)
