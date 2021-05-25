include("cli_arguments.jl")
include("RCCP_Facade.jl")

EXPERIMENT_NAME = "test_sim_2"
# If possible, do not modify the lines below
parsed_args["output-folder"] = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

parsed_args = parse_commandline()
#parsed_args["instances"] = "antoine-skew"
#parsed_args["instances"] = "antoine_11"
parsed_args["instances"] = "toy"
# Solve models
#parsed_args["what"] = "solve"
#parsed_args["model"] = "robust-box"
#run_experiment(parsed_args)
#parsed_args["model"] = "deterministic"
#run_experiment(parsed_args)
# Simulate based on model results
parsed_args["what"] = "simulate"
parsed_args["model"] = "robust-box"
parsed_args["verbose"] = true
parsed_args["solver-verbose"] = true
run_experiment(parsed_args)
#parsed_args["model"] = "deterministic"
#run_experiment(parsed_args)
