# =================================
# RCCP_det_Tests.jl
# RCCP Deterministic Tests
# =================================
include("cli_arguments.jl")
include("RCCP_Facade.jl")

EXPERIMENT_NAME = "run_example_3.1_tests"

parsed_args = parse_commandline()
parsed_args["instances"] = "Example_3.1"
parsed_args["what"] = "solve"
parsed_args["model"] = "robust-box"
run_experiment(parsed_args)
parsed_args["model"] = "deterministic"

parsed_args["output-folder"] = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))
run_experiment(parsed_args)
