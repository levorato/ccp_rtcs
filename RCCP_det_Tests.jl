# =================================
# RCCP_det_Tests.jl
# CCP Deterministic Tests
# =================================

include("cli_arguments.jl")
include("RCCP_Facade.jl")

EXPERIMENT_NAME = "run_deterministic_ccp_tests"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

parsed_args = parse_commandline()
#parsed_args["instances"] = "Example_3.1"
parsed_args["instances"] = "antoine"
parsed_args["what"] = "solve"
parsed_args["model"] = "deterministic"
run_experiment(parsed_args)
