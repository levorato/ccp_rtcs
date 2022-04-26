# ====================================================
# CLI Argument Parse Script for CCP
# ====================================================

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
		"--what"
			help = "What to do when running the program: solve or simulate."
			default = "solve"
        "--machine"
            help = "Machine name in the cluster (used when running experiments in the computer cluster)."
		"--instances"
            help = "Name of the instance group to solve."
		"--model"
            help = "Model used to solve the CCP: deterministic, robust-box or robust-budget."
		"--instance-name"
			help = "Filter to process only instance_name in instance_group (optional)."
			default = ""
		"--solver"
			help = "MILP solver to use : 'Gurobi' or 'CPLEX'"
			default = "CPLEX"
		"--solver-verbose"
			help="Enable solver output (used for debugging purposes)."
			arg_type = Bool
			default = false
		"--simulation-verbose"
			help="Enable simulation output (used for debugging purposes)."
			arg_type = Bool
			default = true
		"--simulation-nthreads"
			help="Number of simulation runs in parallel. Default = 4 threads."
			arg_type = Int
			default = 4
		"--time-limit"
			help = "Procedure time limit. Default = 14400 (seconds)."
			arg_type = Float64
			default = 14400.0
		"--solver-time-limit"
			help = "Solver time limit. Default = 1800 (seconds)."
			arg_type = Float64
			default = 1800.0
		"--max-cores"
			help = "Number of processor cores to use when solving the MILP problems. Default=16."
			arg_type = Int
			default = 4
		"--relative-gap-tol"
			help = "Relative gap tolerance used in MILP solving. Default = 0.02 (2 %)."
			arg_type = Float64
			default = 0.02
		"--rtcs-rounding-error"
			help = "Rounding error allowed in RTCS Operations (e.g. battery and drivable storage / retrieval). Default = 0.01."
			arg_type = Float64
			default = 0.01
		#"--num-runs"
		#	help = "Number of replications (used in simulation experiments). Default=100."
		#	arg_type = Int
		#	default = 100
		"--initial-scenario"
			help = "Number of the initial scenario to be simulated. Default = 1."
			arg_type = Int
			default = 1
		"--final-scenario"
			help = "Number of the last scenario to be simulated. Default = 100."
			arg_type = Int
			default = 100
		"--suffix"
			help = "Result filename suffix (used in parallel experiment executions)."
			default = ""
        "--forecast-type"
			help = "Type of RTCS forecast for uncertain consumption/production values. Default = average (average-based)."
			default = "average"
		"--resume"
			help="If true, resume simulations and reuse existing results saved to disk."
			arg_type = Bool
			default = true
		"--gamma-values"
			help="Gamma values for solving the deterministic or robust-budget CCP models."
			default = [0, 20, 30, 40, 50, 60, 80, 100]
		"--model-policy"
			help = "Model policy to be used in the RTCS simulation ('ignore-model' or 'full-model')."
			default = ["ignore_model", "full_model"]
		"--sim-strategy"
			help = "RTCS Strategy to be used in the simulation (['naive', 'conservative', 'audacious', 'cheapest'])."
			default = ["naive", "conservative", "audacious", "cheapest"]
    end
	parsed_args = parse_args(s)
	println("Parsed command-line args:")
	for (arg,val) in parsed_args
		println("  $arg  =>  $val")
	end
    return parsed_args
end

function get_machine_name()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    return parsed_args["machine"]
end

function setup_gurobi_license()
	println("Setup Gurobi License...")
	if MACHINE != "laptop"
		println("Cluster environment detected, setting up Gurobi license...")
	    machine_name = get_machine_name()
		ENV["GRB_LICENSE_FILE"] = joinpath(home_prefix, "$(machine_name)", "gurobi.lic")
		println("GRB_LICENSE_FILE = $(ENV["GRB_LICENSE_FILE"])")
	else
		println("Local environment detected, ignoring Gurobi license...")
	end
	flush(stdout)
end
