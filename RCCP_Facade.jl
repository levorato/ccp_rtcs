include("config.jl")
include("cli_arguments.jl")
# Include RCCP Deterministic model
include("RCCP_det.jl")
include("RCCP_Robust.jl")
include("RCCP_Robust_Simulation.jl")
include("RCCP_SimUtil.jl")
include("RCCP_FileReader.jl")

using CSV
using DataFrames
using Dates
using Random
using UUIDs

function get_instance_list(instance_group, instance_name)
	basedir = ""
	instance_list = []
	println("Processing $(instance_group) instances...")
	if instance_group == "Example_3.1"
		base_folder = joinpath(project_folder, "instances")
		test_set = ["Example_3.1_A.txt", "Example_3.1_B.txt"]
		for instance_name in test_set
			inputfile = joinpath("$(base_folder)", "$(instance_name)")
			push!(instance_list, (instance_name, inputfile))
		end
	elseif instance_group == "antoine_11"
		base_folder = joinpath(project_folder, "instances")
		#datafile = "../notebooks/data/instance5-contratos_restituicao-ultimo_periodo.txt"
		#datafile = "../notebooks/data/antoine/A_instance2_11scen_1NDU.txt"
		#datafile = "../notebooks/data/antoine/A_instance2_11scen.txt"
		for filename in searchdir(base_folder, ".txt")
			inputfile = joinpath("$(base_folder)", "$(filename)")
			if occursin("A_instance2_11", filename)
				push!(instance_list, (filename, inputfile))
			end
		end
	elseif instance_group == "antoine_oc"
		base_folder = joinpath(project_folder, "instances")
		#datafile = "../notebooks/data/instance5-contratos_restituicao-ultimo_periodo.txt"
		#datafile = "../notebooks/data/antoine/A_instance2_11scen_1NDU.txt"
		#datafile = "../notebooks/data/antoine/A_instance2_11scen.txt"
		for filename in searchdir(base_folder, ".txt")
			inputfile = joinpath("$(base_folder)", "$(filename)")
			if occursin("A_instance2_OC", filename)
				push!(instance_list, (filename, inputfile))
			end
		end
	elseif instance_group == "antoine-skew"
		for filename in searchdir(antoine_instances_folder, ".txt")
			inputfile = joinpath("$(antoine_instances_folder)", "$(filename)")
			push!(instance_list, (filename, inputfile))
		end
	elseif instance_group == "toy"  # Taillard instances 20x5
		for filename in searchdir(toy_instances_folder, ".txt")
			inputfile = joinpath("$(toy_instances_folder)", "$(filename)")
			push!(instance_list, (filename, inputfile))
		end
	elseif occursin("japan", lowercase(instance_group))
        base_folder = joinpath(project_folder, "instances", "japan_microgrid", instance_group)
		for filename in searchdir(base_folder, ".txt")
			inputfile = joinpath("$(base_folder)", "$(filename)")
			push!(instance_list, (filename, inputfile))
		end
	elseif occursin("a_instance", lowercase(instance_group))
		inputfile = joinpath("$(antoine_instances_folder)", "$(instance_group)")
		push!(instance_list, (instance_group, inputfile))
	else
		println("WARN: No instances found!")
	end
	if length(instance_name) > 0
		println("Filtering instance_name = $(instance_name)")
		instance_list = [x for x in instance_list if x[1] == instance_name]
	end
	num_instances = size(instance_list, 1)
	println("# instances to process: $(num_instances)")
	# println("- Instances to process: $(instance_list)")
	return instance_list
	flush(stdout)
end

function run_experiment(parsed_args = nothing)
	setup_gurobi_license()
	if isnothing(parsed_args)
		parsed_args = parse_commandline()
	end
	if !haskey(parsed_args, "instances") || isnothing(parsed_args["instances"])
		println("ERROR: instances argument is mandatory!")
		return
	elseif !haskey(parsed_args, "model") || isnothing(parsed_args["model"])
		println("ERROR: model argument is mandatory!")
		return
	end
	instance_group = parsed_args["instances"]
	instance_name = parsed_args["instance-name"]
	model = parsed_args["model"]
	solver = parsed_args["solver"]
	solver_time_limit = parsed_args["solver-time-limit"]
	time_limit = parsed_args["time-limit"]
	max_cores = parsed_args["max-cores"]
	initial_scenario = parsed_args["initial-scenario"]
	final_scenario = parsed_args["final-scenario"]
	resume = parsed_args["resume"]
	model_policy_list = parsed_args["model-policy"]
	sim_strategy_list = parsed_args["sim-strategy"]
	#num_runs = parsed_args["num-runs"]
	println("Resume simulations? $(resume)")
	println("Number of available threads: $(Threads.nthreads())")
	println("Running CCP Experiment for instance_group = $(instance_group) and model = $(model)...")
	println("Model policy list: $(model_policy_list)")
	println("RTCS Strategy list: $(sim_strategy_list)")
	if model == "robust-budget" || model == "deterministic"
		println("[$(model)] Will process Gamma% = $(parsed_args["gamma-values"]).")
	end
	instance_list = get_instance_list(instance_group, instance_name)
	seed = rand(Int, 2)[1]
    rng = MersenneTwister(abs(seed))
	uid = uuid4(rng)
	batchId = Dates.format(Dates.now(), "yyyy_mm_dd-HH_MM_SS") * "-$(uid)"
	file_prefix = "CCP_model_$(model)_instance_$(instance_group)"
	parsed_args["base_output_path"] = normpath(parsed_args["output-folder"])

	if parsed_args["what"] == "solve"  # solve the model
		output_path = create_full_dir(parsed_args["output-folder"], ["ccp_model", "log"])
		model = parsed_args["model"]
		println("Request to solve CCP model $(model) on instance_group $(instance_group)")
		prev_inputfile = Nothing
		output_file_log = joinpath(normpath(output_path), "$(file_prefix).log")
		output_file_inf_log = joinpath(normpath(output_path), "$(file_prefix)_Infeasible.log")
		scenario_logger = open(output_file_log, "a+")
		for (filename, inputfile) in instance_list
			datetimenow = Dates.now()
			println("Date: $(datetimenow)")
			if inputfile != prev_inputfile
				println("***** Processing instance = $(filename) *****")
				instance_as_dict = read_tabulated_data(inputfile, instance_group, filename)
			end
			println("Calculating $(model) model solution (instance = $(filename))...")
			flush(stdout)
			if model == "robust-box"
				println("Solving for $(model)...")
				obtain_robust_optimization_model_results(parsed_args["base_output_path"], instance_group, model, filename,
					instance_as_dict, solver_time_limit, parsed_args, scenario_logger, [0])
			else
				for Gamma_perc in parsed_args["gamma-values"]
					println("Solving for Gamma_perc = $(Gamma_perc)...")
					obtain_robust_optimization_model_results(parsed_args["base_output_path"], instance_group, model, filename,
						instance_as_dict, solver_time_limit, parsed_args, scenario_logger, Gamma_perc)
				end
			end
			flush(stdout)
			prev_inputfile = inputfile
			flush(scenario_logger)
	    end
		close(scenario_logger)
	elseif parsed_args["what"] == "simulate"  # simulate RTCS using existing model results (from t = 1)
		prev_inputfile = Nothing
		#println("Request $(num_runs) simulation runs.")
		for (filename, inputfile) in instance_list
			datetimenow = Dates.now()
			println("Date: $(datetimenow)")
			if inputfile != prev_inputfile
				println("***** Processing instance = $(filename) *****")
				instance_as_dict = read_tabulated_data(inputfile, instance_group, filename)
			end
			println("Simulating RTCS for model $(model) and instance = $(filename))...")
			flush(stdout)
			simulate_instance(model, instance_group, filename, instance_as_dict, parsed_args;
								model_policy_list=model_policy_list, sim_strategy_list=sim_strategy_list)
			flush(stdout)
			prev_inputfile = inputfile
		end
	end
    println("DONE.\n")
	flush(stdout)
end
