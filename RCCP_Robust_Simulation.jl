# ===================================================================================
# RCCP_Robust_Simulation.jl
# For a given set of CCP model results, perform RTCS Simulation based on scenarios.
# ===================================================================================

using JuMP
using CSV
using DataFrames
using JLD2, FileIO
using Arrow
using Dates
import MathProgBase
using CodecZlib

# Include file reader and util functions
include("RCCP_FileReader.jl")
# Include RCCP Deterministic model
include("RCCP_det.jl")
# Include RCCP Robust model
include("RCCP_Robust.jl")
# Include RCCP Simulation utility functions
include("RCCP_SimUtil.jl")
# Include RCCP RTCS functions
include("RCCP_RTCS.jl")
# Include RCCP Plot functions
#include("RCCP_SimPlot.jl")

function process_individual_scenario(scenario_id, scenario, general_simulation_parameters, instance_simulation_parameters)
    forecast_type = general_simulation_parameters["forecast-type"]
    instance_as_dict = instance_simulation_parameters["instance_as_dict"]
    test_name = instance_simulation_parameters["test_name"]
    sim_strategy = instance_simulation_parameters["sim_strategy"]
    model_policy = instance_simulation_parameters["model_policy"]
    reoptimize = instance_simulation_parameters["reoptimize"]
    scenario_filter = instance_simulation_parameters["scenario_filter"]
    instance_name = instance_simulation_parameters["instance-name"]
    time_limit = instance_simulation_parameters["time_limit"]
    model = instance_simulation_parameters["model"]
    verbose = instance_simulation_parameters["verbose"]
    base_simulation_dir = instance_simulation_parameters["base_simulation_dir"]
    simulation_log_dir = instance_simulation_parameters["simulation-log-dir"]
    period_info = instance_as_dict["period"]
    instance = instance_as_dict["instance"]
    storage = instance_as_dict["storage"]
    drivable = instance_as_dict["drivable"]
    n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbNDU = size(n_drivable_uncertain, 1)
    T, S, C, D, ND, ST = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt)
    num_contracts = instance_as_dict["num_contracts"]
    scenario_list = instance_as_dict["scenarios"]
    nbScenarios = size(scenario_list, 1)
    tid = Threads.threadid()
    # Gamma parameter used by robust-budget model
    Gamma_perc = instance_simulation_parameters["Gamma_perc"]
    Gamma = (Gamma_perc * nbNDU * nbT) / 100.0
    # Open a logfile for writing
    simulation_log_dir = instance_simulation_parameters["simulation-log-dir"]
    output_file_log = joinpath(normpath(simulation_log_dir), model * "_" * test_name * "_" * instance_name * "_thread_$(tid).log")
    general_logger = open(output_file_log, "a+")
    println(general_logger, "[Thread $(tid)] Processing scenario #$(scenario_id), # Scenarios = $(size(scenario_filter, 1))")
    # Define a tracefile for writing
    simulation_trace_dir = get_simulation_trace_dir(instance_simulation_parameters)
    scenario_subpath = create_trace_scenario_filename(
            model, Gamma_perc, test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id)
    output_file_base = joinpath(normpath(simulation_trace_dir), scenario_subpath)
    # Define a logfile for writing scenario trace data (textual format / debug)
    scenario_file_log = output_file_base * ".log"
    scenario_logger = open(scenario_file_log, "w+")
    # Save other trace data (column format)
    output_file_df = output_file_base * ".csv"
    output_file_trace_arrow = output_file_base * ".arrow"
    output_file_var_arrow = output_file_base * "_var.arrow"
    # Create a new dataframe for this scenario
    trace_df, var_df = create_empty_trace_dataframe()
    y_model = nothing
    opt_df = DataFrame()
    sum_time_spent = 0
    sum_reopt_time_spent = 0
    total_processing_time = 0
    trace_df_list = Array[]
    # Run optimization models or obtain solution values if already stored in file
    opt_df_t1 = obtain_robust_optimization_model_results(base_simulation_dir, test_name, model, instance_name,
        instance_as_dict, time_limit, general_simulation_parameters, general_logger)
    if model == "robust-budget"  # Retrieve model result for a specific value of budget Gamma
        opt_df_t1 = opt_df_t1[(opt_df_t1[!, :GammaPerc] .== Gamma_perc), :]
        if nrow(opt_df_t1) == 0
            println(general_logger, "[Thread $(tid)] ERROR: No existing optimal model solution found for budget model Gamma% = $(Gamma_perc).")
            println(scenario_logger, "[Thread $(tid)] ERROR: No existing optimal model solution found for budget model Gamma% = $(Gamma_perc).")
            println("[Thread $(tid)] ERROR: No existing optimal model solution found for budget model Gamma% = $(Gamma_perc).")
        end
    elseif nrow(opt_df_t1) == 0
        println(general_logger, "[Thread $(tid)] ERROR: No existing optimal model solution found for the requested model.")
        println(scenario_logger, "[Thread $(tid)] ERROR: No existing optimal model solution found for the requested model.")
        println("[Thread $(tid)] ERROR: No existing optimal model solution found for the requested model.")
    end
    flush(stdout)
	flush(scenario_logger)
	flush(general_logger)
    try
        if verbose println("===========================================================================") end
        if verbose println(scenario_logger, "===========================================================================") end
        if verbose println("[Thread $(tid)] Scenario #$(scenario_id), # Scenarios = $(size(scenario_filter, 1))") end
        if verbose println(scenario_logger, "Scenario #$(scenario_id), # Scenarios = $(size(scenario_filter, 1))") end
        if verbose println(general_logger, "Scenario #$(scenario_id), # Scenarios = $(size(scenario_filter, 1))") end

        # Calculate the initial battery charge for each battery in ST
        initial_r_t = zeros(Float64, nbSt)
        for s in ST
            initial_r_t[s] = storage[s,:uInit]
        end
        previous_r_td_model = deepcopy(initial_r_t)
        acc_cost_model = 0.0
        previous_batt_levels_model = zeros(Float64, nbSt)
        # Generate the array of drivable devices previous accumulated charge value
        previous_drivable_charge_model = zeros(Float64, nbD)
        scenario_processing_start_time = time_ns()
        for t in 1:nbT   # For each time period t
            Gamma = (Gamma_perc * nbNDU * (nbT - t + 1)) / 100.0
            if verbose println(scenario_logger, "===========================================================================") end
            if verbose
                println(scenario_logger, "(t) Period #$(t) of $(nbT)")
                flush(stdout)
            end
            period_size = period_info[!, :size][t]
            #det_value = 0
            elapsed_time_model = 0.0
            if (!reoptimize) || (t == 1)  # Use model solution starting with t = 1
                if hasproperty(opt_df_t1, :ObjValue)
                    model_value_for_period = opt_df_t1[1, :ObjValue]
                    status_opt = opt_df_t1[1, :Optimal]
                    gap_opt = opt_df_t1[1, :RelGap]
                else
                    model_value_for_period = opt_df_t1[1, :RobValue]
                    status_opt = opt_df_t1[1, :RobOptimal]
                    gap_opt = 0.02
                end
                if t == 1
                    if hasproperty(opt_df_t1, :ModelTime)
                        model_opt_time_spent = opt_df_t1[1, :ModelTime]
                    else
                        model_opt_time_spent = opt_df_t1[1, :RobTime]
                    end
                else
                    model_opt_time_spent = 0.0
                end
                model_solution = read_variables_from_solution_df(model, opt_df_t1, 1, instance, num_contracts, nbNDU)
                # Store a copy of initial model results (t = 1) in this scenario dataframe
                opt_df = store_optimization_results_in_dataframe(model, instance_name, instance_as_dict, scenario_id, t, model_value_for_period,
                                                        model_solution, model_opt_time_spent, status_opt, gap_opt, general_logger, Gamma_perc, Gamma)
                # Fix the contracts from initial robust solution for t = 1 (y)
                y_model = model_solution["y"]
            else   # Use model solution starting with t = t'
                if model == "deterministic"
                    # Reoptimize the Deterministic RCCP Model
                    t1 = time_ns()  # ******************** measure start time
                    model_value_reopt, model_solution_reopt, solve_time_reopt, status_reopt, gap_reopt =
                            re_optimize_deterministic(instance_as_dict, time_limit, general_simulation_parameters, t, y_model,
                                previous_batt_levels_model, previous_drivable_charge_model, general_logger)
                    model_value_for_period = model_value_reopt
                    model_opt_time_spent = solve_time_reopt
                    model_solution = model_solution_reopt
                    elapsed_time_det = time_ns() - t1   # ****************** stop time
                else
                    # TODO Try to use another strategy for budget of uncertainty when re-optimizing the robust model
                    # e.g., considering the number of deviations observed so far.
                    # Reoptimize the Robust RCCP Model
                    t1 = time_ns()  # ******************** measure start time
                    model_value_reopt, model_solution_reopt, solve_time_reopt, status_reopt, gap_reopt =
                            re_optimize_robust(model, instance_as_dict, time_limit, general_simulation_parameters, t, y_model,
                            previous_batt_levels_model, previous_drivable_charge_model, general_logger, Gamma)
                    model_value_for_period = model_value_reopt
                    model_opt_time_spent = solve_time_reopt
                    model_solution = model_solution_reopt
                    elapsed_time_rob = time_ns() - t1  # ******************* stop time
                end
                if status_reopt == false && t > 1   # FIXME workaround for infeasible reopt model result
                    # Read solutions from time period t = 1
                    model_solution = read_variables_from_solution_df(model, opt_df_t1, 1, instance, num_contracts, nbNDU)
                    if hasproperty(opt_df_t1, :ObjValue)
                        model_value_reopt = opt_df_t1[1, :ObjValue]
                    else
                        model_value_reopt = opt_df_t1[1, :RobValue]
                    end
                    status_reopt = true
                    println("WARN: reusing rob solution from first period (t = 1), since current solution is INFEASIBLE.")
                    flush(stdout)
                end
                opt_df = vcat(opt_df, store_optimization_results_in_dataframe(model, instance_name, instance_as_dict, scenario_id, t, model_value_reopt,
                                                        model_solution, solve_time_reopt, status_reopt, gap_reopt, general_logger, Gamma_perc, Gamma))
                if hasproperty(opt_df, :ObjValue)
                    model_value_for_period = opt_df[t, :ObjValue]
                else
                    model_value_for_period = opt_df[t, :RobValue]
                end
                model_opt_time_spent = opt_df[t, :ModelTime]
                model_solution = read_variables_from_solution_df(model, opt_df, t, instance, num_contracts, nbNDU)
            end
            P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t
            q_t_model = zeros(Float64, num_contracts[t])

            # v2 : use P_hat_t = avg_P_hat_t for uncertain load values
            # Obtain average values for uncertain non-drivable devices load, such that P = (Pmin + Pmax)/2
            # TODO Instead of using average values, apply prediction for uncertain devices (Machine Learning)
            P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t
            if !isempty(n_drivable_uncertain)
                if occursin("ml", forecast_type)
                    # TODO Implement Machine Learning RTCS forecast, e.g., forecast_type == "ml-sliding-window"
                    for i in 1:size(n_drivable_uncertain, 1)
                        power = (n_drivable_uncertain[i, :Pmin][t] + n_drivable_uncertain[i, :Pmax][t]) / 2.0
                        P_hat_t[i] = power
                    end
                else
                    if forecast_type != "average"
                        println(general_logger, "\nWARN: Calculating average load values for uncertain non-drivable devices.")
                        println("[CCP Simulator] ERROR: Unknown forecast_type. Assuming average forecast type!")
                        flush(stdout)
                    end
                    for i in 1:size(n_drivable_uncertain, 1)
                        power = (n_drivable_uncertain[i, :Pmin][t] + n_drivable_uncertain[i, :Pmax][t]) / 2.0
                        P_hat_t[i] = power
                    end
                end
            end

            for d in 1:period_size  # For each d in period dt
                if verbose println(scenario_logger, "===========================================================================") end
                if verbose println(scenario_logger, "(d) Microperiod #$(d) of $(period_size) (t = $(t))") end
                t1 = time_ns()  # **************** measure start time
                P_hat_td = Float64[]
                for unds_matrix in scenario
                    push!(P_hat_td, unds_matrix[t][d])  # P_hat_td[s] for s in NDU
                end
                # v2 : use P_hat_t = avg_P_hat_t for uncertain load values
                #for s in 1:nbNDU  # Obtain the P_hat matrix with the sum for the whole period t
                #    P_hat_t[s] += P_hat_td[s]
                #end
                unds_id = 1
                time_spent = 0.0
                # run RTCS and obtain variable values for the current time unit d, period t
                if verbose println(scenario_logger, "\n  $(model)   R T C S") end
                t1 = time_ns()  # *********************** measure start time
                previous_q_t_model = deepcopy(q_t_model)
                q_td, x_td, e_td, r_td, g_td, h_td, gap, cost = rtcs(t, d, period_size, sim_strategy,
                                                            model_policy, instance_as_dict, P_hat_t, P_hat_td,
                                                            y_model, model, model_solution, previous_r_td_model,
                                                            previous_drivable_charge_model, previous_q_t_model,
                                                            scenario_logger, verbose)
                previous_r_td_model = deepcopy(r_td)
                for c in 1:num_contracts[t]
                    q_t_model[c] += q_td[c]
                end
                # For a strange reason, cost is being returned as a one-element array, fix this
                cost = fix_cost_variable_type(cost)
                acc_cost_model += cost
                # IMPORTANT : process drivable devices previous accumulated charge value
                for s in D
                    previous_drivable_charge_model[s] += (drivable[s,:pORc][t] / period_size) * x_td[s]
                end
                if verbose println(scenario_logger, "    ACCUMULATED COST                    : $(acc_cost_model)") end
                elapsed_time_model += (time_ns() - t1)  # ********************** stop time
                push!(trace_df, [model, Gamma_perc, Gamma, sim_strategy, reoptimize, model_policy, forecast_type, scenario_id,
                    t, d, model_value_for_period, (d == 1) ? model_opt_time_spent : 0.0,
                    e_td, gap, cost, elapsed_time_model/1.0e9])
                push!(var_df, [model, Gamma_perc, Gamma, sim_strategy, reoptimize, model_policy, forecast_type, scenario_id, t, d,
                    q_td, x_td, e_td, r_td, g_td, h_td])

                sum_time_spent += model_opt_time_spent
                # FIXME investigate reopt time calculation
                sum_reopt_time_spent += model_opt_time_spent
            end
            previous_batt_levels_model = deepcopy(previous_r_td_model)
        end
        scenario_end_processing_time = time_ns()
        scenario_processing_time = (scenario_end_processing_time - scenario_processing_start_time) / 1.0e9
        total_processing_time += scenario_processing_time
        println(scenario_logger, "\n===========================================================================")
        println(scenario_logger, " Total processing time : $(scenario_processing_time)")
        # TRACE FILE: Create an output file to write the dataframe to disk
        println(general_logger, "\n[Thread $(tid)] Saving trace CSV file to $(output_file_df)...")
        CSV.write(output_file_df, trace_df)
        Arrow.write(output_file_trace_arrow, trace_df; compress=:zstd, ntasks=2)
        Arrow.write(output_file_var_arrow, var_df; compress=:zstd, ntasks=2)
    catch e
        bt = catch_backtrace()
        msg = sprint(showerror, e, bt)
        println(msg)
        println(general_logger, "[Thread $(tid)] Exception thrown: \n" * msg)
    finally
        # Close the scenario log file
        close(scenario_logger)
        # Move the scenario log text file and scenario csv trace file to a zip file to save disk space
        simulation_zip_dir = get_simulation_zip_dir(general_simulation_parameters, instance_name)
        output_file_base = joinpath(normpath(simulation_zip_dir), create_trace_scenario_filename(
                model, Gamma_perc, test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id))
        output_file_zip = output_file_base * ".zip"
        move_files_to_zip_archive(output_file_zip, [scenario_file_log, output_file_df, output_file_var_arrow, output_file_trace_arrow], general_logger)
        println(general_logger, "[Thread $(tid)] Done processing scenario #$(scenario_id), # Scenarios = $(size(scenario_filter, 1)).")
        flush(general_logger)
        if verbose
            println("[Thread $(tid)] Done processing scenario #$(scenario_id), # Scenarios = $(size(scenario_filter, 1)).")
			flush(stdout)
        end
        flush(stdout)
        flush(general_logger)
    end
    close(general_logger)
    return opt_df
end

# Based on a list of scenarios to be simulated, identifies for which scenarios simulation can be skipped, for
# having already being executed.
function process_existing_scenario_runs(simulation_args, run_scenario_list, model, Gamma_perc, test_name, instance_name, sim_strategy,
        model_policy, reoptimize)
    # Define the path where simulation trace is written on disk
    simulation_log_dir = get_simulation_zip_dir(simulation_args, instance_name)
    println("Looking for existing simulations in folder $(simulation_log_dir)...")
    filtered_scenario_id_list = Int64[]
    reuse_scenario_id_list = Int64[]
    for (scenario_id, scenario) in run_scenario_list  # For each scenario
        scenario_subpath = create_trace_scenario_filename(
                model, Gamma_perc, test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id)
        output_file_base = joinpath(normpath(simulation_log_dir), scenario_subpath)
        output_file_trace_zip = output_file_base * ".zip"
	output_file_trace_arrow = scenario_subpath * ".arrow"
        # If both arrow files exist, skips the execution of this scenario
        if isfile(output_file_trace_zip)
	    found = false
            try
	    	r = ZipFile.Reader(output_file_trace_zip)
		for f in r.files
			if GetFileExtension(f.name) == ".arrow" && (f.name == output_file_trace_arrow)
				temp_df = DataFrame(Arrow.Table(f))
				if nrow(temp_df) > 0
					push!(reuse_scenario_id_list, scenario_id)
					found = true
				end
				break
			end
		end
            catch exc
	    	println("Error reading zip trace file: $(exc).")
            end
	    if !found
		push!(filtered_scenario_id_list, scenario_id)
	    end
        else
            push!(filtered_scenario_id_list, scenario_id)
        end
    end
    println("Reusing simulation for scenarios: $(reuse_scenario_id_list)")
	flush(stdout)
    if size(filtered_scenario_id_list, 1) > 0
        new_run_scenario_list = [x for x in run_scenario_list if (x[1] in filtered_scenario_id_list)]
        new_run_scenario_ids = [x[1] for x in run_scenario_list if (x[1] in filtered_scenario_id_list)]
        return new_run_scenario_list, new_run_scenario_ids
    else
        return [], []
    end
end

# Heurística de tempo real. Opcoes:
# a) Olhando apenas o y: Dentro de cada subperíodo de tempo, olhar apenas os contratos que serão utilizados (disponíveis para o
#    período) e as demandas / ofertas de energia dos dispositivos. Com base nisso, compra-se ou vende-se energia por meio e por fora
#    dos contratos. Não se olha para as variáveis que dependem da incerteza.
# b) Olhando o y e as demais variáveis do robusto atreladas ao período t' : usa informações dos contratos utilizados e dos
#    percentuais e quantidades usados da bateria e dos demais dispositivos (variáveis q, r, g, etc.), fornecidos pela multiplicação dos
#    valores das variáveis retornados pelo robusto vezes o valor realizado de P para o período (obtido nos dados do cenário: soma dos
#    valores de P de cada microperíodo). Com esses dados, roda-se a heurística para decidir o quanto usar de cada um.
# c) Olhando o y, as demais variáveis do robusto e o conhecimento das situações de instantes passados.
# Reoptimize: Opções de execução da simulação:
# - Sem re-otimização: roda-se o robusto apenas na primeira vez e a cada instante t' utiliza-se os valores das variáveis por ele
#   retornadas, multiplicando-se os valores pelas realizações do parâmetro incerto P_hat.
# - Com re-otimização: "Solve the Robust Model starting from t' "roda-se novamente o robusto a cada instante t' como se esse fosse
#   o primeiro período de tempo do problema, fixando as variáveis y. As informações correspondentes aos períodos que já passaram
#   são ignoradas.
# Main RTCS simulation function
function simulate(general_simulation_parameters, instance_simulation_parameters, general_logger)
    println("\nNew Simulation")
    instance_as_dict = instance_simulation_parameters["instance_as_dict"]
    test_name = instance_simulation_parameters["test_name"]
    sim_strategy = instance_simulation_parameters["sim_strategy"]
    model_policy = instance_simulation_parameters["model_policy"]
    reoptimize = instance_simulation_parameters["reoptimize"]
    scenario_filter = instance_simulation_parameters["scenario_filter"]
    instance_name = instance_simulation_parameters["instance-name"]
    base_simulation_dir = instance_simulation_parameters["base_simulation_dir"]
    model = instance_simulation_parameters["model"]
    simulation_nthreads = general_simulation_parameters["simulation-nthreads"]
    forecast_type = general_simulation_parameters["forecast-type"]
    resume = general_simulation_parameters["resume"]
    Gamma_perc = instance_simulation_parameters["Gamma_perc"]
    println("Resume simulations? $(resume)")
    println("Simulation strategy is $(sim_strategy)")
    println("Re-optimization is enabled? $(reoptimize)")
    println("RTCS forecast type is: $(forecast_type)")
    println("Number of simulation threads: $(simulation_nthreads)")
    # Obtain instance data
    println("Obtaining instance data...")
    flush(stdout)
    period_info = instance_as_dict["period"]
    instance = instance_as_dict["instance"]
    storage = instance_as_dict["storage"]
    drivable = instance_as_dict["drivable"]
    n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbNDU = size(n_drivable_uncertain, 1)
    T, S, C, D, ND, ST = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt)
    num_contracts = instance_as_dict["num_contracts"]
    scenario_list = instance_as_dict["scenarios"]
    nbScenarios = size(scenario_list, 1)
    if nbScenarios == 0
        println("No scenarios found! Aborting simulation.")
        println(general_logger, "No scenarios found! Aborting simulation.\n")
        flush(stdout)
        return [0.0, 0.0, zeros(Int64, nbC)]
    end
    if size(scenario_filter, 1) > 0
        println(general_logger, "Scenario filter is ENABLED. # Scenarios = $(size(scenario_filter, 1))")
        println("Scenario filter is ENABLED. # Scenarios = $(size(scenario_filter, 1))")
        flush(stdout)
    end

    scenario_start_processing_time = time_ns()
    # TODO Check for existing simulations of scenario_list (reuse existing results)
    run_scenario_list = []
    run_scenario_ids = []
    if size(scenario_filter, 1) > 0
        run_scenario_list = [x for x in enumerate(scenario_list) if (x[1] in scenario_filter)]
        run_scenario_ids = [x[1] for x in enumerate(scenario_list) if (x[1] in scenario_filter)]
    else
        run_scenario_list = [x for x in enumerate(scenario_list)]
        run_scenario_ids = [x[1] for x in enumerate(scenario_list)]
    end
    if resume
        println("Resume simulations is ENABLED. Reusing existing simulations.")
        run_scenario_list, run_scenario_ids = process_existing_scenario_runs(general_simulation_parameters, run_scenario_list, model, Gamma_perc, test_name,
            instance_name, sim_strategy, model_policy, reoptimize)
    else
        println("Resume simulations is DISABLED. Re-running all simulations from the beginning.")
    end
	flush(stdout)
    if size(run_scenario_list, 1) > 0
        df_vector = [DataFrame() for _ in 1:maximum(run_scenario_ids)]
        if simulation_nthreads > 1
            # Create a separated dict for each thread: https://discourse.julialang.org/t/thread-local-dict-for-each-thread/51441/3
            general_parametersn_dict_vector = [deepcopy(general_simulation_parameters) for _ in 1:Threads.nthreads()]
            instance_parameters_dict_vector = [deepcopy(instance_simulation_parameters) for _ in 1:Threads.nthreads()]
            Threads.@threads for (scenario_id, scenario) in collect(run_scenario_list)  # For each scenario in instance file
                tid = Threads.threadid()
                df_vector[scenario_id] = process_individual_scenario(scenario_id, scenario, general_parametersn_dict_vector[tid], instance_parameters_dict_vector[tid])
            end
        else  # no parallelization
            for (scenario_id, scenario) in run_scenario_list  # For each scenario in instance file
                df_vector[scenario_id] = process_individual_scenario(scenario_id, scenario, general_simulation_parameters, instance_simulation_parameters)
            end
        end
		# Concatenate all generated opt_results dataframes and save to disk => DISABLED
        ### opt_df = vcat(df_vector...)
        ### save_optimization_results_dataframe_to_file(base_simulation_dir, model, instance_name, opt_df, general_logger; reopt=true)
        scenario_end_processing_time = time_ns()
        total_processing_time = (scenario_end_processing_time - scenario_start_processing_time) / 1.0e9
        # Dados a serem armazenados nas simulações para análise futura (e.g. box plot):
        #    - Custo da solução (determinístico x robusto)
        #    - Quantidade de energia comprada fora do contrato
        #    - Nível na bateria em cada instante de tempo
        #    - Tempo gasto para rodar o PPL com a re-otimização.
        return total_processing_time
    else
        println("No scenarios to be simulated. All necessary results already exist on disk.")
		flush(stdout)
        return 0.0
    end
end

# Simulate the RTCS for a given microgrid 'instance_name', based on a given CCP 'model'.
function simulate_instance(model, test_name, instance_name, instance_as_dict, general_simulation_parameters;
        scenario_filter = Int64[], save_trace_info = false, model_policy_list = ["ignore_model", "full_model"], # "batteries", "batteries_and_drivable"
        sim_strategy_list = ["conservative", "audacious", "cheapest"], reoptimize_values = [false, true])

    df = DataFrame(Instance = String[], TestName = String[], ModelName = String[],
        Strategy = String[], Reoptimize = Bool[], TimeSpent = Float64[],
        ReOptTime = Float64[], Contracts = Array{Int64, 2}[])
    println("\n===================\n  CCP Simulator\n===================\n")
    time_limit = general_simulation_parameters["time-limit"]
    trace_df, var_df = create_empty_trace_dataframe()
    variable_values_after_uncertainty_df = create_empty_dataframe_variable_values_after_uncertainty()
    # Setup temporary dir for files
    temp_dir = mktempdir(tempdir(); prefix="RCCP_", cleanup=true)
    # Open a logfile for writing
    simulation_log_dir = get_simulation_log_dir(general_simulation_parameters, instance_name * "_" * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
    base_result_filename = model * "_" * test_name * "_" * instance_name
    output_file_log = joinpath(normpath(simulation_log_dir), base_result_filename * ".log")
    scenario_logger = open(output_file_log, "w+")
    simulation_nthreads = general_simulation_parameters["simulation-nthreads"]
    forecast_type = general_simulation_parameters["forecast-type"]
    try
        result_filepath = joinpath(normpath(simulation_log_dir), base_result_filename * ".csv")
        result_file_exists = isfile(result_filepath)
        result_file = open(result_filepath, "a+")
        println("\nSaving consolidated simulation results file to $(result_filepath)...")
        if !result_file_exists
            println(result_file, "model,Gamma_perc,test_name,instance_name,model_policy,sim_strategy,reopt,forecast_type,num_threads,time_spent")
        end
        println(scenario_logger, "===========================================================================")
        println("Simulation for instance $(instance_name)")
        println(scenario_logger, "Simulation for instance $(instance_name)")
        scenario_list = instance_as_dict["scenarios"]
        nbScenarios = size(scenario_list, 1)
        # Run RTCS Simulation
        total_p_time = 0.0
        gamma_list = [0]
        if model == "robust-budget"
            gamma_list = general_simulation_parameters["gamma-values"]
        end
        # Setup initial and final scenario numbers in simulation
        initial_scenario = general_simulation_parameters["initial-scenario"]
    	final_scenario = general_simulation_parameters["final-scenario"]
        println("Initial scenario: $(initial_scenario), Final scenario: $(final_scenario)")
        scenario_filter = [x for x in initial_scenario:final_scenario]
        for Gamma_perc in gamma_list
            for model_policy in model_policy_list
                for sim_strategy in sim_strategy_list
                    reopt_values = deepcopy(reoptimize_values)
                    if model_policy == "ignore_model" # 'ignore_model' uses only y variable info, so no reoptimization is needed
                        reopt_values = [false]
                    end
                    for reopt in reopt_values
                        println("\nSimulation for Gamma=$(Gamma_perc) x $(sim_strategy) x $(model_policy) x ReOpt=$(string(reopt))")
                        println(scenario_logger, "Simulation for Gamma=$(Gamma_perc) x $(sim_strategy) x $(model_policy) x ReOpt=$(string(reopt))")
						flush(stdout)
						flush(scenario_logger)
                        instance_simulation_parameters = Dict("sim_strategy"=>sim_strategy, "reoptimize"=>reopt, "model_policy"=>model_policy,
                            "scenario_filter"=>scenario_filter, "time_limit"=>time_limit, "model"=>model, "instance_as_dict"=>instance_as_dict,
                            "test_name"=>test_name, "instance-name"=>instance_name, "base_simulation_dir"=>normpath(general_simulation_parameters["base_output_path"]),
                            "verbose"=>general_simulation_parameters["simulation-verbose"], "Gamma_perc"=>Gamma_perc, "temp-dir"=>temp_dir,
                            "scenario_filter"=>scenario_filter, "simulation-log-dir"=>simulation_log_dir, "output-folder"=>general_simulation_parameters["output-folder"])
                        p_time = simulate(general_simulation_parameters, instance_simulation_parameters, scenario_logger)
                        println(scenario_logger, "Done. Processing time : $(string(p_time))")
                        println("Done. Processing time : $(string(p_time))")
                        flush(scenario_logger)
                        total_p_time += p_time
                        # Export simulation summary to CSV file
                        println(result_file, "$model,$Gamma_perc,$test_name,$instance_name,$model_policy,$sim_strategy,$reopt,$forecast_type,$(simulation_nthreads),$p_time")
                        flush(result_file)
                    end
                end
            end
        end
        println(scenario_logger, "Simulation done. Processing time : $(string(total_p_time))\n")
        close(result_file)
    	concatenate_trace_df_list(general_simulation_parameters, model, test_name, instance_name, scenario_logger)
    	println("\n[CCP Simulator] DONE.")
        println(scenario_logger, "\n[CCP Simulator] DONE.")
    catch e
        bt = catch_backtrace()
        msg = sprint(showerror, e, bt)
        println(msg)
        println(scenario_logger, "[CCP Simulator] Exception thrown: \n" * msg)
    finally
        flush(stdout)
        flush(scenario_logger)
        # Close the log file
        close(scenario_logger)
    end
end

function fix_cost_variable_type(cost)
    if isa(cost, Array)
        if size(cost, 1) == 1
            return cost[1]
        else
            println("ERROR: Problem converting variable cost !")
            return [0.0]
        end
    else
        return cost
    end
end
