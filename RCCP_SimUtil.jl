# =================================
# RCCP_SimUtil.jl
# RTCS Simulation utility functions.
# =================================

using JLD2
using Arrow
using CSV
using CodecZlib

function get_infeasible_log_path(solver_parameters)
	tid = Threads.threadid()
	return joinpath(normpath(solver_parameters["base_output_path"]), "Infeasible_RCCP_thread_$(tid).log")
end

# Re-otimização a cada período t' :
#   A cada período t', re-otimiza-se, i.e., roda-se novamente o robusto, partindo do instante t'
#   (com as variáveis inteiras y fixas com o valor da solução inicial do robusto), que se transforma
#   num PPL. Ou seja, os parâmetros incertos (P_hat) relativos aos instantes de tempo passados (t < t')
#   passam a ser fixos, assumindo os valores contidos no cenário, uma vez que já se revelaram na simulação.
function re_optimize_robust(model, instance_as_dict, time_limit, solver_parameters, t, y,
        initial_battery, previous_drivable_charge, general_logger, Gamma = 0.0)
    # Solve the Robust Model starting from t', fixing the contracts variable y
    output_file_inf_log = get_infeasible_log_path(solver_parameters)
    infeasible_logger = open(output_file_inf_log, "a+")
    solver_parameters_2 = deepcopy(solver_parameters)
    solver_parameters_2["relative-gap-tol"] = 1e-04  # Set relative gap tolerance to default solver value
	# Since simulation runs are executed in parallel (multithreaded), we can
	# allow only 2 reoptimization threads per simulation run
	solver_parameters_2["max-cores"] = 2
	if model == "robust-box"
		println("Reoptimizing robust box model.")
    	rob_value, rob_solution, solve_time, status, gap = solve_robust_model_box(instance_as_dict, general_logger,
            infeasible_logger, solver_parameters_2, t, y, initial_battery, previous_drivable_charge)
	else  # robust-budget
		println("Reoptimizing robust budget model for Gamma = $(Gamma).")
        rob_value, rob_solution, solve_time, status, gap = solve_robust_model_budget(instance_as_dict, general_logger,
			infeasible_logger, solver_parameters_2, Gamma, t, y, initial_battery, previous_drivable_charge)
	end
    close(infeasible_logger)
    return rob_value, rob_solution, solve_time, status, gap
end

function re_optimize_deterministic(instance_as_dict, time_limit, solver_parameters, t, y,
        initial_battery, previous_drivable_charge, general_logger)
    # Solve the Deterministic Model starting from t', fixing the contracts variable y
    output_file_inf_log = get_infeasible_log_path(solver_parameters)
    infeasible_logger = open(output_file_inf_log, "a+")
    solver_parameters_2 = deepcopy(solver_parameters)
    solver_parameters_2["relative-gap-tol"] = 1e-04  # Set relative gap tolerance to default solver value
	# Since simulation runs are executed in parallel (multithreaded), we can
	# allow only 2 reoptimization threads per simulation run
	solver_parameters_2["max-cores"] = 2
    det_value, det_solution, det_solve_time, det_status, det_gap = solve_deterministic_model_with_t(instance_as_dict,
            general_logger, infeasible_logger, solver_parameters_2, t, y, initial_battery, previous_drivable_charge)
    close(infeasible_logger)
    println("Resolvido model det para t = $(t), opt = $(det_status)")
    return det_value, det_solution, det_solve_time, det_status, det_gap
end

function get_list_of_missing_results(instance_name, result_df, gamma_values)
	# TODO implement existing result recovery
	return []
end

function read_model_results_from_csv(result_filename)
	println("Reading result file: $(result_filename)")
	if isfile(result_filename)
	    println("Result file exists.")
	    df = DataFrame(CSV.File(result_filename, delim=','))
		println("First 6 results: ")
		println(first(df, 6))
		flush(stdout)
	    return df
	else
		println("No existing result file found!")
		flush(stdout)
		return DataFrame()
	end
end

function create_empty_opt_dataframe(model)
	if model == "robust-box"  # solve the robust model under budgeted uncertainty
		return DataFrame(Instance = String[], ScenarioId = Int64[], PeriodId = Int64[],
			Model = String[],
			ObjValue = Float64[], Solution = Array{Int64, 1}[], ModelTime = Float64[],
			Optimal = Bool[], q0 = Array{Float64, 1}[], q_ = Array{Float64, 1}[],
			x0 = Array{Float64, 1}[], x_ = Array{Float64, 1}[], e0 = Array{Float64, 1}[],
			e_ = Array{Float64, 1}[], r0 = Array{Float64, 1}[], r_ = Array{Float64, 1}[],
			g0 = Array{Float64, 1}[], g_ = Array{Float64, 1}[], RelGap = Float64[])
	elseif model == "robust-budget"  # solve the robust model under budgeted uncertainty
        return DataFrame(Instance = String[], ScenarioId = Int64[], PeriodId = Int64[],
            Model = String[], GammaPerc = Float64[], Gamma = Float64[],
            ObjValue = Float64[], Solution = Array{Int64, 1}[], ModelTime = Float64[],
            Optimal = Bool[], q0 = Array{Float64, 1}[], q_ = Array{Float64, 1}[],
            x0 = Array{Float64, 1}[], x_ = Array{Float64, 1}[], e0 = Array{Float64, 1}[],
            e_ = Array{Float64, 1}[], r0 = Array{Float64, 1}[], r_ = Array{Float64, 1}[],
            g0 = Array{Float64, 1}[], g_ = Array{Float64, 1}[], RelGap = Float64[])
	else  #if model == "deterministic"  # solve the deterministic model
        return DataFrame(Instance = String[], ScenarioId = Int64[], PeriodId = Int64[], Model = String[], ObjValue = Float64[],
            Solution = Array{Int64, 1}[], ModelTime = Float64[], Optimal = Bool[],
            q = Array{Float64, 1}[], x = Array{Float64, 1}[], e = Array{Float64, 1}[],
            r = Array{Float64, 1}[], g = Array{Float64, 1}[], h = Array{Float64, 1}[], RelGap = Float64[])
	end
end

# Optimize / run the robust RCCP model (deterministic and robust) for all time intervals between 1 and nbT
# (if reoptimize param is true) or only for t = 1 (otherwise).
# Returns a dataframe with all results obtained.
function run_robust_optimization_models(model, instance_name, instance_as_dict, solver_parameters, general_logger, infeasible_logger)
    # Open output CSV file for results
    file_prefix = "CCP_model_$(model)"
	tid = Threads.threadid()
    full_filename = "$(file_prefix)_thread_$(tid).csv"
	opt_output_path = solver_parameters["base_output_path"]
    result_file_path = joinpath(normpath(opt_output_path), full_filename)
    result_file_exists = isfile(result_file_path)
    result_df = read_model_results_from_csv(result_file_path)
    result_file = open(result_file_path, "a+")
    println("\nSaving to results file to $(result_file_path)...")
    if !result_file_exists
        println(result_file, "model,instance_name,model_parameters,solution_value,time_spent,is_optimal,solution_y,rel_gap")
    end
	if model == "robust-box"  # solve the robust model under budgeted uncertainty
        # Solve the robust model starting from t = 1
        rob_value_t1, rob_solution_t1, rob_solve_time_t1, rob_status_t1, gap_t1 = solve_robust_model_box(instance_as_dict, general_logger, infeasible_logger, parsed_args, 1)
        y = rob_solution_t1["y"]
        println("Solution: y = $(y)")
        println(result_file, "$(model),$(instance_name),none,$(rob_value_t1),$(rob_solve_time_t1),$(rob_status_t1),$(y),$(gap_t1)")
        flush(result_file)
        close(result_file)
        # Save data for t = 1
		opt_df = store_optimization_results_in_dataframe(model, instance_name, instance_as_dict, -1, 1, rob_value_t1, rob_solution_t1,
		                                        rob_solve_time_t1, rob_status_t1, gap_t1, general_logger)
        save_optimization_results_dataframe_to_file(opt_output_path, "robust-box", instance_name, opt_df, general_logger)
        return opt_df
    elseif model == "robust-budget"  # solve the robust model under budgeted uncertainty
        gamma_values = solver_parameters["gamma-values"]
	    run_list = get_list_of_missing_results(instance_name, result_df, gamma_values)
        println("Processing experiment for gamma_values = $(gamma_values).")
        nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance_as_dict["instance"])
        nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
        println("Max number of allowed deviations is $(nbNDU * nbT).")
		df_list = []
        for Gamma_perc in gamma_values
            Gamma = (Gamma_perc * nbNDU * nbT) / 100.0
            println("Gamma% = $(Gamma_perc) => Gamma = $(Gamma)")
            flush(stdout)
            # Solve the robust model starting from t = 1
            rob_value_t1, rob_solution_t1, rob_solve_time_t1, rob_status_t1, gap_t1 = solve_robust_model_budget(instance_as_dict,
				general_logger, infeasible_logger, parsed_args, Gamma, 1)
            y = rob_solution_t1["y"]
            println("Solution: y = $(y)")
            println(result_file, "$(model),$(instance_name),$(Gamma),$(rob_value_t1),$(rob_solve_time_t1),$(rob_status_t1),$(y),$(gap_t1)")
            flush(result_file)
            # Save data for t = 1
			opt_df = store_optimization_results_in_dataframe(model, instance_name, instance_as_dict, -1, 1, rob_value_t1, rob_solution_t1,
			                                        rob_solve_time_t1, rob_status_t1, gap_t1, general_logger, Gamma_perc, Gamma)
			push!(df_list, opt_df)
			save_optimization_results_dataframe_to_file(opt_output_path, "robust-budget", instance_name, opt_df, general_logger)
        end
        close(result_file)
        return vcat(df_list...)
    else  #if model == "deterministic"  # solve the deterministic model
        # Solve the robust model only from t = 1 (no reoptimization)
        # Solve the Deterministic model starting from t = 1
        det_value_t1, det_solution_t1, det_solve_time_t1, det_status_t1, gap_t1 = solve_deterministic_model_with_t(instance_as_dict, general_logger,
                    infeasible_logger, parsed_args, 1)
        y = det_solution_t1["y"]
        println("Solution: y = $(y)")
        println(result_file, "$(model),$(instance_name),none,$(det_value_t1),$(det_solve_time_t1),$(det_status_t1),$(y),$(gap_t1)")
        flush(result_file)
        close(result_file)
        # Save data for t = 1
		opt_df = store_optimization_results_in_dataframe(model, instance_name, instance_as_dict, -1, 1, det_value_t1, det_solution_t1,
		                                        det_solve_time_t1, det_status_t1, gap_t1, general_logger)
        save_optimization_results_dataframe_to_file(opt_output_path, "deterministic", instance_name, opt_df, general_logger)
        return opt_df
    end
end

function flatten_solution(x)
	return collect(Iterators.flatten(x))
end

function store_optimization_results_in_dataframe(model, instance_name, instance_as_dict, scenario_id, period_id, obj_value, model_solution,
                                        model_solve_time, status_opt, rel_gap, general_logger, Gamma_perc = 0.0, Gamma = 0.0)
    instance = instance_as_dict["instance"]
    num_contracts = instance_as_dict["num_contracts"]
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    SetSigma = 1:nbNDU
    T, S, C, D, ND, ST = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt)
	opt_df = create_empty_opt_dataframe(model)
    # only store results if solution if optimal
    if haskey(model_solution, "y") && (rel_gap != Inf)
		if model == "deterministic"  # Instance;ScenarioId;PeriodId;Model;ObjValue;Solution;ModelTime;Optimal;q;x;e;r;g;h;RelGap
	        push!(opt_df, [instance_name, scenario_id, period_id, model,
				obj_value, flatten_solution(model_solution["y"]),
	            model_solve_time, status_opt,
	            flatten_solution(model_solution["q"]), flatten_solution(model_solution["x"]), flatten_solution(model_solution["e"]),
	            flatten_solution(model_solution["r"]), flatten_solution(model_solution["g"]), flatten_solution(model_solution["h"]), rel_gap])
		elseif model == "robust-budget"
			push!(opt_df, [instance_name, scenario_id, period_id, model, Gamma_perc, Gamma,
				obj_value, flatten_solution(model_solution["y"]),
	            model_solve_time, status_opt,
	            flatten_solution(model_solution["q0"]), flatten_solution(model_solution["q_"]),
	            flatten_solution(model_solution["x0"]), flatten_solution(model_solution["x_"]),
	            flatten_solution(model_solution["e0"]), flatten_solution(model_solution["e_"]),
	            flatten_solution(model_solution["r0"]), flatten_solution(model_solution["r_"]),
	            flatten_solution(model_solution["g0"]), flatten_solution(model_solution["g_"]), rel_gap])
		else
			push!(opt_df, [instance_name, scenario_id, period_id, model,
            	obj_value, flatten_solution(model_solution["y"]),
	            model_solve_time, status_opt,
	            flatten_solution(model_solution["q0"]), flatten_solution(model_solution["q_"]),
	            flatten_solution(model_solution["x0"]), flatten_solution(model_solution["x_"]),
	            flatten_solution(model_solution["e0"]), flatten_solution(model_solution["e_"]),
	            flatten_solution(model_solution["r0"]), flatten_solution(model_solution["r_"]),
	            flatten_solution(model_solution["g0"]), flatten_solution(model_solution["g_"]), rel_gap])
		end
    else
        println("ERROR : CCP model $(model) did not obtain optimal solution. Returning empty solution !")
        println(general_logger, "ERROR : CCP model $(model) did not obtain optimal solution. Returning empty solution !")
        println("Model status = $(status_opt)")
        println(general_logger, "Model status = $(status_opt)")
        #error("ERROR : one of the RCCP models did not obtain optimal solution.")
        max_c = 0
        for t in 1:nbT
            max_c = max(max_c, num_contracts[t])
        end
        sol_y = zeros(Int64, nbT, max_c)
		if model == "deterministic"
	        push!(opt_df, [instance_name, scenario_id, period_id, model, obj_value, flatten_solution(sol_y),
	            model_solve_time, status_opt,
				flatten_solution([Float64(0.0) for x=1:nbT, y=1:max_c]), flatten_solution([Float64(0.0) for x=1:nbT, y=D]), flatten_solution([Float64(0.0) for x=1:nbT]),
	            flatten_solution([Float64(0.0) for x=1:nbT+1, y=ST]), flatten_solution([Float64(0.0) for x=1:nbT, y=ST]), flatten_solution([Float64(0.0) for x=1:nbT, y=ST]), Inf])
		elseif model == "robust-budget"
			push!(opt_df, [instance_name, scenario_id, period_id, model, Gamma_perc, Gamma,
				obj_value, flatten_solution(sol_y),
				model_solve_time, status_opt,
	            flatten_solution([Float64(0.0) for x=1:nbT, y=1:max_c]), flatten_solution([Float64(0.0) for x=1:nbT, y=1:max_c, z=1:nbT, k=SetSigma]),
	            flatten_solution([Float64(0.0) for x=1:nbT, y=D]), flatten_solution([Float64(0.0) for x=1:nbT, y=D, z=1:nbT, k=SetSigma]),
	            flatten_solution([Float64(0.0) for x=1:nbT]),
	            flatten_solution([Float64(0.0) for x=1:nbT, y=1:nbT, z=SetSigma]), flatten_solution([Float64(0.0) for x=1:nbT+1, y=ST]),
	            flatten_solution([Float64(0.0) for x=1:nbT+1, y=ST, z=1:nbT, k=SetSigma]),
	            flatten_solution([Float64(0.0) for x=1:nbT, y=ST]), flatten_solution([Float64(0.0) for x=1:nbT, y=ST, z=1:nbT, k=SetSigma]), Inf])
		else
	        push!(opt_df, [instance_name, scenario_id, period_id, model, obj_value, flatten_solution(sol_y),
				model_solve_time, status_opt,
	            flatten_solution([Float64(0.0) for x=1:nbT, y=1:max_c]), flatten_solution([Float64(0.0) for x=1:nbT, y=1:max_c, z=1:nbT, k=SetSigma]),
	            flatten_solution([Float64(0.0) for x=1:nbT, y=D]), flatten_solution([Float64(0.0) for x=1:nbT, y=D, z=1:nbT, k=SetSigma]),
	            flatten_solution([Float64(0.0) for x=1:nbT]),
	            flatten_solution([Float64(0.0) for x=1:nbT, y=1:nbT, z=SetSigma]), flatten_solution([Float64(0.0) for x=1:nbT+1, y=ST]),
	            flatten_solution([Float64(0.0) for x=1:nbT+1, y=ST, z=1:nbT, k=SetSigma]),
	            flatten_solution([Float64(0.0) for x=1:nbT, y=ST]), flatten_solution([Float64(0.0) for x=1:nbT, y=ST, z=1:nbT, k=SetSigma]), Inf])
		end
    end
	return opt_df
end

function read_variables_from_solution_dict(solution_dict)
    q0 = solution_dict["q0"]
    q_ = solution_dict["q_"]
    r0 = solution_dict["r0"]
    r_ = solution_dict["r_"]
    g0 = solution_dict["g0"]
    g_ = solution_dict["g_"]
    e0 = solution_dict["e0"]
    e_ = solution_dict["e_"]
    x0 = solution_dict["x0"]
    x_ = solution_dict["x_"]
    return q0, q_, x0, x_, e0, e_, r0, r_, g0, g_
end

function deflatten_array(y, dims...)
	return reshape(y, dims)
end

function read_variables_from_solution_df(model, solution_df, period_id, instance, num_contracts, nbNDU)
	nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
	max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
	if hasproperty(solution_df, :ObjValue)
		if model == "deterministic"
		    det_solution = Dict()
			det_solution["y"] = deflatten_array(solution_df[period_id, :Solution], nbT, max_contracts_per_period)
		    det_solution["q"] = deflatten_array(solution_df[period_id, :q], nbT, max_contracts_per_period)
		    det_solution["r"] = deflatten_array(solution_df[period_id, :r], nbT+1, nbSt)
		    det_solution["g"] = deflatten_array(solution_df[period_id, :g], nbT, nbSt)
		    det_solution["h"] = deflatten_array(solution_df[period_id, :h], nbT, nbSt)
		    det_solution["e"] = deflatten_array(solution_df[period_id, :e], nbT)
		    det_solution["x"] = deflatten_array(solution_df[period_id, :x], nbT, nbD)
			return det_solution
		else  # Robust model with Linear Decision Rules
		    rob_solution = Dict()
			rob_solution["y"] = deflatten_array(solution_df[period_id, :Solution], nbT, max_contracts_per_period)
		    rob_solution["q0"] = deflatten_array(solution_df[period_id, :q0], nbT, max_contracts_per_period)
		    rob_solution["q_"] = deflatten_array(solution_df[period_id, :q_], nbT, max_contracts_per_period, nbT, nbNDU)
		    rob_solution["x0"] = deflatten_array(solution_df[period_id, :x0], nbT, nbD)
		    rob_solution["x_"] = deflatten_array(solution_df[period_id, :x_], nbT, nbD, nbT, nbNDU)
		    rob_solution["e0"] = deflatten_array(solution_df[period_id, :e0], nbT)
		    rob_solution["e_"] = deflatten_array(solution_df[period_id, :e_], nbT, nbT, nbNDU)
		    rob_solution["r0"] = deflatten_array(solution_df[period_id, :r0], nbT+1, nbSt)
		    rob_solution["r_"] = deflatten_array(solution_df[period_id, :r_], nbT+1, nbSt, nbT+1, nbNDU)
		    rob_solution["g0"] = deflatten_array(solution_df[period_id, :g0], nbT, nbSt)
		    rob_solution["g_"] = deflatten_array(solution_df[period_id, :g_], nbT, nbSt, nbT, nbNDU)
		    return rob_solution
		end
	else  # old dataframe version
		if model == "deterministic"
			det_solution = Dict()
			det_solution["y"] = solution_df[period_id, :DetSolution]
			det_solution["q"] = solution_df[period_id, :q]
			det_solution["r"] = solution_df[period_id, :r]
			det_solution["g"] = solution_df[period_id, :g]
			det_solution["h"] = solution_df[period_id, :h]
			det_solution["e"] = solution_df[period_id, :e]
			det_solution["x"] = solution_df[period_id, :x]
			return det_solution
		else  # Robust model with Linear Decision Rules
			rob_solution = Dict()
			rob_solution["y"] = solution_df[period_id, :RobSolution]
			rob_solution["q0"] = solution_df[period_id, :q0]
			rob_solution["q_"] = solution_df[period_id, :q_]
			rob_solution["x0"] = solution_df[period_id, :x0]
			rob_solution["x_"] = solution_df[period_id, :x_]
			rob_solution["e0"] = solution_df[period_id, :e0]
			rob_solution["e_"] = solution_df[period_id, :e_]
			rob_solution["r0"] = solution_df[period_id, :r0]
			rob_solution["r_"] = solution_df[period_id, :r_]
			rob_solution["g0"] = solution_df[period_id, :g0]
			rob_solution["g_"] = solution_df[period_id, :g_]
			return rob_solution
		end
	end
end

function create_empty_trace_dataframe()
    return DataFrame(Model = String[], GammaPerc = Float64[], Gamma = Float64[],
        Strategy = String[], Reoptimize = Bool[], ModelPolicy = String[],
		ForecastType = String[],
        ScenarioId = Int64[], t = Int64[], d = Int64[], ObjValue = Float64[],
        OptTimeSpent = Float64[],
        e_td = Float64[], gap = Float64[],
        cost = Float64[], RealProcTime = Float64[]),
            DataFrame(Model = String[], GammaPerc = Float64[], Gamma = Float64[],
                Strategy = String[], Reoptimize = Bool[], ModelPolicy = String[],
				ForecastType = String[],
                ScenarioId = Int64[], t = Int64[], d = Int64[],
                # values calculated by RTCS heuristic for period (t, d)
                q_td = Array{Float64, 1}[],
                x_td = Array{Float64, 1}[], e_td = Float64[], r_td = Array{Float64, 1}[],
                g_td = Array{Float64, 1}[], h_td = Array{Float64, 1}[])
end

function create_empty_dataframe_variable_values_after_uncertainty()
    return DataFrame(ModelType = String[],
        Reoptimize = Bool[], ScenarioId = Int64[], t = Int64[], d = Int64[],
        ObjValue = Float64[], OptTimeSpent = Float64[], d0_Method = String[],
        # values calculated by deterministic variables or robust LDRs for time period t
        q_t = Array{Float64, 1}[],
        x_t = Array{Float64, 1}[], e_t = Float64[], r_t = Array{Float64, 1}[],
        g_t = Array{Float64, 1}[], h_t = Array{Float64, 1}[],
        gap = Float64[], cost = Float64[], acc_cost_t = Float64[])
end

function calculate_pure_model_variables_after_uncertainty_is_known(instance_filename, test_name,
                model, instance_name, reoptimize, time_limit, solver_parameters, general_logger)
    variable_values_after_uncertainty_df = create_empty_dataframe_variable_values_after_uncertainty()
    opt_output_path = create_full_dir(normpath(solver_parameters["base_output_path"]))
    # Antoine tabulated format
    println("Read instance file...")
    instance_as_dict = read_tabulated_data(instance_filename)
    period_info = instance_as_dict["period"]
    instance = instance_as_dict["instance"]
    drivable = instance_as_dict["drivable"]
    n_drivable = instance_as_dict["n_drivable"]
    n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
    storage = instance_as_dict["storage"]
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    NDU = 1:nbNDU
    SetSigma = NDU
    T, S, C, D, ND, ST = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt)
    num_contracts = instance_as_dict["num_contracts"]
    scenario_list = instance_as_dict["scenarios"]
    contract = instance_as_dict["contract"]
    period = instance_as_dict["period"]
    nbScenarios = size(scenario_list, 1)
    if nbScenarios == 0
        println("No scenarios found! Aborting.")
        return variable_values_after_uncertainty_df
    end
    lambda = zeros(Float64, nbSt)
    for s in ST
        lambda[s] = storage[s,:lostCoef]
    end
    # Run optimization models or obtain solution values if already stored in file
	# FIXME TODO Retornar apenas as variaveis puras do modelo desejado, talvez filtrando pelo gamma do dataframe
    opt_df = obtain_robust_optimization_model_results(opt_output_path, test_name, model, instance_name,
                                                        instance_as_dict, time_limit, solver_parameters, general_logger)
    # Fix the contracts from initial robust solution for t = 1 (y)
    y_model = opt_df[1, :Solution]
    @assert opt_df[1, :ObjValue] > 0

    scenario_id = 1
    for scenario in scenario_list  # For each scenario in instance file
        #println("Scenario #$(scenario_id) of $(nbScenarios)")
        # Calculate the initial battery charge for each battery in ST
        previous_r_t = zeros(Float64, nbSt)
        for s in ST
            previous_r_t[s] = storage[s,:uInit]
        end
        for t in 1:nbT   # For each time period t
            period_size = period_info[!, :size][t]
            if (!reoptimize) || (t == 1)  # Use model solution starting with t = 1
                obj_value_for_period = opt_df[1, :ObjValue]
                model_opt_time_spent = opt_df[1, :ModelTime]
                model_solution = read_variables_from_solution_df(model, opt_df, 1, instance, num_contracts, nbNDU)
            else   # Use model solution starting with t = t'
                obj_value_for_period = opt_df[t, :ObjValue]
                model_opt_time_spent = opt_df[t, :ModelTime]
                model_solution = read_variables_from_solution_df(model, opt_df, t, instance, num_contracts, nbNDU)
            end
            gap_model = 0.0
            P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t
            for d in 1:period_size
                P_hat_td = Float64[]
                for unds_matrix in scenario
                    push!(P_hat_td, unds_matrix[t][d])  # P_hat_td[s] for s in NDU
                end
                gap_model += sum(P_hat_td[s] for s in 1:nbNDU)  # Ŝ_ND : uncertain non-drivable power (comes from scenario data)
                for s in 1:nbNDU
                    P_hat_t[s] += P_hat_td[s]
                end
            end
			if model == "deterministic"
	            q_t, x_t, e_t, r_t, g_t, h_t = calculate_variables_from_deterministic(instance_as_dict, model_solution, t)
	            # sum(Ŝ_ND) is already stored in variable gap_det (see above)
	            gap_model += sum(n_drivable[s,:pORc][t] for s in ND)  # sum(S_ND)
	            gap_model += sum(q_t[c] for c in 1:num_contracts[t])
	            gap_model += sum(drivable[s,:pORc][t] * x_t[s] for s in D)
	            gap_model += sum(h_t[s] for s in ST) - sum(g_t[s] for s in ST)
	            gap_model += e_t
	            cost = sum( contract[contract[!, :period] .== t, :cost_fix][c] * y_model[t,c] for c in 1:num_contracts[t])
	            cost += sum(contract[contract[!, :period] .== t, :cost_var][c] * q_t[c] for c in 1:num_contracts[t])
	            cost += sum(drivable[s,:cost] * drivable[s,:pORc][t] * x_t[s] for s in D)
	            cost += sum(storage[s,:cost] * (g_t + h_t) for s in ST)
	            cost += (period[t,:cost_out] * e_t)
	            cost += sum(n_drivable_uncertain[sigma, :cost] * P_hat_t[sigma] for sigma in SetSigma)
	            if isa(cost, Array)  # Gambiarra para bug no tipo de dados
	                cost = cost[1]
	            end
	            push!(variable_values_after_uncertainty_df, ["Deterministic", reoptimize, scenario_id, t, 1, obj_value_for_period, model_opt_time_spent,
	                "N/A", q_t, x_t, e_t, r_t, g_t, h_t, gap_model, cost, cost])
			else  # Robust models
	            # Obtain q_t'(.), x_t'(.), e_t'(.), r_t'(.), g_t'(.), based on model solution
	            q_t, x_t, e_t, r_t, g_t, h_t = calculate_variables_from_robust_ldr(instance_as_dict, model_solution, t, P_hat_t, period_size, y_model)
	            delta_r_t = zeros(Float64, nbSt)
	            initial_r_t_d = zeros(Float64, nbSt)
	            for s in ST  # Calculate the difference between current battery charge and charge from the previous period t
	                delta_r_t[s] = r_t[s] - previous_r_t[s]
	                initial_r_t_d[s] = previous_r_t[s]
	            end
	            previous_r_t = deepcopy(r_t)
	            for calc_method in ["average", "all_beginning"]
	                acc_cost = 0.0
	                r_td = zeros(Float64, nbSt)
	                for s in ST
	                    r_td[s] = initial_r_t_d[s]
	                end
	                for d in 1:period_size  # For each d in period dt
	                    P_hat_td = Float64[]
	                    for unds_matrix in scenario
	                        push!(P_hat_td, unds_matrix[t][d])  # P_hat_td[s] for s in NDU
	                    end
	                    previous_r_td_rob = deepcopy(r_td)
	                    for s in ST
	                        r_td[s] += (delta_r_t[s] / period_size)
	                    end
	                    q_td, x_td, e_td, g_td = calculate_model_variable_delta(q_t, x_t, e_t, g_t, period_size, d, calc_method)
	                    #for s in ST   TODO calcular h[s] para o robusto
	                    #    r_td[s] = r_td[s] - h_td[s] + storage[s,:lostCoef] * g_td[s]
	                    #end
	                    h_td = zeros(Float64, nbSt)
	                    gap_model = sum(P_hat_td[s] for s in NDU)  # Ŝ_ND : uncertain non-drivable power (comes from scenario data)
	                    gap_model += sum(n_drivable[s,:pORc][t] / period_size for s in ND)   # S_ND  FIXME Confirmar se eh pra dividir por delta
	                    gap_model += sum(q_td[c] for c in 1:num_contracts[t])
	                    gap_model += sum((drivable[s,:pORc][t] / period_size) * x_td[s] for s in D)
	                    gap_model += sum(((lambda[s]-1)*g_td[s] + r_td[s] - previous_r_td_rob[s]) for s in ST)  #  + sum(h_td[s] for s in ST) - sum(g_td[s] for s in ST)
	                    gap_model += e_td
	                    #println("e_td = $(e_td)")
	                    # only sum contract fixed cost if d == 1 (first micro period)
	                    cost = d == 1 ? sum( contract[contract[!, :period] .== t, :cost_fix][c] * y[t,c] for c in 1:num_contracts[t]) : 0.0
	                    cost += sum(contract[contract[!, :period] .== t, :cost_var][c] * q_td[c] for c in 1:num_contracts[t])
	                    cost += sum(drivable[s,:cost] * (drivable[s,:pORc][t] / period_size) * x_td[s] for s in D)  # sum(S_D)
	                    cost += sum(storage[s,:cost]*(1+lambda[s]) * g_td[s] for s in ST)
	                    cost += sum(n_drivable_uncertain[sigma, :cost] * P_hat_td[sigma] for sigma in SetSigma)  # sum(Ŝ_ND)
	                        # FIXME Consertar o codigo abaixo de soma dos custos da bateria no robusto (substituicao de h_td[s])
	                    cost += sum(storage[s,:cost] * ( g_td[s] + (lambda[s]*g_td[s] + previous_r_td_rob[s] - r_td[s]) ) for s in ST)
	                        # sum(storage[s,:cost] * (g_td[s] + h_td[s]) for s in ST) +
	                    cost += (period[t,:cost_out] * e_td)
	                    acc_cost += cost
	                    push!(variable_values_after_uncertainty_df, ["Robust", reoptimize, scenario_id, t, d, model_value_for_period, model_opt_time_spent,
	                        calc_method, q_td, x_td, e_td, r_td, g_td, h_td, gap_model, cost, acc_cost])
	                end
	            end
			end
        end
        scenario_id += 1
    end
    return variable_values_after_uncertainty_df
end

# Calculate the value of each model variable according to microperiod d
function calculate_model_variable_delta(q, x, e, g, period_size, d, calc_method = "all_beginning")
    if calc_method == "average"  # Divide term (e, q, etc.) by the period size
        return deepcopy(q) ./ period_size, deepcopy(x) ./ period_size, deepcopy(e) ./ period_size, deepcopy(g) ./ period_size
    else  # "full value in the first microperiod d"
        if d == 1
            return deepcopy(q), deepcopy(x), deepcopy(e), deepcopy(g)
        else
            return q .* 0.0, x .* 0.0, e .* 0.0, g .* 0.0
        end
    end
end

function get_optimization_result_base_filename(opt_output_path, model, instance_name)
    return joinpath(normpath(opt_output_path), "RCCP_Sim_OptData_$(model)_$(instance_name)")
end

function GetFileExtension(filename)
    return filename[findlast(isequal('.'),filename):end]
end

function obtain_robust_optimization_model_results(opt_output_path, test_name, model, instance_name,
		instance_as_dict, time_limit, solver_parameters, general_logger)
    base_filename = get_optimization_result_base_filename(opt_output_path, model, instance_name)
	tid = Threads.threadid()
    opt_var_file = base_filename * ".jld"
    opt_csv_file = base_filename * ".csv"
    opt_file_zip = base_filename * ".zip"
    opt_file_arrow = base_filename * ".arrow"
    tmp_file_path = joinpath(tempdir(), base_filename * "_opt_df_thread_$(tid).jld")
    #println(general_logger, "Trying to obtain robust model results from:")
	#println(general_logger, " - Option 1: $(opt_var_file)")
	#println(general_logger, " - Option 2: $(opt_file_zip)")
	#println(general_logger, " - Option 3: $(opt_file_arrow)")
	opt_df = nothing
	l = Threads.ReentrantLock()
	Threads.lock(l)
	n = 0
	try
	    if isfile(opt_file_arrow)
		try
	                println("Skipping optimization runs. Arrow File already exists: $(opt_file_arrow)")
        	        println(general_logger, "Skipping optimization runs. Arrow File already exists: $(opt_file_arrow)")
                	opt_df = DataFrame(Arrow.Table(opt_file_arrow))
	                println(general_logger, "Optimization data read successfully from arrow file.\n")
        	        flush(general_logger)
			n = nrow(opt_df)
		catch exc1
			println(general_logger, "Error reading arrow opt_df file: $(exc1)")
			showerror(general_logger, exc1, catch_backtrace())
		end
	    end
	    if (n == 0) && isfile(opt_file_zip)  # v2.0 feature
	        try
	        	println("Skipping optimization runs. ZIP File already exists: $(opt_file_zip)")
		        println(general_logger, "Skipping optimization runs. ZIP File already exists: $(opt_file_zip)")
		        r = ZipFile.Reader(opt_file_zip)
	        	for f in r.files
	      		      if GetFileExtension(f.name) == ".jld"
	                	println("Reading file from zip: $(f.name)")
		                println(general_logger, "Reading file from zip: $(f.name)")
		                s = read(f, String)
	        	        tmp_file_df = open(tmp_file_path, "w+")
	                	write(tmp_file_df, s)
		                close(tmp_file_df)
		                # Read execution data back to dataframe
	        	        @load tmp_file_path opt_var
	                	opt_df = opt_var["sol"]
		                break
		              end
	    	        end
	        	close(r)
	        	println(general_logger, "Optimization data read successfully from zip file.\n")
			n = nrow(opt_df)
		catch exc2
			println(general_logger, "Error reading JLD/ZIP opt_df file: $(exc2)")
                        showerror(general_logger, exc2, catch_backtrace())
		end
	    end
	    if (n == 0) && isfile(opt_var_file)   # Check if optimization results are cached in file
                println("Skipping optimization runs. File already exists: $(opt_var_file)")
                println(general_logger, "Skipping optimization runs. File already exists: $(opt_var_file)")
                # Read execution data back to dataframe
                #opt_var = load(opt_var_file)
                @load opt_var_file opt_var
                key_list = collect(keys(opt_var))
                opt_df = opt_var[key_list[1]]
                if isa(opt_df, Dict)
                        key_list = collect(keys(opt_df))
                        opt_df = opt_df[key_list[1]]
                end
                println(general_logger, "Optimization data read successfully.\n")
                flush(general_logger)
		n = nrow(opt_df)
	    end
	catch exc
		println(general_logger, "ERROR Unable to read existing model solution. Cause: $(exc)")
		showerror(general_logger, exc, catch_backtrace())
	finally
		if isnothing(opt_df) || (n == 0)
			println(general_logger, "Existing robust model results NOT found. Invoking solver...")
	        	output_file_inf_log = get_infeasible_log_path(solver_parameters)
	        	infeasible_logger = open(output_file_inf_log, "a+")
	        	opt_df = run_robust_optimization_models(model, instance_name, instance_as_dict, solver_parameters, general_logger, infeasible_logger)
	        	close(infeasible_logger)
	    	end
		Threads.unlock(l)
		flush(general_logger)
	end
    return opt_df
end

function save_optimization_results_dataframe_to_file(opt_output_path, model, instance_name, opt_df, general_logger; reopt=false)
    base_filename = get_optimization_result_base_filename(opt_output_path, model, instance_name)
	if reopt
		base_filename *= "_reopt"
	end
    opt_var_file = base_filename * ".arrow"
    opt_csv_file = base_filename * ".csv.gz"
    #println("\nSaving optimization results data file to $(opt_csv_file)...")
    println(general_logger, "\nSaving new optimization results to disk...")
    try
		# Get existing results from existing arrow file
		existing_df = nothing
		if isfile(opt_var_file)
			try
				println(general_logger, "Reading existing optimization runs data. Arrow File already exists: $(opt_var_file)")
				io = open(opt_var_file, "r")
				existing_df = DataFrame(Arrow.Table(io))
				close(io)
				println(general_logger, "Existing optimization data read successfully from arrow file.\n")
			catch exc
				println(general_logger, "ERROR Unable to read existing optimization runs datafile. Cause: $(exc)")
			end
		else
			println(general_logger, "Existing optimization runs data not found on path: $(opt_var_file)")
		end
		# Append new results to existing results
		if !isnothing(existing_df)
			opt_df = vcat(existing_df, opt_df)
		end
		# Write final results
		open(GzipCompressorStream, opt_csv_file, "w") do stream
		   CSV.write(stream, opt_df; delim=';')
	    end
        # Serialize dataframe to file
        Arrow.write(opt_var_file, opt_df; compress=:zstd, ntasks=2)
        println(general_logger, "\nSaved optimization results data files to $(opt_csv_file) and $(opt_var_file).")
        println(general_logger, "Done.")
    catch y2
        println("ERROR Writing optimization output files $(opt_csv_file) and $(opt_var_file). Cause : $(y2).")
        println(general_logger, "ERROR Writing optimization output files $(opt_csv_file) and $(opt_var_file). Cause : $(y2).")
    end
	flush(general_logger)
	flush(stdout)
end
