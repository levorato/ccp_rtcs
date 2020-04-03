# =================================
# RCCP_SimUtil.jl
# RTCS Simulation utility functions.
# =================================

# Re-otimização a cada período t' :
#   A cada período t', re-otimiza-se, i.e., roda-se novamente o robusto, partindo do instante t'
#   (com as variáveis inteiras y fixas com o valor da solução inicial do robusto), que se transforma
#   num PPL. Ou seja, os parâmetros incertos (P_hat) relativos aos instantes de tempo passados (t < t')
#   passam a ser fixos, assumindo os valores contidos no cenário, uma vez que já se revelaram na simulação.
function re_optimize_robust(instance_as_dict, time_limit, cplex_old, t, y,
        initial_battery, previous_drivable_charge, general_logger)
    # Solve the Robust Model starting from t', fixing the contracts variable y
    output_file_inf_log = joinpath(normpath(pwd()), "Infeasible_RCCP.log")
    infeasible_logger = open(output_file_inf_log, "a+")
    rob_value, rob_solution, solve_time, status = solve_robust_model(instance_as_dict, general_logger,
            infeasible_logger, time_limit, cplex_old, t, -100.0, y, initial_battery, previous_drivable_charge)
    close(infeasible_logger)
    return rob_value, rob_solution, solve_time, status
end

function re_optimize_deterministic(instance_as_dict, time_limit, cplex_old, t, y,
        initial_battery, previous_drivable_charge, general_logger)
    # Solve the Deterministic Model starting from t', fixing the contracts variable y
    output_file_inf_log = joinpath(normpath(pwd()), "Infeasible_RCCP.log")
    infeasible_logger = open(output_file_inf_log, "a+")
    det_value, det_solution, det_solve_time, det_status = solve_deterministic_model_with_t(instance_as_dict,
            general_logger, infeasible_logger, time_limit, cplex_old, t, y, initial_battery, previous_drivable_charge)
    close(infeasible_logger)
    println("Resolvido model det para t = $(t), opt = $(det_status)")
    return det_value, det_solution, det_solve_time, det_status
end

# Optimize / run the robust RCCP model (deterministic and robust) for all time intervals between 1 and nbT
# (if reoptimize param is true) or only for t = 1 (otherwise).
# Returns a dataframe with all results obtained.
function run_robust_optimization_models(instance_name, instance_as_dict, time_limit, cplex_old, general_logger, infeasible_logger)
    opt_df = DataFrame(Instance = String[], ScenarioId = Int64[], PeriodId = Int64[], DetValue = Float64[],
        DetSolution = Array{Int64, 2}[], DetTime = Float64[], DetOptimal = Bool[],
        q = Array{Float64, 2}[], x = Array{Float64, 2}[], e = Array{Float64, 1}[],
        r = Array{Float64, 2}[], g = Array{Float64, 2}[], h = Array{Float64, 2}[],
        RobValue = Float64[], RobSolution = Array{Int64, 2}[], RobTime = Float64[],
        RobOptimal = Bool[], q0 = Array{Float64, 2}[], q_ = Array{Float64, 4}[],
        x0 = Array{Float64, 2}[], x_ = Array{Float64, 4}[], e0 = Array{Float64, 1}[],
        e_ = Array{Float64, 3}[], r0 = Array{Float64, 2}[], r_ = Array{Float64, 4}[],
        g0 = Array{Float64, 2}[], g_ = Array{Float64, 4}[])

    # Solve the robust model only from t = 1 (no reoptimization)
    # Solve the Deterministic model starting from t = 1
    det_value_t1, det_solution_t1, det_solve_time_t1, det_status_t1 = solve_deterministic_model_with_t(instance_as_dict,
            general_logger, infeasible_logger, time_limit, cplex_old)
    # Solve the robust model starting from t = 1
    rob_value_t1, rob_solution_t1, rob_solve_time_t1, rob_status_t1 = solve_robust_model(instance_as_dict, general_logger, infeasible_logger, time_limit, cplex_old, 1)
    # Fix the contracts: os valores de y (contratos utilizados) são SEMPRE FIXOS, obtidos pela solução inicial do MIP.
    y_det = det_solution_t1["y"]
    y_rob = rob_solution_t1["y"]
    println("xxx Robust solution value: $(rob_value_t1) ; $(y_rob)")
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance_as_dict["instance"])
    num_contracts = instance_as_dict["num_contracts"]
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    D = 1:nbD
    ND = 1:nbND
    ST = 1:nbSt
    SetSigma = 1:nbNDU

    # Save data for t = 1
    push!(opt_df, [instance_name, -1, 1, det_value_t1, det_solution_t1["y"], det_solve_time_t1, det_status_t1,
        det_solution_t1["q"], det_solution_t1["x"], det_solution_t1["e"],
        det_solution_t1["r"], det_solution_t1["g"], det_solution_t1["h"],
        rob_value_t1, rob_solution_t1["y"], rob_solve_time_t1, rob_status_t1,
        rob_solution_t1["q0"], rob_solution_t1["q_"],
        rob_solution_t1["x0"], rob_solution_t1["x_"],
        rob_solution_t1["e0"], rob_solution_t1["e_"],
        rob_solution_t1["r0"], rob_solution_t1["r_"],
        rob_solution_t1["g0"], rob_solution_t1["g_"]])

    # NOTE: Reoptimization of models is now done inside the simulation procedure, since it depends on the value of
    #       battery levels and, therefore, it depends on the simulation results of the previous period
    reoptimize = false
    if reoptimize  # Solve the Robust Model starting from t' (period_id)
        for period_id in 2:nbT   # For each time period t = {2, ..., nbT}
            # Reoptimize the Deterministic RCCP Model
            det_value_reopt, det_solution_reopt, det_solve_time_reopt, status_det_reopt =
                    re_optimize_deterministic(instance_as_dict, time_limit, cplex_old, period_id, y_det, previous_batt_levels,
                        general_logger)
            # Reoptimize the Robust RCCP Model
            rob_value_reopt, rob_solution_reopt, rob_solve_time_reopt, status_rob_reopt =
                    re_optimize_robust(instance_as_dict, time_limit, cplex_old, period_id, y_rob, previous_batt_levels, general_logger)
            store_optimization_results_in_dataframe(instance_name, instance_as_dict, opt_df, -1, period_id, det_value_reopt,
                                                    det_solution_reopt, det_solve_time_reopt, status_det_reopt,
                                                    rob_value_reopt, rob_solution_reopt, rob_solve_time_reopt, status_rob_reopt,
                                                    general_logger)
        end
    end
    return opt_df
end

function store_optimization_results_in_dataframe(instance_name, instance_as_dict, opt_df, scenario_id, period_id, det_value_reopt, det_solution_reopt,
                                        det_solve_time_reopt, status_det_reopt, rob_value_reopt, rob_solution_reopt, rob_solve_time_reopt,
                                        status_rob_reopt, general_logger)
    instance = instance_as_dict["instance"]
    num_contracts = instance_as_dict["num_contracts"]
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    SetSigma = 1:nbNDU
    T, S, C, D, ND, ST = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt)
    if status_rob_reopt == true && status_det_reopt == true  # only store results if solution if optimal
        push!(opt_df, [instance_name, scenario_id, period_id, det_value_reopt, det_solution_reopt["y"],
            det_solve_time_reopt, status_det_reopt,
            det_solution_reopt["q"], det_solution_reopt["x"], det_solution_reopt["e"],
            det_solution_reopt["r"], det_solution_reopt["g"], det_solution_reopt["h"],
            rob_value_reopt,
            rob_solution_reopt["y"], rob_solve_time_reopt, status_rob_reopt,
            rob_solution_reopt["q0"], rob_solution_reopt["q_"],
            rob_solution_reopt["x0"], rob_solution_reopt["x_"],
            rob_solution_reopt["e0"], rob_solution_reopt["e_"],
            rob_solution_reopt["r0"], rob_solution_reopt["r_"],
            rob_solution_reopt["g0"], rob_solution_reopt["g_"]])
    else
        println("ERROR : one of the RCCP models did not obtain optimal solution. Returning empty solution !")
        println(general_logger, "ERROR : one of the RCCP models did not obtain optimal solution. Returning empty solution !")
        println("Det status = $(status_det_reopt) ; Rob status = $(status_rob_reopt)")
        println(general_logger, "Det status = $(status_det_reopt) ; Rob status = $(status_rob_reopt)")
        #error("ERROR : one of the RCCP models did not obtain optimal solution.")
        max_c = 0
        for t in 1:nbT
            max_c = max(max_c, num_contracts[t])
        end
        push!(opt_df, [instance_name, scenario_id, period_id, det_value_reopt, det_solution_reopt["y"], det_solve_time_reopt, status_det_reopt,
            [Float64(0.0) for x=1:nbT, y=1:max_c], [Float64(0.0) for x=1:nbT, y=D], [Float64(0.0) for x=1:nbT],
            [Float64(0.0) for x=1:nbT+1, y=ST], [Float64(0.0) for x=1:nbT, y=ST], [Float64(0.0) for x=1:nbT, y=ST],
            rob_value_reopt, rob_solution_reopt["y"], rob_solve_time_reopt, status_rob_reopt,
            [Float64(0.0) for x=1:nbT, y=1:max_c], [Float64(0.0) for x=1:nbT, y=1:max_c, z=1:nbT, k=SetSigma],
            [Float64(0.0) for x=1:nbT, y=D], [Float64(0.0) for x=1:nbT, y=D, z=1:nbT, k=SetSigma],
            [Float64(0.0) for x=1:nbT],
            [Float64(0.0) for x=1:nbT, y=1:nbT, z=SetSigma], [Float64(0.0) for x=1:nbT, y=ST],
            [Float64(0.0) for x=1:nbT, y=ST, z=1:nbT, k=SetSigma],
            [Float64(0.0) for x=1:nbT, y=ST], [Float64(0.0) for x=1:nbT, y=ST, z=1:nbT, k=SetSigma]])
    end
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

function read_variables_from_solution_df(solution_df, period_id)
    det_solution = Dict()
    det_solution["q"] = solution_df[period_id, :q]
    det_solution["r"] = solution_df[period_id, :r]
    det_solution["g"] = solution_df[period_id, :g]
    det_solution["h"] = solution_df[period_id, :h]
    det_solution["e"] = solution_df[period_id, :e]
    det_solution["x"] = solution_df[period_id, :x]
    rob_solution = Dict()
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
    return det_solution, rob_solution
end

function create_empty_trace_dataframe()
    return DataFrame(RTCS_Type = String[],
        Strategy = String[], Reoptimize = Bool[], ModelPolicy = String[],
        ScenarioId = Int64[], t = Int64[], d = Int64[], DetValue = Float64[],
        RobValue = Float64[], OptTimeSpent = Float64[],
        e_td = Float64[], gap = Float64[],
        cost = Float64[], RealProcTime = Float64[]),
            DataFrame(RTCS_Type = String[],
                Strategy = String[], Reoptimize = Bool[], ModelPolicy = String[],
                ScenarioId = Int64[], t = Int64[], d = Int64[],
                # values calculated by RTCS heuristic for period (t, d)
                q_td = Array{Float64, 1}[],
                x_td = Array{Float64, 1}[], e_td = Float64[], r_td = Array{Float64, 1}[],
                g_td = Array{Float64, 1}[], h_td = Array{Float64, 1}[])
end

# Read a trace dataframe stored in JLD file
function read_trace_dataframe(test_name, instance_name)
    output_path = joinpath(normpath(pwd()), "output", "simulation", "trace")
    output_file = joinpath(normpath(output_path), test_name *
                    "_RCCP_Sim_TRACE_" * instance_name)
    println("\nReading trace file from $(output_file)...")
    # Read execution data back to dataframe
    trace_var = load(output_file * ".jld")
    trace_df = trace_var["trace"]
    println("\nTrace data read successfully.\n")
    return trace_df
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

function calculate_pure_model_variables_after_uncertainty_is_known(instance_filename,
                test_name, instance_name, reoptimize, time_limit, cplex_old, general_logger)
    variable_values_after_uncertainty_df = create_empty_dataframe_variable_values_after_uncertainty()
    output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "trace"])
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
    opt_df = obtain_robust_optimization_model_results(output_path, test_name, instance_name,
                                                        instance_as_dict, time_limit, cplex_old, general_logger)
    # Fix the contracts from initial robust solution for t = 1 (y)
    y = opt_df[1, :RobSolution]
    y_det = opt_df[1, :DetSolution]
    @assert opt_df[1, :DetValue] > 0
    @assert opt_df[1, :RobValue] > 0

    scenario_id = 0
    for scenario in scenario_list  # For each scenario in instance file
        #println("Scenario #$(scenario_id) of $(nbScenarios)")
        # Calculate the initial battery charge for each battery in ST
        previous_r_t = zeros(Float64, nbSt)
        for s in ST
            previous_r_t[s] = storage[s,:uInit]
        end
        for t in 1:nbT   # For each time period t
            period_size = period_info[:size][t]
            if (!reoptimize) || (t == 1)  # Use model solution starting with t = 1
                det_value_for_period = opt_df[1, :DetValue]
                rob_value_for_period = opt_df[1, :RobValue]
                det_opt_time_spent = opt_df[1, :DetTime]
                rob_opt_time_spent = opt_df[1, :RobTime]
                det_solution, rob_solution = read_variables_from_solution_df(opt_df, 1)
            else   # Use model solution starting with t = t'
                det_value_for_period = opt_df[t, :DetValue]
                rob_value_for_period = opt_df[t, :RobValue]
                det_opt_time_spent = opt_df[t, :DetTime]
                rob_opt_time_spent = opt_df[t, :RobTime]
                det_solution, rob_solution = read_variables_from_solution_df(opt_df, t)
            end
            gap_det = 0.0
            P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t
            for d in 1:period_size
                P_hat_td = Float64[]
                for unds_matrix in scenario
                    push!(P_hat_td, unds_matrix[t][d])  # P_hat_td[s] for s in NDU
                end
                gap_det += sum(P_hat_td[s] for s in 1:nbNDU)  # Ŝ_ND : uncertain non-drivable power (comes from scenario data)
                for s in 1:nbNDU
                    P_hat_t[s] += P_hat_td[s]
                end
            end
            # "Robust"
            # Obtain q_t'(.), x_t'(.), e_t'(.), r_t'(.), g_t'(.), based on model solution
            q_t, x_t, e_t, r_t, g_t, h_t = calculate_variables_from_robust_ldr(instance_as_dict, rob_solution, t, P_hat_t, period_size, y)
            delta_r_t = zeros(Float64, nbSt)
            initial_r_t_d = zeros(Float64, nbSt)
            for s in ST  # Calculate the difference between current battery charge and charge from the previous period t
                delta_r_t[s] = r_t[s] - previous_r_t[s]
                initial_r_t_d[s] = previous_r_t[s]
            end
            previous_r_t = copy(r_t)
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
                    previous_r_td_rob = copy(r_td)
                    for s in ST
                        r_td[s] += (delta_r_t[s] / period_size)
                    end
                    q_td, x_td, e_td, g_td = calculate_model_variable_delta(q_t, x_t, e_t, g_t, period_size, d, calc_method)
                    #for s in ST   TODO calcular h[s] para o robusto
                    #    r_td[s] = r_td[s] - h_td[s] + storage[s,:lostCoef] * g_td[s]
                    #end
                    h_td = zeros(Float64, nbSt)
                    gap_rob = sum(P_hat_td[s] for s in NDU)  # Ŝ_ND : uncertain non-drivable power (comes from scenario data)
                    gap_rob += sum(n_drivable[s,:pORc][t] / period_size for s in ND)   # S_ND  FIXME Confirmar se eh pra dividir por delta
                    gap_rob += sum(q_td[c] for c in 1:num_contracts[t])
                    gap_rob += sum((drivable[s,:pORc][t] / period_size) * x_td[s] for s in D)
                    gap_rob += sum(((lambda[s]-1)*g_td[s] + r_td[s] - previous_r_td_rob[s]) for s in ST)  #  + sum(h_td[s] for s in ST) - sum(g_td[s] for s in ST)
                    gap_rob += e_td
                    #println("e_td = $(e_td)")
                    # only sum contract fixed cost if d == 1 (first micro period)
                    cost = d == 1 ? sum( contract[contract[:period] .== t, :cost_fix][c] * y[t,c] for c in 1:num_contracts[t]) : 0.0
                    cost += sum(contract[contract[:period] .== t, :cost_var][c] * q_td[c] for c in 1:num_contracts[t])
                    cost += sum(drivable[s,:cost] * (drivable[s,:pORc][t] / period_size) * x_td[s] for s in D)  # sum(S_D)
                    cost += sum(storage[s,:cost]*(1+lambda[s]) * g_td[s] for s in ST)
                    cost += sum(n_drivable_uncertain[sigma, :cost] * P_hat_td[sigma] for sigma in SetSigma)  # sum(Ŝ_ND)
                        # FIXME Consertar o codigo abaixo de soma dos custos da bateria no robusto (substituicao de h_td[s])
                    cost += sum(storage[s,:cost] * ( g_td[s] + (lambda[s]*g_td[s] + previous_r_td_rob[s] - r_td[s]) ) for s in ST)
                        # sum(storage[s,:cost] * (g_td[s] + h_td[s]) for s in ST) +
                    cost += (period[t,:cost_out] * e_td)
                    acc_cost += cost
                    push!(variable_values_after_uncertainty_df, ["Robust", reoptimize, scenario_id, t, d, rob_value_for_period, rob_opt_time_spent,
                        calc_method, q_td, x_td, e_td, r_td, g_td, h_td, gap_rob, cost, acc_cost])
                end
            end
            # "deterministic"
            q_t, x_t, e_t, r_t, g_t, h_t = calculate_variables_from_deterministic(instance_as_dict, det_solution, t)
            # sum(Ŝ_ND) is already stored in variable gap_det (see above)
            gap_det += sum(n_drivable[s,:pORc][t] for s in ND)  # sum(S_ND)
            gap_det += sum(q_t[c] for c in 1:num_contracts[t])
            gap_det += sum(drivable[s,:pORc][t] * x_t[s] for s in D)
            gap_det += sum(h_t[s] for s in ST) - sum(g_t[s] for s in ST)
            gap_det += e_t
            cost = sum( contract[contract[:period] .== t, :cost_fix][c] * y_det[t,c] for c in 1:num_contracts[t])
            cost += sum(contract[contract[:period] .== t, :cost_var][c] * q_t[c] for c in 1:num_contracts[t])
            cost += sum(drivable[s,:cost] * drivable[s,:pORc][t] * x_t[s] for s in D)
            cost += sum(storage[s,:cost] * (g_t + h_t) for s in ST)
            cost += (period[t,:cost_out] * e_t)
            cost += sum(n_drivable_uncertain[sigma, :cost] * P_hat_t[sigma] for sigma in SetSigma)
            if isa(cost, Array)  # Gambiarra para bug no tipo de dados
                cost = cost[1]
            end
            push!(variable_values_after_uncertainty_df, ["Deterministic", reoptimize, scenario_id, t, 1, det_value_for_period, det_opt_time_spent,
                "N/A", q_t, x_t, e_t, r_t, g_t, h_t, gap_det, cost, cost])
        end
        scenario_id += 1
    end
    return variable_values_after_uncertainty_df
end

# Calculate the value of each model variable according to microperiod d
function calculate_model_variable_delta(q, x, e, g, period_size, d, calc_method = "all_beginning")
    if calc_method == "average"  # Divide term (e, q, etc.) by the period size
        return copy(q) ./ period_size, copy(x) ./ period_size, copy(e) ./ period_size, copy(g) ./ period_size
    else  # "full value in the first microperiod d"
        if d == 1
            return copy(q), copy(x), copy(e), copy(g)
        else
            return q .* 0.0, x .* 0.0, e .* 0.0, g .* 0.0
        end
    end
end

function get_optimization_result_base_filename(output_path, test_name, instance_name)
    return joinpath(normpath(output_path), test_name * "_RCCP_Sim_OptData_" * instance_name)
end

function obtain_robust_optimization_model_results(output_path, test_name, instance_name, instance_as_dict, time_limit, cplex_old, general_logger)
    base_filename = get_optimization_result_base_filename(output_path, test_name, instance_name)
    opt_var_file = base_filename * ".jld"
    opt_csv_file = base_filename * ".csv"
    opt_file_zip = base_filename * ".zip"
    tmp_file_path = "/tmp/rccp_temp_opt_df.jld"
    println(general_logger, "Trying to obtain robust model results...")
    if isfile(opt_var_file)   # Check if optimization results are cached in file
        println("Skipping optimization runs. File already exists: $(opt_var_file)")
        println(general_logger, "Skipping optimization runs. File already exists: $(opt_var_file)")
        # Read execution data back to dataframe
        opt_var = load(opt_var_file)
        opt_df = opt_var["sol"]
        println(general_logger, "Optimization data read successfully.\n")
    elseif isfile(opt_file_zip)  # v2.0 feature
        println("Skipping optimization runs. ZIP File already exists: $(opt_file_zip)")
        println(general_logger, "Skipping optimization runs. ZIP File already exists: $(opt_file_zip)")
        r = ZipFile.Reader(opt_file_zip)
        for f in r.files
            if extension(f.name) == ".jld"
                println("Reading file from zip: $(f.name)")
                println(general_logger, "Reading file from zip: $(f.name)")
                s = readstring(f)
                tmp_file_df = open(tmp_file_path, "w+")
                write(tmp_file_df, s)
                close(tmp_file_df)
                # Read execution data back to dataframe
                opt_var = load(tmp_file_path)
                opt_df = opt_var["sol"]
                break
            end
        end
        close(r)
        println(general_logger, "Optimization data read successfully from zip file.\n")
    else
        output_file_inf_log = joinpath(normpath(pwd()), "Infeasible_RCCP.log")
        infeasible_logger = open(output_file_inf_log, "a+")
        opt_df = run_robust_optimization_models(instance_name, instance_as_dict, time_limit, cplex_old, general_logger, infeasible_logger)
        close(infeasible_logger)
        save_optimization_results_dataframe_to_file(output_path, test_name, instance_name, opt_df, general_logger)
    end
    return opt_df
end

function save_optimization_results_dataframe_to_file(output_path, test_name, instance_name, opt_df, general_logger)
    base_filename = get_optimization_result_base_filename(output_path, test_name, instance_name)
    opt_var_file = base_filename * ".jld"
    opt_csv_file = base_filename * ".csv"
    #println("\nSaving optimization results data file to $(opt_csv_file)...")
    println(general_logger, "\nSaving optimization results data file to $(opt_csv_file)...")
    try
        CSV.write(opt_csv_file, opt_df; delim=';')
        # Serialize dataframe to file
        save(opt_var_file, "sol", opt_df)
        output_file_zip = base_filename * ".zip"
        move_files_to_zip_archive(output_file_zip, [opt_var_file, opt_csv_file])
        println(general_logger, "Done.")
    catch y2
        println("ERROR Writing output file $(opt_csv_file). Cause : $(y2).")
        println(general_logger, "ERROR Writing output file $(opt_csv_file). Cause : $(y2).")
    end
end
