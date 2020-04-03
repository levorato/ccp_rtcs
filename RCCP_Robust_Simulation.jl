# =================================
# RCCP_Robust_Simulation.jl
# RTCS Simulation based on scenarios.
# =================================

using JuMP
using Cbc
using CSV
using DataFrames
#using DataArrays
using CPLEX
using JLD2, FileIO
import MathProgBase

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

function sim_heuristic()
    # Heurística de tempo real. Opcoes:
    # a) Olhando apenas o y: Dentro de cada subperíodo de tempo, olhar apenas os contratos que serão utilizados (disponíveis para o
    #    período) e as demandas / ofertas de energia dos dispositivos. Com base nisso, compra-se ou vende-se energia por meio e por fora
    #    dos contratos. Não se olha para as variáveis que dependem da incerteza.
    # b) Olhando o y e as demais variáveis do robusto atreladas ao período t' : usa informações dos contratos utilizados e dos
    #    percentuais e quantidades usados da bateria e dos demais dispositivos (variáveis q, r, g, etc.), fornecidos pela multiplicação dos
    #    valores das variáveis retornados pelo robusto vezes o valor realizado de P para o período (obtido nos dados do cenário: soma dos
    #    valores de P de cada microperíodo). Com esses dados, roda-se a heurística para decidir o quanto usar de cada um.
    # c) Olhando o y, as demais variáveis do robusto e o conhecimento das situações de instantes passados.

end

# Reoptimize: Opções de execução da simulação:
# - Sem re-otimização: roda-se o robusto apenas na primeira vez e a cada instante t' utiliza-se os valores das variáveis por ele
#   retornadas, multiplicando-se os valores pelas realizações do parâmetro incerto P_hat.
# - Com re-otimização: "Solve the Robust Model starting from t' "roda-se novamente o robusto a cada instante t' como se esse fosse
#   o primeiro período de tempo do problema, fixando as variáveis y. As informações correspondentes aos períodos que já passaram
#   são ignoradas.
function simulate(instance_as_dict, test_name, instance_name, output_path, sim_strategy, model_policy, time_limit, cplex_old,
        general_logger, scenario_filter = Int64[], reoptimize = false, verbose = true)
    println("\nNew Simulation")
    println("Simulation strategy is $(sim_strategy)")
    println("Re-optimization is enabled? $(reoptimize)")
    output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "trace"])
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
        return [0.0, 0.0, zeros(Int64, nbC)]
    end
    if size(scenario_filter, 1) > 0
        println(general_logger, "Scenario filter is ENABLED. Filter = $(scenario_filter)")
        println("Scenario filter is ENABLED. Filter = $(scenario_filter)")
    end
    sum_time_spent = 0
    sum_reopt_time_spent = 0
    trace_df_list = Array[]

    # Run optimization models or obtain solution values if already stored in file
    opt_df = obtain_robust_optimization_model_results(output_path, test_name,
        instance_name, instance_as_dict, time_limit, cplex_old, general_logger)
    # Fix the contracts from initial robust solution for t = 1 (y)
    y_rob = opt_df[1, :RobSolution]
    y_det = opt_df[1, :DetSolution]

    scenario_id = -1
    total_processing_time = 0.0
    for scenario in scenario_list  # For each scenario in instance file
        scenario_id += 1
        println(general_logger, "Processing scenario #$(scenario_id) of $(nbScenarios)...")
        if size(scenario_filter, 1) > 0 && !(scenario_id in scenario_filter)
            println(general_logger, "Skipping scenario #$(scenario_id) of $(nbScenarios) due to filter.")
            continue
        end
        # Define a logfile for writing
        output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "log", instance_name])
        scenario_subpath = create_trace_scenario_filename(
                test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id)
        output_file_base = joinpath(normpath(output_path), scenario_subpath)
        output_file_log = output_file_base * ".log"
        scenario_logger = open(output_file_log, "w+")
        # Define a tracefile for writing
        output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "trace", instance_name])
        output_file_base = joinpath(normpath(output_path), scenario_subpath)
        output_file_df = output_file_base * ".csv"
        output_file_jld2 = output_file_base * ".jld2"
        # Create a new dataframe for this scenario
        trace_df, var_df = create_empty_trace_dataframe()
        try
            if verbose println("===========================================================================") end
            if verbose println(scenario_logger, "===========================================================================") end
            if verbose println("Scenario #$(scenario_id) of $(nbScenarios)") end
            if verbose println(scenario_logger, "Scenario #$(scenario_id) of $(nbScenarios)") end
            if verbose println(general_logger, "Scenario #$(scenario_id) of $(nbScenarios)") end

            # Calculate the initial battery charge for each battery in ST
            initial_r_t = zeros(Float64, nbSt)
            for s in ST
                initial_r_t[s] = storage[s,:uInit]
            end
            previous_r_td_det = copy(initial_r_t)
            previous_r_td_rob = copy(initial_r_t)
            acc_cost_det = 0.0
            acc_cost_rob = 0.0
            previous_batt_levels_det = zeros(Float64, nbSt)
            previous_batt_levels_rob = zeros(Float64, nbSt)
            # Generate the array of drivable devices previous accumulated charge value
            previous_drivable_charge_det = zeros(Float64, nbD)
            previous_drivable_charge_rob = zeros(Float64, nbD)
            scenario_processing_start_time = time_ns()
            for t in 1:nbT   # For each time period t
                if verbose println(scenario_logger, "===========================================================================") end
                if verbose println(scenario_logger, "(t) Period #$(t) of $(nbT)") end
                period_size = period_info[:size][t]
                #det_value = 0
                elapsed_time_det = 0.0
                elapsed_time_rob = 0.0
                if (!reoptimize) || (t == 1)  # Use model solution starting with t = 1
                    det_value_for_period = opt_df[1, :DetValue]
                    rob_value_for_period = opt_df[1, :RobValue]
                    if t == 1
                        det_opt_time_spent = opt_df[1, :DetTime]
                        rob_opt_time_spent = opt_df[1, :RobTime]
                    else
                        det_opt_time_spent = 0.0
                        rob_opt_time_spent = 0.0
                    end
                    det_solution, rob_solution = read_variables_from_solution_df(opt_df, 1)
                else   # Use model solution starting with t = t'
                    # Reoptimize the Deterministic RCCP Model
                    t1 = time_ns()  # ******************** measure start time
                    det_value_reopt, det_solution_reopt, det_solve_time_reopt, status_det_reopt =
                            re_optimize_deterministic(instance_as_dict, time_limit, cplex_old, t, y_det,
                                previous_batt_levels_det, previous_drivable_charge_det, general_logger)
                    det_value_for_period = det_value_reopt
                    det_opt_time_spent = det_solve_time_reopt
                    det_solution = det_solution_reopt
                    elapsed_time_det = time_ns() - t1   # ****************** stop time
                    # Read solutions from previous time period (t - 1)
                    det_solution_tminus1, rob_solution_tminus1 = read_variables_from_solution_df(opt_df, t - 1)
                    if status_det_reopt == false && t > 1   # FIXME workaround for infeasible reopt model result
                        det_value_reopt = opt_df[t - 1, :DetValue]
                        det_solution_reopt = det_solution_tminus1
                        status_det_reopt = true
                        println("WARN: reusing det solution from previous period (t - 1), since current solution is INFEASIBLE.")
                    end
                    # Reoptimize the Robust RCCP Model
                    t1 = time_ns()  # ******************** measure start time
                    rob_value_reopt, rob_solution_reopt, rob_solve_time_reopt, status_rob_reopt =
                            re_optimize_robust(instance_as_dict, time_limit, cplex_old, t, y_rob,
                            previous_batt_levels_rob, previous_drivable_charge_rob, general_logger)
                    rob_value_for_period = rob_value_reopt
                    rob_opt_time_spent = rob_solve_time_reopt
                    rob_solution = rob_solution_reopt
                    elapsed_time_rob = time_ns() - t1  # ******************* stop time
                    if status_rob_reopt == false && t > 1   # FIXME workaround for infeasible reopt model result
                        rob_value_reopt = opt_df[t - 1, :RobValue]
                        rob_solution_reopt = rob_solution_tminus1
                        status_rob_reopt = true
                        println("WARN: reusing rob solution from previous period (t - 1), since current solution is INFEASIBLE.")
                    end
                    store_optimization_results_in_dataframe(instance_name, instance_as_dict, opt_df, scenario_id, t, det_value_reopt,
                                                            det_solution_reopt, det_solve_time_reopt, status_det_reopt,
                                                            rob_value_reopt, rob_solution_reopt, rob_solve_time_reopt, status_rob_reopt,
                                                            general_logger)
                    det_value_for_period = opt_df[t, :DetValue]
                    rob_value_for_period = opt_df[t, :RobValue]
                    det_opt_time_spent = opt_df[t, :DetTime]
                    rob_opt_time_spent = opt_df[t, :RobTime]
                    det_solution, rob_solution = read_variables_from_solution_df(opt_df, t)
                end
                P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t
                q_t_det = zeros(Float64, num_contracts[t])
                q_t_rob = zeros(Float64, num_contracts[t])

                # v2 : use P_hat_t = avg_P_hat_t for uncertain load values
                # Obtain average values for uncertain non-drivable devices load, such that P = (Pmin + Pmax)/2
                P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t
                if !isempty(n_drivable_uncertain)
                    println(general_logger, "\nNOTE: Calculating average load values for uncertain non-drivable devices.")
                    for i in 1:size(n_drivable_uncertain, 1)
                        power = (n_drivable_uncertain[i, :Pmin][t] + n_drivable_uncertain[i, :Pmax][t]) / 2.0
                        P_hat_t[i] = power
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
                    if verbose println(scenario_logger, "\n  DETERMINISTIC   R T C S") end
                    previous_q_t_det = copy(q_t_det)
                    q_td, x_td, e_td, r_td, g_td, h_td, gap, cost = rtcs(t, d, period_size, sim_strategy,
                                                                model_policy, instance_as_dict, P_hat_t, P_hat_td,
                                                                y_det, "deterministic", det_solution, previous_r_td_det,
                                                                previous_drivable_charge_det, previous_q_t_det,
                                                                scenario_logger, verbose)
                    previous_r_td_det = copy(r_td)
                    for c in 1:num_contracts[t]
                        q_t_det[c] += q_td[c]
                    end
                    # For a strange reason, cost is being returned as a one-element array, fix this
                    cost = fix_cost_variable_type(cost)
                    acc_cost_det += cost
                    # IMPORTANT : process drivable devices previous accumulated charge value
                    for s in D
                        previous_drivable_charge_det[s] += (drivable[s,:pORc][t] / period_size) * x_td[s]
                    end
                    if verbose println(scenario_logger, "    ACCUMULATED COST                    : $(acc_cost_det)") end
                    elapsed_time_det += (time_ns() - t1)  # ********************** stop time
                    push!(trace_df, ["Deterministic", sim_strategy, reoptimize, model_policy, scenario_id,
                        t, d, det_value_for_period, rob_value_for_period, (d == 1) ? det_opt_time_spent : 0.0,
                        e_td, gap, cost, elapsed_time_det/1.0e9])
                    push!(var_df, ["Deterministic", sim_strategy, reoptimize, model_policy, scenario_id,
                        t, d,
                        q_td, x_td, e_td, r_td, g_td, h_td])


                    if verbose println(scenario_logger, "\n  ROBUST   R T C S") end
                    t1 = time_ns()  # *********************** measure start time
                    previous_q_t_rob = copy(q_t_rob)
                    q_td, x_td, e_td, r_td, g_td, h_td, gap, cost = rtcs(t, d, period_size, sim_strategy,
                                                            model_policy, instance_as_dict, P_hat_t, P_hat_td,
                                                            y_rob, "robust", rob_solution, previous_r_td_rob,
                                                            previous_drivable_charge_rob, previous_q_t_rob,
                                                            scenario_logger, verbose)
                    previous_r_td_rob = copy(r_td)
                    for c in 1:num_contracts[t]
                        q_t_rob[c] += q_td[c]
                    end

                    # For a strange reason, cost is being returned as a one-element array, fix this
                    cost = fix_cost_variable_type(cost)
                    acc_cost_rob += cost
                    # IMPORTANT : process drivable devices previous accumulated charge value
                    for s in D
                        previous_drivable_charge_rob[s] += (drivable[s,:pORc][t] / period_size) * x_td[s]
                    end
                    if verbose println(scenario_logger, "    ACCUMULATED COST                    : $(acc_cost_rob)") end
                    elapsed_time_rob += (time_ns() - t1)  # ********************** stop time
                    push!(trace_df, ["Robust", sim_strategy, reoptimize, model_policy, scenario_id,
                        t, d, det_value_for_period, rob_value_for_period, (d == 1) ? rob_opt_time_spent : 0.0,
                        e_td,
                        gap, cost, elapsed_time_rob/1.0e9])
                    push!(var_df, ["Robust", sim_strategy, reoptimize, model_policy, scenario_id,
                        t, d,
                        q_td, x_td, e_td, r_td, g_td,
                        h_td])
                    sum_time_spent += det_opt_time_spent + rob_opt_time_spent
                    sum_reopt_time_spent += det_opt_time_spent + rob_opt_time_spent
                end
                previous_batt_levels_det = copy(previous_r_td_det)
                previous_batt_levels_rob = copy(previous_r_td_rob)
            end
            scenario_end_processing_time = time_ns()
            scenario_processing_time = (scenario_end_processing_time - scenario_processing_start_time) / 1.0e9
            total_processing_time += scenario_processing_time
            println(scenario_logger, "\n===========================================================================")
            println(scenario_logger, " Total processing time : $(scenario_processing_time)")
            # TRACE FILE: Create an output file to write the dataframe to disk
            println(general_logger, "\nSaving trace CSV file to $(output_file_df)...")
            CSV.write(output_file_df, trace_df)
            @save output_file_jld2 trace_df var_df
            #push!(trace_df_list, trace_df)
        catch e
            bt = catch_backtrace()
            msg = sprint(showerror, e, bt)
            println(msg)
            println(general_logger, "Exception thrown: \n" * msg)
        finally
            # Close the log file
            close(scenario_logger)
            flush(general_logger)
            # Move the scenario log text file and scenario csv trace file to a zip file to save disk space
            output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "zip", instance_name])
            output_file_base = joinpath(normpath(output_path), create_trace_scenario_filename(
                    test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id))
            output_file_zip = output_file_base * ".zip"
            move_files_to_zip_archive(output_file_zip, [output_file_log, output_file_df, output_file_jld2])
        end
    end
    save_optimization_results_dataframe_to_file(output_path, test_name, instance_name, opt_df, general_logger)
    # Dados a serem armazenados nas simulações para análise futura (e.g. box plot):
    #    - Custo da solução (determinístico x robusto)
    #    - Quantidade de energia comprada fora do contrato
    #    - Nível na bateria em cada instante de tempo
    #    - Tempo gasto para rodar o PPL com a re-otimização.
    return trace_df_list, total_processing_time
end

function test_solve_robust_from_t_2()
    datafile = pwd() * "/../notebooks/data/antoine/A_instance2_11scen_1NDU_shift_test.txt"
    instance_as_dict = read_tabulated_data(datafile)
    output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "log"])
    output_file_log = joinpath(normpath(output_path), "test_solve_robust_from_t_2" * "_A_instance2_11scen_1NDU_shift_test" * ".log")
    output_file_inf_log = joinpath(normpath(output_path), "Infeasible_solve_robust_from_t_2" * "_A_instance2_11scen_1NDU_shift_test" * ".log")
    scenario_logger = open(output_file_log, "w+")
    infeasible_logger = open(output_file_inf_log, "w+")
    rob_value, rob_solution, solve_time, status = solve_robust_model(instance_as_dict, scenario_logger, infeasible_logger, 1800.0, false, 2)
    close(scenario_logger)
    close(infeasible_logger)
end

function simulate_antoine_instance(base_folder, test_set, cplex_old, scenario_filter = Int64[])
    test_name = "Antoine"
    df = DataFrame(Instance = String[], TestName = String[], ModelName = String[],
        Strategy = String[], Reoptimize = Bool[], TimeSpent = Float64[],
        ReOptTime = Float64[], Contracts = Array{Int64, 2}[])
    println("\n===================\n  RCCP Simulator\n===================\n")
    time_limit = 14000.0
    trace_df, var_df = create_empty_trace_dataframe()
    variable_values_after_uncertainty_df = create_empty_dataframe_variable_values_after_uncertainty()
    for instance_name in test_set
        instance_filename = joinpath(normpath(base_folder), instance_name)
        # Open a logfile for writing
        output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "log"])
        output_file_log = joinpath(normpath(output_path), test_name * "_" * instance_name * ".log")
        scenario_logger = open(output_file_log, "w+")
        try
            println(scenario_logger, "===========================================================================")
            println("Simulation for instance $(instance_name)")
            println(scenario_logger, "Simulation for instance $(instance_name)")
            # Read instance file
            instance_as_dict = read_tabulated_data(instance_filename)
            scenario_list = instance_as_dict["scenarios"]
            nbScenarios = size(scenario_list, 1)
            # For each period t, obtain "pure" model variables values after uncertainty is revealed
            #println(scenario_logger, "Calculating pure model variables...")
            #for reoptimize in [false]  # false
            #    variable_values_after_uncertainty_df_2 = calculate_pure_model_variables_after_uncertainty_is_known(instance_filename,
            #                                                test_name, instance_name, reoptimize, time_limit, cplex_old, scenario_logger)
            #    variable_values_after_uncertainty_df = vcat(variable_values_after_uncertainty_df, variable_values_after_uncertainty_df_2)
            #end
            output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "trace"])
            # Run RTCS Simulation
            total_p_time = 0.0
            for model_policy in ["ignore_model", "full_model"] # "batteries", "batteries_and_drivable"
                if model_policy == "full_model"
                    for sim_strategy in ["conservative", "audacious", "cheapest"]
                        for reopt in [false, true]
                            println("\nSimulation for $(sim_strategy) x $(model_policy) x ReOpt=$(string(reopt))")
                            println(scenario_logger, "Simulation for $(sim_strategy) x $(model_policy) x ReOpt=$(string(reopt))")
                            trace_df_2_list, p_time = simulate(instance_as_dict, test_name, instance_name, output_path, sim_strategy,
                                model_policy, time_limit, cplex_old, scenario_logger, scenario_filter, reopt)
                            println(scenario_logger, "Done. Processing time : $(string(p_time))")
                            println("Done. Processing time : $(string(p_time))")
                            flush(scenario_logger)
                            total_p_time += p_time
                            #trace_df = vcat(trace_df, trace_df_2)
                        end
                    end
                else  # 'ignore_model' uses only y variable info, so no reoptimization is needed
                    reopt = false
                    for sim_strategy in ["conservative", "audacious", "cheapest"]
                        println(scenario_logger, "Simulation for $(sim_strategy) x $(model_policy) x ReOpt=$(string(reopt))")
                        trace_df_2_list, p_time = simulate(instance_as_dict, test_name, instance_name, output_path, sim_strategy,
                            model_policy, time_limit, cplex_old, scenario_logger, scenario_filter, reopt)
                        println(scenario_logger, "Done. Processing time : $(string(p_time))")
                        flush(scenario_logger)
                        total_p_time += p_time
                        #trace_df = vcat(trace_df, trace_df_2)
                    end
                end
            end
            println(scenario_logger, "Simulation done. Processing time : $(string(total_p_time))")
            # Save variable_values_after_uncertainty_df to CSV and JLD files
            #output_file_vv = joinpath(normpath(output_path), test_name * "_RCCP_VariableValues_" * instance_name)
            #println("\nSaving variable values CSV file to $(output_file_vv)...")
            #println(scenario_logger, "\nSaving variable values CSV file to $(output_file_vv)...")
            #CSV.write(output_file_vv * ".csv", variable_values_after_uncertainty_df)
            # Serialize dataframe to file
            #save(output_file_vv * ".jld", "variable_values", variable_values_after_uncertainty_df)

            # Serialize dataframe to file
            #save(output_file * ".jld", "trace", trace_df)
            # push!(df, [[instance_name, test_name, "Robust_RCCP_Sim", sim_strategy, reopt] ; result_as_array])
            # Create an output file to write the dataframe to disk
            #output_file = joinpath(normpath(output_path), test_name * "_Robust_RCCP_Sim_Summary.csv")
            #println("\nSaving results file to $(output_file)...")
            #CSV.write(output_file, df)
            concatenate_trace_df_list(test_name, instance_name, scenario_logger)
            println("\nDone.")
            println(scenario_logger, "\nDone.")
        catch e
            bt = catch_backtrace()
            msg = sprint(showerror, e, bt)
            println(msg)
            println(scenario_logger, "Exception thrown: \n" * msg)
        finally
            # Close the log file
            close(scenario_logger)
        end
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

function run_simulation()
    # run simulation tests
    test_set = ["A_instance2_1NDU_Cons_1000s_skewed-left.txt", "A_instance2_1NDU_Cons_1000s_skewed-right.txt", "A_instance2_1NDU_Cons_1000s_uniform.txt"]
    
    base_folder = pwd() * "/instances/full_instances/"
    cplex_old = true
    #scenario_filter = [0, 1]
    simulate_antoine_instance(base_folder, test_set, cplex_old)  #, scenario_filter)
end


#run_simulation()

