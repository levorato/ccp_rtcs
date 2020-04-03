include("RCCP_Robust_Simulation.jl")

function simulate_antoine_instance2(base_folder, test_set, cplex_old, scenario_filter = Int64[])
    test_name = "Antoine"
    strategies_to_simulate = ["conservative", "audacious"]  # , "cheapest"]
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
            println(scenario_logger, "Calculating pure model variables...")
            for reoptimize in [false]  # false
                variable_values_after_uncertainty_df_2 = calculate_pure_model_variables_after_uncertainty_is_known(instance_filename,
                                                            test_name, instance_name, reoptimize, time_limit, cplex_old, scenario_logger)
                variable_values_after_uncertainty_df = vcat(variable_values_after_uncertainty_df, variable_values_after_uncertainty_df_2)
            end
            output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "trace"])
            # Run RTCS Simulation
            total_p_time = 0.0
            for model_policy in ["ignore_model", "full_model"] # "batteries", "batteries_and_drivable"
                if model_policy == "full_model"
                    for sim_strategy in strategies_to_simulate
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
                    for sim_strategy in strategies_to_simulate
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
            output_file_vv = joinpath(normpath(output_path), test_name * "_RCCP_VariableValues_" * instance_name)
            println("\nSaving variable values CSV file to $(output_file_vv)...")
            println(scenario_logger, "\nSaving variable values CSV file to $(output_file_vv)...")
            CSV.write(output_file_vv * ".csv", variable_values_after_uncertainty_df)
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


function test_sim()
    # run simulation tests
    # test_solve_robust_from_t_2()
    #base_folder = pwd() * "/../notebooks/data/antoine/"
    #test_set = ["A_instance2_11scen.txt"]  #  "A_instance2_11scen_1NDU.txt"]

    #base_folder = pwd() * "/instances/"
    #test_set = ["A_instance2_HighOC.txt"] # ["A_instance2_11scen_noNDU.txt", "A_instance2_11scen_noNDU.txt"]
    #test_set = ["A_instance2_1000s_1NDU_x1000COST.txt", "A_instance2_1000s_x1000COST.txt"]  # "A_instance2_11scen_x1000COST.txt"
    #test_set = ["A_instance3_100s_OC0.1406.txt", "A_instance3_100s_OC0.2812.txt", "A_instance3_100s_OC0.5625.txt", "A_instance3_100s_OC1.125.txt"]

    #base_folder = pwd() * "/instance_gen/output/scenarios/"
    #test_set = ["A_instance2_100s.txt"]
    #test_set = ["A_instance6_extreme_100s.txt", "A_instance4_100s.txt"]

    base_folder = pwd() * "/instances/toy"
    already_done = ["OC19_Ct_default_DS3.5_ST3_NDU_2x.txt",
"OC19_Ct_default_DS3.5_ST3_NDU_4x.txt",
"OC19_Ct_default_DS3.5_ST3_NDU_default.txt",
"OC19_Ct_default_DS3.5_ST6_NDU_2x.txt",
"OC19_Ct_default_DS3.5_ST6_NDU_4x.txt",
"OC19_Ct_default_DS3.5_ST6_NDU_default.txt",
"OC19_Ct_default_DS3.5_ST12_NDU_2x.txt",
"OC19_Ct_default_DS3.5_ST12_NDU_4x.txt",
"OC19_Ct_default_DS3.5_ST12_NDU_default.txt",
"OC19_Ct_default_DS3.5_ST24_NDU_2x.txt",
"OC19_Ct_default_DS3.5_ST24_NDU_4x.txt",
"OC19_Ct_default_DS3.5_ST24_NDU_default.txt",
"OC19_Ct_default_DS7.5_ST3_NDU_4x.txt",
"OC19_Ct_default_DS7.5_ST3_NDU_default.txt",
"OC19_Ct_default_DS7.5_ST6_NDU_4x.txt",
"OC19_Ct_default_DS7.5_ST6_NDU_default.txt",
"OC19_Ct_default_DS7.5_ST12_NDU_2x.txt",
"OC19_Ct_default_DS7.5_ST12_NDU_4x.txt",
"OC19_Ct_default_DS7.5_ST12_NDU_default.txt",
"OC19_Ct_default_DS7.5_ST24_NDU_2x.txt",
"OC19_Ct_default_DS7.5_ST24_NDU_4x.txt",
"OC19_Ct_default_DS7.5_ST24_NDU_default.txt",
"OC19_Ct_default_DS13.0_ST3_NDU_4x.txt",
"OC19_Ct_default_DS13.0_ST3_NDU_default.txt",
"OC19_Ct_default_DS13.0_ST6_NDU_4x.txt",
"OC19_Ct_default_DS13.0_ST6_NDU_default.txt",
"OC19_Ct_default_DS13.0_ST12_NDU_2x.txt",
"OC19_Ct_default_DS13.0_ST12_NDU_4x.txt",
"OC19_Ct_default_DS13.0_ST12_NDU_default.txt",
"OC19_Ct_default_DS13.0_ST24_NDU_2x.txt",
"OC19_Ct_default_DS13.0_ST24_NDU_4x.txt",
"OC19_Ct_default_DS13.0_ST24_NDU_default.txt",
"OC19_Ct_default_DS26.0_ST3_NDU_4x.txt",
"OC19_Ct_default_DS26.0_ST3_NDU_default.txt",
"OC19_Ct_default_DS26.0_ST6_NDU_4x.txt",
"OC19_Ct_default_DS26.0_ST6_NDU_default.txt",
"OC19_Ct_default_DS26.0_ST12_NDU_2x.txt",
"OC19_Ct_default_DS26.0_ST12_NDU_4x.txt",
"OC19_Ct_default_DS26.0_ST12_NDU_default.txt",
"OC19_Ct_default_DS26.0_ST24_NDU_2x.txt",
"OC19_Ct_default_DS26.0_ST24_NDU_4x.txt",
"OC19_Ct_default_DS26.0_ST24_NDU_default.txt",
"OC19_Ct_high_DS3.5_ST3_NDU_2x.txt",
"OC19_Ct_high_DS3.5_ST3_NDU_4x.txt",
"OC19_Ct_high_DS3.5_ST3_NDU_default.txt",
"OC19_Ct_high_DS3.5_ST6_NDU_2x.txt",
"OC19_Ct_high_DS3.5_ST6_NDU_4x.txt",
"OC19_Ct_high_DS3.5_ST6_NDU_default.txt",
"OC19_Ct_high_DS3.5_ST12_NDU_2x.txt",
"OC19_Ct_high_DS3.5_ST12_NDU_4x.txt",
"OC19_Ct_high_DS3.5_ST12_NDU_default.txt",
"OC19_Ct_high_DS3.5_ST24_NDU_2x.txt",
"OC19_Ct_high_DS3.5_ST24_NDU_4x.txt",
"OC19_Ct_high_DS3.5_ST24_NDU_default.txt",
"OC19_Ct_high_DS7.5_ST3_NDU_4x.txt",
"OC19_Ct_high_DS7.5_ST3_NDU_default.txt",
"OC19_Ct_high_DS7.5_ST6_NDU_2x.txt",
"OC19_Ct_high_DS7.5_ST6_NDU_4x.txt",
"OC19_Ct_high_DS7.5_ST6_NDU_default.txt",
"OC19_Ct_high_DS7.5_ST12_NDU_2x.txt",
"OC19_Ct_high_DS7.5_ST12_NDU_4x.txt",
"OC19_Ct_high_DS7.5_ST12_NDU_default.txt",
"OC19_Ct_high_DS7.5_ST24_NDU_2x.txt",
"OC19_Ct_high_DS7.5_ST24_NDU_4x.txt",
"OC19_Ct_high_DS7.5_ST24_NDU_default.txt",
"OC19_Ct_high_DS13.0_ST3_NDU_4x.txt",
"OC19_Ct_high_DS13.0_ST3_NDU_default.txt",
"OC19_Ct_high_DS13.0_ST6_NDU_4x.txt",
"OC19_Ct_high_DS13.0_ST6_NDU_default.txt",
"OC19_Ct_high_DS13.0_ST12_NDU_2x.txt",
"OC19_Ct_high_DS13.0_ST12_NDU_4x.txt",
"OC19_Ct_high_DS13.0_ST12_NDU_default.txt",
"OC19_Ct_high_DS13.0_ST24_NDU_2x.txt",
"OC19_Ct_high_DS13.0_ST24_NDU_4x.txt",
"OC19_Ct_high_DS13.0_ST24_NDU_default.txt",
"OC19_Ct_high_DS26.0_ST3_NDU_4x.txt",
"OC19_Ct_high_DS26.0_ST3_NDU_default.txt",
"OC19_Ct_high_DS26.0_ST6_NDU_4x.txt",
"OC19_Ct_high_DS26.0_ST6_NDU_default.txt",
"OC19_Ct_high_DS26.0_ST12_NDU_2x.txt",
"OC19_Ct_high_DS26.0_ST12_NDU_4x.txt",
"OC19_Ct_high_DS26.0_ST12_NDU_default.txt",
"OC19_Ct_high_DS26.0_ST24_NDU_2x.txt",
"OC19_Ct_high_DS26.0_ST24_NDU_4x.txt",
"OC19_Ct_high_DS26.0_ST24_NDU_default.txt",
"OC19_Ct_low_DS3.5_ST3_NDU_2x.txt",
"OC19_Ct_low_DS3.5_ST3_NDU_4x.txt",
"OC19_Ct_low_DS3.5_ST3_NDU_default.txt",
"OC19_Ct_low_DS3.5_ST6_NDU_2x.txt",
"OC19_Ct_low_DS3.5_ST6_NDU_4x.txt",
"OC19_Ct_low_DS3.5_ST6_NDU_default.txt",
"OC19_Ct_low_DS3.5_ST12_NDU_2x.txt",
"OC19_Ct_low_DS3.5_ST12_NDU_4x.txt",
"OC19_Ct_low_DS3.5_ST12_NDU_default.txt",
"OC19_Ct_low_DS3.5_ST24_NDU_2x.txt",
"OC19_Ct_low_DS3.5_ST24_NDU_4x.txt",
"OC19_Ct_low_DS3.5_ST24_NDU_default.txt",
"OC19_Ct_low_DS7.5_ST3_NDU_2x.txt",
"OC19_Ct_low_DS7.5_ST3_NDU_4x.txt",
"OC19_Ct_low_DS7.5_ST3_NDU_default.txt",
"OC19_Ct_low_DS7.5_ST6_NDU_2x.txt",
"OC19_Ct_low_DS7.5_ST6_NDU_4x.txt",
"OC19_Ct_low_DS7.5_ST6_NDU_default.txt",
"OC19_Ct_low_DS7.5_ST12_NDU_2x.txt",
"OC19_Ct_low_DS7.5_ST12_NDU_4x.txt",
"OC19_Ct_low_DS7.5_ST12_NDU_default.txt",
"OC19_Ct_low_DS7.5_ST24_NDU_2x.txt",
"OC19_Ct_low_DS7.5_ST24_NDU_4x.txt",
"OC19_Ct_low_DS7.5_ST24_NDU_default.txt",
"OC19_Ct_low_DS13.0_ST3_NDU_2x.txt",
"OC19_Ct_low_DS13.0_ST3_NDU_4x.txt",
"OC19_Ct_low_DS13.0_ST3_NDU_default.txt",
"OC19_Ct_low_DS13.0_ST6_NDU_2x.txt",
"OC19_Ct_low_DS13.0_ST6_NDU_4x.txt",
"OC19_Ct_low_DS13.0_ST6_NDU_default.txt",
"OC19_Ct_low_DS13.0_ST12_NDU_2x.txt",
"OC19_Ct_low_DS13.0_ST12_NDU_4x.txt",
"OC19_Ct_low_DS13.0_ST12_NDU_default.txt",
"OC19_Ct_low_DS13.0_ST24_NDU_2x.txt",
"OC19_Ct_low_DS13.0_ST24_NDU_4x.txt",
"OC19_Ct_low_DS13.0_ST24_NDU_default.txt",
"OC19_Ct_low_DS26.0_ST3_NDU_2x.txt",
"OC19_Ct_low_DS26.0_ST3_NDU_4x.txt",
"OC19_Ct_low_DS26.0_ST3_NDU_default.txt",
"OC19_Ct_low_DS26.0_ST6_NDU_2x.txt",
"OC19_Ct_low_DS26.0_ST6_NDU_4x.txt",
"OC19_Ct_low_DS26.0_ST6_NDU_default.txt",
"OC19_Ct_low_DS26.0_ST12_NDU_2x.txt",
"OC19_Ct_low_DS26.0_ST12_NDU_4x.txt",
"OC19_Ct_low_DS26.0_ST12_NDU_default.txt",
"OC19_Ct_low_DS26.0_ST24_NDU_2x.txt",
"OC19_Ct_low_DS26.0_ST24_NDU_4x.txt",
"OC19_Ct_low_DS26.0_ST24_NDU_default.txt",
"OC19_Ct_med_DS3.5_ST3_NDU_2x.txt",
"OC19_Ct_med_DS3.5_ST3_NDU_4x.txt",
"OC19_Ct_med_DS3.5_ST12_NDU_2x.txt",
"OC19_Ct_med_DS3.5_ST12_NDU_4x.txt",
"OC19_Ct_med_DS3.5_ST12_NDU_default.txt",
"OC19_Ct_med_DS3.5_ST24_NDU_2x.txt",
"OC19_Ct_med_DS3.5_ST24_NDU_4x.txt",
"OC19_Ct_med_DS3.5_ST24_NDU_default.txt",
"OC19_Ct_med_DS13.0_ST3_NDU_4x.txt",
"OC19_Ct_med_DS13.0_ST3_NDU_default.txt",
"OC19_Ct_med_DS13.0_ST6_NDU_4x.txt",
"OC19_Ct_med_DS13.0_ST6_NDU_default.txt",
"OC19_Ct_med_DS13.0_ST12_NDU_2x.txt",
"OC19_Ct_med_DS13.0_ST12_NDU_4x.txt",
"OC19_Ct_med_DS13.0_ST12_NDU_default.txt",
"OC19_Ct_med_DS13.0_ST24_NDU_2x.txt",
"OC19_Ct_med_DS13.0_ST24_NDU_4x.txt",
"OC19_Ct_med_DS13.0_ST24_NDU_default.txt",
"OC19_Ct_med_DS26.0_ST3_NDU_4x.txt",
"OC19_Ct_med_DS26.0_ST3_NDU_default.txt",
"OC19_Ct_med_DS26.0_ST6_NDU_default.txt",
"OC19_Ct_med_DS26.0_ST12_NDU_2x.txt",
"OC19_Ct_med_DS26.0_ST12_NDU_4x.txt",
"OC19_Ct_med_DS26.0_ST12_NDU_default.txt",
"OC19_Ct_med_DS26.0_ST24_NDU_2x.txt",
"OC19_Ct_med_DS26.0_ST24_NDU_4x.txt",
"OC19_Ct_med_DS26.0_ST24_NDU_default.txt"]
    all_files = readdir(base_folder)
    unprocessed_files = [x for x in all_files if !(x in already_done)]
    files_to_process = [x for x in unprocessed_files if startswith(x, "OC19_Ct_")]
    tamanho = size(files_to_process)
    test_set = files_to_process
    println("unprocessed_files = $(tamanho)")
    #cd("/projetos/RCCP")

    cplex_old = true
    scenario_filter = [x for x in 0:100]
    simulate_antoine_instance(base_folder, test_set, cplex_old) #, scenario_filter)
end

test_sim()
