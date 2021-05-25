include("RCCP_Robust_Simulation.jl")

function compute_statistics(base_folder, test_set)
    df = DataFrame(Instance = String[], NDUId = Int64[], ScenarioId = Int64[],
        avg_perc_obs = Float64[], avg_perc_min = Float64[],
        avg_perc_max = Float64[], perc_energy = Float64[])
    summary_df = DataFrame(Instance = String[], NDUId = Int64[],
        avg_perc_obs = Float64[], avg_perc_min = Float64[],
        avg_perc_max = Float64[], avg_perc_energy = Float64[])
    verbose = false
    println("\n===================\n  RCCP UNDS Statistics\n===================\n")
    output_path = create_full_dir(normpath(EXPERIMENT_OUTPUT_FOLDER), ["output", "stats"])
    output_file_df = joinpath(normpath(output_path), "stats_full.csv")
    output_file_df_summary = joinpath(normpath(output_path), "stats_summary.csv")
    for instance_name in test_set
        instance_filename = joinpath(normpath(base_folder), instance_name)
        try
            println("Analysis for instance $(instance_name)")
            flush(stdout)
            # Read instance file
            instance_as_dict = read_tabulated_data(instance_filename)
            # Run Analysis
            nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance_as_dict["instance"])
            nbTau = nbT
            n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
            nbNDU = size(n_drivable_uncertain, 1)
            T, S, C, D, ND, ST, NDU = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt, nbNDU)
            SetTau = 1:nbTau
            SetSigma = 1:nbNDU
            println("Number of periods is $(nbT)")
            println("Number of uncertain devices is $(nbNDU)")
            period_info = instance_as_dict["period"]
            scenario_list = instance_as_dict["scenarios"]
            nbScenarios = size(scenario_list, 1)
            println("Number of scenarios is $(nbScenarios)")
            if nbScenarios == 0
                println("No scenarios found! Aborting analysis.")
                continue
            end
            v_bar = zeros(Float64, nbTau, nbNDU)
            v_hat = zeros(Float64, nbTau, nbNDU)
            p_min = zeros(Float64, nbTau, nbNDU)
            p_max = zeros(Float64, nbTau, nbNDU)
            for sigma in SetSigma
                for tau in SetTau
                    p_min[tau,sigma] = n_drivable_uncertain[sigma, :Pmin][tau]
                    pdt_min = n_drivable_uncertain[sigma, :Pdt_min][tau]
                    p_max[tau,sigma] = n_drivable_uncertain[sigma, :Pmax][tau]
                    pdt_max = n_drivable_uncertain[sigma, :Pdt_max][tau]
                    v_hat[tau,sigma] = Float64(abs(p_max[tau,sigma] - p_min[tau,sigma]) / 2.0)   # peak value interval length
                    v_bar[tau,sigma] = Float64((p_max[tau,sigma] + p_min[tau,sigma]) / 2.0)   # nominal / mean value
                    @assert v_hat[tau,sigma] >= 0
                    if p_max[tau,sigma] > 0 || p_min[tau,sigma] > 0
                        @assert (p_min[tau,sigma] + v_hat[tau,sigma] == v_bar[tau,sigma])
                    elseif p_min[tau,sigma] < 0 || p_max[tau,sigma] < 0
                        #println("sigma=$(sigma), tau=$(tau) : v_bar=$(v_bar[tau,sigma]), v_hat=$(v_hat[tau,sigma]), p_max - v_hat=$(p_max - v_hat[tau,sigma])")
                        @assert abs((p_max[tau,sigma] - v_hat[tau,sigma]) - v_bar[tau,sigma]) <= EPS
                    end
                end
            end
            scenario_id = -1
            avg_perc_min_allsc = zeros(Float64, nbNDU)
            avg_perc_max_allsc = zeros(Float64, nbNDU)
            avg_perc_obs_allsc = zeros(Float64, nbNDU)
            sum_energy_allsc = zeros(Float64, nbNDU)
            for scenario in scenario_list  # For each scenario in instance file
                scenario_id += 1
                perc_min = zeros(Float64, nbT, nbNDU)
                perc_max = zeros(Float64, nbT, nbNDU)
                perc_obs = zeros(Float64, nbT, nbNDU)
                avg_perc_min = zeros(Float64, nbNDU)
                avg_perc_max = zeros(Float64, nbNDU)
                avg_perc_obs = zeros(Float64, nbNDU)
                sum_energy = zeros(Float64, nbNDU)
                sum_v_bar = zeros(Float64, nbNDU)
                for t in 1:nbT   # For each time period t
                    period_size = period_info[!, :size][t]
                    P_hat_t = zeros(Float64, nbNDU)  # P_hat for a whole time period t for each uncertain device
                    for d in 1:period_size  # For each d in period dt
                        P_hat_td = Float64[]
                        for unds_matrix in scenario
                            push!(P_hat_td, unds_matrix[t][d])  # P_hat_td[s] for s in NDU
                        end
                        for s in 1:nbNDU  # Obtain the P_hat matrix with the sum for the whole period t
                            P_hat_t[s] += P_hat_td[s]
                            sum_energy[s] += P_hat_td[s]
                        end
                    end
                    for s in 1:nbNDU  # Obtain the P_hat matrix with the sum for the whole period t
                        perc_obs[t,s] = P_hat_t[s] / v_bar[t,s]  # v_bar : nominal values (average between Pmin and Pmax) for each uncertain device
                        avg_perc_obs[s] += perc_obs[t,s]
                        perc_min[t,s] = p_min[t,s] / v_bar[t,s]
                        avg_perc_min[s] += perc_min[t,s]
                        perc_max[t,s] = p_max[t,s] / v_bar[t,s]
                        avg_perc_max[s] += perc_max[t,s]
                        sum_v_bar[s] += v_bar[t,s]
                    end
                end
                perc_energy = zeros(Float64, nbNDU)
                for s in 1:nbNDU  # Calculate mean values for the percentage of min, max and observed values for P_hat
                    avg_perc_obs[s] /= nbT
                    avg_perc_min[s] /= nbT
                    avg_perc_max[s] /= nbT
                    perc_energy[s] = sum_energy[s] / sum_v_bar[s]
                    push!(df, [instance_name, s, scenario_id,
                        avg_perc_obs[s], avg_perc_min[s],
                        avg_perc_max[s], perc_energy[s]])
                    avg_perc_min_allsc[s] += avg_perc_min[s]
                    avg_perc_max_allsc[s] += avg_perc_max[s]
                    avg_perc_obs_allsc[s] += avg_perc_obs[s]
                    sum_energy_allsc[s] += perc_energy[s]
                end
            end
            for s in 1:nbNDU
                avg_perc_min_allsc[s] /= nbScenarios
                avg_perc_max_allsc[s] /= nbScenarios
                avg_perc_obs_allsc[s] /= nbScenarios
                sum_energy_allsc[s] /= nbScenarios
                push!(summary_df, [instance_name, s,
                    avg_perc_obs_allsc[s], avg_perc_min_allsc[s],
                    avg_perc_max_allsc[s], sum_energy_allsc[s]])
            end
        catch e
            bt = catch_backtrace()
            msg = sprint(showerror, e, bt)
            println(msg)
        end
    end
    println("Analysis done.")
    println("\nSaving full CSV file to $(output_file_df)...")
    CSV.write(output_file_df, df)
    println("\nSaving summary CSV file to $(output_file_df_summary)...")
    CSV.write(output_file_df_summary, summary_df)
    println("\nDone.")
    flush(stdout)
end

function test()
    #base_folder = toy_instances_folder
    base_folder = antoine_instances_folder
    files_to_process = readdir(base_folder)
    tamanho = size(files_to_process)
    test_set = files_to_process
    compute_statistics(base_folder, test_set)
end

test()
