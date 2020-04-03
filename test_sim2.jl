include("RCCP_Robust_Simulation.jl")

function test_simulation()
    # run simulation tests
    # test_solve_robust_from_t_2()
    #test_set = ["A_instance2_11scen_1contract.txt"]  # "A_instance2_11scen_1NDU.txt" , "A_instance2_11scen.txt"]  # "A_instance2_11scen_1contract.txt"
    #test_set = ["A_instance2_11scen_1NDU.txt"] # ["A_instance2_11scen_noNDU.txt", "A_instance2_11scen_noNDU.txt"]
    test_set = ["A_instance2_100s_skewed-left.txt", "A_instance2_100s_uniform.txt", "A_instance2_100s_skewed-right.txt"]
    #test_set = ["A_instance2_100s.txt"]
    #test_set = ["A_instance2_1NDU_100s.txt"]

    base_folder = pwd() * "/instances/"
    #base_folder = pwd() * "/../notebooks/data/antoine/"
    #base_folder = pwd() * "/instance_gen/output/scenarios/"
    cplex_old = true
    scenario_filter = [x for x in 0:10]
    simulate_antoine_instance(base_folder, test_set, cplex_old, scenario_filter)
end

test_simulation()
