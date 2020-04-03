================================================
   README for CCP / RTCS Simulation in Julia
================================================

* System requirements:
    - Julia 0.6.4
    - CPLEX 12.6.8 or above
    - The following Julia packages must also be installed: JuMP, Cbc, CSV, DataFrames, CPLEX, JLD2, FileIO, MathProgBase, ZipFile, DataStructures

* Entry point:
    1. In order to select the problem instance to be used in the simulation, please edit the file 'RCCP_Robust_Simulation.jl' and change the function 'run_simulation()', 
        adjusting the 'base_folder' to the path to the instances folder and the 'test_set' array to the list of instance filenames. 
    2. Please run the following Julia file: 'run_sim.jl' (the working folder must be the folder containing the Julia files, i.e. the Julia program must be run from the folder 'ccp_rtcs'). 
       e.g. of command line invocation :~/ccp_rtcs$  julia run_sim.jl
    3. The simulation logs and trace files will be available in the folders 'output/simulation/log' and 'output/simulation/trace', respectively. 
    4. The trace simulation result file is a ZIP file containing the term 'VariableValues' (and the suffix of the instance file: e.g. Antoine_RCCP_VariableValues_A_instance2_1000s_uniform.txt.zip) 
        and will be located in the folder 'output/simulation/trace'. Inside this ZIP file, there is a CSV file containing the simulation trace for each scenario, time period t and time slot d.
    5. In order to process all ZIP trace files and obtain a human-readable simulation results summary, please run the python script 'python/compare_multiple_parameters.py', passing, as argument (--folders), 
        the list of folders containing the trace ZIP files generated in the previous step. The python script will automatically expand the zip files and process their content. 
    6. In the folder 'python/output', the python script will produce 5 CSV files (Det_vs_Rob_Wins_t1... until Det_vs_Rob_Wins_t5...), with different comparisons for each simulated instance file
       (filename column) and simulation parameters (columns Strategy, Reoptimize and ModelPolicy -- see simulation parameters description below).
    7. These produced CSV files include scenario-wise comparisons, such as how many times the robust RTCS obtained smaller solution values (cost), when compared to the deterministic RTCS, given the same scenario. 


* Simulation Parameters: 

RTCS without Reoptimization:
  - Fix the engaged contracts array [y] to the solution provided by the MILP model
    with t_0 = 1. This is the list of the engaged contracts to be used in all time
    periods t in {1, ..., t_bar}.
  - Use the solution values for model variables [q, x, e, r, g, h] provided by
    the MILP model with t_0 = 1.

RTCS with Reoptimization:
  - Fix the engaged contracts array [y] to the solution provided by the MILP model
    with t_0 = 1. This is the list of the engaged contracts to be used in all time
    periods t in {1, ..., t_bar}.
  - Use the solution values for model variables [q, x, e, r, g, h] provided by
    the MILP model with t_0 = t_now.
  - The RTCP with model_policy == "full_model" is the only one which uses the
    reoptimization feature.

Model policies:
- use_model (called 'full' in the paper) : start with model-suggested variable values for q (contract usage) and x (drivable usage) (FULLY TRUST THE MODEL PREDICTIONS)
- ignore_model (called 'simple' in the paper) : start with the minimum quota for each engaged contract ($q[c] = pi_minus_d[c]$) and power required drivable systems,
                    by satisfying their minimum power requirement for the period.

Simulation strategies (called 'gap policy' in the paper):
- Cheapest : Prioritize RTCS operations (contracts, batteries, drivable) according to the cost.

# Simulation strategies for positive gap:
# - cheapest : always execute the cheapest operation first (using either batteries, drivable or contracts);
# - conservative : use batteries first (to store positive gap), then use sell contracts (i.e. try to preserve batteries to worst-case periods);
# - audacious : use sell contracts first, then store positive gap in batteries.

# Simulation strategies for negative gap:
# - cheapest : always execute the cheapest operation first (using either batteries, drivable or contracts);
# - audacious : use batteries first (to compensate for neg gap), then use contracts;
# - conservative : use contracts first, then use batteries (i.e. try to preserve batteries to worst-case periods).



General RTCS Algorithm (for each time period t):
  # Sum energy consumption/production from non-drivable (certain and uncertain)
  gap = sum_NDU + sum_ND

  If model_policy == "full_model":
    - Buy/Sell from contracts according to model solution [q(t)]     => gap += sum_q
    - Send energy to drivable devices according to solution [x(t)]   => gap += sum_x
    - Send/retrieve to/from batteries according to solution [g,h(t)] => gap += sum_gh
    - Update current battery levels
  Else:
    - Buy/Sell the minimum quota for each engaged contract [y(t,c)]  => gap += sum_q
    - TURN ON and POWER required DRIVABLE consuming devices          => gap -= sum_x_cons
  End If

  If Gap > 0:      # energy left over
      If sim_strategy == "cheapest":
          - Execute operations ({STORE in batteries, SELL to contracts, SEND to drivable})
                according to smallest cost first
      Else:
          If sim_strategy == "conservative":
              - STORE energy gap in BATTERIES                          => gap -= sum_g
              - SELL  energy gap to engaged sell contracts             => gap -= sum_q_sell
          Else:  # "audacious"
              - SELL  energy gap to engaged sell contracts             => gap -= sum_q_sell
              - STORE energy gap in batteries                          => gap -= sum_g
          End If
          - SEND  energy gap to DRIVABLE consuming systems             => gap -= sum_x
      End If
      - Throw remaining energy away                                => gap = 0
  Else (Gap < 0):  # lacking energy
      If sim_strategy == "cheapest":
          - Execute operations ({Turn on drivable, RETRIEVE from batteries, BUY from contracts, BUY OUT of contract})
                    according to smallest cost first
      Else:
          - TURN ON  DRIVABLE producing devices                        => gap += sum_x_prod
          If sim_strategy == "audacious":
              - RETRIEVE energy gap from batteries_ST                  => gap += sum_h
              - BUY energy gap from engaged buy contracts              => gap += sum_q_buy
          Else:  "conservative"
              - BUY energy gap from engaged buy contracts              => gap += sum_q_buy
              - RETRIEVE energy gap from batteries_ST                  => gap += sum_h
          End If
      End If
      - BUY remaining gap out of contracts if needed                   => e_t = gap
  End If

