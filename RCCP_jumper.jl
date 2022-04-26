#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# example/inventory.jl
# Implements the inventory management example from
#   A. Ben-Tal, A. Goryashko, E. Guslitzer, A. Nemirovski
#   "Adjustable Robust Solutions of Uncertain Linear Programs"
# Requires a linear optimization solver.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Gurobi
using CSV
using DataFrames
#using DataArrays

# Include file reader and util functions
include("RCCP_FileReader.jl")

function solve_robust_model_with_t(instance_as_dict::Dict, general_logger, infeasible_logger, solver_params,
                                            t0 = 1, fixed_contracts = Int64[], initial_battery = Float64[],
                                            previous_drivable_charge = Float64[])
    verbose = false
    instance = instance_as_dict["instance"]
    period = instance_as_dict["period"]
    contract = instance_as_dict["contract"]
    drivable = instance_as_dict["drivable"]
    n_drivable = instance_as_dict["n_drivable"]
    n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
    storage = instance_as_dict["storage"]
    scenarios = instance_as_dict["scenarios"]
    filepath = instance_as_dict["filepath"]
    # Rounding error allowed in RTCS Operations (e.g. battery and drivable storage / retrieval)
    RTCS_ROUNDING_ERROR = solver_params["rtcs-rounding-error"]

        # Instantiating the model indices
    println("Indices: $(instance)")
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbND = size(n_drivable, 1)
    T = t0:nbT
    C = 1:nbC
    D = 1:nbD
    ND = 1:nbND
    ST = 1:nbSt
    # T = {0, . . . , t̄ − 1} => 1:nbT
    # t̄ (t_bar) follows the last period index
    t_bar = nbT + 1
    # Uncertain indices: TODO revisar
    nbNDU = size(n_drivable_uncertain, 1)
    nbTau = nbT # TODO Validar isso!
    SetTau = t0:nbTau
    SetSigma = 1:nbNDU
    println("Number of periods is $(nbT)")
    println("Solve starting at period $(t0)")
    println("Number of uncertain devices is $(nbNDU)")

    # Setup robust model
    mdet = RobustModel(solver=GurobiSolver(MIPGap=0.005))

    v_min = zeros(Float64, nbTau, nbNDU)
    v_max = zeros(Float64, nbTau, nbNDU)
    for sigma in SetSigma
        for tau in SetTau
            p_min = n_drivable_uncertain[sigma, :Pmin][tau]
            p_max = n_drivable_uncertain[sigma, :Pmax][tau]
            v_min[tau,sigma] = p_min
            v_max[tau,sigma] = p_max
        end
    end
    # Uncertain parameter: power at each time stage lies in a interval
    @uncertain(mdet, v_min[tau,sigma] <= P_NDU[tau=SetTau, sigma=SetSigma] <= v_max[tau,sigma])


    # find out the number of contracts in each time period
    num_contracts = zeros(Int64, nbT)
    for c in C
        t = contract[c,:period]
        num_contracts[t] += 1
    end
    max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
    pi_minus = zeros(Float64, nbT, max_contracts_per_period)
    pi_plus = zeros(Float64, nbT, max_contracts_per_period)
    for t in 1:nbT
        contracts_in_period_t = contract[contract[!, :period] .== t, :]
        for c in 1:num_contracts[t]
            pi_minus[t, c] = contracts_in_period_t[!, :min_period][c]
            pi_plus[t, c] = contracts_in_period_t[!, :max_period][c]
        end
    end

    # (12) 1 if engages contract, 0 otherwise
    @variable(mdet, y[t=T, c=1:num_contracts[t]], Bin)
    # Contract fixation : if parameter fixed_contracts is given
    if size(fixed_contracts, 1) > 0
        engaged_v = Int64[]
        for t in 1:t0-1, c in 1:num_contracts[t]
            push!(engaged_v, fixed_contracts[t, c])
        end
        @constraint(mdet, fix_c[t=t0:nbT,c=1:num_contracts[t]], y[t, c] == fixed_contracts[t, c])
        for t in t0:nbT
            for c in 1:num_contracts[t]
                #@constraint(mdet, y[t, c] == fixed_contracts[t, c])
                push!(engaged_v, fixed_contracts[t, c])
            end
        end
        println("Fixing contract engagement variables to $(engaged_v)")
    end

    # amount consumed (< 0) / provided (> 0) by the client in contract i
    #@variable(mdet, q[t=T, c=1:num_contracts[t]])
    @adaptive(mdet, q[t=T, c=1:max_contracts_per_period], policy=Affine, depends_on=P_NDU[SetTau,SetSigma])
    #num_contracts[t]

    # (15) percentage of time period p in which drivable system s is turned on
    @adaptive(mdet, 0 <= x[T, D] <= 1, policy=Affine, depends_on=P_NDU[SetTau,SetSigma])
    # (13) amount stored in system s at the beginning of time period t
    @adaptive(mdet, r[t0:nbT+1, ST] >= 0, policy=Affine, depends_on=P_NDU[t0:nbT-1,SetSigma])

    # Initial battery storage fixation : if parameter initial_battery is given
    if size(initial_battery, 1) > 0
        for s in ST
            # treat rounding errors of battery levels
            if (initial_battery[s] - storage[s,:uMin] < -EPS)
                println("WARN: assert initial_battery[s] >= storage[s,:uMin] : $(initial_battery[s]) < $(storage[s,:uMin])")
                println(general_logger, "[Det Model, t0 = $(t0)] WARN: assert initial_battery[s] >= storage[s,:uMin] : $(initial_battery[s]) < $(storage[s,:uMin])")
                if abs(initial_battery[s] - storage[s,:uMin]) / storage[s,:uMin] <= RTCS_ROUNDING_ERROR
                    initial_battery[s] = storage[s,:uMin]
                else
                    println(general_logger, "[Det Model, t0 = $(t0)] initial_battery[s] must be >= storage[s,:uMin]. Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                    flush(general_logger)
                    error("initial_battery[s] must be >= storage[s,:uMin]. Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                end
            end
            if (initial_battery[s] - storage[s,:uMax] > EPS)
                println("WARN: assert initial_battery[s] <= storage[s,:uMax] : $(initial_battery[s]) > $(storage[s,:uMax])")
                println(general_logger, "[Det Model, t0 = $(t0)] WARN: assert initial_battery[s] <= storage[s,:uMax] : $(initial_battery[s]) > $(storage[s,:uMax])")
                if abs(initial_battery[s] - storage[s,:uMax]) / storage[s,:uMax] <= RTCS_ROUNDING_ERROR
                    initial_battery[s] = storage[s,:uMax]
                else
                    println(general_logger, "[Det Model, t0 = $(t0)] initial_battery[s] must be <= storage[s,:uMax] Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                    flush(general_logger)
                    error("initial_battery[s] must be <= storage[s,:uMax] Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                end
            end
            @constraint(mdet, r[t0, s] == initial_battery[s])
        end
        println("Fixing battery level variables to $(initial_battery)")
    else
        @constraint(mdet, [s = ST], r[t0, s] == storage[s,:uInit])
    end

    # Initial drivable charge fixation
    if size(previous_drivable_charge, 1) > 0
        println("Fixing drivable initial charge to $(previous_drivable_charge)")
    end
    # (13) amount absorbed by system s during the time period t and that will be already stored in system s
    # at the end of this time period
    @adaptive(mdet, g[T, ST] >= 0, policy=Affine, depends_on=P_NDU[SetTau,SetSigma])
    # (13) amount refunded by system s during time period t and that was stored in system s
    # at the beginning of this time period
    @adaptive(mdet, h[T, ST] >= 0, policy=Affine, depends_on=P_NDU[SetTau,SetSigma])
    # (14) extra amount of electricity requested by the client to the parter
    # (out of any engaged contract) in order to satisfy his needs at time period t
    @adaptive(mdet, e[T] >= 0, policy=Affine, depends_on=P_NDU[SetTau,SetSigma])

    # Problem constraints
    # (2) Electricity balance at each time period
    @constraint(mdet, c2[t = T], sum(q[t, c] for c in 1:num_contracts[t])
                        + sum(n_drivable[s,:pORc][t] for s in ND) + sum(P_NDU[t, sigma] for sigma in SetSigma)
                        + sum(drivable[s,:pORc][t] * x[t, s] for s in D)
                        + sum(h[t, s] for s in ST) - sum(g[t, s] for s in ST) + e[t] >= 0)

    # (3) The amount consumed from system s during a time period t must be at most the amount stored
    @constraint(mdet, c3[t = T, s = ST], h[t, s] - r[t, s] <= 0 )
    # (4) and (5) Storage capacity of systems must be respected
    # NOTE : Initial battery storage is defined above! See variable declaration section.
    @constraint(mdet, c4[s = ST, t = (t0):(nbT+1)], r[t, s] >= storage[s,:uMin])
    @constraint(mdet, c5[s = ST, t = (t0):(nbT+1)], r[t, s] <= storage[s,:uMax])
    # (6) Electricity stored on system s at the next time period (taking into account the loss coefficient of the system)
    @constraint(mdet, c6[t = T, s = ST], r[t+1, s] - r[t, s] + h[t, s] - storage[s,:lostCoef] * g[t, s] == 0)
    # (7) Maximum quantity of energy that can be absorbed by a drivable system during a time period
    @constraint(mdet, c7[t = T, s = ST], g[t, s] - storage[s,:maxAbsorption] * period[t,:size] <= 0)
    # (8) Maximum quantity of energy that can be refunded by a drivable system during a time period
    @constraint(mdet, c8[t = T, s = ST], h[t, s] - storage[s,:maxRefund] * period[t,:size] <= 0)
    # (9) Consumption/production restrictions of contracts must be respected:
    buy_ct = Array[]
    sell_ct = Array[]
    for t in T, c in 1:num_contracts[t]
        if pi_plus[t, c] > 0
            push!(buy_ct, [t, c])
        elseif pi_plus[t, c] < 0
            push!(sell_ct, [t, c])
        end
    end
    # Buy contracts constraints
    # Sell contracts constraints
    min_period = zeros(Float64, nbT, max_contracts_per_period)
    max_period = zeros(Float64, nbT, max_contracts_per_period)
    for t in T
        c_t = 1
        for c in C
            if contract[c,:period] == t
                min_period[t, c_t] = contract[c,:min_period]
                max_period[t, c_t] = contract[c,:max_period]
                c_t += 1
            end
        end
    end
    for t in T, c in 1:num_contracts[t]
        if pi_plus[t, c] > 0
            @constraint(mdet, q[t, c] <= max_period[t, c] * y[t, c])
            @constraint(mdet, q[t, c] >= min_period[t, c] * y[t, c])
        else
            @constraint(mdet, q[t, c] >= max_period[t, c] * y[t, c])
            @constraint(mdet, q[t, c] <= min_period[t, c] * y[t, c])
        end
    end
    # (10) and (11) Minimum usage of drivable system s:
    for t in T, s in D
        #println("t = $(t), s = $(s): $(drivable[s,Symbol("pORc_min")][t])")
        previous_charge = 0.0
        if size(previous_drivable_charge, 1) > 0
            previous_charge = previous_drivable_charge[s]
        end
        if drivable[s,:pORc][t] > 0 || drivable[s,:pORc_min][t] > 0 ||  drivable[s,:pORc_max][t] > 0   # P_D[t,s] > 0
            if drivable[s,:pORc_min][t] > drivable[s,:pORc_max][t]
                println(general_logger, "[Det Model, t0 = $(t0)] Assert failed : [drivable D] P > 0 => Pmin <= Pmax")
            end
            @assert drivable[s,:pORc_min][t] <= drivable[s,:pORc_max][t] "[drivable D] P > 0 => Pmin <= Pmax"
            @constraint(mdet, previous_charge + sum(drivable[s,:pORc][t1] * x[t1, s] for t1 in T if t1 <= t) - drivable[s,:pORc_min][t] >= 0 )
        elseif drivable[s,:pORc][t] < 0 || drivable[s,:pORc_min][t] < 0 || drivable[s,:pORc_max][t] < 0
            if drivable[s,:pORc_min][t] < drivable[s,:pORc_max][t]
                println(general_logger, "[Det Model, t0 = $(t0)] Assert failed : [drivable D] P < 0 => Pmin >= Pmax")
            end
            @assert drivable[s,:pORc_min][t] >= drivable[s,:pORc_max][t] "[drivable D] P < 0 => Pmin >= Pmax"
            @constraint(mdet, previous_charge + sum(drivable[s,:pORc][t1] * x[t1, s] for t1 in T if t1 <= t) - drivable[s,:pORc_min][t] <= 0 )
        end
    end

    # Objective function
    cost_fix = zeros(Float64, nbT, max_contracts_per_period)
    cost_var = zeros(Float64, nbT, max_contracts_per_period)
    for t in T
        c_t = 1
        for c in C
            if contract[c,:period] == t
                cost_fix[t, c_t] = contract[c,:cost_fix]
                cost_var[t, c_t] = contract[c,:cost_var]
                c_t += 1
            end
        end
    end
    # Objective: minimize total cost of production
    v_NDU = zeros(Float64, nbNDU)
    for sigma in SetSigma
       v_NDU[sigma] = n_drivable_uncertain[sigma, :cost]
    end
    @variable(mdet, F)  # Overall cost
    @objective(mdet, Min, F)
    @constraint(mdet, F >= sum( cost_fix[t, c] * y[t, c] + cost_var[t, c] * q[t, c] for t in T, c in 1:num_contracts[t] ) +
        sum(drivable[s,:cost] * drivable[s,:pORc][t] * x[t, s] for t in T, s in D) +
        sum(storage[s,:cost] * (g[t, s] + h[t, s]) for t in T, s in ST) +
        sum(period[t,:cost_out] * e[t] for t in T) + sum(v_NDU[s] * P_NDU[t, s] for t in T, s in SetSigma))

    # Solve the problem
    println("===================================================================\nSolving deterministic model for t0 = $(t0)...")
    status = solve(mdet)
    # Print variables
    solution = Dict()
    sol_y = zeros(Int64, nbT, max_contracts_per_period)
    max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
    if status == :Optimal
        println("\n===========  R O B U S T    S O L U T I O N  ( t0 = $(t0) ) ===========\n")
        println("Optimal Objective Function value: ", getobjectivevalue(mdet))
        println("Solve time : ", getsolvetime(mdet))
        var_y = getvalue(y)
        if t0 > 1  # if starting from period > 1
            for t in 1:(t0-1), c in 1:num_contracts[t]
                if size(fixed_contracts, 1) > 0
                    sol_y[t,c] = fixed_contracts[t,c]
                else
                    sol_y[t,c] = -1
                    println("WARN: Returning invalid contract solution y. Check RCCP Deterministic Model procedure !")
                    println(general_logger, "[Det Model, t0 = $(t0)] WARN: Returning invalid contract solution y. Check RCCP Deterministic Model procedure !")
                end
            end
        end
        for t in T, c in 1:num_contracts[t]
            sol_y[t,c] = trunc(Int, getvalue(y[t, c]))
        end
        solution["y"] = sol_y
        return getobjectivevalue(mdet), solution, getsolvetime(mdet), (status == :Optimal ? true : false)
    else
        #env = CPLEX.Env()
        #CPLEX.set_logfile(env, "cplex.log")
        #println("$(CPLEX.getIIS())")
        #CPLEX.close_CPLEX(env)
        println(infeasible_logger, "t0 = $(t0); Infeasible det model : $(filepath)")
        println("Optimal solution not found! Best bound: ", getobjectivebound(mdet))
        println(general_logger, "[Det Model, t0 = $(t0)] WARN: Optimal solution not found! Best bound: ", getobjectivebound(mdet))
        solution["y"] = sol_y
        return getobjectivebound(mdet), solution, getsolvetime(mdet), (status == :Optimal ? true : false)
    end
end

datafile = "../notebooks/data/antoine/A_instance2_11scen.txt"
instance_as_dict = read_input_data(datafile, false)
general_logger = open(joinpath(tempdir(), "rccp_concat_log.txt"), "w+")
infeasible_logger = open(joinpath(tempdir(), "rccp_concat_log_infeasible.txt"), "w+")
solve_robust_model_with_t(instance_as_dict, general_logger, infeasible_logger)
close(general_logger)
