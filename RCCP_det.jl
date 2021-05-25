# =================================
# RCCP_det.jl
# CCP Deterministic MILP Model
# =================================

using JuMP
using CPLEX
using CSV
using DataFrames
#using DataArrays

# Include file reader and util functions
include("RCCP_FileReader.jl")

function solve_deterministic_model(filename::String, solver_params)
    instance_as_dict = read_input_data(filename, false)
    return solve_deterministic_model(instance_as_dict, solver_params)
end

function solve_deterministic_model(instance_as_dict::Dict, solver_params)
    instance = instance_as_dict["instance"]
    period = instance_as_dict["period"]
    contract = instance_as_dict["contract"]
    drivable = instance_as_dict["drivable"]
    n_drivable = instance_as_dict["n_drivable"]
    n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
    storage = instance_as_dict["storage"]
    scenarios = instance_as_dict["scenarios"]

    println("\n===========  D E T E R M I N I S T I C    M O D E L  ( t0 = 1 ) ===========\n")
    println("\nSolving CCP Deterministic model for t0 = 1...")
    # If there is any uncertain non-drivable device, include these devices to non-drivable list
    #   such that P = (Pmin + Pmax)/2
    if !isempty(n_drivable_uncertain)
        println("\nWARN: Including uncertain non-drivable devices as normal non-drivable using average.")
        # TODO FIXME fix the following line, uncommenting the code of size(...)
        for i in 1:size(n_drivable_uncertain, 1)
            power = (n_drivable_uncertain[i, :Pmin] + n_drivable_uncertain[i, :Pmax]) / 2.0
            push!(n_drivable, [n_drivable_uncertain[i, :name], n_drivable_uncertain[i, :cost], power])
        end
        println("\nNEW Non Drivable dataframe:")
        showall(n_drivable)
    end

    solver_parameters_2 = deepcopy(solver_params)
    solver_parameters_2["relative-gap-tol"] = 1e-04  # Set relative gap tolerance to default solver value
    mdet = create_jump_model(solver_parameters_2)
    # Instantiating the model indices
    println("Indices: $(instance)")
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbND = size(n_drivable, 1)
    T = 1:nbT
    C = 1:nbC
    D = 1:nbD
    ND = 1:nbND
    ST = 1:nbSt

    # Instantiating the decision variables
    # (12) 1 if engages contract, 0 otherwise
    @variable(mdet, y[C], Bin)
    # amount consumed (< 0) / provided (> 0) by the client in contract i
    @variable(mdet, q[C])
    # (15) percentage of time period p in which drivable system s is turned on
    @variable(mdet, 0 <= x[T, D] <= 1)
    # (13) amount stored in system s at the beginning of time period t
    @variable(mdet, r[1:nbT+1, ST] >= 0)
    # (13) amount absorbed by system s during the time period t and that will be already stored in system s
    # at the end of this time period
    @variable(mdet, g[T, ST] >= 0)
    # (13) amount refunded by system s during time period t and that was stored in system s
    # at the beginning of this time period
    @variable(mdet, h[T, ST] >= 0)
    # (14) extra amount of electricity requested by the client to the parter
    # (out of any engaged contract) in order to satisfy his needs at time period t
    @variable(mdet, e[T] >= 0)

    # Problem constraints
    # (2) Electricity balance at each time period
    for t in T
        #println("t is $(t):")
        #for s in D
        #    println("s = $(s): $(n_drivable[s,:pORc][t]) $(drivable[s,:pORc][t])")
        #end
        @constraint(mdet, sum(q[c] for c in C if contract[c,:period] == t)
                        + sum(n_drivable[s,:pORc][t] for s in ND)
                        + sum(drivable[s,:pORc][t] * x[t, s] for s in D)
                        + sum(storage[s,:lostCoef] * h[t, s] for s in ST) - sum(g[t, s] for s in ST) + e[t] >= 0)
    # FIXME: por que o balanco tem que >= 0 e nao exatamente == 0 ? Por conta dos sistemas de dissipacao?
    # FIXME: Isso pode gerar distorcoes do tipo: a cada t, o sistema vai enviar o maximo de energia (q[c])
    # FIXME: para a empresa eletrica, independente do que ele eh capaz de gerar...
    end
    # (3) The amount consumed from system s during a time period t must be at most the amount stored
    for t in T
        for s in ST
            @constraint(mdet, h[t, s] - r[t, s] <= 0 )
        end
    end
    # (4) and (5) Storage capacity of systems must be respected
    for s in ST
        @constraint(mdet, r[1, s] == storage[s,:uInit])
        for t in 1:nbT+1  # FIXME
            @constraint(mdet, r[t, s] >= storage[s,:uMin])
            @constraint(mdet, r[t, s] <= storage[s,:uMax])
        end
    end
    # (6) Electricity stored on system s at the next time period (taking into account the loss coefficient of the system)
    for t in T
        for s in ST
            @constraint(mdet, r[t+1, s] - r[t, s] + h[t, s] - storage[s,:lostCoef] * g[t, s] == 0)
        end
    end
    # (7) Maximum quantity of energy that can be absorbed by a drivable system during a time period
    for t in T
        for s in ST
            @constraint(mdet, g[t, s] - storage[s,:maxAbsorption] * period[t,:size] <= 0)
        end
    end
    # (8) Maximum quantity of energy that can be refunded by a drivable system during a time period
    for t in T
        for s in ST
            @constraint(mdet, h[t, s] - storage[s,:maxRefund] * period[t,:size] <= 0)
        end
    end
    # (9) Consumption/production restrictions of contracts must be respected:
    # FIXME duvida: o que limita o valor produzido de energia q[c] (eg. devolvido a empresa eletrica)?
    # FIXME quando limitamos, por exemplo, o valor q[c] = 1000, por meio do max_period, estabelecemos um limite do
    # FIXME contrato, mas algo limita a quantidade de energia produzida pelo consumidor (e.g. painel eletrico)?
    for t in T
        for c in C
            if contract[c,:period] == t
                if contract[c,:min_period] >= 0
                    @constraint(mdet, q[c] <= contract[c,:max_period] * y[c])
                    @constraint(mdet, q[c] >= contract[c,:min_period] * y[c])
                else
                    @constraint(mdet, q[c] >= contract[c,:max_period] * y[c])
                    @constraint(mdet, q[c] <= contract[c,:min_period] * y[c])
                end
            end
        end
    end
    # (10) and (11) Minimum usage of drivable system s:
    for t in T
        for s in D
            #println("t = $(t), s = $(s): $(drivable[s,Symbol("pORc_min")][t])")
            if drivable[s,:pORc][t] > 0 || drivable[s,:pORc_min][t] > 0 ||  drivable[s,:pORc_max][t] > 0
                @constraint(mdet, sum(drivable[s,:pORc][t1] * x[t1, s] for t1 in T if t1 <= t) - drivable[s,:pORc_min][t] >= 0 )
            elseif drivable[s,:pORc][t] < 0 || drivable[s,:pORc_min][t] < 0 || drivable[s,:pORc_max][t] < 0
                @constraint(mdet, sum(drivable[s,:pORc][t1] * x[t1, s] for t1 in T if t1 <= t) - drivable[s,:pORc_min][t] <= 0 )
            end
        end
    end

    # Objective function
    @objective(mdet, Min, sum( contract[c,:cost_fix] * y[c] + contract[c,:cost_var] * q[c] for c in C) +
        sum(drivable[s,:cost] * drivable[s,:pORc][t] * x[t, s] for t in T, s in D) +
        sum(storage[s,:cost] * (g[t, s] + h[t, s]) for t in T, s in ST) +
        sum(period[t,:cost_out] * e[t] for t in T) )

    #println(mdet)

    # Solve the problem
    status = optimize!(mdet)
    # Print variables
    solution = Dict()
    if JuMP.termination_status(mdet) == MOI.OPTIMAL
        println("\n===========  S O L U T I O N ===========\n")
        println("Optimal Objective Function value: ", objective_value(mdet))
        println("Solve time : ", solve_time(mdet))
        sol_y = zeros(Int64, nbC)
        for c in C
            sol_y[c] = trunc(Int, value(y[c]))
        end
        solution["y"] = sol_y
        if verbose
            println("Optimal Solutions:")
            println("x = ", value.(x))
            println("y = ", value.(y))
            println("q = ", value.(q))
            println("r = ", value.(r))
            println("g = ", value.(g))
            println("h = ", value.(h))
            println("e = ", value.(e))
        end
        return objective_value(mdet), solution, solve_time(mdet), true
    else
        println("Optimal solution not found! Best bound: ", objective_bound(mdet))
        return objective_bound(mdet), solution, solve_time(mdet), false
    end
end

# Solve the Deterministic RCCP MILP Model starting from period t0
# Fix the contracts y according to fixed_contracts parameter, if parameter is given
# Fix the initial battery levels (parameter "initial_battery"), if parameter is given
function solve_deterministic_model_with_t(instance_as_dict::Dict, general_logger, infeasible_logger, solver_params, t0 = 1,
                                            fixed_contracts = Int64[], initial_battery = Float64[],
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

    # If there is any uncertain non-drivable device, include these devices to non-drivable list
    #   such that P = (Pmin + Pmax)/2
    if !isempty(n_drivable_uncertain)
        println("\nWARN: Including uncertain non-drivable devices as normal non-drivable using average.")
        # TODO FIXME fix the following line, uncommenting the code of size(...)
        for i in 1:size(n_drivable_uncertain, 1)
            power = (n_drivable_uncertain[i, :Pmin] + n_drivable_uncertain[i, :Pmax]) / 2.0
            push!(n_drivable, [n_drivable_uncertain[i, :name], n_drivable_uncertain[i, :cost], power])
        end
        if verbose
            println("\nNEW Non Drivable dataframe:")
            showall(n_drivable)
        end
    end

    # Instantiating the model indices
    println("Indices: $(instance)")
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
    nbND = size(n_drivable, 1)
    T = t0:nbT
    C = 1:nbC
    D = 1:nbD
    ND = 1:nbND
    ST = 1:nbSt
    println("\n===========  D E T E R M I N I S T I C    M O D E L  ( t0 = $(t0) ) ===========\n")
    println("Total number of periods is $(nbT)")
    println("Solve starting at period $(t0)")

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

    solver_parameters_2 = deepcopy(solver_params)
    solver_parameters_2["relative-gap-tol"] = 1e-04  # Set relative gap tolerance to default solver value
    mdet = create_jump_model(solver_parameters_2)

    # Instantiating the decision variables
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
    @variable(mdet, q[t=T, c=1:num_contracts[t]])
    # (15) percentage of time period p in which drivable system s is turned on
    @variable(mdet, 0 <= x[T, D] <= 1)
    # (13) amount stored in system s at the beginning of time period t
    @variable(mdet, r[t0:nbT+1, ST] >= 0)
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
        @constraint(mdet, init_batt[s = ST], r[t0, s] == storage[s,:uInit])
    end
    # Initial drivable charge fixation
    if size(previous_drivable_charge, 1) > 0
        println("Fixing drivable initial charge to $(previous_drivable_charge)")
    end
    # (13) amount absorbed by system s during the time period t and that will be already stored in system s
    # at the end of this time period
    @variable(mdet, g[T, ST] >= 0)
    # (13) amount refunded by system s during time period t and that was stored in system s
    # at the beginning of this time period
    @variable(mdet, h[T, ST] >= 0)
    # (14) extra amount of electricity requested by the client to the parter
    # (out of any engaged contract) in order to satisfy his needs at time period t
    @variable(mdet, e[T] >= 0)

    # Problem constraints
    # (2) Electricity balance at each time period
    @constraint(mdet, c2[t = T], sum(q[t, c] for c in 1:num_contracts[t])
                        + sum(n_drivable[s,:pORc][t] for s in ND)
                        + sum(drivable[s,:pORc][t] * x[t, s] for s in D)
                        + sum(storage[s,:lostCoef] * h[t, s] for s in ST) - sum(g[t, s] for s in ST) + e[t] >= 0)
    # FIXME: por que o balanco tem que >= 0 e nao exatamente == 0 ? Por conta dos sistemas de dissipacao?
    # FIXME: Isso pode gerar distorcoes do tipo: a cada t, o sistema vai enviar o maximo de energia (q[c])
    # FIXME: para a empresa eletrica, independente do que ele eh capaz de gerar...

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
    # FIXME duvida: o que limita o valor produzido de energia q[c] (eg. devolvido a empresa eletrica)?
    # FIXME quando limitamos, por exemplo, o valor q[c] = 1000, por meio do max_period, estabelecemos um limite do
    # FIXME contrato, mas algo limita a quantidade de energia produzida pelo consumidor (e.g. painel eletrico)?
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
#    @constraint(mdet, c9a[x = 1:size(buy_ct, 1)], q[buy_ct[x][1], buy_ct[x][2]] <= contract[contract[!, :period] .== buy_ct[x][1], :max_period][buy_ct[x][2]] * y[buy_ct[x][1], buy_ct[x][2]])
#    @constraint(mdet, c9b[x = 1:size(buy_ct, 1)], q[buy_ct[x][1], buy_ct[x][2]] >= contract[contract[!, :period] .== buy_ct[x][1], :min_period][buy_ct[x][2]] * y[buy_ct[x][1], buy_ct[x][2]])
    # Sell contracts constraints
#    @constraint(mdet, c9c[x = 1:size(sell_ct, 1)], q[buy_ct[x][1], buy_ct[x][2]] >= contract[contract[!, :period] .== buy_ct[x][1], :max_period][buy_ct[x][2]] * y[buy_ct[x][1], buy_ct[x][2]])
#    @constraint(mdet, c9d[x = 1:size(sell_ct, 1)], q[buy_ct[x][1], buy_ct[x][2]] <= contract[contract[!, :period] .== buy_ct[x][1], :min_period][buy_ct[x][2]] * y[buy_ct[x][1], buy_ct[x][2]])
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
    @objective(mdet, Min, sum( cost_fix[t, c] * y[t, c] + cost_var[t, c] * q[t, c] for t in T, c in 1:num_contracts[t] ) +
        sum(drivable[s,:cost] * drivable[s,:pORc][t] * x[t, s] for t in T, s in D) +
        sum(storage[s,:cost] * (g[t, s] + h[t, s]) for t in T, s in ST) +
        sum(period[t,:cost_out] * e[t] for t in T) )

    #println(mdet)

    # Solve the problem
    println("===================================================================\nSolving deterministic model for t0 = $(t0)...")
    # Save model to LP file
    #output_path = joinpath(normpath(EXPERIMENT_OUTPUT_FOLDER), "output", "models", "rccp_det_t" * string(t0) * ".lp")
    #writeLP(mdet, output_path)
    status = optimize!(mdet)
    # Print variables
    solution = Dict()
    sol_y = zeros(Int64, nbT, max_contracts_per_period)
    max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
    if JuMP.termination_status(mdet) == MOI.OPTIMAL
        rel_gap = MOI.get(mdet, MOI.RelativeGap())
        println("\n===========  D E T E R M I N I S T I C    S O L U T I O N  ( t0 = $(t0) ) ===========\n")
        println("Optimal Objective Function value: ", objective_value(mdet))
        println("Relative gap = $(rel_gap)")
        println("Solve time : ", solve_time(mdet))
        var_y = value.(y)
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
            sol_y[t,c] = trunc(Int, value(y[t, c]))
        end
        solution["y"] = sol_y
        # println("\namount consumed (< 0) / provided (> 0) by the client in contract i")
        solution["q"] = capture_variable_solution_as_array_2d_custom_q(q, t0, nbT, num_contracts)
        # println("\namount stored in system s at the beginning of time period t ")
        solution["r"] = capture_variable_solution_as_array_2d_custom_r(r, t0, nbT, ST)
        # println("\namount absorbed by system s during the time period t and that will be already stored in system s at the end of this time period")
        solution["g"] = capture_variable_solution_as_array_2d(g, t0, nbT, ST)
        # amount refunded by system s during time period t and that was stored in system s at the beginning of this time period
        solution["h"] = capture_variable_solution_as_array_2d(h, t0, nbT, ST)
        # println("\nextra amount of electricity requested by the client to the parter (out of any engaged contract) in order to satisfy his needs at time period t")
        solution["e"] = capture_variable_solution_as_array_1d(e, t0, nbT)
        # println("\npercentage of time period p in which drivable system s is turned on")
        solution["x"] = capture_variable_solution_as_array_2d(x, t0, nbT, D)
        if verbose
            println("Optimal Solutions:")
            println("x = ", value.(x))
            println("y = ", value.(y))
            println("q = ", value.(q))
            println("r = ", value.(r))
            println("g = ", value.(g))
            println("h = ", value.(h))
            println("e = ", value.(e))
            flush(stdout)
        end
        return objective_value(mdet), solution, solve_time(mdet), true, rel_gap
    elseif has_values(mdet) && (JuMP.termination_status(mdet) in [MOI.TIME_LIMIT, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.LOCALLY_SOLVED])
        # Recover the best integer solution found (suboptimal) and return it
        rel_gap = MOI.get(mdet, MOI.RelativeGap())
        println("\n===========  D E T E R M I N I S T I C    S O L U T I O N  ( t0 = $(t0) ) ===========\n")
        println("SubOptimal Objective Function value: ", objective_value(mdet))
        println("Relative gap = $(rel_gap)")
        println("Solve time : ", solve_time(mdet))
        println(general_logger, "[Deterministic Model, t0 = $(t0)] WARN: Time limit exceeded. Optimal solution not found! Best bound: ", objective_bound(mdet))
        var_y = value.(y)
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
            sol_y[t,c] = trunc(Int, value(y[t, c]))
        end
        solution["y"] = sol_y
        # println("\namount consumed (< 0) / provided (> 0) by the client in contract i")
        solution["q"] = capture_variable_solution_as_array_2d_custom_q(q, t0, nbT, num_contracts)
        # println("\namount stored in system s at the beginning of time period t ")
        solution["r"] = capture_variable_solution_as_array_2d_custom_r(r, t0, nbT, ST)
        # println("\namount absorbed by system s during the time period t and that will be already stored in system s at the end of this time period")
        solution["g"] = capture_variable_solution_as_array_2d(g, t0, nbT, ST)
        # amount refunded by system s during time period t and that was stored in system s at the beginning of this time period
        solution["h"] = capture_variable_solution_as_array_2d(h, t0, nbT, ST)
        # println("\nextra amount of electricity requested by the client to the parter (out of any engaged contract) in order to satisfy his needs at time period t")
        solution["e"] = capture_variable_solution_as_array_1d(e, t0, nbT)
        # println("\npercentage of time period p in which drivable system s is turned on")
        solution["x"] = capture_variable_solution_as_array_2d(x, t0, nbT, D)
        if verbose
            println("SubOptimal Solutions:")
            println("x = ", value.(x))
            println("y = ", value.(y))
            println("q = ", value.(q))
            println("r = ", value.(r))
            println("g = ", value.(g))
            println("h = ", value.(h))
            println("e = ", value.(e))
            flush(stdout)
        end
        return objective_value(mdet), solution, solve_time(mdet), true, rel_gap
    else  # Infeasible
        #env = CPLEX.Env()
        #CPLEX.set_logfile(env, "cplex.log")
        #println("$(CPLEX.getIIS())")
        #CPLEX.close_CPLEX(env)
        println(infeasible_logger, "t0 = $(t0); Infeasible det model : $(filepath)")
        println("Optimal solution not found! Best bound: ", objective_bound(mdet))
        println(general_logger, "[Det Model, t0 = $(t0)] WARN: Optimal solution not found! Best bound: ", objective_bound(mdet))
        solution["y"] = sol_y
        return objective_bound(mdet), solution, solve_time(mdet), false, Inf
    end
end

function capture_variable_solution_as_array_1d(var_value, t0, nbT)
    var_sol = [Float64(0.0) for x=1:nbT]
    for x in t0:nbT
        var_sol[x] = trunc_if_less_than_eps(value(var_value[x]))
    end
    return var_sol
end

function capture_variable_solution_as_array_2d(var_value, t0, nbT, range2)
    var_sol = [Float64(0.0) for x=1:nbT, y=range2]
    for x in t0:nbT, y in range2
        var_sol[x, y] = trunc_if_less_than_eps(value(var_value[x, y]))
    end
    return var_sol
end

function capture_variable_solution_as_array_2d_custom_q(var_value, t0, nbT, num_contracts)
    max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
    var_sol = [Float64(0.0) for x=1:nbT, y=1:max_contracts_per_period]
    for x in t0:nbT, y in 1:num_contracts[x]
        var_sol[x, y] = trunc_if_less_than_eps(value(var_value[x, y]))
    end
    return var_sol
end

function capture_variable_solution_as_array_2d_custom_r(var_value, t0, nbT, range2)
    var_sol = [Float64(0.0) for x=1:nbT+1, y=range2]
    for t in t0:nbT+1, y in range2
        var_sol[t, y] = trunc_if_less_than_eps(value(var_value[t, y]))
    end
    return var_sol
end

function calculate_r_td_from_deterministic_r(r, t, nbSt)
    r_td = [Float64(0.0) for x=1:nbSt] # zeros(Float64, nbSt)
    for y in 1:nbSt
        #r_td[y] = r[t, y]
        push!(r_td, r[t, y])
    end
    return r_td
end

# Calculate the value of each deterministic variable (q, x etc.) by dividing by period_size (simple average)
function calculate_independent_term_det_delta(q, x, e, r, g, h, period_size)
    return q ./ period_size, x ./ period_size, e ./ period_size, r ./ period_size, g ./ period_size, h ./ period_size
end

function calculate_variables_from_deterministic(instance_as_dict, det_solution, t)
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance_as_dict["instance"])
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    T, S, C, D, ND, ST, NDU = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt, nbNDU)
    SetSigma = NDU
    num_contracts = instance_as_dict["num_contracts"]
    storage = instance_as_dict["storage"]
    q = det_solution["q"]
    r = det_solution["r"]
    g = det_solution["g"]
    h = det_solution["h"]
    e = det_solution["e"]
    x = det_solution["x"]
    q_t = zeros(Float64, nbC)
    x_t = zeros(Float64, nbD)
    e_t = e[t]
    r_t = zeros(Float64, nbSt)
    g_t = zeros(Float64, nbSt)
    h_t = zeros(Float64, nbSt)
    for s in D  # x[t,s] : percentage of time period t in which drivable system s is turned on
        x_t[s] = x[t,s]
    end
    for c in 1:num_contracts[t]  # Rule 3 -> Average value
        q_t[c] = q[t, c]
    end
    for s in ST  # Rule 4 -> Average values
        u_max = storage[s,:uMax]
        r_t[s] = r[t, s]
        if u_max != 0
            g_t[s] = g[t, s]
            h_t[s] = h[t, s]
        end
    end
    return q_t, x_t, e_t, r_t, g_t, h_t
end
