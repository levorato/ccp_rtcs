# ============================================================================
# RCCP_RTCS.jl
# RTCS: Defines a set of actions to be taken by the client system at each time
#       unit d ∈ {1, . . . , δ t } of each time period t ∈ T.
# =============================================================================

# Include file reader and util functions
include("RCCP_FileReader.jl")
# Include RTCS queue and operation functions
include("RCCP_RTCS_Queue.jl")

using ZipFile

# Executes the RTCS for the deterministic model at time period t and time unit d.
# Considering an optimal deterministic solution (y, q, r, g, h, x) for formulation DET-L1,
# we can define the following set of actions to be taken at each time unit.
# Consider a time unit d ∈ {1, . . . , δ t } at a time period t ∈ T . The RTCS follows.
# 1. For each s ∈ SD such that P_ts > 0: s must be (kept) turned on if and only if d ≤ ceil(δ_t * x_ts).
# 2. For each s ∈ SD such that P_ts < 0: s must be (kept) turned on if and only if d ≤ floor(δ_t * x_ts).
# 3. For each c ∈ C_t : set q_tc(d) = y_tc * q_tc / δ_t .
# 4. For each s ∈ B such that u_max != 0: set h_ts(d) = h_ts / δ_t and g_ts(d) = g_ts / δ_t .

# Executes the RTCS for the robust model at time period t and time unit d.
# the possible actions that must be considered at each time unit d ∈ {1, . . . , δ t },
# of each time period t ∈ T , are:
# • turn on/off a production/consumption system s ∈ S D ,
# • define a quantity q c t (d) of energy to be bought/sold under an engaged contract c ∈ C t ,
# • define a quantity g s t (d) to be absorbed by storage system s ∈ B,
# • define a quantity h ts (d) to be refunded by storage system s ∈ B.
# Clearly, a subset of actions taken in a time unit must obey the same constraints considered in
# the last section for each time period.
function rtcs(t, d, period_size, sim_strategy, model_policy, instance_as_dict, P_hat_t, P_hat_td,
                        y, model_type, model_solution, previous_r_td, previous_drivable_charge,
                        previous_q_t, logger, verbose)
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance_as_dict["instance"])
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    T, S, C, D, ND, ST, NDU = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt, nbNDU)
    SetSigma = NDU
    num_contracts = instance_as_dict["num_contracts"]
    drivable = instance_as_dict["drivable"]
    n_drivable = instance_as_dict["n_drivable"]
    n_drivable_uncertain = instance_as_dict["n_drivable_uncertain"]
    storage = instance_as_dict["storage"]
    contract = instance_as_dict["contract"]
    period = instance_as_dict["period"]
    C_t = 1:num_contracts[t]
    contracts_in_period_t = contract[contract[!, :period] .== t, :]

    if model_type == "deterministic"
		q_t, x_t, e_t, r_t, g_t, h_t = calculate_variables_from_deterministic(instance_as_dict, model_solution, t)
    else  # Robust
		# Obtain q_t'(.), x_t'(.), e_t'(.), r_t'(.), g_t'(.), based on model solution
        q_t, x_t, e_t, r_t, g_t, h_t = calculate_variables_from_robust_ldr(instance_as_dict, model_solution, t, P_hat_t, period_size, y)
    end
    # The variables below will be calculated according to the chosen RTCS strategy
    q_td = zeros(Float64, num_contracts[t])  # amount consumed (< 0) / provided (> 0) by the client in contract c -> OK!
    x_td = zeros(Float64, nbD)
    e_td = ZERO  # extra amount of electricity requested by the client to the parter (out of any contract) ->
    g_td = zeros(Float64, nbSt)  # amount absorbed by system s (in ST) during the time period t -> OK!
    h_td = zeros(Float64, nbSt)  # amount refunded by system s (in ST) during time period t -> OK!
    # amount stored in system s (in ST) at the beginning of time period t -> VALUES COME FROM THE RCCP MODEL VARIABLE FROM PREVIOUS ITERATION
    r_td = deepcopy(previous_r_td)

    # print model-suggested variable values in trace log file
    print_model_values_for_variables(logger, t, d, num_contracts, period_size, y, q_t, nbD, D, drivable, x_t, nbSt, ST, storage, g_t, h_t)

    certain_ND_power = zeros(Float64, nbND)
    for s in ND
        certain_ND_power[s] = n_drivable[s,:pORc][t] / period_size
    end
    engaged_contracts = zeros(Int64, num_contracts[t])
    num_engaged_contracts = 0
    for c in C_t
        engaged_contracts[c] = y[t, c]
        num_engaged_contracts += y[t, c]
    end
    sum_D_prod = ZERO
    sum_D_cons = ZERO
    D_power = zeros(Float64, nbD)
    for s in D
        D_power[s] = drivable[s,:pORc][t] / period_size
        if drivable[s,:pORc][t] > EPS
            sum_D_prod += drivable[s,:pORc][t] / period_size
        elseif drivable[s,:pORc][t] < -EPS
            sum_D_cons += drivable[s,:pORc][t] / period_size
        end
    end
    # *** Déficit ou superavit: a cada microperíodo, podemos ter balanço positivo ou negativo.
    # No segundo caso, pode ser necessário, por ex, usar energia da bateria ou mesmo comprar
    # fora de contrato para fechar o período com o balanço positivo.
    # gap = sum (q_t'(.), x_t'(.), g_t'(.), e_t'(.), r_t'(.)) + Ŝ_ND + S_ND
    sum_NDU = (nbNDU > 0) ? sum(P_hat_td[s] for s in NDU) : ZERO  # Ŝ_ND : uncertain non-drivable power (comes from scenario data)
    sum_ND = (nbND > 0) ? sum(n_drivable[s,:pORc][t] / period_size for s in ND) : ZERO  # S_ND  FIXME Confirmar se eh pra dividir por delta

    if verbose println(logger, "    POWER PRODUCTION / CONSUMPTION for microperiod :") end
    if verbose println(logger, "        [Ŝ_ND] Uncertain non-drivable devices : sum = $(sum_NDU) : $(P_hat_td)") end
    if verbose println(logger, "        [S_ND] Certain non-drivable           : sum = $(sum_ND) : $(certain_ND_power)") end
    if verbose println(logger, "        [S_D] Consuming drivable              : sum = $(sum_D_cons) : $(D_power)") end
    if verbose println(logger, "        [S_D] Producing drivable              : sum = $(sum_D_prod) : $(D_power)") end
    if verbose println(logger, "        [q]    Contracts (sell / buy)         : $(q_td)") end
    if verbose println(logger, "        [y]    Engaged contracts              : $(engaged_contracts)") end
    if verbose println(logger, "        [y]    # of Engaged contracts         : $(num_engaged_contracts)") end
    gap = sum_NDU + sum_ND
    if verbose
        print(logger, "    Base gap (Ŝ_ND + S_ND) : ")
        if abs(gap) < EPS  println(logger, "BASE GAP ZERO") else  println(logger, "$(gap)") end
    end
    # Energy balance equation
    # OPTIONS FOR GAP CALCULATION
    # A) USE ALL MODEL SUGGESTED VARIABLE VALUES FOR q (contracts) and x (drivable) (FULLY TRUST THE MODEL)
    if model_policy == "full_model"
        if verbose println(logger, "    POLICY : Full-model policy") end
        # *** Use model values for contract usage
        q_td_model = calculate_model_value_for_q_td(t, num_contracts, period_size, y, q_t)
		sum_q = 0.0
        q_td = zeros(Float64, num_contracts[t])
        for c in C_t
            pi_minus_t = contracts_in_period_t[!, :min_period][c]
            pi_plus_t = contracts_in_period_t[!, :max_period][c]
            pi_minus_d = contracts_in_period_t[!, :min_delta][c]
            pi_plus_d = contracts_in_period_t[!, :max_delta][c]
            if y[t, c] == 1
                if pi_plus_t > 0
                    # buy contract  ==>> buy the minimum amount allowed in microperiod d or the minimum for period t / period_size
		    		min_buy = max(pi_minus_d, (pi_minus_t / period_size))
					min_buy = calculate_min_buy_given_contract_constraints(min_buy, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, (d == period_size), previous_q_t[c])
					# 1. The maximum energy to be bought (max_buy) is limited by contract upper limits (either period or delta plus)
					# (including energy already bought in previous deltas)
					max_buy = min(pi_plus_d, pi_plus_t - previous_q_t[c])
					# 2. The maxium energy to be bought cannot be greater than the maximum_per_period_t / period_size
					max_buy = min(max_buy, pi_plus_t / period_size)
					if q_td_model[c] < min_buy
						q_td[c] = min_buy
						if q_td_model[c] < min_buy - EPS  println(logger, "WARN: q_td_model[c] < min_buy : $(q_td_model[c]) < $(min_buy)")  end
					elseif q_td_model[c] > max_buy
						q_td[c] = max_buy
						if q_td_model[c] > max_buy + EPS  println(logger, "WARN: q_td_model[c] > max_buy : $(q_td_model[c]) > $(max_buy)")  end
					else  # q_td_model[c] is OK
						q_td[c] = q_td_model[c]
					end
					sum_q += q_td[c]
                    if verbose println(logger, "        Buying initial energy amount from contract (c = $(c)) : qtd = $(q_td[c])") end
                else    # sell contract ==>> sell the minimum amount allowed in microperiod d
					min_sell = -max(abs(pi_minus_d), abs((pi_minus_t / period_size)))  # TODO use the floor function to avoid propagating rounding errors
					min_sell = calculate_min_sell_given_contract_constraints(min_sell, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, (d == period_size), previous_q_t[c])
					max_sell = min(pi_plus_d, pi_plus_t - previous_q_t[c])
					max_sell = min(max_buy, pi_plus_t / period_size)
					if abs(q_td_model[c]) < abs(min_buy)
						q_td[c] = min_sell
						if abs(q_td_model[c]) < abs(min_buy) - EPS  println(logger, "WARN: q_td_model[c] < min_sell : $(q_td_model[c]) < $(min_sell)")  end
					elseif abs(q_td_model[c]) > abs(max_buy)
						q_td[c] = max_sell
						if abs(q_td_model[c]) > abs(max_buy) + EPS  println(logger, "WARN: q_td_model[c] > max_sell : $(q_td_model[c]) > $(max_sell)")  end
					else  # q_td_model[c] is OK
						q_td[c] = q_td_model[c]
					end
					sum_q += q_td[c]
                    if verbose println(logger, "        Selling initial energy amount to contract (c = $(c)) : qtd = $(q_td[c])") end
                end
            end
        end
        gap += sum_q
        # *** Use model values for drivable (variables x)
        x_td = calculate_model_value_for_x_td(drivable, nbD, D, period_size, t, d, x_t)
		if nbD > 0
        	gap += sum((drivable[s,:pORc][t] / period_size) * x_td[s] for s in D)
		end
        if verbose println(logger, "    Gap after full model usage (+ C_t + drivable_D) : $(gap)") end
        # BATTERY LEVELS WILL BE CALCULATED ON THE NORMAL RTCS ALGORITHM (NEXT SECTION)
    else
        # Include the minimum quota for each engaged contract
        sum_q = 0.0
        q_td = zeros(Float64, num_contracts[t])
        for c in C_t  # min_use_per_delta
            pi_minus_t = contracts_in_period_t[!, :min_period][c]
            pi_plus_t = contracts_in_period_t[!, :max_period][c]
            pi_minus_d = contracts_in_period_t[!, :min_delta][c]
            pi_plus_d = contracts_in_period_t[!, :max_delta][c]
            if y[t, c] == 1
                if pi_plus_t > 0
                    # buy contract  ==>> buy the minimum amount allowed in microperiod d or the minimum for period t / period_size
		    		# TODO Respect min and max per period and per delta! e.g. Contract	0	0.033	0.009	5000	7000	300	1200
                    # For each delta, must buy at least pi_minus_d, but in the last delta for the period, should have bought at least pi_minus_t !
                    min_buy = max(pi_minus_d, (pi_minus_t / period_size))
					min_buy = calculate_min_buy_given_contract_constraints(min_buy, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, (d == period_size), previous_q_t[c])
					sum_q += min_buy
                    q_td[c] = min_buy
                    if verbose println(logger, "        Buying minimum energy amount from contract (c = $(c)) : qtd = $(min_buy)") end
                else    # sell contract ==>> sell the minimum amount allowed in microperiod d
					min_sell = -max(abs(pi_minus_d), abs((pi_minus_t / period_size)))  # TODO use the floor function to avoid propagating rounding errors
					min_sell = calculate_min_sell_given_contract_constraints(min_sell, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, (d == period_size), previous_q_t[c])
                    sum_q += min_sell
                    q_td[c] = min_sell
                    if verbose println(logger, "        Selling minimum energy amount to contract (c = $(c)) : qtd = $(min_sell)") end
                end
            end
        end
        gap += sum_q
        if verbose println(logger, "        New gap after minimum contract buy    : $(gap)") end
        # Power required drivable systems, due to their minimum power requirement
        gap, x_td = power_required_drivable_consuming_system(drivable, D, nbT, t, gap, x_td, period_size,
                        previous_drivable_charge, logger, verbose)
        if verbose println(logger, "        New gap after minimum drivable supply : $(gap)") end
    end
    # Must end the period with balance >= 0
    if abs(gap) < EPS
        gap = ZERO
    end
    if verbose println(logger, "    Initial gap (+ C_t + drivable_D) = $(gap)") end

    if gap > EPS  # energy left over
        if verbose println(logger, "    STATUS : Energy left over (GAP > 0)") end
        # Simulation strategies:
        # - cheapest : always execute the cheapest operation first (using either batteries, drivable or contracts);
        # - conservative : use batteries first (to store positive gap), then use sell contracts (i.e. try to preserve batteries to worst-case periods);
        # - audacious : use sell contracts first, then store positive gap in batteries.
        if sim_strategy == "cheapest"
            gap, q_td, x_td, r_td, g_td, h_td = send_positive_gap_to_cheapest_operations(contracts_in_period_t, period_size,
                    C_t, y, q_td, previous_q_t, t, d, drivable, D, nbT, x_td, previous_drivable_charge, storage, ST, r_td, g_td, h_td, gap,
                    logger, verbose)
        else
            if sim_strategy == "conservative"
                # try to store the energy gap in the batteries
                gap, r_td, g_td = store_energy_gap_in_batteries(storage, ST, gap, r_td, g_td, logger, verbose)
                # if there is still energy left over, sell via 'Sell' contracts
                gap, q_td = sell_energy_gap_to_engaged_sell_contracts(contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
            else  # "audacious"
                # first, sell via 'Sell' contracts
                gap, q_td = sell_energy_gap_to_engaged_sell_contracts(contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
                # try to store the remaining energy in the batteries
                gap, r_td, g_td = store_energy_gap_in_batteries(storage, ST, gap, r_td, g_td, logger, verbose)
            end
            # if there is still energy left over, try to feed even more the drivable devices
            gap, x_td = send_energy_gap_to_drivable_consuming_system(drivable, D, nbT, t, gap, x_td, period_size,
                                previous_drivable_charge, logger, verbose)
        end
        # if there is still energy left over, throw it away (i.e. leave gap > 0)
        if gap > ZERO
            if verbose println(logger, "    WARN : There is still energy left over, throw it away") end
        end
    elseif gap < -EPS  # Lacking energy
        if verbose println(logger, "    STATUS : Lacking energy (GAP < 0)") end
        # Simulation strategies:
        # - cheapest : always execute the cheapest operation first (using either batteries, drivable or contracts);
        # - audacious : use batteries first (to compensate for neg gap), then use contracts;
        # - conservative : use contracts first, then use batteries (i.e. try to preserve batteries to worst-case periods).
        if sim_strategy == "cheapest"
            gap, q_td, x_td, r_td, g_td, h_td = fill_negative_gap_with_cheapest_operations(contracts_in_period_t,
                    period_size, C_t, y, q_td, previous_q_t, t, d, drivable, D, nbT, x_td, previous_drivable_charge, storage, ST,
                    r_td, g_td, h_td, gap, logger, verbose)
        else
            if sim_strategy == "audacious"  # "BCD", "BDC", "CBD", "CDB", "DBC", "DCB"
                # get extra energy from batteries, if they exist
                gap, r_td, h_td = retrieve_energy_gap_from_batteries(storage, ST, gap, r_td, h_td, logger, verbose)
                # if still needs energy, try to buy from any engaged contract
                gap, q_td_dummy = retrieve_energy_gap_from_engaged_buy_contracts(contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
            else  # "conservative"
                # buy all energy needed from engaged contracts and then out of contract and preserve batteries
                gap, q_td_dummy = retrieve_energy_gap_from_engaged_buy_contracts(contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
                # If it still needs energy, get from batteries, if they exist
                gap, r_td, h_td = retrieve_energy_gap_from_batteries(storage, ST, gap, r_td, h_td, logger, verbose)
            end
			# Try to power on any drivable devices that produce energy
            gap, x_td = turn_on_drivable_producing_devices(drivable, D, t, gap, x_td, period_size,
                                previous_drivable_charge, logger, verbose)
        end
        # if still needs energy, buy the rest out of contract ('e' variable)
        if gap < ZERO
            e_td = abs(gap)
            gap = ZERO
            if verbose println(logger, "    (e) WARN : Buying energy out of contract : qtd = $(e_td)") end
        end
    else
        gap = ZERO
        if verbose println(logger, "    STATUS : GAP == 0  NOTHING TO DO") end
    end

    # Ensure the quantity bought from contract c is under the maximum allowed (pi_plus)
	for c in C_t
		if y[t, c] == 1
			if d == period_size
				pi_plus_t = contracts_in_period_t[!, :max_period][c]
				pi_plus_d = contracts_in_period_t[!, :max_delta][c]
				sum_q_t_c = previous_q_t[c] + q_td[c]
				if pi_plus_t > 0  # buy contract
					if !(sum_q_t_c - pi_plus_t <= EPS)
						println("[RTCS] c = $(c), t = $(t), d = $(d) : sum_q_t[c] > pi_plus_t : $(sum_q_t_c) > $(pi_plus_t) ")
					end
					@assert sum_q_t_c - pi_plus_t <= EPS "[RTCS] sum_q_t[c] <= pi_plus_t"
				else  # sell contract
					if !(abs(sum_q_t_c) - abs(pi_plus_t) <= EPS)
						println("[RTCS] c = $(c), t = $(t), d = $(d) : sum_q_t[c] < pi_plus_t : $(sum_q_t_c) > $(pi_plus_t) ")
					end
					@assert abs(sum_q_t_c) - abs(pi_plus_t) <= EPS "[RTCS] sum_q_t[c] <= pi_plus_t"
				end
			end
		else
			@assert q_td[c] == 0
			@assert previous_q_t[c] == 0
		end
	end

    if verbose println(logger, "    Final gap = $(gap)") end
    if verbose println(logger, "    (r) Initial battery levels          : $(previous_r_td)") end
    if verbose println(logger, "    (r) Final battery levels            : $(r_td)") end
    if verbose println(logger, "    (x) Final usage of drivable devices : $(x_td)") end
    if verbose println(logger, "    (q) Final usage of contracts        : $(q_td)") end
    if verbose println(logger, "    (g) Final usage of batteries        : $(g_td)") end
    if verbose println(logger, "    (h) Final refund of batteries       : $(h_td)") end
    if verbose println(logger, "    (e) Final out of contract usage     : $(e_td)") end
    if verbose && e_td > ZERO println(logger, "    (e) e > 0") end
    # only sum contract fixed cost if d == 1 (first micro period)
    # FIXME Devo somar o custo dos dispositivos NDU no calculo do custo ?
    cost = d == 1 ? sum( contract[contract[!, :period] .== t, :cost_fix][c] * y[t,c] for c in 1:num_contracts[t]) : ZERO
    cost += sum(contract[contract[!, :period] .== t, :cost_var][c] * q_td[c] for c in 1:num_contracts[t])
	if size(drivable, 1) > 0
    	cost += sum(drivable[s,:cost] * (drivable[s,:pORc][t] / period_size) * x_td[s] for s in D)
	end
	if size(storage, 1) > 0
    	cost += sum(storage[s,:cost] * (g_td[s] + h_td[s]) for s in ST)
	end
    cost += (period[t,:cost_out] * e_td)
	if size(n_drivable_uncertain, 1) > 0
    	cost += sum(n_drivable_uncertain[sigma, :cost] * P_hat_td[sigma] for sigma in SetSigma)
	end
	if size(n_drivable, 1) > 0
    	cost += sum(n_drivable[s,:cost] * n_drivable[s,:pORc][t] / period_size for s in ND)
	end
    if verbose println(logger, "    PERIOD COST                         : $(cost)") end
    flush(logger)
    return q_td, x_td, e_td, r_td, g_td, h_td, gap, cost  # values calculated by RTCS heuristic for period (t, d)
end

function calculate_model_value_for_x_td(drivable, nbD, D, period_size, t, d, x_t)
    x_td = zeros(Float64, nbD)  # percentage of time period t in which drivable system s (in D) is turned on -> OK!
    # VARIABLE X : DEFAULT RULE FOR DRIVABLE SYSTEM ==>> AVERAGE VALUE FOR PERIOD T
    # Must use (d - 1) instead of d because the period number in Julia starts with 1, not 0 !
    for s in D  # x[t,s] : percentage of time period t in which drivable system s is turned on
        P_D_ts = drivable[s,:pORc][t]
        if x_t[s] > EPS  # Avoid rounding up period 1 (x_td[s] = 1.0) if x_t[s] == 0
            if P_D_ts > EPS  # Rule 1
                if (d - 1) <= ceil(period_size * x_t[s])
                    x_td[s] = 1.0
                end
            elseif P_D_ts < -EPS # Rule 2
                if (d - 1) <= floor(period_size * x_t[s])
                    x_td[s] = 1.0
                end
            end
        end
    end
    return x_td
end

function calculate_model_value_for_q_td(t, num_contracts, period_size, y, q_t)
    q_td = zeros(Float64, num_contracts[t])  # amount consumed (< 0) / provided (> 0) by the client in contract c -> OK!
    for c in 1:num_contracts[t]  # Rule 3 -> Average value
        q_td[c] = y[t, c] * (q_t[c] / period_size)
    end
    return q_td
end

function print_model_values_for_variables(logger, t, d, num_contracts, period_size, y, q_t, nbD, D, drivable, x_t, nbSt, ST, storage, g_t, h_t)
    println(logger, "    MODEL SUGGESTIONS  :")
    x_td = calculate_model_value_for_x_td(drivable, nbD, D, period_size, t, d, x_t)
    println(logger, "        (x_t)  Model-suggested drivable devices % usage  : $(x_t)")
    println(logger, "        (x_td) Model-suggested drivable devices % usage  : $(x_td)")
    q_td = calculate_model_value_for_q_td(t, num_contracts, period_size, y, q_t)
    q_t = [q_t[i] for i in 1:num_contracts[t]]
    sum_q = sum(q_td[c] for c in 1:num_contracts[t])
    println(logger, "        (q_t)  Model-suggested usage of contracts (avg)  : $(q_t)")
    println(logger, "        (q_td) Model-suggested usage of contracts (avg)  : $(q_td)")
    println(logger, "        [q]    Contracts (sell / buy)         : sum = $(sum_q) : $(q_td)")
    # VARIABLES G AND H : DEFAULT RULE FOR BATTERY SYSTEMS ==>> AVERAGE VALUE FOR PERIOD T
    g_td = zeros(Float64, nbSt)  # amount absorbed by system s (in ST) during the time period t -> OK!
    h_td = zeros(Float64, nbSt)  # amount refunded by system s (in ST) during time period t -> OK!
    for s in ST  # Rule 4 -> Average values
        u_max = storage[s,:uMax]
        if u_max != ZERO
            g_td[s] = g_t[s] / period_size
            if size(h_t, 1) > 0
                h_td[s] = h_t[s] / period_size
            end
        end
    end
    println(logger, "        (g_t)  Model-suggested usage of batteries        : g = $(g_t); h = $(h_t)")
    println(logger, "        (g_td) Model-suggested usage of batteries (avg)  : g = $(g_td); h = $(h_td)")
	flush(logger)
end
