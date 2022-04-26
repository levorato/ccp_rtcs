# ============================================================================
# RCCP_RTCS_Ops.jl
# RTCS basic operation functions.
# =============================================================================

# Execute operation on battery device s
# maxAbsorption : maximum quantity of energy that can be absorbed by system s during ∆ time units (in each microperiod p)
function store_energy_gap_in_battery(s, storage, ST, gap, r_td, g_td, logger, verbose)
	maxAbsorption = storage[s,:maxAbsorption]
    energy = min(abs(gap), abs(storage[s,:uMax] - r_td[s]))
	energy = min(energy, maxAbsorption)
    energy_after_loss = energy * storage[s,:lostCoef]
    r_td[s] += energy * storage[s,:lostCoef]
    g_td[s] += energy  # amount absorbed by system s
    gap -= energy
    cost = storage[s,:cost]
    if verbose println(logger, "        Storing extra energy in battery (s = $(s), unit_cost = $(cost)) : qtd = $(energy), effective = $(energy_after_loss) ; gap = $(gap)") end
    return gap, r_td, g_td
end

# Execute operation on battery device s
# maxRefund : maximum quantity of energy that can be refunded by system s during ∆ time units (in each microperiod p)
function retrieve_energy_gap_from_battery(s, storage, ST, gap, r_td, h_td, logger, verbose)
	# try to get as much energy as possible from the batteries
	maxRefund = storage[s,:maxRefund]
    energy = min(abs(gap), abs(r_td[s] - storage[s,:uMin]))
	energy = min(energy, maxRefund)
    r_td[s] -= energy
    h_td[s] += energy * storage[s,:lostCoef]  # amount refunded by system s
    gap += energy * storage[s,:lostCoef]
    cost = storage[s,:cost]
    if verbose println(logger, "        Retrieving energy from battery (s = $(s), unit_cost = $(cost)) : qtd = $(energy) ; gap = $(gap)") end
    return gap, r_td, h_td
end

# Execute operation on contract c
function retrieve_energy_gap_from_engaged_buy_contract(c, contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
    pi_minus_t = contracts_in_period_t[!, :min_period][c]
	pi_plus_t = contracts_in_period_t[!, :max_period][c]
	pi_minus_d = contracts_in_period_t[!, :min_delta][c]
	pi_plus_d = contracts_in_period_t[!, :max_delta][c]
	# 1. The maximum energy to be bought (max_buy) is limited by contract upper limits (either period or delta plus)
	# (including energy already bought in previous deltas)
	bought_so_far = previous_q_t[c] + q_td[c]
    max_buy = min(pi_plus_d, pi_plus_t - previous_q_t[c])
	#println("max_buy = $(max_buy)")
	# 2. The maxium energy to be bought cannot be greater than the maximum_per_period_t / period_size
	max_buy = min(max_buy, pi_plus_t / period_size)
	#println("max_buy = $(max_buy)")
	# 3. Remove the amount of energy already bought for this contract in the current delta d (q_td[c])
	max_buy = max(0, max_buy - q_td[c])
	#println("max_buy = $(max_buy)")
	# 4. The gap (energy to be bought) is limited by max_buy
	#println("c = $(c), t = $(t), d = $(d), pi_plus_t = $(pi_plus_t), pi_plus_d = $(pi_plus_d)")
	#println("*** previous_q_t[c] = $(previous_q_t[c]), q_td[c] = $(q_td[c]), bought_so_far = $(bought_so_far)")
	#println("*** gap = $(gap), max_buy = $(max_buy)")
	must_buy = min(abs(gap), max_buy)
	# 5. The minium energy to be bought is at least the contract lower limits (either period or delta minus, the max of them)
	#    But the amount already bought in this delta (q_td[c]) must be subtracted !
	min_buy = max(0, max(pi_minus_d, pi_minus_t / period_size) - q_td[c])
	# 6. The minimum energy to be bought is the maximum between min_buy and the gap (to ensure contract lower limits)
	if min_buy > max_buy
		println(logger, "*** ERROR : min_buy > max_buy : $(min_buy) > $(max_buy)")
	end
	@assert min_buy <= max_buy
	buy = max(must_buy, min_buy)
	remaining = pi_plus_t - previous_q_t[c] - q_td[c]
	#println("*** min_buy = $(min_buy), q_td[c] = $(q_td[c]), remaining_to_buy = $(remaining), buy = $(buy)")

	energy = buy  # TODO use the floor function to avoid propagating rounding errors
	energy = calculate_min_buy_given_contract_constraints(energy, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, (d == period_size), previous_q_t[c] + abs(q_td[c]))

    q_td[c] += energy
	#println("=> Will buy $(energy) units from contract, q_td[c] = $(q_td[c])")
	sum_q_t_c = previous_q_t[c] + q_td[c]
	gap += energy
    cost_var = contracts_in_period_t[!, :cost_var][c]
    if verbose
		println(logger, "        Retrieving energy from engaged buy contract (c = $(c), cost_var = $(cost_var)) : qtd = $(energy) ; gap = $(gap)")
	end
    return gap, q_td
end

# Execute operation on contract c
# TODO rewrite rules for this operation !
function sell_energy_gap_to_engaged_sell_contract(c, contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
	pi_minus_t = contracts_in_period_t[!, :min_period][c]
	pi_plus_t = contracts_in_period_t[!, :max_period][c]
	pi_minus_d = contracts_in_period_t[!, :min_delta][c]
	pi_plus_d = contracts_in_period_t[!, :max_delta][c]
    max_sell = -min(abs(pi_plus_d), abs(pi_plus_t) - abs(previous_q_t[c]))
	# 2. The maxium energy to be sold cannot be greater than the maximum_per_period_t / period_size
	max_sell = -min(abs(max_sell), abs(pi_plus_t / period_size))
	# 3. Remove the amount of energy already sold for this contract in the current delta d (q_td[c])
	max_sell = -max(0, abs(max_sell) - abs(q_td[c]))

	must_sell = -min(abs(gap), abs(max_sell))
	# 5. The minium energy to be bought is at least the contract lower limits (either period or delta minus, the max of them)
	#    But the amount already bought in this delta (q_td[c]) must be deducted !
	min_sell = min(0, -max(abs(pi_minus_d), abs(pi_minus_t / period_size)) - abs(q_td[c]))
	# 6. The minimum energy to be bought is the maximum between min_buy and the gap (to ensure contract lower limits)
	if abs(min_sell) > abs(max_sell)
		println(logger, "*** ERROR : min_sell > max_sell : $(min_sell) > $(max_sell)")
	end
	@assert abs(min_sell) <= abs(max_sell)
	buy = min(must_buy, min_buy)

    energy = -min(abs(gap), abs(max_sell) - abs(q_td[c]))
	energy = calculate_min_sell_given_contract_constraints(energy, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, (d == period_size), previous_q_t[c])

    q_td[c] += energy   # In sell contracts, q[t,c] is negative by definition !
    gap += energy
    cost_var = contracts_in_period_t[!, :cost_var][c]
    if verbose println(logger, "        Selling extra energy to engaged sell contract (c = $(c), cost_var = $(cost_var)) : qtd = $(energy) ; gap = $(gap)") end
    return gap, q_td
end

# Execute the operation on drivable system s
function send_energy_gap_to_drivable_consuming(s, drivable, D, nbT, t, gap, x_td, period_size,
    previous_drivable_charge, logger, verbose)
    # ATTENTION ! power values are negative, since devices are consumer devices
    max_capacity = findmin(drivable[s,:pORc_max])[1]
    capacity_left = previous_drivable_charge[s] - max_capacity
    # Take into account what has already been sent to the drivable device (required min power)
    P_D_ts = ((1.0 - x_td[s]) * drivable[s,:pORc][t]) / period_size
    if P_D_ts < ZERO   # Drivable consumes energy
        # turn on device
        energy = min(gap, abs(P_D_ts))
        energy = min(energy, capacity_left)
        gap -= energy
        x_td[s] += energy / abs(P_D_ts)
        cost = drivable[s,:cost]
        if verbose println(logger, "        Sending extra energy to drivable consumer (s = $(s), cost = $(cost)) : qtd = $(energy) ; gap = $(gap) ; x[s] = $(x_td[s])") end
    end
    return gap, x_td
end

# Execute operation on drivable device s
function turn_on_drivable_producing_device(s, drivable, D, t, gap, x_td, period_size,
        previous_drivable_charge, logger, verbose)
    # x[t,s] : percentage of time period t in which drivable system s is turned on
    # Take into account what has already been sent to the drivable device (required min power)
    P_D_ts = ((1.0 - x_td[s]) * drivable[s,:pORc][t]) / period_size
    energy = min(abs(gap), P_D_ts)
    # turn on Drivable
    x_td[s] += energy / P_D_ts
    gap += energy
    cost = drivable[s,:cost]
    if verbose println(logger, "        Turning ON drivable producer (s = $(s), unit_cost = $(cost)) : qtd = $(energy) ; gap = $(gap) ; x_td[s] = $(x_td[s])") end
    return gap, x_td
end

function calculate_min_buy_given_contract_constraints(buy_qty, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, last_delta, previous_q_t)
    min_buy = buy_qty
    min_buy = min(min_buy, pi_plus_d)
    @assert (pi_plus_t - previous_q_t > EPS) "[RTCS] pi_plus_t[c] > q_t[c]"
    min_buy = min(min_buy, pi_plus_t - previous_q_t)
	sum_q_t_c = previous_q_t + min_buy
    if last_delta  # indicates if d == period_size
		# Ensure contract lower and upper limits for time period t are respected
		if sum_q_t_c - pi_minus_t < -EPS
			println("[RTCS] contract_buy : sum_q_t[c] < pi_minus_t : $(sum_q_t_c) < $(pi_minus_t) ")
		end
		@assert sum_q_t_c - pi_minus_t >= -EPS "[RTCS] sum_q_t_c >= pi_minus_t"
    end
	if sum_q_t_c - pi_plus_t > EPS
		println("[RTCS] contract_buy : sum_q_t[c] > pi_plus_t : $(sum_q_t_c) > $(pi_plus_t) ")
	end
	#@assert sum_q_t_c - pi_plus_t <= EPS "[RTCS] contract_buy : sum_q_t[c] <= pi_plus_t"
    return min_buy
end

function calculate_min_sell_given_contract_constraints(sell_qty, pi_minus_d, pi_minus_t, pi_plus_d, pi_plus_t, last_delta, previous_q_t)
    min_sell = sell_qty
    min_sell = -min(abs(min_sell), abs(pi_plus_d))
    @assert (abs(pi_plus_t) - abs(previous_q_t) > EPS) "[RTCS] pi_plus_t[c] > q_t[c]"
    min_sell = -min(abs(min_sell), abs(pi_plus_t) - abs(previous_q_t))
	sum_q_t_c = previous_q_t + min_sell
    if last_delta  # indicates if d == period_size
		# Ensure contract lower and upper limits for time period t are respected
		if abs(sum_q_t_c) - abs(pi_minus_t) < EPS
			println("[RTCS] contract_sell : sum_q_t[c] < pi_minus_t : $(sum_q_t_c) < $(pi_minus_t) ")
		end
        @assert abs(sum_q_t_c) - abs(pi_minus_t) >= EPS "[RTCS] sum_q_t_c >= pi_minus_t"
    end
	if abs(sum_q_t_c) - abs(pi_plus_t) > EPS
		println("[RTCS] contract_sell : sum_q_t[c] > pi_plus_t : $(sum_q_t_c) > $(pi_plus_t) ")
	end
	#@assert abs(sum_q_t_c) - abs(pi_plus_t) <= EPS "[RTCS] contract_sell : sum_q_t[c] <= pi_plus_t"
    return min_sell
end

function find_contract_buy_above_limit(c, pi_plus_t, pi_plus_d, sum_q_t_c, gap)
	if sum_q_t_c - pi_plus_t > EPS
		diff = sum_q_t_c - pi_plus_t
		@assert diff < pi_plus_d
		return diff
		#println("[RTCS] t = $(t), d = $(d) : sum_q_t[c] > pi_plus_t : $(sum_q_t_c) > $(pi_plus_t) ")
	end
	return 0
end

function find_contract_sell_above_limit(c, pi_plus_t, pi_plus_d, sum_q_t_c, gap)
	if abs(sum_q_t_c) - abs(pi_plus_t) > EPS
		diff = abs(sum_q_t_c) - abs(pi_plus_t)
		@assert diff < abs(pi_plus_d)
		return diff
		#println("[RTCS] t = $(t), d = $(d) : sum_q_t[c] > pi_plus_t : $(sum_q_t_c) > $(pi_plus_t) ")
	end
	return 0
end
