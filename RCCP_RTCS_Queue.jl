# ============================================================================
# RCCP_RTCS_Queue.jl
# RTCS priority queue and utility functions.
# =============================================================================

using DataStructures

# Include RTCS util and operation functions
include("RCCP_RTCS_Ops.jl")

# In the periority queue key, the first letter represents the operation type ('b' for batteries,
# 'd' for drivable, 'c' for contracts) and the rest of the string contains the device or contract id.
function extract_id_from_key(key)
    return parse(Int64, key[2:length(key)])
end

function store_energy_gap_in_batteries(storage, ST, gap, r_td, g_td, logger, verbose)
    pq = PriorityQueue()
    store_energy_gap_in_batteries_enqueue(pq, storage, ST, gap, r_td, g_td, logger, verbose)
    while !isempty(pq) && (gap > EPS)
        key = dequeue!(pq)  # Obtain the minimum cost storage device
        s = extract_id_from_key(key)
        gap, r_td, g_td = store_energy_gap_in_battery(s, storage, ST, gap, r_td, g_td, logger, verbose)
    end
    return gap, r_td, g_td
end

function store_energy_gap_in_batteries_enqueue(pq, storage, ST, gap, r_td, g_td, logger, verbose)
    for s in ST
        if r_td[s] < storage[s,:uMax] - EPS
            key = 'b' * string(s)
            pq[key] = storage[s,:cost]  # enqueue sorage elements by their cost
        else
            if verbose println(logger, "        CANNOT store extra energy in battery (s = $(s)). Exceeds storage limit of $(storage[s,:uMax])") end
        end
    end
    return pq
end

function retrieve_energy_gap_from_batteries(storage, ST, gap, r_td, h_td, logger, verbose)
    pq = PriorityQueue()
    retrieve_energy_gap_from_batteries_enqueue(pq, storage, ST, gap, r_td, h_td, logger, verbose)
    while !isempty(pq) && (gap < -EPS)
        key = dequeue!(pq)  # Obtain the minimum cost storage device
        s = extract_id_from_key(key)
        gap, r_td, h_td = retrieve_energy_gap_from_battery(s, storage, ST, gap, r_td, h_td, logger, verbose)
    end
    return gap, r_td, h_td
end

function retrieve_energy_gap_from_batteries_enqueue(pq, storage, ST, gap, r_td, h_td, logger, verbose)
    # TODO change the order of the batteries used, according to lostCoef
    for s in ST
        if r_td[s] > storage[s,:uMin] + EPS
            key = 'b' * string(s)
            pq[key] = storage[s,:cost]  # enqueue sorage elements by their cost
        else
            if verbose println(logger, "        CANNOT retrieve energy from battery (s = $(s)). Reached minimum storage limit of $(storage[s,:uMin])") end
        end
    end
end

function retrieve_energy_gap_from_engaged_buy_contracts(contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
    pq = PriorityQueue()
    retrieve_energy_gap_from_engaged_buy_contracts_enqueue(pq, contracts_in_period_t, period_size, C_t, y, t, gap, q_td, previous_q_t, logger, verbose)
    while !isempty(pq) && (gap < -EPS)
        key = dequeue!(pq)  # Obtain the minimum cost contract
        c = extract_id_from_key(key)
        gap, q_td = retrieve_energy_gap_from_engaged_buy_contract(c, contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
    end
    return gap, q_td
end

function retrieve_energy_gap_from_engaged_buy_contracts_enqueue(pq, contracts_in_period_t, period_size, C_t, y, t, gap, q_td, previous_q_t, logger, verbose)
    for c in C_t
        pi_plus_t = contracts_in_period_t[:max_period][c]
        pi_plus_d = contracts_in_period_t[:max_delta][c]
        if y[t,c] == 1 && pi_plus_t > 0 # if the contract is active and is a 'buy' contract
            max_buy = min(pi_plus_d, pi_plus_t - previous_q_t[c])
            if abs(q_td[c]) < abs(max_buy) - EPS
                # enqueue sorage elements by their cost
                key = 'c' * string(c)
                pq[key] = contracts_in_period_t[:cost_var][c]
            else
                if verbose println(logger, "        CANNOT retrieve energy from contract (c = $(c)). Exceeds limit of $(max_buy)") end
            end
        end
    end
end

function send_energy_gap_to_drivable_consuming_system(drivable, D, nbT, t, gap, x_td, period_size,
        previous_drivable_charge, logger, verbose)
    pq = PriorityQueue()
    send_energy_gap_to_drivable_consuming_system_enqueue(pq, drivable, D, nbT, t, gap, x_td, period_size,
            previous_drivable_charge, logger, verbose)
    while !isempty(pq) && (gap > EPS)
        key = dequeue!(pq)  # Obtain the minimum cost drivable device
        s = extract_id_from_key(key)
        gap, x_td = send_energy_gap_to_drivable_consuming(s, drivable, D, nbT, t, gap, x_td, period_size,
            previous_drivable_charge, logger, verbose)
    end
    return gap, x_td
end

function send_energy_gap_to_drivable_consuming_system_enqueue(pq, drivable, D, nbT, t, gap, x_td, period_size,
        previous_drivable_charge, logger, verbose)
    for s in D  # x[t,s] : percentage of time period t in which drivable system s is turned on
        P_D_ts = drivable[s,:pORc][t] / period_size
        if P_D_ts < ZERO && x_td[s] + EPS < 1.0  # If device consumes energy and can still receive more energy
            # ATTENTION ! power values are negative, since devices are consumer devices
            max_capacity = findmin(drivable[s,:pORc_max])[1]
            capacity_left = previous_drivable_charge[s] - max_capacity
            if capacity_left > ZERO
                key = 'd' * string(s)
                pq[key] = drivable[s,:cost]  # enqueue drivable elements by their cost
            else
                if verbose println(logger, "        CANNOT send more energy to drivable consumer (s = $(s)). Exceeds the max_capacity = $(max_capacity).") end
            end
        end
    end
end

# OBS : No need to respect [Pmin, Pmax] of producer drivable device!
function turn_on_drivable_producing_devices(drivable, D, t, gap, x_td, period_size,
        previous_drivable_charge, logger, verbose)
    pq = PriorityQueue()
    turn_on_drivable_producing_devices_enqueue(pq, drivable, D, t, gap, x_td, period_size,
            previous_drivable_charge, logger, verbose)
    while !isempty(pq) && (gap < -EPS)
        key = dequeue!(pq)  # Obtain the minimum cost drivable device
        s = extract_id_from_key(key)
        gap, x_td = turn_on_drivable_producing_device(s, drivable, D, t, gap, x_td, period_size,
                previous_drivable_charge, logger, verbose)
    end
    return gap, x_td
end

function turn_on_drivable_producing_devices_enqueue(pq, drivable, D, t, gap, x_td, period_size,
        previous_drivable_charge, logger, verbose)
    for s in D
        if drivable[s,:pORc][t] > ZERO  # Drivable produces energy (P > 0)
            # Take into account what has already been sent to the drivable device (required min power)
            P_D_ts = ((1.0 - x_td[s]) * drivable[s,:pORc][t]) / period_size
            if P_D_ts > EPS
                key = 'd' * string(s)
                pq[key] = drivable[s,:cost]  # enqueue drivable elements by their cost
            else
                if verbose println(logger, "        CANNOT get more energy from drivable producer (s = $(s)). Exceeds 100% usage.") end
            end
        end
    end
end

function sell_energy_gap_to_engaged_sell_contracts(contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
    pq = PriorityQueue()
    sell_energy_gap_to_engaged_sell_contracts_enqueue(pq, contracts_in_period_t, period_size, C_t, y, t, gap, q_td, previous_q_t, logger, verbose)
    while !isempty(pq) && (gap > EPS)
        key = dequeue!(pq)  # Obtain the minimum cost contract
        c = extract_id_from_key(key)
        gap, q_td = sell_energy_gap_to_engaged_sell_contract(c, contracts_in_period_t, period_size, C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
    end
    return gap, q_td
end

function sell_energy_gap_to_engaged_sell_contracts_enqueue(pq, contracts_in_period_t, period_size, C_t, y, t, gap, q_td, previous_q_t, logger, verbose)
    for c in C_t
        pi_plus_t = contracts_in_period_t[:max_period][c]
        pi_minus_t = contracts_in_period_t[:min_period][c]
        pi_minus_d = contracts_in_period_t[:min_delta][c]
        pi_plus_d = contracts_in_period_t[:max_delta][c]
        if y[t,c] == 1 && pi_plus_t < 0 # if the contract is active and is a 'sell' contract
            max_sell = -min(abs(pi_plus_d), abs(pi_plus_t) - abs(previous_q_t[c]))
            if abs(q_td[c]) < abs(max_sell) - EPS
                # enqueue sorage elements by their cost
                key = 'c' * string(c)
                pq[key] = contracts_in_period_t[:cost_var][c]
            else
                if verbose println(logger, "        CANNOT sell energy to contract (c = $(c)). Exceeds limit of $(max_sell).") end
            end
        end
    end
end

function send_positive_gap_to_cheapest_operations(contracts_in_period_t, period_size, C_t, y, q_td, previous_q_t,
        t, d, drivable, D, nbT, x_td, previous_drivable_charge, storage, ST, r_td, g_td, h_td, gap, logger, verbose)
    pq = PriorityQueue()
    # store_energy_gap_in_batteries()
    store_energy_gap_in_batteries_enqueue(pq, storage, ST, gap, r_td, g_td, logger, verbose)
    # sell_energy_gap_to_engaged_sell_contracts()
    sell_energy_gap_to_engaged_sell_contracts_enqueue(pq, contracts_in_period_t, period_size, C_t, y, t,
            gap, q_td, previous_q_t, logger, verbose)
    # send_energy_gap_to_drivable_consuming_system()
    send_energy_gap_to_drivable_consuming_system_enqueue(pq, drivable, D, nbT, t, gap, x_td, period_size,
            previous_drivable_charge, logger, verbose)
    while !isempty(pq) && (gap > EPS)
        key = dequeue!(pq)  # Obtain the minimum cost storage device
        s = extract_id_from_key(key)
        if key[1] == 'c'
            gap, q_td = sell_energy_gap_to_engaged_sell_contract(s, contracts_in_period_t, period_size,
                C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
        elseif key[1] == 'd'
            gap, x_td = send_energy_gap_to_drivable_consuming(s, drivable, D, nbT, t, gap, x_td, period_size,
                previous_drivable_charge, logger, verbose)
        elseif key[1] == 'b'
            gap, r_td, g_td = store_energy_gap_in_battery(s, storage, ST, gap, r_td, g_td, logger, verbose)
        end
    end
    return gap, q_td, x_td, r_td, g_td, h_td
end

function fill_negative_gap_with_cheapest_operations(contracts_in_period_t, period_size, C_t, y, q_td, previous_q_t,
        t, d, drivable, D, nbT, x_td, previous_drivable_charge, storage, ST, r_td, g_td, h_td, gap, logger, verbose)
    pq = PriorityQueue()
    # retrieve_energy_gap_from_batteries()
    retrieve_energy_gap_from_batteries_enqueue(pq, storage, ST, gap, r_td, h_td, logger, verbose)
    # retrieve_energy_gap_from_engaged_buy_contracts()
    retrieve_energy_gap_from_engaged_buy_contracts_enqueue(pq, contracts_in_period_t, period_size,
            C_t, y, t, gap, q_td, previous_q_t, logger, verbose)
    # turn_on_drivable_producing_devices()
    turn_on_drivable_producing_devices_enqueue(pq, drivable, D, t, gap, x_td, period_size,
            previous_drivable_charge, logger, verbose)
    while !isempty(pq) && (gap < -EPS)
        key = dequeue!(pq)  # Obtain the minimum cost operation
        s = extract_id_from_key(key)
        if key[1] == 'c'
            gap, q_td = retrieve_energy_gap_from_engaged_buy_contract(s, contracts_in_period_t, period_size,
                            C_t, y, t, d, gap, q_td, previous_q_t, logger, verbose)
        elseif key[1] == 'd'
            gap, x_td = turn_on_drivable_producing_device(s, drivable, D, t, gap, x_td, period_size,
                            previous_drivable_charge, logger, verbose)
        elseif key[1] == 'b'
            gap, r_td, h_td = retrieve_energy_gap_from_battery(s, storage, ST, gap, r_td, h_td, logger, verbose)
        end
    end
    return gap, q_td, x_td, r_td, g_td, h_td
end

# Power all consumer drivable systems which need an amount of energy (obligatory)
function power_required_drivable_consuming_system(drivable, D, nbT, t, gap, x_td, period_size,
        previous_drivable_charge, logger, verbose)
    for s in D  # x[t,s] : percentage of time period t in which drivable system s is turned on
        P_D_ts = drivable[s,:pORc][t] / period_size
        if P_D_ts < ZERO  # Drivable consumes energy
            # ATTENTION ! power values are negative, since devices are consumer devices
            # Lookahead for the next time periods starting with current period t
            capacity_needed = previous_drivable_charge[s] - sum(drivable[s,:pORc_min][t1] for t1 in t:nbT)
            if capacity_needed < ZERO  # minimum energy already achieved, given all upcoming time periods
                x_td[s] = ZERO
            else
                # turn on device
                energy = min(capacity_needed, abs(P_D_ts))
                gap -= energy
                x_td[s] = energy / abs(P_D_ts)
                if verbose println(logger, "        Sending required energy to drivable consumer (s = $(s)) :
                    qtd = $(energy) ; gap = $(gap) ; x[s] = $(x_td[s]) ; capacity_needed = $(capacity_needed)")
                end
            end
        end
    end
    return gap, x_td
end
