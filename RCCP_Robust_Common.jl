# ==================================================================
# RCCP_Robust_Common.jl
# RCCP Robust MILP Model Helper Functions.
# ==================================================================

include("RCCP_FileReader.jl")

function print_detailed_solution()
    println("\namount consumed (< 0) / provided (> 0) by the client in contract i")
    println("q0[t, c] = ", value(q0))
    println("q_[T, c, SetTau, SetSigma] = ", value(q_))

    println("\namount stored in system s at the beginning of time period t ")
    println("r0[1:nbT+1, ST] = ", value(r0))
    println("r_[1:nbT+1, ST, SetTau, SetSigma] = ", value(r_))

    println("\namount absorbed by system s during the time period t and that will be already stored in system s at the end of this time period")
    println("g0[T, ST] = ", value(g0))
    println("g_[T, ST, SetTau, SetSigma] = ", value(g_))

    println("\nextra amount of electricity requested by the client to the parter (out of any engaged contract) in order to satisfy his needs at time period t")
    println("e0[T] = ", value(e0))
    println("e_[T, SetTau, SetSigma] = ", value(e_))

    println("\npercentage of time period p in which drivable system s is turned on")
    println("x0[T, D] = ", value(x0))
    println("x_[T, D, SetTau, SetSigma] = ", value(x_))
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

function capture_variable_solution_as_array_3d(var_value, t0, nbT, range3)
    var_sol = [Float64(0.0) for x=1:nbT, y=1:nbT, z=range3]
    for t in t0:nbT, y in t0:t, z=range3
        var_sol[t, y, z] = trunc_if_less_than_eps(value(var_value[t, y, z]))
    end
    return var_sol
end

function capture_variable_solution_as_array_4d(var_value, t0, nbT, range2, range4)
    var_sol = [Float64(0.0) for x=1:nbT, y=range2, z=1:nbT, k=range4]
    for t in t0:nbT, y in range2, z=t0:t, k=range4
        var_sol[t, y, z, k] = trunc_if_less_than_eps(value(var_value[t, y, z, k]))
    end
    return var_sol
end

function capture_variable_solution_as_array_4d_custom_q(var_value, t0, nbT, num_contracts, range4)
    max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
    var_sol = [Float64(0.0) for x=1:nbT, y=1:max_contracts_per_period, z=1:nbT, k=range4]
    for t in t0:nbT, c in 1:num_contracts[t], z=t0:t, k=range4
        var_sol[t, c, z, k] = trunc_if_less_than_eps(value(var_value[t, c, z, k]))
    end
    return var_sol
end

function capture_variable_solution_as_array_4d_custom_r(var_value, t0, nbT, range2, range4)
    var_sol = [Float64(0.0) for x=1:nbT+1, y=range2, z=1:nbT+1, k=range4]
    for t in t0:nbT+1, y in range2, z=t0:t, k=range4
        var_sol[t, y, z, k] = trunc_if_less_than_eps(value(var_value[t, y, z, k]))
    end
    return var_sol
end

# Calculate the independent term of each affine function (q(.), x(.) etc.)
function calculate_independent_term_rob_delta(q0, x0, e0, r0, g0, period_size, d, calc_method = "all_beginning")
    if calc_method == "average"  # Divide independent term (e0, q0, etc.) by the period size
        return q0 ./ period_size, x0 ./ period_size, e0 ./ period_size, r0 ./ period_size, g0 ./ period_size
    else  # "full value in the first microperiod d"
        if d == 1
            return q0, x0, e0, r0, g0
        else
            return q0 .* 0.0, x0 .* 0.0, e0 .* 0.0, r0 .* 0.0, g0 .* 0.0
        end
    end
end

# Calculate the robust model variable values (q, x, e, etc.) for a specified time period t
function calculate_variables_from_robust_ldr(instance_as_dict, rob_solution, t, P_hat_t, period_size, y)
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance_as_dict["instance"])
    nbNDU = size(instance_as_dict["n_drivable_uncertain"], 1)
    T, S, C, D, ND, ST, NDU = obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt, nbNDU)
    SetSigma = NDU
    num_contracts = instance_as_dict["num_contracts"]
    storage = instance_as_dict["storage"]
    C_t = 1:num_contracts[t]
    q0 = rob_solution["q0"]
    q_ = rob_solution["q_"]
    x0 = rob_solution["x0"]
    x_ = rob_solution["x_"]
    e0 = rob_solution["e0"]
    e_ = rob_solution["e_"]
    r0 = rob_solution["r0"]
    r_ = rob_solution["r_"]
    g0 = rob_solution["g0"]
    g_ = rob_solution["g_"]

    # FIXME Implement alternative method which uses the whole independent term (D_0) in period d == 1
    # For independent terms from each microperiod d (D_0), use either average values
    q0_d, x0_d, e0_d, r0_d, g0_d = calculate_independent_term_rob_delta(q0, x0, e0, r0, g0, period_size, 1, "all_beginning")

    e_t = e0_d[t] + sum(e_[t, t, sigma] * P_hat_t[sigma] for sigma in SetSigma)
    q_t = Float64[]
    for c in C_t
        q_tc = 0.0
        # if the contract is active
        @assert y[t,c] >= 0
        q_tc += y[t,c] * (q0_d[t, c] + sum(q_[t, c, t, sigma] * P_hat_t[sigma] for sigma in SetSigma))
        push!(q_t, q_tc)
    end
    x_t = Float64[]
    for s in D
        temp = x0_d[t, s] + sum(x_[t, s, t, sigma] * P_hat_t[sigma] for sigma in SetSigma)
        push!(x_t, temp)
    end
    r_t = Float64[]
    g_t = Float64[]
    h_t = Float64[]
    for s in ST
        temp_r = r0_d[t, s] + sum(r_[t, s, t, sigma] * P_hat_t[sigma] for sigma in SetSigma)
        temp_g = g0_d[t, s] + sum(g_[t, s, t, sigma] * P_hat_t[sigma] for sigma in SetSigma)
        push!(r_t, temp_r)
        push!(g_t, temp_g)
    end
    return q_t, x_t, e_t, r_t, g_t, h_t
end
