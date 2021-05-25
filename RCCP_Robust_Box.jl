# ==================================================================
# RCCP_Robust_Box.jl
# RCCP Robust MILP Model Based on the Box Uncertainty Set
# ==================================================================
using JuMP
using CSV
using DataFrames
import MathProgBase

include("RCCP_Robust_Common.jl")

# Solve the Robust RCCP MILP Model starting from period t0
# Fix the contracts y according to fixed_contracts parameter, if parameter is given
function solve_robust_model_box(instance_as_dict::Dict, general_logger, infeasible_logger, solver_params, t0 = 1,
                                fixed_contracts = Int64[],
                                initial_battery = Float64[], previous_drivable_charge = Float64[])
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

    println("\n===========  R O B U S T    B O X    M O D E L  ( t0 = $(t0) ) ===========\n")
    mrob = create_jump_model(solver_params)
    #mrob = Model(solver=GurobiSolver(TimeLimit=1800))
    # Instantiating the model indices
    println("Indices: $(instance)")
    nbT, nbS, nbC, nbD, nbND, nbSt = obtain_instance_parameters(instance)
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

    v_bar = zeros(Float64, nbTau, nbNDU)
    v_hat = zeros(Float64, nbTau, nbNDU)
    for sigma in SetSigma
        for tau in SetTau
            p_min = n_drivable_uncertain[sigma, :Pmin][tau]
            pdt_min = n_drivable_uncertain[sigma, :Pdt_min][tau]
            p_max = n_drivable_uncertain[sigma, :Pmax][tau]
            pdt_max = n_drivable_uncertain[sigma, :Pdt_max][tau]
            # Model sanity Check
            #if p_min < 0 || p_max < 0
            #    @assert p_min >= p_max "[uncertain ND] pmax < 0 => pmin >= pmax"
            #    @assert pdt_min >= pdt_max "[uncertain ND] pmax < 0 => pdt_min >= pdt_max"
            #elseif p_min > 0 || p_max > 0
            #    @assert p_min <= p_max "[uncertain ND] pmax > 0 => pmin <= pmax"
            #    @assert pdt_min <= pdt_max "[uncertain ND] pmax > 0 => pdt_min <= pdt_max"
            #end
            v_hat[tau,sigma] = Float64(abs(p_max - p_min) / 2.0)   # peak value interval length
            v_bar[tau,sigma] = Float64((p_max + p_min) / 2.0)   # nominal / mean value
            @assert v_hat[tau,sigma] >= 0
            if p_max > 0 || p_min > 0
                @assert abs(p_min + v_hat[tau,sigma] - v_bar[tau,sigma]) <= EPS
            elseif p_min < 0 || p_max < 0
                #println("sigma=$(sigma), tau=$(tau) : v_bar=$(v_bar[tau,sigma]), v_hat=$(v_hat[tau,sigma]), p_max - v_hat=$(p_max - v_hat[tau,sigma])")
                @assert abs((p_max - v_hat[tau,sigma]) - v_bar[tau,sigma]) <= EPS
            end
        end
    end

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
    alpha = zeros(Float64, nbT, max_contracts_per_period)
    pi_minus = zeros(Float64, nbT, max_contracts_per_period)
    pi_plus = zeros(Float64, nbT, max_contracts_per_period)
    v_C = zeros(Float64, nbT, max_contracts_per_period)
    cost_fix = zeros(Float64, nbT, max_contracts_per_period)
    cost_var = zeros(Float64, nbT, max_contracts_per_period)
    for t in T
        contracts_in_period_t = contract[contract[!, :period] .== t, :]
        for c in 1:num_contracts[t]
            alpha[t, c] = contract[contract[!, :period] .== t, :cost_var][c]
            pi_minus[t, c] = contracts_in_period_t[!, :min_period][c]
            pi_plus[t, c] = contracts_in_period_t[!, :max_period][c]
            v_C[t, c] = contracts_in_period_t[!, :cost_fix][c]
        end
        # validation code
        c_t = 1
        for c in C
            if contract[c,:period] == t
                cost_fix[t, c_t] = contract[c,:cost_fix]
                cost_var[t, c_t] = contract[c,:cost_var]
                @assert cost_fix[t, c_t] == v_C[t, c_t]
                @assert cost_var[t, c_t] == alpha[t, c_t]
                c_t += 1
            end
        end
    end
    ## sum(q[c] for c in C if contract[c,:period] == t)

    P_ND = zeros(Float64, nbT, nbND)
    for s in ND, t in T
        P_ND[t,s] = n_drivable[s,:pORc][t]
    end

    v_D = zeros(Float64, nbD)
    P_D = zeros(Float64, nbT, nbD)
    Pmin = zeros(Float64, nbT, nbD)
    Pmax = zeros(Float64, nbT, nbD)
    for s in D
        v_D[s] = drivable[s,:cost]
        for t in T
            P_D[t,s] = drivable[s,:pORc][t]
            Pmin[t,s] = drivable[s,:pORc_min][t]
            Pmax[t,s] = drivable[s,:pORc_max][t]
        end
    end

    v_ST = zeros(Float64, nbSt)
    theta_ref = zeros(Float64, nbSt)
    theta_abs = zeros(Float64, nbSt)
    lambda = zeros(Float64, nbSt)
    u_min = zeros(Float64, nbSt)
    u_max = zeros(Float64, nbSt)
    u0 = zeros(Float64, nbSt)
    for s in ST
        v_ST[s] = storage[s,:cost]
        theta_ref[s] = storage[s,:maxRefund]
        theta_abs[s] = storage[s,:maxAbsorption]
        lambda[s] = storage[s,:lostCoef]
        u_min[s] = storage[s,:uMin]
        u_max[s] = storage[s,:uMax]
        u0[s] = storage[s,:uInit]
    end

    beta = zeros(Float64, nbT)
    delta = zeros(Float64, nbT)
    for t in T
        beta[t] = period[t,:cost_out]
        delta[t] = period[t,:size]
    end
    v_NDU = zeros(Float64, nbNDU)
    for sigma in SetSigma
       v_NDU[sigma] = n_drivable_uncertain[sigma, :cost]
    end

    # Instantiating the decision variables
    # \varsigma in S_ND_uncertain, \tau in

    # ==== Here and now variables ====
    # (12) 1 if engages contract, 0 otherwise
    # Ineq. 44 (?)
    @variable(mrob, y[t=T, c=1:num_contracts[t]], Bin)
    # Contract fixation : if parameter fixed_contracts is given
    if size(fixed_contracts, 1) > 0
        engaged_v = Int64[]
        for t in 1:t0-1, c in 1:num_contracts[t]
            push!(engaged_v, fixed_contracts[t, c])
        end
        for t in t0:nbT
            for c in 1:num_contracts[t]
                @constraint(mrob, y[t, c] == fixed_contracts[t, c])
                push!(engaged_v, fixed_contracts[t, c])
            end
        end
        println("Fixing contract engagement variables to $(engaged_v)")
    end

    # ==== Wait and see variables ====
    # amount consumed (< 0) / provided (> 0) by the client in contract i
    ### @variable(mrob, q[C, T])
    @variable(mrob, q0[t=T, c=1:num_contracts[t]])
    @variable(mrob, q_[t=T, c=1:num_contracts[t], t0:t, SetSigma])
    # (13) amount stored in system s at the beginning of time period t
    ### @variable(mrob, r[1:nbT+1, ST])
    @variable(mrob, r0[t0:nbT+1, ST])
    @variable(mrob, r_[t=t0:nbT+1, ST, tau=t0:nbT+1, SetSigma])  # [t=t0:nbT+1, ST, tau=t0:t, SetSigma]  # x[i=1:2, j=i:2]
    # Initial battery storage fixation : if parameter initial_battery is given
    if size(initial_battery, 1) > 0
        for s in ST
            # treat rounding errors of battery levels
            if (initial_battery[s] - storage[s,:uMin] < -EPS)
                println("WARN: assert initial_battery[s] >= storage[s,:uMin] : $(initial_battery[s]) < $(storage[s,:uMin])")
                println(general_logger, "[Robust Model, t0 = $(t0)] WARN: assert initial_battery[s] >= storage[s,:uMin] : $(initial_battery[s]) < $(storage[s,:uMin])")
                if abs(initial_battery[s] - storage[s,:uMin]) / storage[s,:uMin] <= RTCS_ROUNDING_ERROR
                    initial_battery[s] = storage[s,:uMin]
                else
                    println(general_logger, "[Robust Model, t0 = $(t0)] initial_battery[s] must be >= storage[s,:uMin]. Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                    flush(general_logger)
                    error("initial_battery[s] must be >= storage[s,:uMin]. Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                end
            end
            if (initial_battery[s] - storage[s,:uMax] > EPS)
                println("WARN: assert initial_battery[s] <= storage[s,:uMax] : $(initial_battery[s]) > $(storage[s,:uMax])")
                println(general_logger, "[Robust Model, t0 = $(t0)] WARN: assert initial_battery[s] <= storage[s,:uMax] : $(initial_battery[s]) > $(storage[s,:uMax])")
                if abs(initial_battery[s] - storage[s,:uMax]) / storage[s,:uMax] <= RTCS_ROUNDING_ERROR
                    initial_battery[s] = storage[s,:uMax]
                else
                    println(general_logger, "[Robust Model, t0 = $(t0)] initial_battery[s] must be <= storage[s,:uMax] Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                    flush(general_logger)
                    error("initial_battery[s] must be <= storage[s,:uMax] Difference exceeds rounding error of $(RTCS_ROUNDING_ERROR).")
                end
            end
            u0[s] = initial_battery[s]
        end
        println("Fixing battery level variables to $(initial_battery)")
    end
    # Initial drivable charge fixation
    if size(previous_drivable_charge, 1) > 0
        println("Fixing drivable initial charge to $(previous_drivable_charge)")
    end
    # (13) amount absorbed by system s during the time period t and that will be already stored in system s
    # at the end of this time period
    ### @variable(mrob, g[T, ST])
    @variable(mrob, g0[T, ST])
    @variable(mrob, g_[t=T, ST, tau=t0:t, SetSigma])
    # WARNING: in the robust model, the variable h(t,s) was eliminated, since the other constraints can
    #          be rewritten by using inequalities (22)
    ### @variable(mrob, h[T, ST] >= 0)
    # (14) extra amount of electricity requested by the client to the parter
    # (out of any engaged contract) in order to satisfy his needs at time period t
    ### @variable(mrob, e[T])
    @variable(mrob, e0[T])
    @variable(mrob, e_[t=T, tau=t0:t, SetSigma])
    # (15) percentage of time period p in which drivable system s is turned on
    ### @variable(mrob, x[T, D])
    @variable(mrob, x0[T, D])
    @variable(mrob, x_[t=T, D, tau=t0:t, SetSigma])

    # ==> Caso em que a desigualdade eh <= 0


    # FIXME Assumindo que r0[s] == r0[1,s]
    @variable(mrob, E)
    @variable(mrob, f0_34)
    @variable(mrob, f_34[SetTau,SetSigma])
    #@variable(mrob, r00[ST])
    @constraint(mrob, f0_34 == -E + sum(v_C[t,c]*y[t,c] + alpha[t,c]*q0[t,c] for t in T, c in 1:num_contracts[t])
            + sum( sum(v_D[s]*P_D[t,s]*x0[t,s] for s in D)
            + sum(v_ST[s]*(1+lambda[s])*g0[t,s] for s in ST) + beta[t]*e0[t] for t in T)
            + sum(v_ST[s]*(u0[s] - r0[t_bar,s]) for s in ST))
    #println(x)
    @constraint(mrob, [tau=SetTau, sigma=SetSigma], f_34[tau, sigma] == sum((  sum(alpha[t,c]*q_[t,c,tau,sigma] for c in 1:num_contracts[t]) +
                sum(v_D[s]*P_D[t,s]*x_[t,s,tau,sigma] for s in D) + sum(v_ST[s]*(1+lambda[s])*g_[t,s,tau,sigma] for s in ST) + beta[t]*e_[t,tau,sigma]
                )  for t in tau:(t_bar-1))
                - sum(v_ST[s]*r_[t_bar,s,tau,sigma] for s in ST) + v_NDU[sigma])
    #println(x)


    @variable(mrob, n_34[SetTau, SetSigma])

    @constraint(mrob, f0_34 + sum(v_bar[tau, sigma] * f_34[tau, sigma]
                + v_hat[tau, sigma] * n_34[tau, sigma]
            for tau in SetTau, sigma in SetSigma) <= 0)
    @constraint(mrob, [tau=SetTau, sigma=SetSigma], f_34[tau, sigma] >= -n_34[tau, sigma])
    @constraint(mrob, [tau=SetTau, sigma=SetSigma], f_34[tau, sigma] <= n_34[tau, sigma])
    #@constraint(mrob, [tau=SetTau, sigma=SetSigma], n_34[tau, sigma] >= 0)

    # Caso em que a desigualdade eh >= 0

    @variable(mrob, f0_35[T])
    @variable(mrob, f_35[t=T, tau=t0:t, SetSigma])
    # TODO FIXME incluir if no uso do indice c para pegar apenas os contratos relativos ao periodo c no tempo t
    @constraint(mrob, [t=T], f0_35[t] == sum(q0[t,c] for c in 1:num_contracts[t])
        + sum(P_D[t,s]*x0[t,s] for s in D)
        # FIXME Adapting to new balance constraint where h[t,s] is now multiplied by lambda[s]
        + sum(((lambda[s]*lambda[s]-1)*g0[t,s] + lambda[s] * (r0[t,s] - r0[t+1,s])) for s in ST) + e0[t] + sum(P_ND[t,s] for s in ND))
        ### + sum(((lambda[s]-1)*g0[t,s] + r0[t,s] - r0[t+1,s]) for s in ST) + e0[t] + sum(P_ND[t,s] for s in ND))
    @constraint(mrob, [t=T, tau=t0:t-1, sigma=SetSigma], f_35[t,tau,sigma] ==
        sum(q_[t,c,tau,sigma] for c in 1:num_contracts[t])
        + sum(P_D[t,s]*x_[t,s,tau,sigma] for s in D)
        # FIXME Adapting to new balance constraint where h[t,s] is now multiplied by lambda[s]
        + sum( ((lambda[s]*lambda[s]-1)*g_[t,s,tau,sigma] + lambda[s] * (r_[t,s,tau,sigma] - r_[t+1,s,tau,sigma]) ) for s in ST)
        ### + sum(((lambda[s]-1)*g_[t,s,tau,sigma] + r_[t,s,tau,sigma] - r_[t+1,s,tau,sigma]) for s in ST)
        + e_[t,tau,sigma])
    # Adaptacao da equacao acima para tau = t
    @constraint(mrob, [t=T, sigma=SetSigma], f_35[t,t,sigma] ==
        sum(q_[t,c,t,sigma] for c in 1:num_contracts[t])
        + 1
        + sum(P_D[t,s]*x_[t,s,t,sigma] for s in D)
        # FIXME Adapting to new balance constraint where h[t,s] is now multiplied by lambda[s]
        + sum(((lambda[s]*lambda[s]-1)*g_[t,s,t,sigma] - lambda[s] * r_[t+1,s,t,sigma]) for s in ST)
        ### + sum(((lambda[s]-1)*g_[t,s,t,sigma] - r_[t+1,s,t,sigma]) for s in ST)
        + e_[t,t,sigma])

    @variable(mrob, n_35[t=T, tau=t0:t, SetSigma])

    @constraint(mrob, [t=T], f0_35[t] + sum(v_bar[tau, sigma] * f_35[t, tau, sigma]
                - v_hat[tau, sigma] * n_35[t, tau, sigma]
            for tau in t0:t, sigma in SetSigma) >= 0)
    @constraint(mrob, [t=T, tau=t0:t, sigma=SetSigma], f_35[t, tau, sigma] >= -n_35[t, tau, sigma])
    @constraint(mrob, [t=T, tau=t0:t, sigma=SetSigma], f_35[t, tau, sigma] <= n_35[t, tau, sigma])

    #@constraint(mrob, [t=T, tau=SetTau, sigma=SetSigma], n_35[t, tau, sigma] >= 0)

    @variable(mrob, f0_36[T, ST])
    @variable(mrob, f_36[t=T, ST, tau=t0:t, SetSigma])
    @variable(mrob, n_36[t=T, ST, tau=t0:t, SetSigma])

    # Writing f_0 and f_ in function of the original problem variables
    @constraint(mrob, [t=T, s=ST], f0_36[t, s] == lambda[s]*g0[t,s] - r0[t+1,s])
    # TODO Fazer tau=t0:t ou tau=t0:t-1 ?
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_36[t,s,tau,sigma] == lambda[s]*g_[t,s,tau,sigma] - r_[t+1,s,tau,sigma])

    # Linear Decision Rules (ARC-L1)
    @constraint(mrob, [t=T, s=ST], 0 >= f0_36[t, s] + sum((v_bar[tau, sigma] * f_36[t, s, tau, sigma]
                + v_hat[tau, sigma] * n_36[t, s, tau, sigma])
            for tau in t0:t, sigma in SetSigma) )
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], -n_36[t, s, tau, sigma] <= f_36[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_36[t, s, tau, sigma] <= n_36[t, s, tau, sigma])

    # TODO o que fazer aqui ?
    #@variable(mrob, f0_37)
    #@constraint(mrob, [s=ST], r00[s] == u0[s])
    @constraint(mrob, [s=ST], r0[t0, s] == u0[s])  #  r00[s] == u0[s]
    @constraint(mrob, [s=ST, tau=t0:t0, sigma=SetSigma], r_[t0, s, tau, sigma] == 0)

    @variable(mrob, f0_38a[t0:nbT+1, ST])
    @variable(mrob, f_38a[t=t0:nbT+1, ST, tau=t0:t-1, SetSigma])
    @variable(mrob, n_38a[t=t0:nbT+1, ST, tau=t0:t-1, SetSigma])

    # ∀ t ∈ T ∪ { t̄}
    @constraint(mrob, [t=t0:nbT+1, s=ST], f0_38a[t, s] == u_min[s] - r0[t, s])
    # TODO Fazer tau=t0:t-1 ou tau=t0:t ?
    @constraint(mrob, [t=t0:nbT+1, s=ST, tau=t0:t-1, sigma=SetSigma], f_38a[t, s, tau, sigma] == -r_[t, s, tau, sigma])

    # TODO verificar o tau in 1:(nbTau-1)
    @constraint(mrob, [t=t0:nbT+1, s=ST], f0_38a[t, s] + sum((v_bar[tau, sigma] * f_38a[t, s, tau, sigma]
                + v_hat[tau, sigma] * n_38a[t, s, tau, sigma]) for tau in t0:t-1, sigma in SetSigma) <= 0 )
    @constraint(mrob, [t=t0:nbT+1, s=ST, tau=t0:t-1, sigma=SetSigma], -n_38a[t, s, tau, sigma] <= f_38a[t, s, tau, sigma])
    @constraint(mrob, [t=t0:nbT+1, s=ST, tau=t0:t-1, sigma=SetSigma], f_38a[t, s, tau, sigma] <= n_38a[t, s, tau, sigma])

    # TODO: Rever a questao da somatoria de Tau ser de 0 a t (como resolver o indice 0 em Julia)
    # TODO Resolver o T = T_{bar}!
    @variable(mrob, f0_38b[t0:nbT+1, ST])
    @variable(mrob, f_38b[t=t0:nbT+1, ST, tau=t0:t-1, SetSigma])
    @variable(mrob, n_38b[t=t0:nbT+1, ST, tau=t0:t-1, SetSigma])

    # ∀ t ∈ T ∪ { t̄}
    @constraint(mrob, [t=t0:nbT+1, s=ST], f0_38b[t, s] == r0[t,s] - u_max[s])
    # TODO Fazer tau=t0:t-1 ou tau=t0:t ?
    @constraint(mrob, [t=t0:nbT+1, s=ST, tau=t0:t-1, sigma=SetSigma], f_38b[t, s, tau, sigma] == r_[t, s, tau, sigma])

    # TODO verificar o tau in 1:(nbTau-1)
    @constraint(mrob, [t=t0:nbT+1, s=ST], f0_38b[t, s]
        + sum((v_bar[tau, sigma] * f_38b[t, s, tau, sigma] + v_hat[tau, sigma] * n_38b[t, s, tau, sigma])
            for tau in t0:t-1, sigma in SetSigma) <= 0 )
    @constraint(mrob, [t=t0:nbT+1, s=ST, tau=t0:t-1, sigma=SetSigma], -n_38b[t, s, tau, sigma] <= f_38b[t, s, tau, sigma])
    @constraint(mrob, [t=t0:nbT+1, s=ST, tau=t0:t-1, sigma=SetSigma], f_38b[t, s, tau, sigma] <= n_38b[t, s, tau, sigma])

    # TODO Corrigir a terceira linha abaixo no artigo em Latex

    @variable(mrob, f0_39[T, ST])
    @variable(mrob, f_39[t=T, ST, tau=t0:t, SetSigma])
    @variable(mrob, n_39[t=T, ST, tau=t0:t, SetSigma])

    # TODO Fazer tau=t0:t-1 ou tau=t0:t ?
    @constraint(mrob, [t=T, s=ST], f0_39[t, s] == lambda[s]*g0[t,s] + r0[t,s] - r0[t+1,s])
    @constraint(mrob, [t=T, s=ST, tau=t0:t-1, sigma=SetSigma], f_39[t, s, tau, sigma] == lambda[s]*g_[t,s,tau,sigma]
        + r_[t,s,tau,sigma] - r_[t+1,s,tau,sigma])
    @constraint(mrob, [t=T, s=ST, sigma=SetSigma], f_39[t, s, t, sigma] == lambda[s]*g_[t, s, t, sigma]
        - r_[t+1,s,t,sigma])

    @constraint(mrob, [t=T, s=ST], f0_39[t, s] + sum((v_bar[tau, sigma] * f_39[t, s, tau, sigma]
                - v_hat[tau, sigma] * n_39[t, s, tau, sigma])
            for tau in t0:t, sigma in SetSigma) >= 0 )
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], -n_39[t, s, tau, sigma] <= f_39[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_39[t, s, tau, sigma] <= n_39[t, s, tau, sigma])

    @variable(mrob, f0_40[T, ST])
    @variable(mrob, f_40[t=T, ST, tau=t0:t, SetSigma])
    @variable(mrob, n_40[t=T, ST, tau=t0:t, SetSigma])

    @constraint(mrob, [t=T, s=ST], f0_40[t, s] == g0[t,s] - theta_abs[s]*delta[t])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_40[t, s, tau, sigma] == g_[t,s,tau,sigma])

    @constraint(mrob, [t=T, s=ST], f0_40[t, s] + sum((v_bar[tau, sigma] * f_40[t, s, tau, sigma]
                + v_hat[tau, sigma] * n_40[t, s, tau, sigma])
            for tau in t0:t, sigma in SetSigma) <= 0 )
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], -n_40[t, s, tau, sigma] <= f_40[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_40[t, s, tau, sigma] <= n_40[t, s, tau, sigma])

    @variable(mrob, f0_41[T, ST])
    @variable(mrob, f_41[t=T, ST, tau=t0:t, SetSigma])
    @variable(mrob, n_41[t=T, ST, tau=t0:t, SetSigma])

    # TODO Fazer tau=t0:t-1 ou tau=t0:t ?
    @constraint(mrob, [t=T, s=ST], f0_41[t, s] == lambda[s]*g0[t,s] + r0[t,s] - r0[t+1,s] - theta_ref[s]*delta[t])
    @constraint(mrob, [t=T, s=ST, tau=t0:t-1, sigma=SetSigma], f_41[t, s, tau, sigma] == g_[t, s, tau, sigma] + r_[t, s, tau, sigma]
        - r_[t+1,s,tau,sigma])
    @constraint(mrob, [t=T, s=ST, sigma=SetSigma], f_41[t, s, t,sigma] == g_[t,s,t,sigma] - r_[t+1,s,t,sigma])

    @constraint(mrob, [t=T, s=ST], f0_41[t, s] + sum((v_bar[tau, sigma] * f_41[t, s, tau, sigma]
                + v_hat[tau, sigma] * n_41[t, s, tau, sigma])
            for tau in t0:t, sigma in SetSigma) <= 0)
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], -n_41[t, s, tau, sigma] <= f_41[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_41[t, s, tau, sigma] <= n_41[t, s, tau, sigma])

    # TODO nao eh possivel reduzir o numero de variaveis usando o indice c in C diretamente! Modificar para usar c e t
    # f_43a e f_43b foram removidos para economizar variaveis => reaproveitando a f_42a e f_42b
    @variable(mrob, f0_42a[t=T, c=1:num_contracts[t]])
    @variable(mrob, f_42a[t=T, c=1:num_contracts[t], tau=t0:t, SetSigma])
    @variable(mrob, n_42a[t=T, c=1:num_contracts[t], tau=t0:t, SetSigma])

    for t in T, c in 1:num_contracts[t]
        if pi_plus[t, c] > 0
            @constraint(mrob, f0_42a[t,c] == -pi_minus[t,c]*y[t,c] + q0[t,c])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42a[t, c, tau, sigma] == q_[t, c, tau, sigma])

            # Caso em que a desigualdade eh >= 0
            @constraint(mrob, f0_42a[t, c] + sum(v_bar[tau, sigma] * f_42a[t, c, tau, sigma]
                        - v_hat[tau, sigma] * n_42a[t, c, tau, sigma]
                    for tau in t0:t, sigma in SetSigma) >= 0)
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42a[t, c, tau, sigma] >= -n_42a[t, c, tau, sigma])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42a[t, c, tau, sigma] <= n_42a[t, c, tau, sigma])
        else  # if contract[c,:max_period] < 0 || contract[c,:min_period] < 0
            @constraint(mrob, f0_42a[t, c] == -pi_plus[t, c]*y[t, c] + q0[t, c])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42a[t, c, tau, sigma] == q_[t, c, tau, sigma])

            # Caso em que a desigualdade eh >= 0
            @constraint(mrob, f0_42a[t, c] + sum((v_bar[tau, sigma] * f_42a[t, c, tau, sigma]
                        - v_hat[tau, sigma] * n_42a[t, c, tau, sigma])
                    for tau in t0:t, sigma in SetSigma) >= 0)
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], -n_42a[t, c, tau, sigma] <= f_42a[t, c, tau, sigma])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42a[t, c, tau, sigma] <= n_42a[t, c, tau, sigma])
        end
    end

    @variable(mrob, f0_42b[t=T, c=1:num_contracts[t]])
    @variable(mrob, f_42b[t=T, c=1:num_contracts[t], tau=t0:t, SetSigma])
    @variable(mrob, n_42b[t=T, c=1:num_contracts[t], tau=t0:t, SetSigma])

    for t in T, c in 1:num_contracts[t]
        if pi_plus[t, c] > 0
            @constraint(mrob, f0_42b[t, c] == -pi_plus[t, c]*y[t, c] + q0[t, c])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42b[t, c, tau, sigma] == q_[t, c, tau, sigma])
        else  # if contract[c,:max_period] < 0 || contract[c,:min_period] < 0
            @constraint(mrob, f0_42b[t, c] == -pi_minus[t, c]*y[t, c] + q0[t, c])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42b[t, c, tau, sigma] == q_[t, c, tau, sigma])
        end
        # Caso em que a desigualdade eh <= 0
        @constraint(mrob, f0_42b[t, c] + sum(v_bar[tau, sigma] * f_42b[t, c, tau, sigma]
                    + v_hat[tau, sigma] * n_42b[t, c, tau, sigma]
                for tau in t0:t, sigma in SetSigma) <= 0)
        @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42b[t, c, tau, sigma] >= -n_42b[t, c, tau, sigma])
        @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_42b[t, c, tau, sigma] <= n_42b[t, c, tau, sigma])
    end

    @variable(mrob, f0_44[T, D])
    @variable(mrob, f_44[t=T, D, tau=t0:t, SetSigma])
    @variable(mrob, n_44[t=T, D, tau=t0:t, SetSigma])

    # f_45 foi removida para economizar variaveis => reaproveitando na f_44
    for t in T, s in D
        previous_charge = 0.0
        if size(previous_drivable_charge, 1) > 0
            previous_charge = previous_drivable_charge[s]
        end
        if P_D[t,s] > 0
            if Pmin[t,s] > Pmax[t,s]
                println(general_logger, "[Robust Model, t0 = $(t0)] Assertion error : [drivable D] P[t,s] > 0 => Pmin <= Pmax")
            end
            @assert Pmin[t,s] <= Pmax[t,s] "[drivable D] P[t,s] > 0 => Pmin <= Pmax"
            @constraint(mrob, f0_44[t, s] == previous_charge + sum(P_D[t1,s]*x0[t1,s] for t1 in t0:t) - Pmin[t,s])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_44[t, s, tau, sigma] ==
                sum(P_D[t1,s] * x_[t1,s,tau,sigma] for t1 in tau:t))

            @constraint(mrob, f0_44[t, s] + sum((v_bar[tau, sigma] * f_44[t, s, tau, sigma]
                        - v_hat[tau, sigma] * n_44[t, s, tau, sigma])
                    for tau in t0:t, sigma in SetSigma) >= 0)
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], -n_44[t, s, tau, sigma] <= f_44[t, s, tau, sigma])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_44[t, s, tau, sigma] <= n_44[t, s, tau, sigma])
        else
            if Pmin[t,s] < Pmax[t,s]
                println(general_logger, "[Robust Model, t0 = $(t0)] Assertion error : [drivable D] P[t,s] < 0 => Pmin >= Pmax")
            end
            @assert Pmin[t,s] >= Pmax[t,s] "[drivable D] P[t,s] < 0 => Pmin >= Pmax"
            @constraint(mrob, f0_44[t, s] == previous_charge + sum(P_D[t1,s]*x0[t1,s] for t1 in t0:t) - Pmin[t,s])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_44[t, s, tau,sigma] ==
                sum(P_D[t1,s] * x_[t1,s,tau,sigma] for t1 in tau:t))

            @constraint(mrob, f0_44[t, s] + sum((v_bar[tau, sigma] * f_44[t, s, tau, sigma]
                        + v_hat[tau, sigma] * n_44[t, s, tau, sigma])
                    for tau in t0:t, sigma in SetSigma) <= 0)
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], -n_44[t, s, tau, sigma] <= f_44[t, s, tau, sigma])
            @constraint(mrob, [tau=t0:t, sigma=SetSigma], f_44[t, s, tau, sigma] <= n_44[t, s, tau, sigma])
        end
    end
    # Already defined in previous variable definiton

    # The following constraint must only be valid from the 2nd time period on
    @variable(mrob, f0_47a[t=(t0+1):(nbT+1), s=ST])
    @variable(mrob, f_47a[t=(t0+1):(nbT+1), s=ST, tau=t0:t-1, SetSigma])
    @variable(mrob, n_47a[t=(t0+1):(nbT+1), s=ST, tau=t0:t-1, SetSigma])

    @constraint(mrob, [t=(t0+1):(nbT+1), s=ST], f0_47a[t,s] == r0[t,s])
    @constraint(mrob, [t=(t0+1):(nbT+1), s=ST, tau=t0:t-1, sigma=SetSigma], f_47a[t,s,tau,sigma] == r_[t,s,tau,sigma])

    @constraint(mrob, [t=(t0+1):(nbT+1), s=ST], f0_47a[t, s] + sum((v_bar[tau, sigma] * f_47a[t, s, tau, sigma]
                - v_hat[tau, sigma] * n_47a[t, s, tau, sigma])
            for tau in t0:t-1, sigma in SetSigma) >= 0)
    @constraint(mrob, [t=(t0+1):(nbT+1), s=ST, tau=t0:t-1, sigma=SetSigma], -n_47a[t, s, tau, sigma] <= f_47a[t, s, tau, sigma])
    @constraint(mrob, [t=(t0+1):(nbT+1), s=ST, tau=t0:t-1, sigma=SetSigma], f_47a[t, s, tau, sigma] <= n_47a[t, s, tau, sigma])

    @variable(mrob, f0_47b[t=T, s=ST])
    @variable(mrob, f_47b[t=T, s=ST, tau=t0:t, SetSigma])
    @variable(mrob, n_47b[t=T, s=ST, tau=t0:t, SetSigma])

    @constraint(mrob, [t=T, s=ST], f0_47b[t,s] == g0[t,s])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_47b[t,s,tau,sigma] == g_[t,s,tau,sigma])

    @constraint(mrob, [t=T, s=ST], f0_47b[t, s] + sum((v_bar[tau, sigma] * f_47b[t, s, tau, sigma]
                - v_hat[tau, sigma] * n_47b[t, s, tau, sigma])
            for tau in t0:t, sigma in SetSigma) >= 0)
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], -n_47b[t, s, tau, sigma] <= f_47b[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=ST, tau=t0:t, sigma=SetSigma], f_47b[t, s, tau, sigma] <= n_47b[t, s, tau, sigma])

    @variable(mrob, f0_48[T])
    @variable(mrob, f_48[t=T, tau=t0:t, SetSigma])
    @variable(mrob, n_48[t=T, tau=t0:t, SetSigma])

    @constraint(mrob, [t=T], f0_48[t] == e0[t])
    @constraint(mrob, [t=T, tau=tau=t0:t, sigma=SetSigma], f_48[t,tau,sigma] == e_[t,tau,sigma])

    @constraint(mrob, [t=T], f0_48[t] + sum(v_bar[tau, sigma] * f_48[t, tau, sigma]
                - v_hat[tau, sigma] * n_48[t, tau, sigma]
            for tau in t0:t, sigma in SetSigma) >= 0)
    @constraint(mrob, [t=T, tau=t0:t, sigma=SetSigma], f_48[t, tau, sigma] >= -n_48[t, tau, sigma])
    @constraint(mrob, [t=T, tau=t0:t, sigma=SetSigma], f_48[t, tau, sigma] <= n_48[t, tau, sigma])
    #@constraint(mrob, [t=T, tau=SetTau, sigma=SetSigma], n_48[t, tau, sigma] >= 0)

    #x = @constraint(mrob, [t=T], e0[t] >= 0)
    #println(x)
    #x = @constraint(mrob, [t=T, tau=t0:t, sigma=SetSigma], e_[t,tau,sigma] >= 0)
    #println(x)

    @variable(mrob, f0_49a[T, D])
    @variable(mrob, f_49a[t=T, D, tau=t0:t, SetSigma])
    @variable(mrob, n_49a[t=T, D, tau=t0:t, SetSigma])

    @constraint(mrob, [t=T, s=D], f0_49a[t, s] == x0[t,s] - 1)
    @constraint(mrob, [t=T, s=D, tau=t0:t, sigma=SetSigma], f_49a[t,s,tau,sigma] == x_[t,s,tau,sigma])

    @constraint(mrob, [t=T, s=D], f0_49a[t, s] + sum((v_bar[tau,sigma] * f_49a[t, s, tau, sigma]
                + v_hat[tau,sigma] * n_49a[t, s, tau, sigma])
                    for tau in t0:t, sigma in SetSigma) <= 0)
    @constraint(mrob, [t=T, s=D, tau=t0:t, sigma=SetSigma], -n_49a[t, s, tau, sigma] <= f_49a[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=D, tau=t0:t, sigma=SetSigma], f_49a[t, s, tau, sigma] <= n_49a[t, s, tau, sigma])

    @variable(mrob, f0_49b[T, D])
    @variable(mrob, f_49b[t=T, D, tau=t0:t, SetSigma])
    @variable(mrob, n_49b[t=T, D, tau=t0:t, SetSigma])

    @constraint(mrob, [t=T, s=D], f0_49b[t, s] == x0[t,s])
    @constraint(mrob, [t=T, s=D, tau=t0:t, sigma=SetSigma], f_49b[t,s,tau,sigma] == x_[t,s,tau,sigma])

    @constraint(mrob, [t=T, s=D], f0_49b[t, s] + sum((v_bar[tau,sigma] * f_49b[t, s, tau, sigma]
                - v_hat[tau,sigma] * n_49b[t, s, tau, sigma])
            for tau in t0:t, sigma in SetSigma) >= 0)
    @constraint(mrob, [t=T, s=D, tau=t0:t, sigma=SetSigma], -n_49b[t, s, tau, sigma] <= f_49b[t, s, tau, sigma])
    @constraint(mrob, [t=T, s=D, tau=t0:t, sigma=SetSigma], f_49b[t, s, tau, sigma] <= n_49b[t, s, tau, sigma])

    @objective(mrob, Min, E)
    println("===================================================================\nSolving robust model for t0 = $(t0)...")
    status = optimize!(mrob)
    solution = Dict()
    sol_y = zeros(Int64, nbT, max_contracts_per_period)
    max_contracts_per_period = 0
    for t in 1:nbT
        max_contracts_per_period = max(max_contracts_per_period, num_contracts[t])
    end
    if JuMP.termination_status(mrob) == MOI.OPTIMAL
        rel_gap = MOI.get(mrob, MOI.RelativeGap())
        println("\n===========  R O B U S T    B O X    S O L U T I O N  ( t0 = $(t0) ) ===========\n")
        println("Optimal Objective Function value : ", objective_value(mrob))
        println("Relative gap = $(rel_gap)")
        println("Solve time : ", solve_time(mrob))
        var_y = value.(y)
        if t0 > 1  # if starting from period > 1
            for t in 1:(t0-1), c in 1:num_contracts[t]
                if size(fixed_contracts, 1) > 0
                    sol_y[t,c] = fixed_contracts[t,c]
                else
                    sol_y[t,c] = -1
                    println("[Robust Model, t0 = $(t0)] WARN: Returning invalid contract solution y. Check RCCP Robust Model procedure !")
                    println(general_logger, "WARN: Returning invalid contract solution y. Check RCCP Robust Model procedure !")
                end
            end
        end
        for t in T, c in 1:num_contracts[t]
            sol_y[t,c] = trunc(Int, value(y[t, c]))
        end
        solution["y"] = sol_y
        # println("\namount consumed (< 0) / provided (> 0) by the client in contract i")
        solution["q0"] = capture_variable_solution_as_array_2d_custom_q(q0, t0, nbT, num_contracts)
        solution["q_"] = capture_variable_solution_as_array_4d_custom_q(q_, t0, nbT, num_contracts, SetSigma)
        # println("\namount stored in system s at the beginning of time period t ")
        solution["r0"] = capture_variable_solution_as_array_2d_custom_r(r0, t0, nbT, ST)
        solution["r_"] = capture_variable_solution_as_array_4d_custom_r(r_, t0, nbT, ST, SetSigma)
        # println("\namount absorbed by system s during the time period t and that will be already stored in system s at the end of this time period")
        solution["g0"] = capture_variable_solution_as_array_2d(g0, t0, nbT, ST)
        solution["g_"] = capture_variable_solution_as_array_4d(g_, t0, nbT, ST, SetSigma)
        # println("\nextra amount of electricity requested by the client to the parter (out of any engaged contract) in order to satisfy his needs at time period t")
        solution["e0"] = capture_variable_solution_as_array_1d(e0, t0, nbT)
        solution["e_"] = capture_variable_solution_as_array_3d(e_, t0, nbT, SetSigma)
        # println("\npercentage of time period p in which drivable system s is turned on")
        solution["x0"] = capture_variable_solution_as_array_2d(x0, t0, nbT, D)
        solution["x_"] = capture_variable_solution_as_array_4d(x_, t0, nbT, D, SetSigma)
        return objective_value(mrob), solution, solve_time(mrob), true, rel_gap
    elseif has_values(mrob) && (JuMP.termination_status(mrob) in [MOI.TIME_LIMIT, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.LOCALLY_SOLVED])
        # Recover the best integer solution found (suboptimal) and return it
        rel_gap = MOI.get(mrob, MOI.RelativeGap())
        println("\n===========  R O B U S T    B O X    S O L U T I O N  ( t0 = $(t0) ) ===========\n")
        println("Optimal Objective Function value : ", objective_value(mrob))
        println("Relative gap = $(rel_gap)")
        println("Solve time : ", solve_time(mrob))
        println(general_logger, "[Robust Box Model, t0 = $(t0)] WARN: Time limit exceeded. Optimal solution not found! Best bound: ", objective_bound(mrob))
        var_y = value.(y)
        if t0 > 1  # if starting from period > 1
            for t in 1:(t0-1), c in 1:num_contracts[t]
                if size(fixed_contracts, 1) > 0
                    sol_y[t,c] = fixed_contracts[t,c]
                else
                    sol_y[t,c] = -1
                    println("[Robust Model, t0 = $(t0)] WARN: Returning invalid contract solution y. Check RCCP Robust Model procedure !")
                    println(general_logger, "WARN: Returning invalid contract solution y. Check RCCP Robust Model procedure !")
                end
            end
        end
        for t in T, c in 1:num_contracts[t]
            sol_y[t,c] = trunc(Int, value(y[t, c]))
        end
        solution["y"] = sol_y
        # println("\namount consumed (< 0) / provided (> 0) by the client in contract i")
        solution["q0"] = capture_variable_solution_as_array_2d_custom_q(q0, t0, nbT, num_contracts)
        solution["q_"] = capture_variable_solution_as_array_4d_custom_q(q_, t0, nbT, num_contracts, SetSigma)
        # println("\namount stored in system s at the beginning of time period t ")
        solution["r0"] = capture_variable_solution_as_array_2d_custom_r(r0, t0, nbT, ST)
        solution["r_"] = capture_variable_solution_as_array_4d_custom_r(r_, t0, nbT, ST, SetSigma)
        # println("\namount absorbed by system s during the time period t and that will be already stored in system s at the end of this time period")
        solution["g0"] = capture_variable_solution_as_array_2d(g0, t0, nbT, ST)
        solution["g_"] = capture_variable_solution_as_array_4d(g_, t0, nbT, ST, SetSigma)
        # println("\nextra amount of electricity requested by the client to the parter (out of any engaged contract) in order to satisfy his needs at time period t")
        solution["e0"] = capture_variable_solution_as_array_1d(e0, t0, nbT)
        solution["e_"] = capture_variable_solution_as_array_3d(e_, t0, nbT, SetSigma)
        # println("\npercentage of time period p in which drivable system s is turned on")
        solution["x0"] = capture_variable_solution_as_array_2d(x0, t0, nbT, D)
        solution["x_"] = capture_variable_solution_as_array_4d(x_, t0, nbT, D, SetSigma)
        return objective_value(mrob), solution, solve_time(mrob), true, rel_gap
    else
        println(infeasible_logger, "t0 = $(t0); Infeasible robust model : $(filepath)")
        println("Optimal solution not found! Best bound: ", objective_bound(mrob))
        println(general_logger, "WARN: Optimal solution not found! Best bound: ", objective_bound(mrob))
        solution["y"] = sol_y
        return objective_bound(mrob), solution, solve_time(mrob), false, Inf
    end
end
