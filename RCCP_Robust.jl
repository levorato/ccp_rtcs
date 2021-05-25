# =================================
# RCCP_Robust.jl
# RCCP Robust MILP Model Wrapper
# =================================

include("RCCP_Robust_Box.jl")
include("RCCP_Robust_Budget.jl")

function solve_robust_model(instance_as_dict::Dict, general_logger, infeasible_logger, solver_params, t0 = 1,
                                fixed_contracts = Int64[],
                                initial_battery = Float64[], previous_drivable_charge = Float64[])
    if solver_params["model"] == "robust-box"
        return solve_robust_model_box(instance_as_dict, general_logger, infeasible_logger, solver_params, t0,
                                        fixed_contracts, initial_battery, previous_drivable_charge)
    elseif solver_params["model"] == "robust-budget"
        return solve_robust_model_budget(instance_as_dict, general_logger, infeasible_logger, solver_params, t0,
                                        fixed_contracts, initial_battery, previous_drivable_charge)
    else
        println("ERROR: Invalid robust model type! => $(solver_params["model"])")
        flush(stdout)
    end
end
