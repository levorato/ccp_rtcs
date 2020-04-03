# =================================
# RCCP_det_Tests.jl
# RCCP Deterministic Tests
# =================================

# Include RCCP Deterministic model
include("RCCP_det.jl")

#datafile = "../notebooks/data/instance5-contratos_restituicao-ultimo_periodo.txt"
#solve_deterministic_model(datafile)

#datafile = "../notebooks/data/antoine/A_instance2_11scen_1NDU.txt"
#datafile = "../notebooks/data/antoine/A_instance2_11scen.txt"
datafile = "./instances/Example_3.1_A.txt"
#solve_deterministic_model(datafile, 1800, false)
general_logger = open("/tmp/rccp_concat_log.txt", "w+")
infeasible_logger = open("/tmp/rccp_concat_log_infeasible.txt", "w+")
instance_as_dict = read_input_data(datafile, false)
solve_deterministic_model_with_t(instance_as_dict, general_logger, infeasible_logger)
close(general_logger)
