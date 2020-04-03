#  RCCP_FileReader.jl
using CSV
using DataFrames
#using DataArrays
using CPLEX
using JLD2, FileIO
using ZipFile

# Epsilon used in floating point comparisons
EPS = 1e-5
ZERO = zero(1.0)
# Rounding error allowed in RTCS Operations (e.g. battery and drivable storage / retrieval)
RTCS_ROUNDING_ERROR = 0.01
CPLEX_OPT_GAP = 0.02  # gap of cplex adjusted to 2%

function create_full_dir(basepath, folders)
    fullpath = basepath
    for folder in folders
        fullpath = joinpath(fullpath, folder)
        if ! isdir(fullpath)
            mkdir(fullpath)
        end
    end
    return fullpath
end

# Method to read Xpress input files
function read_xpress_file(filepath)
    # Read instance file
    file = open(filepath)
    text = ""
    count = 0
    for ln in eachline(file)
        if length(ln) > 0 && (! contains(ln, "!"))  # ignore comments and empty lines
            # remove trailing (num) containing line number of array
            ln = replace(ln, r"\((\d+)\)", s"")
            # remove character delimiters of array (left and right brackets)
            #ln = replace(ln, r"\[", "")
            ln = replace(ln, r"\]", s"")
            ln = lstrip(ln)
            #println("$(count): $(ln)")
            if length(ln) > 0
                text = text * ln
                #if ! contains(ln, r"(\w+)\:")
                text = text * "\n"
                #end
            end
            #println("$(ln)")
        end
        count += 1
    end
    #println("$(text)")
    input_structs = split(text, r"(\w+)\:")
    # Read each struct from file (period, contract, drivable, ndrivable and storage)
    count = 0
    instance = []
    period = DataFrame()
    contract = DataFrame()
    drivable = DataFrame()
    n_drivable = DataFrame()
    storage = DataFrame()
    for structs in input_structs
       #println("$(count):")
       if length(structs) > 0
            #println("$(structs)")
            # Read each record from the struct (array, integer, etc.)
            records = split(structs, r"\[")
            #println("$(records)")
            if count != 4 && count != 5
                record_as_string = ""
                for record in records
                    item = strip(record)
                    if length(item) > 0
                        record_as_string *= item * "\n"
                        #println("$(item)")
                    end
                end
                #println("Record: $(record_as_string)")
            end
            if count == 1
                instance = readdlm(IOBuffer(record_as_string), Int64)
                println("Instance: $(instance)")
            elseif count == 2
                arr = readdlm(IOBuffer(record_as_string))
                period = convert(DataFrame, arr)
                rename!(period, :x1 => :size, :x2 => :cost_out)
                println("Period: $(period)")
            elseif count == 3
                arr = readdlm(IOBuffer(record_as_string))
                contract = convert(DataFrame, arr)
                rename!(contract, :x1 => :period, :x2 => :cost_fix, :x3 => :cost_var,
                    :x4 => :min_period, :x5 => :max_period, :x6 => :min_delta, :x7 => :max_delta)
                # Convert period to Integer
                #contract[:period] = convert(DataVector{Int64}, contract[:period])
                contract[:period] = map(x -> convert(Int64, x), contract[:period])
                println("Contract: $(contract)")
            elseif count == 4
                record_as_string = String["", "", ""]
                item_count = 0
                for record in records
                    item = strip(record)
                    if length(item) > 0
                        record_as_string[item_count % 3 + 1] *= item * "\n"
                        item_count += 1
                    end
                end
                #println("Record 4: $(record_as_string[0])")
                drivable_cost = readdlm(IOBuffer(record_as_string[1]))
                drivable_pORc = readdlm(IOBuffer(record_as_string[2]))
                drivable_pORc_min = readdlm(IOBuffer(record_as_string[3]))
                drivable = convert(DataFrame, drivable_cost)
                rename!(drivable, :x1 => :cost)
                drivable[:pORc] = [ drivable_pORc[i, :] for i in 1:size(drivable_pORc, 1) ]
                drivable[:pORc_min] = [ drivable_pORc_min[i, :] for i in 1:size(drivable_pORc_min, 1) ]
                println("Drivable: $(drivable)")
            elseif count == 5
                record_as_string = ["", ""]
                item_count = 0
                for record in records
                    item = strip(record)
                    if length(item) > 0
                        if item_count % 2 == 0
                           record_as_string[1] *= item * "\n"
                        else
                           record_as_string[2] *= item * "\n"
                        end
                        #println("$(item)")
                        item_count += 1
                    end
                end
                #println("Record 5: $(record_as_string)")
                n_drivable_cost = readdlm(IOBuffer(record_as_string[1]))
                n_drivable_pORc = readdlm(IOBuffer(record_as_string[2]))
                n_drivable = convert(DataFrame, n_drivable_cost)
                rename!(n_drivable, :x1 => :cost)
                n_drivable[:pORc] = [ n_drivable_pORc[i, :] for i in 1:size(n_drivable_pORc, 1) ]
                println("Non Drivable: $(n_drivable)")
                #println("Non Drivable Cost: $(n_drivable_cost)")
                #println("Non Drivable pORc: $(n_drivable_pORc)")
            elseif count == 6
                arr = readdlm(IOBuffer(record_as_string))
                storage = convert(DataFrame, arr)
                rename!(storage, :x1 => :cost, :x2 => :uMin, :x3 => :uMax,
                    :x4 => :uInit, :x5 => :lostCoef, :x6 => :maxAbsorption, :x7 => :maxRefund)
                println("Storage: $(storage)")
            end
       end
       count += 1
    end
    println("\n\nInstance read successfully.\n\n")

    instance_as_dict = Dict()
    instance_as_dict["instance"] = instance
    instance_as_dict["period"] = period
    instance_as_dict["contract"] = contract
    instance_as_dict["drivable"] = drivable
    instance_as_dict["n_drivable"] = n_drivable
    instance_as_dict["n_drivable_uncertain"] = Array[]
    instance_as_dict["storage"] = storage
    instance_as_dict["scenarios"] = Array[]
    return instance_as_dict
end

function read_tabulated_data(filepath, verbose = false)
    # Read instance file
    println("Processing input file: $(filepath)...")
    file = open(filepath)
    text = ""
    bag = Dict()
    bag["Period"] = Array[]
    bag["Contract"] = Array[]
    bag["DrivableS"] = Array[]
    bag["NDrivableS"] = Array[]
    bag["Storage"] = Array[]
    bag["UNDS"] = Array[]
    bag["Scenarios"] = Array[]
    lines = String[]
    for ln in eachline(file)
        if length(ln) > 0 && (! contains(ln, "//"))  # ignore comments and empty lines
            ln = strip(ln)
            #println("$(count): $(ln)")
            if length(ln) > 0
                push!(lines, ln)
            end
        end
    end
    count = 1
    reading_scenario = false
    numberOfPeriods = 0
    while count <= length(lines)
        ln = lines[count]
        ln = replace(ln, r" +" => " ")  # remove duplicate spaces
        input = split(ln, r" |\t")
        filter!(x->x!="", input)
        #println("Split: $(input)")
        if contains(ln, "Instance")
            numberOfPeriods = parse(Int, input[2])
            nbOfSystems = parse(Int, input[3])
            bag["Instance"] = [numberOfPeriods, nbOfSystems]
            println("Number of periods is $(numberOfPeriods)")
            println("Number of systems is $(nbOfSystems)")
        elseif contains(ln, "Horizon")
            bag["Horizon"] = ln
        elseif contains(ln, "Period") && (! reading_scenario)
            duration = parse(Float64, input[3])
            costOutOfContract = parse(Float64, input[4])
            push!(bag["Period"], [duration, costOutOfContract])
        elseif contains(ln, "Contract")
            C = input[3:end]
            array2 = [parse(Float64, c) for c in C]
            contracts = [parse(Float64, input[2]); array2]
            #println("Contracts: $(contracts)")
            push!(bag["Contract"], contracts)
        elseif contains(ln, "NDrivableS")
            dev_name_and_cost = split(strip(ln), r" |\t")
            filter!(x->x!="", dev_name_and_cost)
            dev_name = dev_name_and_cost[2]
            cost = parse(Float64, dev_name_and_cost[3])
            # Capture the next line containing P
            P = split(lines[count + 1], r" |\t")[2:end]
            filter!(x->x!="", P)
            count += 1
            push!(bag["NDrivableS"], [dev_name, cost, [parse(Int, x) for x in P]])
        elseif contains(ln, "DrivableS")
            ln = strip(ln)
            dev_name_and_cost = split(strip(ln), r" |\t")
            filter!(x->x!="", dev_name_and_cost)
            dev_name = dev_name_and_cost[2]
            cost = parse(Float64, dev_name_and_cost[3])
            # Capture the next 3 lines containing P, Pmin and Pmax
            P = split(strip(lines[count + 1]), r" |\t")[2:end]
            filter!(x->x!="", P)
            Pmin = split(strip(lines[count + 2]), r" |\t")[2:end]
            filter!(x->x!="", Pmin)
            Pmax = split(strip(lines[count + 3]), r" |\t")[2:end]
            filter!(x->x!="", Pmax)
            count += 3
            push!(bag["DrivableS"], [dev_name, cost, [parse(Float64, x) for x in P], [parse(Float64, x) for x in Pmin],
                [parse(Float64, x) for x in Pmax]])
        elseif contains(ln, "Storage")
            S = input[3:end]
            array2 = [input[2]; [parse(Float64, s) for s in S]]
            push!(bag["Storage"], array2)
        elseif contains(ln, "UNDS") && (! reading_scenario )
            dev_name_and_cost = split(strip(ln), r" |\t")
            filter!(x->x!="", dev_name_and_cost)
            #println("$(ln) : $(dev_name_and_cost)")
            dev_name = dev_name_and_cost[2]
            cost = parse(Float64, dev_name_and_cost[3])
            # Capture the next 4 lines containing Pmin, Pmax, Pdt_min and Pdt_max
            ln = replace(strip(lines[count + 1]), r" +|\t+" => " ")
            Pmin = split(ln, r" |\t")[2:end]
            filter!(x->x!="", Pmin)
            ln = replace(strip(lines[count + 2]), r" +|\t+" => " ")
            Pmax = split(ln, r" |\t")[2:end]
            filter!(x->x!="", Pmax)
            ln = replace(strip(lines[count + 3]), r" +|\t+" => " ")
            Pdt_min = split(ln, r" |\t")[2:end]
            filter!(x->x!="", Pdt_min)
            ln = replace(strip(lines[count + 4]), r" +|\t+" => " ")
            Pdt_max = split(ln, r" |\t")[2:end]
            filter!(x->x!="", Pdt_max)
            count += 4
            push!(bag["UNDS"], [dev_name, cost, [parse(Float64, x) for x in Pmin], [parse(Float64, x) for x in Pmax],
                [parse(Float64, x) for x in Pdt_min], [parse(Float64, x) for x in Pdt_max]])
        elseif contains(ln, "NBScenarios")
            num_scenarios = parse(Int, input[2])
            println("The number of scenarios is $(num_scenarios)")
            unds_count = length(bag["UNDS"])
            println("The number of uncertain devices is $(unds_count)")
            reading_scenario = true
            local_count = 1
            # For each UNDS, capture the next num_period lines containing scenario data
            for i in 1:num_scenarios
                local_count += 1
                if verbose
                    println("Scenario #$(i)")
                end
                scenario_vect = Array[]
                for j in 1:unds_count
                    if verbose
                        println("UNDS #$(j)")
                    end
                    local_count += 1
                    matrix = Array[]
                    for k in 1:numberOfPeriods
                        ln = lines[count + local_count]
                        input = split(ln, r" |\t")
                        filter!(x->x!="", input)
                        local_count += 1
                        matrix_row = [input[i] for i in 3:length(input)]
                        martrix_row_as_float = [parse(Float64, elem) for elem in matrix_row]
                        push!(matrix, martrix_row_as_float)
                    end
                    push!(scenario_vect, matrix)
                    if verbose
                        println("Scenario $(i) x UNDS $(j): $(matrix)")
                    end
                end
                push!(bag["Scenarios"], scenario_vect)
            end
            count += local_count
        end
        count += 1
    end
    #println("\nPeriod data: $(bag["Period"])\n")
    #println("DrivableS: $(bag["DrivableS"])\n")
    #println("NDrivableS: $(bag["NDrivableS"])\n")
    if verbose
        println("UNDS: $(bag["UNDS"])\n")
    end
    #println("Storage: $(bag["Storage"])\n")
    #println("Contract: $(bag["Contract"])\n")

    # Read each struct from Dicts (period, contract, drivable, ndrivable and storage)
    instance = []
    period = DataFrame(size = Float64[], cost_out = Float64[])
    contract = DataFrame(period = Float64[], cost_fix = Float64[], cost_var = Float64[],
                    min_period = Float64[], max_period = Float64[], min_delta = Float64[], max_delta = Float64[])
    drivable = DataFrame(name = String[], cost = Float64[], pORc = Array[], pORc_min = Array[], pORc_max = Array[])
    n_drivable = DataFrame(name = String[], cost = Float64[], pORc = Array[])
    n_drivable_uncertain = DataFrame(name = String[], cost = Float64[], Pmin = Array[], Pmax = Array[],
                Pdt_min = Array[], Pdt_max = Array[])
    storage = DataFrame(name = String[], cost = Float64[], uMin = Float64[], uMax = Float64[],
                uInit = Float64[], lostCoef = Float64[], maxAbsorption = Float64[], maxRefund = Float64[])

    # Instance
    # time period T
    nbT = length(bag["Period"])
    # desconhecido?
    nbS = -1
    # number of contracts
    nbC = length(bag["Contract"])
    # number of Drivable devices
    nbD = length(bag["DrivableS"])
    # number of Non Drivable devices
    nbND = length(bag["NDrivableS"])
    # number of storage systems
    nbSt = length(bag["Storage"])
    instance = [nbT, nbS, nbC, nbD, nbND, nbSt]
    println("\nInstance: $(instance)")
    nbNDU = size(bag["UNDS"])
    println("\nnbNDU: $(nbNDU)")

    # Period
    for i in 1:nbT
        push!(period, bag["Period"][i]')
    end
    # Convert period size to Integer
    #period[:size] = convert(DataVector{Int64}, period[:size])
    period[:size] = map(x -> convert(Int64, x), period[:size])
    if verbose
        println("\nPeriod:")
        showall(period)
    end

    # Contract
    for i in 1:nbC
        push!(contract, bag["Contract"][i]')
        # Assertions on contract
        @assert contract[i, :min_period] <= contract[i, :max_period]
        @assert contract[i, :min_delta] <= contract[i, :max_delta]
        # Add 1 to all contract period numbers (they must begin at 1, not 0!)
        contract[i, :period] += 1
    end
    # Convert period to Integer
    #contract[:period] = convert(DataVector{Int64}, contract[:period])
    contract[:period] = map(x -> convert(Int64, x), contract[:period])
    if verbose
        println("\nContract:")
        showall(contract)
    end

    # Drivable
    for i in 1:nbD
        x = bag["DrivableS"][i]
        reshape(x, (1, :))
        push!(drivable, x)
    end
    if verbose
        println("\nDrivable:")
        showall(drivable)
    end

    # Non Drivable
    for i in 1:nbND
        x = bag["NDrivableS"][i]
        reshape(x, (1, :))
        push!(n_drivable, x)
    end
    if verbose
        println("\nNon Drivable:")
        showall(n_drivable)
    end

    # Uncertain non drivable
    for i in 1:length(bag["UNDS"])
        x = bag["UNDS"][i]
        reshape(x, (1, :))
        #println("x: $(x)")
        push!(n_drivable_uncertain, x)
    end
    if verbose
        println("\nUncertain Non Drivable:")
        showall(n_drivable_uncertain)
    end

    # Storage
    for i in 1:nbSt
        x = bag["Storage"][i]
        reshape(x, (1, :))
        push!(storage, x)
    end
    # FIXME Correct storage loss coefficient by changing its value x to (1 - x)
    for s in 1:nbSt
        storage[s,:lostCoef] = 1 - storage[s,:lostCoef]
    end
    if verbose
        println("\nStorage:")
        showall(storage)
    end
    # find out the number of contracts in each time period
    num_contracts = zeros(Int64, nbT)
    for c in 1:nbC
        t = contract[c,:period]
        num_contracts[t] += 1
    end
    println("\n\nInstance read successfully.\n\n")

    instance_as_dict = Dict()
    instance_as_dict["instance"] = instance
    instance_as_dict["period"] = period
    instance_as_dict["contract"] = contract
    instance_as_dict["drivable"] = drivable
    instance_as_dict["n_drivable"] = n_drivable
    instance_as_dict["n_drivable_uncertain"] = n_drivable_uncertain
    instance_as_dict["storage"] = storage
    instance_as_dict["scenarios"] = bag["Scenarios"]
    instance_as_dict["num_contracts"] = num_contracts
    instance_as_dict["filepath"] = filepath
    return instance_as_dict
end

function read_input_data(datafile, verbose = false)
    s = open(datafile) do file
        read(file, String)
    end
    if contains(s, "!")
        # Xpress format
        return read_xpress_file(datafile)
    elseif contains(s, "//")
        # Antoine tabulated format
        return read_tabulated_data(datafile, verbose)
    end
end

function obtain_instance_parameters(instance)
    return [instance[i] for i in 1:length(instance)]
end

function obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt)
    return 1:nbT, 1:nbS, 1:nbC, 1:nbD, 1:nbND, 1:nbSt
end

function obtain_instance_ranges(nbT, nbS, nbC, nbD, nbND, nbSt, nbNDU)
    return 1:nbT, 1:nbS, 1:nbC, 1:nbD, 1:nbND, 1:nbSt, 1:nbNDU
end

function create_cplex_model(time_limit, cplex_old = false, robust = false, relative_gap_tol = CPLEX_OPT_GAP)
    if cplex_old
        if robust  # for robust model, set relative gap tolerance to get solution in viable time
            return Model(solver=CplexSolver(CPX_PARAM_EPGAP=relative_gap_tol,CPX_PARAM_SCRIND=0))
        end
        return Model(solver=CplexSolver(CPX_PARAM_TILIM=time_limit,CPX_PARAM_SCRIND=0))
    else
        return Model(solver=CplexSolver(CPXPARAM_TimeLimit=time_limit,CPX_PARAM_SCRIND=0))
    end
end

function trunc_if_less_than_eps(value::Float64)
    if abs(value) < EPS
        return 0.0
    end
    return value
end

function create_trace_scenario_filename(test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id)
    return test_name * "_" * instance_name * "_" * sim_strategy * "_" * model_policy * "_ReOpt-" * string(reoptimize) * "_scen" * string(scenario_id)
end

function generate_all_model_parameter_combination(model_policy_params = ["full_model", "ignore_model"],
            sim_strategy_params = ["conservative", "audacious", "cheapest"])
    parameters_to_process = Array[]
    for model_policy in model_policy_params
        if model_policy == "full_model"
            for sim_strategy in sim_strategy_params
                for reoptimize in [true, false]
                    push!(parameters_to_process, (model_policy, sim_strategy, reoptimize))
                end
            end
        else  # 'ignore_model' uses only y variable info, so no reoptimization is needed
            reoptimize = false
            for sim_strategy in sim_strategy_params
                push!(parameters_to_process, (model_policy, sim_strategy, reoptimize))
            end
        end
    end
    return parameters_to_process
end

extension(url::String) = try match(r"\.[A-Za-z0-9]+$", url).match catch "" end

function concatenate_trace_df_list(test_name, instance_name, general_logger, df_filepath = "")
    trace_files_path = create_full_dir(normpath(pwd()), ["output", "simulation", "zip", instance_name])
    if length(df_filepath) > 0
        trace_files_path = normpath(df_filepath, instance_name)
    end
    output_path = create_full_dir(normpath(pwd()), ["output", "simulation", "trace"])
    output_file = joinpath(normpath(output_path), test_name * "_RCCP_Sim_TRACE_" * instance_name)
    println(general_logger, "Concatenating individual trace_df dataframes...")
    trace_df_full, var_df_full = create_empty_trace_dataframe()
    tmp_file_path = "/tmp/rcc_temp_df.csv"
    println(general_logger, "Processing dir : $(trace_files_path)")
    for file in readdir(trace_files_path)
        if extension(file) == ".zip"
            trace_file_zip = joinpath(trace_files_path, file)
            #println("Reading zip file : $(file)...")
            try
                r = ZipFile.Reader(trace_file_zip)
                for f in r.files
                    if extension(f.name) == ".jld2"
                        try
                            println("Reading JLD2 file : $(f.name)...")
                            s = readstring(f)
                            tmp_file_df = open(tmp_file_path, "w+")
                            write(tmp_file_df, s)
                            close(tmp_file_df)
                            @load tmp_file_path trace_df # var_df
                            trace_df_full = vcat(trace_df_full, trace_df)
                            # var_df_full = vcat(var_df_full, var_df)
                        catch y1
                            println("ERROR Reading JLD2 file $(f.name). Cause : $(y1).")
                            println(general_logger, "ERROR Reading JLD2 file $(f.name). Cause : $(y1).")
                        end
                    end
                end
                close(r)
            catch y2
                println("ERROR Reading ZIP file $(file). Cause : $(y2).")
                println(general_logger, "ERROR Reading ZIP file $(file). Cause : $(y2).")
            end
        end
    end
    #rm(tmp_file_path)
    output_file_csv = output_file * ".csv"
    println("\nSaving full trace CSV file to $(output_file_csv)...")
    println(general_logger, "\nSaving full trace CSV file to $(output_file_csv)...")
    CSV.write(output_file_csv, trace_df_full)
    output_file_zip = output_file * ".zip"
    move_files_to_zip_archive(output_file_zip, [output_file_csv])
    println(general_logger, "Concatenation done.")
    println("Concatenation done.")
end

# Move the files in the array 'files_to_compress' to the zip file 'output_file_zip'
function move_files_to_zip_archive(output_file_zip, files_to_compress)
    w = ZipFile.Writer(output_file_zip)
    for filepath in files_to_compress
        filename = basename(filepath)
        f = ZipFile.addfile(w, filename, method=ZipFile.Deflate)
        s = read(filepath, String)
        write(f, s)
    end
    close(w)
    # Delete the files in files_to_compress
    for filepath in files_to_compress
        rm(filepath)
    end
end

function read_file_from_zip_archive(zip_file_path, filename_to_read)
    r = ZipFile.Reader(zip_file_path)
    for f in r.files
        if filename_to_read == f.name
            s = readstring(f)
        end
    end
    close(r)
    return s
end
