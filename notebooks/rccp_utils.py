import glob
import os
import pathlib
import pyarrow as pa
import pyarrow.feather as feather
import zipfile
import pandas as pd
import tempfile
import shutil


def get_instance_list(project_folder, antoine_instances_folder, toy_instances_folder, japan_instances_folder, instance_group):
    basedir = ""
    instance_list = []
    print("Processing, ", instance_group, " instances...")
    if instance_group == "Example_3.1":
        base_folder = os.path.join(project_folder, "instances")
        test_set = ["Example_3.1_A.txt", "Example_3.1_B.txt"]
        for instance_name in test_set:
            inputfile = os.path.join(base_folder, instance_name)
            instance_list.append((instance_name, inputfile))
        #end
    elif instance_group == "antoine_11":
        base_folder = os.path.join(project_folder, "instances")
        #datafile = "../notebooks/data/instance5-contratos_restituicao-ultimo_periodo.txt"
        #datafile = "../notebooks/data/antoine/A_instance2_11scen_1NDU.txt"
        #datafile = "../notebooks/data/antoine/A_instance2_11scen.txt"
        for filename in glob.glob(os.path.join(base_folder, "*.txt")):
            inputfile = os.path.join(base_folder, filename)
            if "A_instance2_11" in filename:
                instance_list.append((filename, inputfile))
            #end
        #end
    elif instance_group == "antoine_oc":
        base_folder = os.path.join(project_folder, "instances")
        #datafile = "../notebooks/data/instance5-contratos_restituicao-ultimo_periodo.txt"
        #datafile = "../notebooks/data/antoine/A_instance2_11scen_1NDU.txt"
        #datafile = "../notebooks/data/antoine/A_instance2_11scen.txt"
        for filename in glob.glob(os.path.join(base_folder, "*.txt")):
            inputfile = os.path.join(base_folder, filename)
            if "A_instance2_OC" in filename:
                instance_list.append((filename, inputfile))
            #end
        #end
    elif instance_group == "antoine-skew":
        for filename in glob.glob(os.path.join(antoine_instances_folder, "*.txt")):
            inputfile = os.path.join(antoine_instances_folder, filename)
            instance_list.append((filename, inputfile))
        #end
    elif instance_group == "toy":  # Taillard instances 20x5
        for filename in glob.glob(os.path.join(toy_instances_folder, "*.txt")):
            inputfile = os.path.join(toy_instances_folder, filename)
            instance_list.append((filename, inputfile))
        #end
    elif instance_group == "japan-10":
        for filename in glob.glob(os.path.join(japan_instances_folder, instance_group, "*.txt")):
            inputfile = os.path.join(japan_instances_folder, instance_group, filename)
            instance_list.append((filename, inputfile))
        #end
    else:
        print("WARN: No instances found!")
    # end
    num_instances = len(instance_list)
    print("# instances to process: ", num_instances)
    # print("- Instances to process: ", instance_list)
    return instance_list
# end function

def create_trace_scenario_filename(model, Gamma_perc, test_name, instance_name, sim_strategy, model_policy, reoptimize, scenario_id):
    if model == "robust-budget":
        return model + string(int(ceil(Gamma_perc))) + "_" + test_name + "_" + instance_name + "_" + sim_strategy + "_" + model_policy + "_ReOpt-" + string(reoptimize) + "_scen" + string(scenario_id)
    else:
        return model + "_" + test_name + "_" + instance_name + "_" * sim_strategy + "_" + model_policy + "_ReOpt-" + string(reoptimize) + "_scen" + string(scenario_id)
    #end
# end function

def read_concatenated_trace_df(result_path):
    print("Looking for existing simulations in folder", result_path, "...")
    print("Concatenating individual trace_df dataframes...")
    df_list = []
    tempdir = tempfile.mkdtemp()
    for file in glob.glob(os.path.join(result_path.strip(), "*.zip")):
        if (pathlib.Path(file).suffix == ".zip"):
            trace_file_zip = file  # os.path.join(result_path.strip(), file)
            ##print("Reading zip file :", trace_file_zip, "...")
            output_file_trace_arrow = pathlib.Path(trace_file_zip).stem + ".arrow"
            try:
                r = zipfile.ZipFile(trace_file_zip, 'r')
                for f in r.namelist():
                    if pathlib.Path(f).suffix == ".arrow" and output_file_trace_arrow == f:
                        try:
                            #arrow_file_obj = r.open(f)
                            tmppath = os.path.join(tempdir, f)
                            r.extract(f, tempdir)
                            ##print("Reading arrow file :", f, "...")
                            trace_df = feather.read_feather(tmppath)
                            df_list.append(trace_df)
                        except Exception as y1:
                            print("ERROR Reading arrow file", f, ". Cause :", y1)
                            #arrow_file_obj.close()
                            continue
                        #end
                #end
                r.close()
            except Exception as y2:
                print("ERROR Reading ZIP trace file", file, ". Cause :", y2)
                r.close()
                continue
            #end
        #end
    #end
    shutil.rmtree(tempdir)
    print("Concatenation done.")
    return pd.concat(df_list)
#end

