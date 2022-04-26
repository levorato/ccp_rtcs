#!/usr/bin/env python
# coding: utf-8

# # Process CCP simulation results - Japan microgrid

run_from_collab = False

import pandas as pd
import numpy as np
import os, fnmatch
import glob
import gzip
from sys import platform
import pyarrow as pa
import pathlib
import zipfile
import tempfile
import shutil


def get_instance_list(project_folder, antoine_instances_folder, toy_instances_folder, japan_instances_folder, instance_group):
    basedir = ""
    instance_list = []
    print("Processing " + instance_group + " instances...")
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
        print('antoine path: ', antoine_instances_folder)
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
        print('Japan path: ', os.path.join(japan_instances_folder, instance_group))
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

def read_concatenated_trace_df(result_path, model_name):
    print("Looking for existing simulations in folder", result_path, "...")
    print("Concatenating individual trace_df dataframes...")
    df_list = []
    tempdir = tempfile.mkdtemp()
    for file in glob.glob(os.path.join(result_path.strip(), "*.zip")):
        if (pathlib.Path(file).suffix == ".zip") and (model_name in file):
            trace_file_zip = os.path.join(result_path.strip(), file)
            print("Reading zip file :", trace_file_zip, "...")
            output_file_trace_arrow = pathlib.Path(trace_file_zip).stem + "_var.arrow"
            print('Arrow file: ', output_file_trace_arrow)
            try:
                r = zipfile.ZipFile(trace_file_zip, 'r')
                for f in r.namelist():
                    if pathlib.Path(f).suffix == ".arrow" and output_file_trace_arrow == f:
                        try:
                            #arrow_file_obj = r.open(f)
                            tmppath = os.path.join(tempdir, f)
                            r.extract(f, tempdir)
                            ##print("Reading arrow file :", f, "...")
                            trace_df = pd.read_feather(tmppath)  #feather.read_feather(tmppath)
                            df_list.append(trace_df)
                            # print(trace_df.head())
                        except Exception as y1:
                            print("ERROR Reading arrow file", f, ". Cause :", y1)
                            #arrow_file_obj.close()
                            return pd.DataFrame()
                        #end
                #end
                r.close()
            except Exception as y2:
                print("ERROR Reading ZIP trace file", file, ". Cause :", y2)
                r.close()
                return pd.DataFrame()
            #end
        #end
    #end
    shutil.rmtree(tempdir)
    print("Concatenation done.")
    if len(df_list) > 0:
        return pd.concat(df_list)
    else:
        return pd.DataFrame()
#end


# Import custom python file from github repo: https://changhsinlee.com/colab-import-python/
if run_from_collab:
    get_ipython().system('pip install requests')
    import requests
    # Save python as file to colab working directory
    # If you are using GitHub, make sure you get the "Raw" version of the code
    url = 'https://raw.githubusercontent.com/levorato/ccp_rtcs/master/notebooks/rccp_utils.py'
    r = requests.get(url)
    # make sure your filename is the same as how you want to import 
    with open('rccp_utils.py', 'w') as f:
        f.write(r.text)
    # now we can import


# ## 1. Process result files

# ### 1.1. Setup project folders
if run_from_collab:
    from google.colab import drive
    drive.mount('/content/gdrive/')
    gdrive_folder = '/content/gdrive/MyDrive'
else:
    gdrive_folder = '/users/mlevorato'
print('gdrive_folder=', gdrive_folder)

project_folder = os.path.join(gdrive_folder, 'doutorado', 'robusto', 'RCCP')
utc_instances_folder = os.path.join(project_folder, "instances", "antoine_skew")
toy_instances_folder = os.path.join(project_folder, "instances", "toy")
instances_folder = os.path.join(project_folder, "instances")
japan_instances_folder = os.path.join(project_folder, "instances", "japan_microgrid")
output_folder = os.path.join(gdrive_folder, "rccp_experiments")
results_folder = os.path.join(gdrive_folder, "rccp_results")
print("*** Project folder is", project_folder)
print("*** Instances folder is",  instances_folder)
print("*** Output folder is", output_folder)

pd.set_option('max_columns', None)

experiment_list = ["run_sim_japan_forecast_avg", "run_sim"]
instance_group_list = ["japan-10", "antoine-skew"]
experiment_folder_list = [os.path.join(output_folder, exp) for exp in experiment_list]
simulated_model_list = ["robust-budget", "robust-box", "robust-budget"]
forecast_type_list = ["average"]  # average-based RTCS forecast

for instance_group, experiment_folder, experiment in zip(instance_group_list, experiment_folder_list, experiment_list):
    print("Processing instance_group = ", instance_group)
    print('Experiment_folder = ', experiment_folder)
    instances_to_process = get_instance_list(project_folder, utc_instances_folder, toy_instances_folder,
                                             japan_instances_folder, instance_group)
    print("instances_to_process: ", instances_to_process)
    robust_model_names = ['robust-budget_{}'.format(_) for _ in [0, 20, 40, 60, 80, 100]]
    det_model_names = ['deterministic_{}'.format(_) for _ in [0, 50, 100]]
    for model_name in ['robust-box'] + robust_model_names + det_model_names:
        print('Procesing model ', model_name)
        for instance_path in instances_to_process:
            instance_name = os.path.basename(instance_path[1])
            print('Processing instance ', instance_name)
            result_path = os.path.join(experiment_folder, "output", "simulation", "zip", instance_name)
            df = read_concatenated_trace_df(result_path, model_name)
            df['InstanceName'] = instance_name
            # ### Create the output folders for processed results
            reportfolder = os.path.join(output_folder, 'consolidated_results')
            reportfolder_df = os.path.join(reportfolder, 'df')
            if not os.path.exists(reportfolder_df):
                os.makedirs(reportfolder_df)
            print('Saving files on folder: ' + reportfolder)

            # ### Save the consolidated dataframe to file
            if len(df.index) > 0:
                outfilepath = os.path.join(reportfolder_df, experiment + '.' + instance_group + '.' + instance_name + '.' + model_name + '.var-results.pkl.gz')
                df.to_pickle(outfilepath)
                print("Saved consolidated pickle file to: ", outfilepath)
            else:
                print('Skipping empty dataframe.')

