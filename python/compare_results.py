#!/usr/bin/env python
# coding: utf-8

# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas
# python-sympy python-nose python-tk
# OR pip install pandas <etc.>

import sys, getopt
import csv
import os
import os.path
import argparse
import pandas as pd
import scipy.stats as stats
import numpy as np
import matplotlib as mpl
from statsmodels.stats.multicomp import MultiComparison
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm
from statsmodels.formula.api import ols
#import researchpy as rp
import fnmatch
import copy

## agg backend is used to create plot as a .png file
mpl.use('agg')

import matplotlib.pyplot as plt


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main(argv):

    csv.field_size_limit(100000)
    pd.set_option('display.max_columns', None)
    mpl.rcParams['figure.figsize'] = (16.0, 12.0)
    mpl.style.use('ggplot')

    parser = argparse.ArgumentParser(description='Compare RCCP result files, generating plots.')
    parser.add_argument("--plots", type=str2bool, nargs='?',
                        const=True, default=True,
                        help="Generate plots.")
    parser.add_argument("--tukey", type=str2bool, nargs='?',
                        const=True, default=True,
                        help="Execute ANOVA / Tukey stat tests.")
    parser.add_argument('--folders', nargs='+',
                        help='the folders containing the result files (one for each experiment)')
    parser.add_argument('--filefilter', default='*_TRACE_*.csv', required=False,
                        help='the file extension for result files (default: *.csv)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    plots = args.plots
    tukey = args.tukey
    args = parser.parse_args()

    print('Input folders are ', folders)
    print('File filter is ', filter)
    processResult(folders, filter, plots, tukey)


def processResult(folders, filter, gen_plots, gen_tukey):
    # for each algorithm, determines the best solution value
    for folder in folders:
        print("Processing folder " + ''.join(folder))

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            files.sort()
            for filename in fnmatch.filter(files, filter):
                processTraceFile(folder, filename, gen_plots, gen_tukey)
            # end for
        # end for
    # end process all result files of an instance / filename


def aggregate_simulation_samples_by_scenario(df):
    filter_df = copy.deepcopy(df)
    filter_df = filter_df.filter(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy",
            "ScenarioId", "DetValue", "RobValue", "gap", "cost", "e_td"])
    # Group by scenario_id only
    group_df = filter_df.groupby(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy", "ScenarioId"])\
        .agg({'gap' : ['sum'], 'cost' : ['sum'], 'e_td' : ['sum'], 'DetValue' : ['mean'], 'RobValue' : ['mean']})
    group_df.columns = ["_".join(x) for x in group_df.columns.ravel()]
    group_df_each_scenario_and_period = group_df.reset_index()
    #print_df_head(group_df_each_scenario_and_period, 'Group by scenario_id only (sum_by_scenario_df)')
    return group_df_each_scenario_and_period.copy()


def aggregate_simulation_samples_by_scenario_and_period(df):
    filter_df = copy.deepcopy(df)
    filter_df = filter_df.filter(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy",
            "ScenarioId", "t", "DetValue", "RobValue", "gap", "cost", "e_td"])
    # Group by scenario_id and t
    group_df = filter_df.groupby(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy", "ScenarioId", "t",
                                  "DetValue", "RobValue"]).agg(['sum'])
    group_df.columns = ["_".join(x) for x in group_df.columns.ravel()]
    group_df_each_scenario_and_period = group_df.reset_index()
    #print_df_head(group_df_each_scenario_and_period, 'Group by scenario_id and t (samples_df)')
    return group_df_each_scenario_and_period.copy()


# For each period t, calculates the average cost considering all scenarios
def aggregate_simulation_samples_by_period(df):
    filter_df = copy.deepcopy(df)
    filter_df = filter_df.filter(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy",
                          "ScenarioId", "t", "gap", "cost", "e_td"])
    # 1. For each period t and scenario_id, calculate the sum of costs over all microperiods d
    group_df = filter_df.groupby(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy",
                                  "ScenarioId", "t"]).agg(['sum'])
    group_df.columns = ["_".join(x) for x in group_df.columns.ravel()]
    group_df = group_df.reset_index()


    # 2. For each period t, calculate the average cost given all scenarios
    group_df = group_df.filter(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy",
                                 "t", "gap_sum", "cost_sum", "e_td_sum"])
    group_df = group_df.groupby(["RTCS_Type", "Strategy", "Reoptimize", "ModelPolicy",
                                 "t"]).agg(['mean'])
    group_df.columns = ["_".join(x) for x in group_df.columns.ravel()]
    group_df_each_period = group_df.reset_index()
    #print_df_head(group_df_each_period, 'Group only by t (mean_sum_df)')
    return group_df_each_period.copy()


def print_df_head(df, name):
    print('\n=========================================\n' + name + '\n=========================================\n')
    print("Grouped dataframe deterministic:\n" + str(df[df['RTCS_Type'] == "Deterministic"].head()))
    print("Grouped dataframe robust:\n" + str(df[df['RTCS_Type'] == "Robust"].head()))

def print_full_df(df, name):
    print('\n=========================================\n' + name + '\n=========================================\n')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df)


# TODO Reescrever utilizando o valor da FO do Robusto obtido a partir do arquivo de log de otimizacao :
# Antoine_RCCP_Sim_OptData_A_instance2_11scen.txt_.csv
def find_problematic_simulation_values(df):
    #problematic_df_det = df[df['RTCS_Type'] == "Deterministic"]
    #problematic_df_det = problematic_df_det[problematic_df_det['DetValue_mean'] > problematic_df_det['cost_sum']]
    #print_df_head(problematic_df_det, 'Problematic Deterministic simulation cost values')

    problematic_df_rob = df[df['RTCS_Type'] == "Robust"]
    problematic_df_rob = problematic_df_rob[problematic_df_rob['RobValue_mean'] < problematic_df_rob['cost_sum']]
    problematic_df_rob['Diff_Cost_Values'] = np.subtract(problematic_df_rob['cost_sum'], problematic_df_rob['RobValue_mean'])
    #print_full_df(problematic_df_rob, 'Problematic Robust simulation cost values')
    return problematic_df_rob


def processTraceFile(folder, filename, gen_plots, gen_tukey):
    print("\nProcessing trace file : " + filename + "...\n")

    base_df = pd.read_csv(os.path.join(folder, filename), quotechar="\"")  # header=0,
#                     dtype = {'RTCS_Type': object, 'Strategy': object, 'Reoptimize': object, 'ModelPolicy' : object,
#                              'ScenarioId' : np.int64, 't' : np.int64, 'd' : np.int64, 'DetValue' : np.float64,
#                              'RobValue' : np.float64, 'OptTimeSpent' : np.float64, 'q_td' : object,
#                              'x_td' : object, 'e_td' : object, 'r_td' : object, 'g_td' : object,
#                              'h_td' : object, 'gap' : np.float64, 'cost' : np.float64 })
    #print("Dataframe data types:\n" + str(base_df.describe()))
    #print("Read to dataframe:\n" + str(base_df.head()))

    output_path = os.path.join(os.getcwd(), "output", filename)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Sum over all microperiods d, group by scenario_id and t
    samples_df = aggregate_simulation_samples_by_scenario_and_period(base_df)
    # Sum over all microperiods d and time periods t, group by scenario
    sum_by_scenario_df = aggregate_simulation_samples_by_scenario(base_df)
    # For each strategy, totalize the number of times the Det RTCS won over the Rob RTCS, given all scenarios
    totalize_n_samples_det_vs_rob_wins(sum_by_scenario_df, output_path)
    totalize_n_samples_which_scenarios_each_strategy_wins(sum_by_scenario_df, output_path)

    mean_sum_df = aggregate_simulation_samples_by_period(base_df)
    probl_df = find_problematic_simulation_values(sum_by_scenario_df)
    probl_df.to_csv(os.path.join(output_path, 'Problematic_Robust_Sim_Values.csv'))

    if gen_plots:
        generate_plots(mean_sum_df, samples_df, sum_by_scenario_df, output_path)
        generate_comparison_plot_reopt_yes_no(mean_sum_df, output_path)
        generate_comparison_plot_each_strategy(mean_sum_df, output_path)
    # end if gen_plots

    if gen_tukey:
        full_output_path = os.path.join(output_path, "acc_cost")
        if not os.path.exists(full_output_path):
            os.makedirs(full_output_path)

        #run_paired_statistical_tests(series, output_path)
        for reopt in mean_sum_df['Reoptimize'].unique():
            print("==================================\nStatistical tests for ReOpt=", str(reopt), "\n==================================\n")
            acc_cost_samples_dict, acc_cost_samples_df = generate_acc_cost_samples_as_dict_and_df(
                mean_sum_df, sum_by_scenario_df, reopt)
            run_n_samples_stats_tests(acc_cost_samples_dict, acc_cost_samples_df, full_output_path, str(reopt))

            # Generate box plots for the accumulated cost measure
            boxplot = acc_cost_samples_df.boxplot(column='cost_sum', by='combination', rot=20)

            plt.savefig(os.path.join(full_output_path, 'BoxPlot_acc_cost_ReOpt_' + str(reopt) + '.png'))
            plt.savefig(os.path.join(full_output_path, 'BoxPlot_acc_cost_ReOpt_' + str(reopt) + '.pdf'))
        # end for reopt
    # end if gen_tukey


def generate_acc_cost_samples_as_dict_and_df(df, sum_by_scenario_df, reopt):
    acc_cost_samples_dict = dict()
    y = 'cost_sum'
    which_df = copy.deepcopy(sum_by_scenario_df)
    # For each strategy, policy and reoptimization flag
    for strategy in df['Strategy'].unique():
        for model_policy in df['ModelPolicy'].unique():
            #for reopt in df['Reoptimize'].unique():
                combination_name = strategy + " x " + model_policy + " x ReOpt=" + str(reopt)
                boxplot_df2 = which_df[(which_df['Strategy'] == strategy) &
                                       (which_df['ModelPolicy'] == model_policy) & (
                                               which_df['Reoptimize'] == reopt)]
                boxplot_df2 = boxplot_df2.filter(["RTCS_Type", y])
                #boxplot_df2[y] = pd.to_numeric(boxplot_df2[y], errors='coerce')
                # If resulting dataframe is empty, ignore this combination (does not exist)
                if boxplot_df2.empty:
                    continue
                for rtcs_type in df['RTCS_Type'].unique():
                    # Store the series in a dictionary for future ANOVA and Tukey tests
                    acc_cost_samples_dict[combination_name + " x " + rtcs_type + "RTCS"] = \
                        np.array(boxplot_df2[boxplot_df2['RTCS_Type'] == rtcs_type][y].values, dtype=float)
                # end for rtcs_type
            # end for
        # end for
    # end for
    which_df = which_df[which_df['Reoptimize'] == reopt]
    which_df["combination"] = which_df['Strategy'] + '_' + which_df['ModelPolicy'] + '_ReOpt=' + which_df['Reoptimize'].map(str) + '_' + which_df['RTCS_Type']
    sub_df = which_df.filter(["ScenarioId", "combination", "cost_sum"])
    sub_df_table = pd.pivot_table(sub_df, values=['cost_sum'], index=['ScenarioId'], columns=['combination'], aggfunc=np.mean)
    # https://stackoverflow.com/questions/38951345/how-to-get-rid-of-multilevel-index-after-using-pivot-table-pandas
    sub_df_table.index.name = None
    #print("========================\nacc_cost_samples_df\n======================\n", sub_df)
    acc_cost_samples_df = sub_df

    return acc_cost_samples_dict, acc_cost_samples_df  # , which_df["combination"].unique().tolist()


def generate_plots(df, samples_df, sum_by_scenario_df, output_path):
    print("Starting plots generation...\n")
    # Determ vs Robust cost, gap and OC plots
    acc_series = dict()
    acc_y_labels = dict()
    for y in ['cost_sum', 'gap_sum', 'e_td_sum']:
        acc_series[y] = dict()
        acc_y_labels[y] = []
    # end for y
    # For each strategy, policy and reoptimization flag
    for strategy in df['Strategy'].unique():
        for model_policy in df['ModelPolicy'].unique():
            for reopt in df['Reoptimize'].unique():
                combination_name = strategy + " x " + model_policy + " x ReOpt=" + str(reopt)
                print("\nProcessing plots for ", combination_name)
                for y in ['cost_sum_mean', 'gap_sum_mean', 'e_td_sum_mean']:
                    print("    * Bar and line plots for measure: ", y)
                    plot_df = df[ (df['Strategy'] == strategy) &
                                  (df['ModelPolicy'] == model_policy) & (df['Reoptimize'] == reopt) ]
                    plot_df = plot_df.filter(["t", "RTCS_Type", y])
                    # If resulting dataframe is empty, ignore this combination (does not exist)
                    if plot_df.empty:
                        continue
                    p_table = pd.pivot_table(plot_df, values=y, index=['t'], columns=['RTCS_Type'], aggfunc=np.sum)

                    title = strategy + "_" + model_policy + "\n_ReOpt" + str(reopt)

                    #print("\n PTable for " + title + ": \n\n" + str(p_table))
                    p_table.plot(kind='bar', title=title)
                    #plt.show()
                    full_output_path = os.path.join(output_path, y)
                    if not os.path.exists(full_output_path):
                        os.makedirs(full_output_path)
                    plt.savefig(os.path.join(full_output_path, 'BarPlot_' + title + ".png"))
                    plt.savefig(os.path.join(full_output_path, 'BarPlot_' + title + ".pdf"))

                    p_table.plot(kind='line', title=title)
                    plt.savefig(os.path.join(full_output_path, 'LinePlot_' + title + ".png"))
                    plt.savefig(os.path.join(full_output_path, 'LinePlot_' + title + ".pdf"))
                # end for y
                # ===================================== BOX PLOT GENERATION
                for y in ['cost_sum', 'gap_sum']:
                    print("    * Box plots for measure: ", y)
                    # ======================== Generate box plots for accumulated measures ( sum over all periods (t, d) )
                    which_df = copy.deepcopy(sum_by_scenario_df)
                    boxplot_df2 = which_df[(which_df['Strategy'] == strategy) &
                                           (which_df['ModelPolicy'] == model_policy) & (
                                                   which_df['Reoptimize'] == reopt)]
                    #print("*** " + str(boxplot_df2))
                    boxplot_df2 = boxplot_df2.filter(["RTCS_Type", y])
                    # If resulting dataframe is empty, ignore this combination (does not exist)
                    if boxplot_df2.empty:
                        continue

                    acc_series[y][title] = dict()
                    acc_y_labels[y].append(title)
                    labels = df['RTCS_Type'].unique()
                    for rtcs_type in labels:
                        acc_series[y][title][rtcs_type] = boxplot_df2[boxplot_df2['RTCS_Type'] == rtcs_type][y].values
                        print('(acc) y = ' + str(y) + ', RTCS = ' + str(rtcs_type) + ' : ' + str(acc_series[y][title][rtcs_type]))
                    # end for rtcs_type

                    boxplot_df = samples_df[(samples_df['Strategy'] == strategy) &
                                 (samples_df['ModelPolicy'] == model_policy) & (samples_df['Reoptimize'] == reopt)]
                    boxplot_df = boxplot_df.filter(["RTCS_Type", "t", y])
                    #p_table = pd.pivot_table(plot_df, values=y, index=['t'], columns=['RTCS_Type'], aggfunc=np.sum)
                    #boxplot_df.reset_index()
                    #boxplot_df.set_index(['t', 'RTCS_Type'])
                    #boxplot_df.boxplot(by=['t'])
                    series = dict()
                    ylabels = df['t'].unique()
                    labels = df['RTCS_Type'].unique()
                    for t in ylabels:
                        series[t] = dict()
                        for rtcs_type in labels:
                            series[t][rtcs_type] = boxplot_df[(boxplot_df['RTCS_Type'] == rtcs_type)
                                                              & (boxplot_df['t'] == t)][y].values
                            #print('t = ' + str(t) + ', RTCS = ' + str(rtcs_type) + ' : ' + str(series[t][rtcs_type]))
                        # end for rtcs_type
                    # end for t
                    full_output_path = os.path.join(output_path, y)
                    if not os.path.exists(full_output_path):
                        os.makedirs(full_output_path)
                    result_filepath = os.path.join(full_output_path, 'BoxPlot_EachPeriod_' + title + ".png")
                    #generate_box_plot_horizontal(series, labels, ylabels, result_filepath)
                    # ======================== Generate box plots for each period t
                    generate_box_plot(series, labels, ylabels, result_filepath, y)

                # end for y
            # end for reopt
        # end for model_policy
    # end for strategy
    #pp = PdfPages('foo.pdf')
    #pp.savefig(plot1)
    #pp.savefig(plot2)
    #pp.savefig(plot3)
    #pp.close()
    # ======================== Generate box plots for accumulated measures ( sum over all periods (t, d) )
    for y in ['cost_sum', 'gap_sum']:
        result_filepath = os.path.join(output_path, "BoxPlot_AccAllStrategies_" + str(y) + ".png")
        generate_box_plot_horizontal(acc_series[y], labels, acc_y_labels[y], result_filepath, y)
    # end for y


def run_paired_statistical_tests(series, output_path):
    Anova = []
    Tukey = []
    labels = []

    # 1. Mann - Whitney - Wilcoxon(MWW) Rank Sum test  ================================
    # The MWW RankSum test is a useful test to determine if two distributions are
    # significantly different or not. Unlike the t - test, the RankSum test does
    # not assume that the data are normally distributed, potentially providing a
    # more accurate assessment of the data sets.
    # =================================================================================
    # A RankSum test will provide a P value indicating whether or not the two distributions are the same.
    # To run Wilcoxon signed rank test, sample sizes must be equal!
    if len(series[0]) != len(series[1]):
        raise Exception('Sample sizes must be equal!')
    z_stat, p_val = stats.ranksums(series[0], series[1])
    print("MWW RankSum P for series 1 and 2 =", p_val)
    # If P <= 0.05, we are highly confident that the distributions significantly differ, and can claim that
    # the samples had a significant impact on the measured value. If the samples do not significantly differ,
    # we could expect a result such as P == 0.99.
    # With P > 0.05, we must say that the distributions do not significantly differ.

    # 2. Independent T-test  ==========================================================
    # https://pythonfordatascience.org/independent-t-test-python/
    # TODO Create python notebook with statistical analysis of data
    # Like every test, this inferential statistic test has assumptions.
    # The assumptions that the data must meet in order for the test results to be valid are
    #  a. The samples are independently and randomly drawn
    #  b. The distribution of the residuals between the two groups should follow the normal distribution
    #  c. The variances between the two groups are equal
    # If any of these assumptions are violated then another test should be used. The dependent variable
    # (outcome being measured) should be continuous which is measured on an interval or ratio scale.
    stats.ttest_ind(series[0], series[1])


    # Wilcoxon.append(stats.mannwhitneyu(series[0], series[1]))
    # Applies ANOVA to verify that the samples differ
    f, p = stats.f_oneway(series[0], series[1])
    if p > 0.05:
        print("WARN: p-value above 0.05! Samples may not differ!")


    # runs tukey's test to compare all the samples
    for item in series[0]:
        Tukey.append((labels[0], item))
    for item in series[1]:
        Tukey.append((labels[1], item))

    names = ['id', 'ttt']
    formats = ['U30', 'f8']
    dtype = dict(names=names, formats=formats)
    # array = np.rec.array(ttt_tukey, dtype=dtype)

    # res = pairwise_tukeyhsd(array['ttt'], array['id'], alpha=0.05)
    # res.plot_simultaneous()
    # plt.show()
    # print res
    # print res.meandiffs[0]
    # df.to_csv(os.path.join(output_path, 'paired_stats_tests.csv'))


# http://web.archive.org/web/20161011044754/http://statsmodels.sourceforge.net:80/devel/examples/generated/example_interactions.html
def run_n_samples_stats_tests(n_samples_dict, n_samples_df, output_path, out_file_suffix):
    print("\nStarting ANOVA / Tukey tests...")
    file = open(os.path.join(output_path, "stats_tests_ReOpt_" + str(out_file_suffix) + ".txt"), "w+")
    # 1. One-way analysis of variance (ANOVA)  ==========================================================
    # If you need to compare more than two data sets at a time, an ANOVA is your best bet.
    # Another way of putting this is that we can use ANOVA for testing hypothesis. A common approach is to assume
    # that the data sets are samples of the same distribution (i.e. the null hypothesis is that their means are equal).
    # Rejecting the null hypothesis would imply that at least one of the means is different.
    # In other words our null hypothesis is that the means of all populations are equal.
    # If p <= 0.05, we reject the null hypothesis and we conclude that at least one of the means is different from at
    # least one other population mean.
    #f, p = stats.f_oneway(n_samples_dict['conservative x full_model x ReOpt=False x RobustRTCS'],
    #                      n_samples_dict['conservative x ignore_model x ReOpt=False x RobustRTCS'])
    f, p = stats.f_oneway(*[value for key, value in n_samples_dict.items()])
    # value for key, value in n_samples_dict.items())
    file.write("\n==========================================\n")
    file.write('One-way ANOVA')
    file.write("\n==========================================\n")
    file.write('F value: ' + str(f) + '\n')
    file.write('P value: ' + str(p) + '\n\n')
    # If P > 0.05, we can claim with high confidence that the means of the results of all n experiments
    # are not significantly different.
    # *** The thing with one-way ANOVA is that although we now know that there is difference in the performance of
    # each group, we do not know know exactly who performs best or worst. This is why the analysis of variance is
    # often followed by a post hoc analysis.

    # https://pythonfordatascience.org/anova-python/
    # Getting summary statistics
    #table1 = rp.summary_cont(n_samples_df['cost_sum'])
    #print(table1)
    #table2 = rp.summary_cont(n_samples_df['cost_sum'].groupby(n_samples_df['combination']))
    #print(table2)
    # ANOVA with statsmodels
    results = ols('cost_sum ~ C(combination)', data=n_samples_df).fit()
    #results = ols('libido ~ C(dose)', data=n_samples_df).fit()
    file.write("\n==========================================\n")
    file.write("ANOVA")
    file.write("\n==========================================\n")
    file.write(str(results.summary()) + '\n\n')
    aov_table = sm.stats.anova_lm(results, typ=2)
    file.write(str(aov_table) + '\n')
    file.write("\n==========================================\n")
    file.write("Shapiro test for normality (residuals)")
    file.write("\n==========================================\n")
    file.write(str(stats.shapiro(results.resid)) + '\n\n')

    # 2. Tukey's range test   ===========================================================================
    # Tukey's range test, named after the American mathematician John Tukey, is a common method used as post hoc
    # analysis after one-way ANOVA. This test compares all possible pairs and we can use it to precisely identify
    # difference between two means that's greater than the expected standard error.
    # http://cleverowl.uk/2015/07/01/using-one-way-anova-and-tukeys-test-to-compare-data-sets/
    # https://pt.coursera.org/lecture/data-analysis-tools/python-lesson-9-post-hoc-tests-for-anova-uMNjz
    mc = MultiComparison(n_samples_df['cost_sum'], n_samples_df['combination'])
    result = mc.tukeyhsd()
    file.write("\n==========================================\n")
    file.write("Tukey HSD post-hoc comparison")
    file.write("\n==========================================\n")
    file.write(str(result) + '\n')
    file.write(str(mc.groupsunique) + '\n\n')
    file.close()
    print("Stat tests done.\n")


def det_wins(det, rob):
    if det < rob:
        return 1
    return 0


def rob_wins(det, rob):
    if rob < det:
        return 1
    return 0


def totalize_n_samples_det_vs_rob_wins(df, output_path):
    print("\nTotalizing Det vs Rob wins given all scenarios...")

    table = pd.pivot_table(df, values='cost_sum', index=['Strategy', 'Reoptimize', 'ModelPolicy', 'ScenarioId'],
                           columns=['RTCS_Type'], aggfunc=np.sum)

    table['det_wins'] = np.vectorize(det_wins)(table['Deterministic'], table['Robust'])
    table['rob_wins'] = np.vectorize(rob_wins)(table['Deterministic'], table['Robust'])
    filename = os.path.join(output_path, "Det_vs_Rob_Wins_Per_Scenario.csv")
    table.to_csv(filename)

    df2 = table.groupby(level=[0, 1, 2])['det_wins', 'rob_wins'].sum()
    filename = os.path.join(output_path, "Det_vs_Rob_Wins_Summary.csv")
    df2.to_csv(filename)

    print("Det vs Rob totalization done.\n")


def totalize_n_samples_which_scenarios_each_strategy_wins(df, output_path):
    print("\nTotalizing how much scenarios each strategy wins all the others...")

    table = pd.pivot_table(df, values='cost_sum', index=['ScenarioId'],
                           columns=['Strategy', 'Reoptimize', 'ModelPolicy', 'RTCS_Type'], aggfunc=np.sum)
    filename = os.path.join(output_path, "All_Strategies_Costs_Per_Scenario.csv")
    table.to_csv(filename)
    table_0_1 = table.apply(lambda x: x == x.min(), axis=1).astype(int)

    #table['det_wins'] = np.vectorize(det_wins)(table['Deterministic'], table['Robust'])
    #table['rob_wins'] = np.vectorize(rob_wins)(table['Deterministic'], table['Robust'])
    filename = os.path.join(output_path, "All_Strategies_Wins_Per_Scenario.csv")
    table_0_1.to_csv(filename)

    df2 = table.groupby(level=[0, 1, 2])['det_wins', 'rob_wins'].sum()
    filename = os.path.join(output_path, "All_Strategies_Wins_Summary.csv")
    df2.to_csv(filename)

    print("Det vs Rob totalization done.\n")


def generate_box_plot_horizontal(series, labels, ylabels, result_filepath, y_axis_name):
    # group boxplots - http://stackoverflow.com/questions/20365122/how-to-make-a-grouped-boxplot-graph-in-matplotlib
    height = 10 * len(ylabels)
    padding = 60
    bbox_to_anchor = (0.98, 0.965)
    fig, axes = plt.subplots(nrows=len(ylabels), sharex=True, figsize=(12, height))
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)

    # draw temporary colored lines and use them to create a legend
    hB, = plt.plot([1, 1], '#800080')
    hR, = plt.plot([1, 1], '#DAA520')
    fig.legend((hB, hR), ('Det', 'Rob'), loc='upper right', shadow=True,
               bbox_to_anchor=bbox_to_anchor, ncol=2, borderaxespad=0.)
    # bbox_to_anchor=(0., 1.02, 1., .102))  # bbox_to_anchor=(0, 1))
    hB.set_visible(False)
    hR.set_visible(False)
    # Create an axes instance
    #ax = fig.add_subplot(111)

    print(ylabels)
    axis_count = -1
    for ax, t in zip(axes, ylabels):
        #print("Processing t = " + str(t))
        axis_count += 1

        bp = ax.boxplot([ series[t][item] for item in labels ], patch_artist=True, vert=False) #, showfliers=False)
        #ax.set(yticklabels=labels) #, ylabel=name)
        ax.set(yticklabels=['Det', 'Rob'])  # , ylabel=name)
        ax.set_ylabel(str(t), rotation = 0, labelpad = padding)
        #if axis_count == 0:
        #    ax.set_xlabel("Time period (t)", labelpad=20)
        #    ax.xaxis.set_label_position('top')
        ax.tick_params(axis='y', which='major', pad=5)

        #ax.margins(0.05)  # Optional

        ## change outline color, fill color and linewidth of the boxes
        count = 0
        for box in bp['boxes']:
            # change outline color
            if count % 2 == 0:
                box.set(color='#DAA520', linewidth=2)
            else:
                box.set(color='#7570b3', linewidth=2)
            # change fill color
            if count % 2 == 0:
                #box.set(facecolor='#1b9e77')
                box.set(facecolor='#FFFF66')
            else:
                box.set(facecolor='#800080')
            count += 1

        ## change color and linewidth of the whiskers
        for whisker in bp['whiskers']:
            whisker.set(color='#7570b3', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='#7570b3', linewidth=2)

        ## change color and linewidth of the medians
        count = 0
        for median in bp['medians']:
            if count % 2 == 0:
                median.set(color='#654321', linewidth=2)
            else:
                median.set(color='#b2df8a', linewidth=2)
            count += 1

        ## change the style of fliers and their fill
        for flier in bp['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.5)

        ## Remove top axes and right axes ticks
        #ax.get_xaxis().tick_bottom()
        if axis_count == 0:
            ax.get_xaxis().tick_top()
        #if axis_count == 1:


        ax.get_yaxis().tick_left()

    ## Custom y-axis labels
    #ax.set_yticklabels(instance_names)

    # Save the figure
    plt.tight_layout()
    plt.subplots_adjust(hspace=.01)

    #fig.show()
    print("           Saving box plot png file to " + result_filepath + '-hor_box_plot.png')
    fig.savefig(result_filepath + '-hor_box_plot.png') #, bbox_inches='tight')
    pp = PdfPages(result_filepath + '-hor_box_plot.pdf')
    pp.savefig(plt.gcf())
    pp.close()


# Vertical box plots
def generate_box_plot(series, labels, ylabels, result_filepath, y_axis_name):
    # group boxplots - http://stackoverflow.com/questions/20365122/how-to-make-a-grouped-boxplot-graph-in-matplotlib
    width = 20
    height = 6
    padding = 5
    bbox_to_anchor = (0.98, 0.965)
    rotation = 90
    #if instance_names[0].rfind('file_') >= 0:
    #    width = 9
    #    height = 6 #20
    #    padding = 5
    #    bbox_to_anchor = (0.98, 0.955)
    #    rotation = 60
    fig, axes = plt.subplots(ncols=len(ylabels), sharey=True, figsize=(width, height)) #ncols=len(labels)) #, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)

    # Create a figure instance
    #fig = plt.figure(1, figsize=(30, 20))

    # draw temporary colored lines and use them to create a legend
    hB, = plt.plot([1, 1], '#800080')
    hR, = plt.plot([1, 1], '#DAA520')
    fig.legend((hR, hB), ('Det', 'Rob'), loc='upper right', shadow=True,
               bbox_to_anchor=bbox_to_anchor, ncol=2, borderaxespad=0.)
    # bbox_to_anchor=(0., 1.02, 1., .102))  # bbox_to_anchor=(0, 1))
    hB.set_visible(False)
    hR.set_visible(False)
    # Create an axes instance
    #ax = fig.add_subplot(111)

    axis_count = -1
    for ax, t in zip(axes, ylabels):
        #print("Processing instance t = " + str(t))
        axis_count += 1

        bp = ax.boxplot([ series[t][item] for item in labels ], patch_artist=True, vert=True) #, showfliers=False)
        #ax.set(yticklabels=labels) #, ylabel=name)
        ax.set(xticklabels=['(D)', '(R)'])  # , ylabel=name)
        # remove the .g file extension from instance name
        name = ylabels[axis_count]
        ax.set_xlabel(str(t), rotation = rotation, labelpad = padding)
        if axis_count == 0:
            ax.set_ylabel(y_axis_name, labelpad=5)
            #ax.yaxis.set_label_position('top')
        ax.tick_params(axis='y', which='major', pad=1)

        #ax.margins(0.05)  # Optional

        ## change outline color, fill color and linewidth of the boxes
        count = 0
        for box in bp['boxes']:
            # change outline color
            if count % 2 == 0:
                box.set(color='#DAA520', linewidth=2)
            else:
                box.set(color='#7570b3', linewidth=2)
            # change fill color
            if count % 2 == 0:
                #box.set(facecolor='#1b9e77')
                box.set(facecolor='#FFFF66')
            else:
                box.set(facecolor='#800080')
            count += 1

        ## change color and linewidth of the whiskers
        for whisker in bp['whiskers']:
            whisker.set(color='#7570b3', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='#7570b3', linewidth=2)

        ## change color and linewidth of the medians
        count = 0
        for median in bp['medians']:
            if count % 2 == 0:
                median.set(color='#654321', linewidth=2)
            else:
                median.set(color='#b2df8a', linewidth=2)
            count += 1

        ## change the style of fliers and their fill
        for flier in bp['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.5)

        ## Remove top axes and right axes ticks
        #ax.get_xaxis().tick_bottom()
        #if axis_count == 0:
        #    ax.get_xaxis().tick_top()
        #if axis_count == 1:

        ax.get_xaxis().tick_bottom()

    ## Custom y-axis labels
    #ax.set_yticklabels(instance_names)

    # Save the figure
    plt.tight_layout()
    plt.subplots_adjust(hspace=.01)

    #fig.show()
    print("           Saving box plot png file to " + result_filepath + '-ver_box_plot.png')
    fig.savefig(result_filepath + '-ver_box_plot.png') #, bbox_inches='tight')
    pp = PdfPages(result_filepath + '-ver_box_plot.pdf')
    pp.savefig(plt.gcf())
    pp.close()


def generate_comparison_plot_reopt_yes_no(df, result_filepath):
    figs = []
    if not df.empty:
        print('           Generating Deterministic vs Robust comparison plots...')
        for y in ['cost_sum_mean', 'gap_sum_mean']:
            titles = []
            df_list = []
            nb_combinations = 0
            for strategy in df['Strategy'].unique():  # For each strategy, policy and reoptimization flag
                for model_policy in df['ModelPolicy'].unique():
                    for reopt in df['Reoptimize'].unique():
                        plot_df = df[(df['Strategy'] == strategy) &
                                     (df['ModelPolicy'] == model_policy) & (df['Reoptimize'] == reopt)]
                        plot_df = plot_df.filter(["t", "RTCS_Type", y])
                        df_list.append(plot_df)
                        titles.append(strategy + "_" + model_policy + "_ReOpt=" + str(reopt))
                        nb_combinations += 1
                    # end for reopt
                # end for model_policy
            # end for strategy
            height = 6 * int(nb_combinations / 2)
            padding = 40
            bbox_to_anchor = (0.98, 0.965)
            fig, axes = plt.subplots(nrows=int(nb_combinations / 2), ncols=2, sharex=False, sharey=False, figsize=(20, height))
            for ax, i in zip(axes, range(0, int(nb_combinations / 2))):
                # left column : reopt = no
                if not df_list[2 * i].empty:
                    p_table = pd.pivot_table(df_list[2 * i], values=y, index=['t'], columns=['RTCS_Type'], aggfunc=np.sum)
                    p_table.plot(ax = ax[0], kind='line', title=titles[2 * i])
                # right column : reopt = yes
                if not df_list[2 * i + 1].empty:
                    p_table = pd.pivot_table(df_list[2 * i + 1], values=y, index=['t'], columns=['RTCS_Type'], aggfunc=np.sum)
                    p_table.plot(ax = ax[1], kind='line', title=titles[2 * i + 1])
            # end for ax, i
            # ajuste da margem entre os subplots
            plt.subplots_adjust(hspace=0.6)
            fig.tight_layout()
            #fig.show()
            figs.append(fig)
            # Export plots do PDF file
            full_output_path = os.path.join(result_filepath, y)
            if not os.path.exists(full_output_path):
                os.makedirs(full_output_path)
            with PdfPages(os.path.join(full_output_path, 'Plot_' + str(y) + '_ReOpt_Yes_vs_No.pdf')) as pdf:
                for fig in figs:
                    plt.figure(fig.number)
                    pdf.savefig()
                    # pdf.savefig(plt.gcf())
                plt.close()
        # end for y
        print('           Done.')
    # end if


# TODO trocar os graficos desta funcao para um grafico que acumule os valores anteriores
def generate_comparison_plot_each_strategy(df, result_filepath):
    figs = []
    if not df.empty:
        print('           Generating Deterministic/Robust plots comparing all strategies...')
        print(str(df.head()))
        for y in ['cost_sum_mean', 'gap_sum_mean']:
            full_output_path = os.path.join(result_filepath, y)
            if not os.path.exists(full_output_path):
                os.makedirs(full_output_path)
            titles = []
            df_list = []
            nb_combinations = 0
            # Line Plots with one line for each strategy, policy and reoptimization flag
            plot_df = df
            source_col_loc = plot_df.columns.get_loc('RTCS_Type')  # column position starts from 0
            plot_df['Params'] = plot_df.iloc[:, source_col_loc :source_col_loc + 4].apply(
                lambda x: "_".join(x.astype(str)), axis=1)

            agg_dfs = []
            for param in plot_df['Params'].unique():
                sub_df = plot_df[(plot_df['Params'] == param)].copy()
                new = [sub_df[y].values[0]]
                for i in range(1, len(sub_df.index)):
                    new.append(new[i - 1] + sub_df[y].values[i])
                acc_value_col_name = 'acc_' + y
                sub_df[acc_value_col_name] = new
                agg_dfs.append(sub_df)
            # end for
            agg_df = pd.concat(agg_dfs)
            agg_df.to_csv(os.path.join(full_output_path, 'DF_' + acc_value_col_name + '_CompareStrategy_Det_Rob.csv'))
            plot_df = agg_df
            for rtcs in plot_df['RTCS_Type'].unique():
                plot_df2 = plot_df[(plot_df['RTCS_Type'] == rtcs)]
                plot_df2 = plot_df2.filter(["t", "Params", acc_value_col_name])
                df_list.append(plot_df2)
                titles.append(rtcs)
                nb_combinations += 1
            for strategy in plot_df['Strategy'].unique():
                plot_df2 = plot_df[(plot_df['Strategy'] == strategy)]
                plot_df2 = plot_df2.filter(["t", "Params", acc_value_col_name])
                df_list.append(plot_df2)
                titles.append(strategy)
                nb_combinations += 1
            height = 9 * int(nb_combinations / 2)
            padding = 40
            bbox_to_anchor = (0.98, 0.965)
            fig, axes = plt.subplots(nrows=int(nb_combinations / 2), ncols=2, sharex=False, sharey=False, figsize=(20, height))
            for ax, i in zip(axes, range(0, int(nb_combinations / 2))):
                # left column
                p_table = pd.pivot_table(df_list[2 * i], values=acc_value_col_name, index=['t'], columns=['Params'], aggfunc=np.sum)
                p_table.plot(ax = ax[0], kind='line', title=titles[2 * i])
                # right column
                p_table = pd.pivot_table(df_list[2 * i + 1], values=acc_value_col_name, index=['t'], columns=['Params'], aggfunc=np.sum)
                p_table.plot(ax = ax[1], kind='line', title=titles[2 * i + 1])
            # end for ax, i
            # ajuste da margem entre os subplots
            plt.subplots_adjust(hspace=0.6)
            fig.tight_layout()
            #fig.show()
            figs.append(fig)
            # Export plots do PDF file
            with PdfPages(os.path.join(full_output_path, 'Plot_' + acc_value_col_name + '_CompareStrategy_Det_Rob.pdf')) as pdf:
                for fig in figs:
                    plt.figure(fig.number)
                    pdf.savefig()
                    # pdf.savefig(plt.gcf())
                plt.close()
        # end for y
        print('           Done.')
    # end if


# Originally in http://stackoverflow.com/questions/11882393/matplotlib-disregard-outliers-when-plotting
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


if __name__ == "__main__":
    main(sys.argv[1:])
