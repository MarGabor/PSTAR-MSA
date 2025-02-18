import os
import subprocess
from datetime import datetime
import traceback
import matplotlib.pyplot
import pandas
import argparse
import matplotlib.pyplot
import numpy as np
from scipy import stats
#from sklearn.linear_model import LinearRegression
#from sklearn.metrics import r2_score

def set_script_path():
    
    global SCRIPT_DIR
    SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
    os.chdir(SCRIPT_DIR)

#writing error messages to log
#void (str)
def errorFct(errorMsg):
    try:
        print(errorMsg)
        log_path = os.path.join(SCRIPT_DIR, 'job_plot.log')
        with open(log_path, 'a') as errFile:
            logEntry = "\n [%s] %s" % ((datetime.now()).strftime("%d/%m/%Y %H:%M:%S"), errorMsg)
            errFile.write(logEntry)
    except:
        print(traceback.format_exc())
        open_file_err_msg = "Error log file could not be opened."
        if errorMsg == 'mp':
            return open_file_err_msg
        print(open_file_err_msg)
        exit(1)

#create dir in specified path
#void str        
def create_dir_safely(path): 
    
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    except:
        errMsg = "Failed to create \"%s\" directory. Insufficient writing priviliges in specified output path?" % os.path.split(path)[1]
        errorFct(errMsg)
        exit(1)

def get_job_name_list(jobs_path):

    job_name_list = []
    for job_name in os.listdir(jobs_path):
        if os.path.isdir(os.path.join(jobs_path, job_name)):
            job_name_list.append(job_name)

    return job_name_list

def get_red_index_path_dict(jobs_path, job_name_list):

    aln_data_path_dict = {}
    for job_name in job_name_list:
        DATA_path = os.path.join(jobs_path, job_name, 'DATA')
        aln_data_path_list = []
        for aln_data_name in os.listdir(DATA_path):
            if os.path.isdir(os.path.join(DATA_path, aln_data_name)):
                aln_data_path_list.append(os.path.join(DATA_path, aln_data_name))
        aln_data_path_dict[job_name] = aln_data_path_list.copy()
    
    red_index_path_dict = {}
    for job_name, aln_data_path_list in aln_data_path_dict.items():
        red_index_path_list = []
        for aln_data_path in aln_data_path_list:
            aln_data_name = os.path.split(aln_data_path)[1]
            red_index_name = 'red_'+aln_data_name+'_exact_match.index'
            red_index_path = os.path.join(aln_data_path, red_index_name)
            red_index_path_list.append(red_index_path)
        red_index_path_dict[job_name] = red_index_path_list.copy()

    return red_index_path_dict

def get_score_file_path_dict(jobs_path, job_name_list):

    score_file_path_dict = {}
    for job_name in job_name_list:
        EVAL_path = os.path.join(jobs_path, job_name, 'EVAL')
        score_file_path = os.path.join(EVAL_path, 'scores_per_alignment.tsv')
        score_file_path_dict[job_name] = score_file_path

    return score_file_path_dict

def get_sample_num_dict(red_index_path_dict):

    sample_num_dict = {}
    for job_name, red_index_path_list in red_index_path_dict.items():
        seq_num_set = set()
        for red_index_path in red_index_path_list:
            red_index_df = pandas.read_csv(red_index_path, sep='\t')
            seq_num = len(red_index_df)
            seq_num_set.add(seq_num)
        if len(seq_num_set) == 1:
            seq_num = seq_num_set.pop()
            sample_num = ((seq_num * seq_num) - seq_num)/2
            if sample_num <= 1000:
                sample_num_dict[job_name] = sample_num
            else:
                sample_num_dict[job_name] = 1000
        else:
            err_msg = 'Number of sequences differs between reduced alignment index files in job %s.' % (job_name)
            errorFct(err_msg)
            get_sample_num_dict[job_name] = 'NA'

    return sample_num_dict

def write_aggregated_score_file(score_file_path_dict, sample_num_dict, outputs_path, aln_len_dict):

    score_df_dict = {}
    aligner_list = ['muscle','learnMSA','learnMSA_deprecated','learnMSA2','homfamref','mafft_sparsecore',
                    'clustalo','famsa','t_coffee','baseline']
    for job_name, score_file_path in score_file_path_dict.items():
        job_score_list = []
        try:
            job_score_df = pandas.read_csv(score_file_path, sep='\t')
        except FileNotFoundError:
            for aligner in aligner_list:
                job_score_list.append('NA')
            job_score_list.append(sample_num_dict[job_name])
            job_score_list.append(aln_len_dict[job_name])
            continue
        for aligner in aligner_list:
            for row in job_score_df.itertuples():
                if aligner in row._1:
                    job_score_list.append(row.score)
                    break
            else:
                job_score_list.append('NA')
        job_score_list.append(sample_num_dict[job_name])
        job_score_list.append(aln_len_dict[job_name])
        score_df_dict[job_name] = job_score_list.copy()

    column_list = aligner_list.copy()
    column_list.append('num_of_samples')
    column_list.append('num_of_aln_seq')
    total_score_df = pandas.DataFrame.from_dict(score_df_dict, orient='index', columns=column_list)
    csv_path = os.path.join(outputs_path, 'jacc_aggregated_scores_ref.tsv')
    total_score_df.to_csv(csv_path, sep='\t')

    return total_score_df

def get_aln_lengths(alns_path, job_name_list):

    aln_path_dict = {}
    aln_name_list = os.listdir(alns_path)
    for job_name in job_name_list:
        for aln_dir_name in aln_name_list:
            if aln_dir_name.startswith(job_name):
                aln_path = os.path.join(alns_path, aln_dir_name)
                aln_path_dict[job_name] = aln_path
                break

    aln_len_dict = {}
    for job_name in job_name_list:
        aln_len_counter = 0
        with open(os.path.join(aln_path_dict[job_name],'learnMSA_'+os.path.split(aln_path_dict[job_name])[1]), 'r') as aln_file:
            for line in aln_file:
                if line.startswith('>'):
                    aln_len_counter += 1
        aln_len_dict[job_name] = aln_len_counter

    return aln_len_dict

def remove_outliers(agg_data_df):
    #new_agg_data_df = pandas.DataFrame()
    column_list = ['index','alignment','muscle','learnMSA','learnMSA_deprecated','learnMSA2','homfamref','mafft_sparsecore',
                   'clustalo','famsa','t_coffee','baseline','num_of_samples','num_of_aln_seq']
    row_list = []
    for row in agg_data_df.itertuples():
        if int(row.num_of_samples) > 1:
            row_list.append(row)

    new_agg_data_df = pandas.DataFrame.from_records(row_list, columns=column_list)

    return new_agg_data_df

def calc_weighted_mean(two_columns, aln_name):
    summand_list = []
    denom = 0
    for row in two_columns.itertuples():
        summand_list.append(float(row.aln_name)*float(row.num_of_samples))
        denom += int(row.num_of_samples)
    numerator = 0
    for summand in summand_list:
        numerator += summand
    weighted_mean = numerator/denom

    return weighted_mean

def sort_row_descending(row, agg_data_df):
    score_aligner_dict = {}
    float_list = []
    score_aligner_dict[row.muscle] = 'muscle'
    float_list.append(row.muscle)
    score_aligner_dict[row.learnMSA] = 'learnMSA'
    float_list.append(row.learnMSA)
    score_aligner_dict[row.learnMSA_deprecated] = 'learnMSA_deprecated'
    float_list.append(row.learnMSA_deprecated)
    score_aligner_dict[row.learnMSA2] = 'learnMSA2'
    float_list.append(row.learnMSA2)
    score_aligner_dict[row.homfamref] = 'homfamref'
    float_list.append(row.homfamref)
    score_aligner_dict[row.mafft_sparsecore] = 'mafft_sparsecore'
    float_list.append(row.mafft_sparsecore)
    score_aligner_dict[row.clustalo] = 'clustalo'
    float_list.append(row.clustalo)
    score_aligner_dict[row.famsa] = 'famsa'
    float_list.append(row.famsa)
    score_aligner_dict[row.t_coffee] = 't_coffee'
    float_list.append(row.t_coffee)
    score_aligner_dict[row.baseline] = 'baseline'
    float_list.append(row.baseline)

    sorted_float_list = sorted(float_list, key=float, reverse=True)
    col_plot_order_list = []
    for score in sorted_float_list:
        col_plot_order_list.append(agg_data_df.columns.get_loc(score_aligner_dict[score]))

    return col_plot_order_list

def estimate_coef(x, y):
  # number of observations/points
  n = np.size(x)

  # mean of x and y vector
  m_x = np.mean(x)
  m_y = np.mean(y)

  # calculating cross-deviation and deviation about x
  SS_xy = np.sum(y*x) - n*m_y*m_x
  SS_xx = np.sum(x*x) - n*m_x*m_x

  # calculating regression coefficients
  b_1 = SS_xy / SS_xx
  b_0 = m_y - b_1*m_x

  return (b_0, b_1)

def plot_aggregated_data_stem(agg_data_full_path):
    agg_data_df_with_outliers = pandas.read_csv(agg_data_full_path, sep='\t')

    agg_data_df_no_outliers = remove_outliers(agg_data_df_with_outliers)
    #agg_data_df_no_outliers = agg_data_df_with_outliers

    #calculating means and writing legend labels
    muscle_legend_lab = 'muscle, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['muscle'].mean(),4), round(agg_data_df_no_outliers['muscle'].median(),4))
    learnMSA_legend_lab = 'learnMSA, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['learnMSA'].mean(),4), round(agg_data_df_no_outliers['learnMSA'].median(),4))
    learnMSA_dep_legend_lab = 'learnMSA_deprecated, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['learnMSA_deprecated'].mean(),4), round(agg_data_df_no_outliers['learnMSA_deprecated'].median(),4))
    learnMSA2_legend_lab = 'learnMSA2, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['learnMSA2'].mean(),4), round(agg_data_df_no_outliers['learnMSA2'].median(),4))
    homfamref_legend_lab = 'homfamref, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['homfamref'].mean(),4), round(agg_data_df_no_outliers['homfamref'].median(),4))
    mafft_legend_lab = 'mafft_sparsecore, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['mafft_sparsecore'].mean(),4), round(agg_data_df_no_outliers['mafft_sparsecore'].median(),4))
    clustalo_legend_lab = 'clustalo, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['clustalo'].mean(),4), round(agg_data_df_no_outliers['clustalo'].median(),4))
    famsa_legend_lab = 'famsa, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['famsa'].mean(),4), round(agg_data_df_no_outliers['famsa'].median(),4))
    t_coffee_legend_lab = 't_coffee, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['t_coffee'].mean(),4), round(agg_data_df_no_outliers['t_coffee'].median(),4))
    baseline_legend_lab = 'baseline, M/Mdn: %s/%s' % (round(agg_data_df_no_outliers['baseline'].mean(),4), round(agg_data_df_no_outliers['baseline'].median(),4))

    column_ind_label_dict = {}
    column_ind_markerfmt_dict = {}
    column_ind_linefmt_dict = {}
    column_ind_aligner_name_dict = {}
    for column in agg_data_df_no_outliers.columns.to_list():
        name = column
        if column == 'muscle':
            label = muscle_legend_lab
            markerfmt = 'C0.'
            linefmt = 'C0-'
        elif column == 'learnMSA':
            label = learnMSA_legend_lab
            markerfmt = 'C1.'
            linefmt = 'C1-'
        elif column == 'learnMSA_deprecated':
            label = learnMSA_dep_legend_lab
            markerfmt = 'C2.'
            linefmt = 'C2-'
        elif column == 'learnMSA2':
            label = learnMSA2_legend_lab
            markerfmt = 'C3.'
            linefmt = 'C3-'
        elif column == 'homfamref':
            label = homfamref_legend_lab
            markerfmt = 'C4.'
            linefmt = 'C4-'
        elif column == 'mafft_sparsecore':
            label = mafft_legend_lab
            markerfmt = 'C5.'
            linefmt = 'C5-'
        elif column == 'clustalo':
            label = clustalo_legend_lab
            markerfmt = 'C6.'
            linefmt = 'C6-'
        elif column == 'famsa':
            label = famsa_legend_lab
            markerfmt = 'C7.'
            linefmt = 'C7-'
        elif column == 't_coffee':
            label = t_coffee_legend_lab
            markerfmt = 'C8.'
            linefmt = 'C8-'
        elif column == 'baseline':
            label = baseline_legend_lab
            markerfmt = 'C9.'
            linefmt = 'C9-'
        else:
            continue

        column_ind_aligner_name_dict[agg_data_df_no_outliers.columns.get_loc(column)] = name
        column_ind_label_dict[agg_data_df_no_outliers.columns.get_loc(column)] = label
        column_ind_markerfmt_dict[agg_data_df_no_outliers.columns.get_loc(column)] = markerfmt
        column_ind_linefmt_dict[agg_data_df_no_outliers.columns.get_loc(column)] = linefmt

    best_score_in_row_list = []
    for row in agg_data_df_no_outliers.itertuples():
        col_order_list = sort_row_descending(row, agg_data_df_no_outliers)
        best_scoring_col_num = col_order_list[0]
        best_score_in_row_list.append(getattr(row, column_ind_aligner_name_dict[best_scoring_col_num]))

    #agg_data_df_no_outliers['best_score'] = best_score_in_row_list

    agg_data_df = agg_data_df_no_outliers.sort_values(by='homfamref', ascending=False, ignore_index=True)

    print(agg_data_df)

    #muscle_weight_mean = calc_weighted_mean(agg_data_df['muscle','num_of_samples'], 'muscle')

    custom_xtick_labels = []
    main_fig = matplotlib.pyplot.figure()
    main_ax = main_fig.add_subplot()
    for row in agg_data_df.itertuples():
        col_plot_order_list = sort_row_descending(row, agg_data_df)
        #get col_num of second best after homfamref, or best if better than homfamref
        if column_ind_aligner_name_dict[col_plot_order_list[0]] != 'homfamref':
            marked_best = col_plot_order_list[0]
            del col_plot_order_list[0]
            col_plot_order_list.append(marked_best)
        elif column_ind_aligner_name_dict[col_plot_order_list[0]] == 'homfamref':
            marked_best = col_plot_order_list[1]
            del col_plot_order_list[1]
            col_plot_order_list.append(marked_best)
        #counter = 0
        for col_num in col_plot_order_list:
            #if counter == len(col_plot_order_list)-1:
                #column_ind_linefmt_dict[col_num] = column_ind_linefmt_dict[col_num][0:2]+':'

            main_ax.stem([row.Index], [getattr(row, column_ind_aligner_name_dict[col_num])], markerfmt=column_ind_markerfmt_dict[col_num],
                         linefmt=column_ind_linefmt_dict[col_num], label=column_ind_label_dict[col_num])
            #counter += 1

        
        
        #revert tmp changes to linefmt dict
        #column_ind_linefmt_dict[col_num] = column_ind_linefmt_dict[col_num][0:2]+'-'

        custom_xtick_labels.append(str(row.alignment)+'/'+str(int(row.num_of_samples))+'/'+str(row.num_of_aln_seq))

    patch_list = []
    label_list = []
    for key, name in column_ind_aligner_name_dict.items():

        patch = matplotlib.patches.Patch(color=column_ind_markerfmt_dict[key][0:2])
        patch_list.append(patch)
        
        label_list.append(column_ind_label_dict[key])

    matplotlib.pyplot.legend(loc='lower left', bbox_to_anchor=(-0.15,1.002), fontsize=6, ncol=3, labels=label_list, handles=patch_list)
    
    xtick_range = range(0, len(agg_data_df['alignment']))

    #get linear regression line for baseline
    #linear_regressor = LinearRegression()
    #linear_regressor.fit(xtick_range, agg_data_df['baseline'])
    #baseline_pred = linear_regressor.predict(xtick_range)
    #rsq = r2_score(agg_data_df['baseline'], baseline_pred)
    #m = linear_regressor.coef_
    #c = linear_regressor.intercept_

    #linear regression without scikit
    Y = agg_data_df['baseline'].tolist()
    X = np.array(list(xtick_range))
    slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
    #intercept, slope = estimate_coef(list(xtick_range), Y)
    rsq = r_value**2
    baseline_pred = []
    baseline_pred.clear()
    for x in xtick_range:
        baseline_pred.append((slope*(x))+intercept)
    out_msg= "y = %sx + %s; r-squared: %s" % (str(slope),str(intercept), str(rsq))
    print(out_msg)

    #add baseline prediction to plot
    main_ax.plot(xtick_range, baseline_pred, color='C9')

    main_ax.set_xticks(xtick_range, custom_xtick_labels, rotation=90, ha='center', fontsize=2)
    main_fig.savefig('data_plot_ref.png', dpi=1200)

def plot_hist(agg_data_full_path):
    agg_data_df = pandas.read_csv(agg_data_full_path, sep='\t')

    muscle_diff_list = []
    learnMSA_diff_list = []
    learnMSA_dep_diff_list = []
    learnMSA2_diff_list = []
    mafft_sparsecore_diff_list = []
    clustalo_diff_list = []
    famsa_diff_list = []
    t_coffee_diff_list = []
    for row in agg_data_df.itertuples():
        muscle_diff_list.append(float(row.homfamref)-float(row.muscle))
        learnMSA_diff_list.append(float(row.homfamref)-float(row.learnMSA))
        learnMSA_dep_diff_list.append(float(row.homfamref)-float(row.learnMSA_deprecated))
        learnMSA2_diff_list.append(float(row.homfamref)-float(row.learnMSA2))
        mafft_sparsecore_diff_list.append(float(row.homfamref)-float(row.mafft_sparsecore))
        clustalo_diff_list.append(float(row.homfamref)-float(row.clustalo))
        famsa_diff_list.append(float(row.homfamref)-float(row.famsa))
        t_coffee_diff_list.append(float(row.homfamref)-float(row.t_coffee))

    agg_data_df['muscle_diff'] = muscle_diff_list
    agg_data_df['learnMSA_diff'] = learnMSA_diff_list
    agg_data_df['learnMSA_deprecated_diff'] = learnMSA_dep_diff_list
    agg_data_df['learnMSA2_diff'] = learnMSA2_diff_list
    agg_data_df['mafft_sparsecore_diff'] = mafft_sparsecore_diff_list
    agg_data_df['clustalo_diff'] = clustalo_diff_list
    agg_data_df['famsa_diff'] = famsa_diff_list
    agg_data_df['t_coffee_diff'] = t_coffee_diff_list

    #define xticks
    custom_xtick_labels = agg_data_df['alignment']
    xtick_range = range(0, len(agg_data_df['alignment']))

    #initialize figures and axes
    muscle_fig = matplotlib.pyplot.figure()
    muscle_ax = muscle_fig.add_subplot()
    learnMSA_fig = matplotlib.pyplot.figure()
    learnMSA_ax = learnMSA_fig.add_subplot()
    learnMSA_dep_fig = matplotlib.pyplot.figure()
    learnMSA_dep_ax = learnMSA_dep_fig.add_subplot()
    learnMSA2_fig = matplotlib.pyplot.figure()
    learnMSA2_ax = learnMSA2_fig.add_subplot()
    mafft_sparsecore_fig = matplotlib.pyplot.figure()
    mafft_sparsecore_ax = mafft_sparsecore_fig.add_subplot()
    clustalo_fig = matplotlib.pyplot.figure()
    clustalo_ax = clustalo_fig.add_subplot()
    famsa_fig = matplotlib.pyplot.figure()
    famsa_ax = famsa_fig.add_subplot()
    t_coffee_fig = matplotlib.pyplot.figure()
    t_coffee_ax = t_coffee_fig.add_subplot()

    #provide y-data and plot histograms
    muscle_ax.bar(xtick_range, agg_data_df['muscle_diff'])
    muscle_ax.set_xticks(xtick_range, custom_xtick_labels, rotation=90, ha='center', fontsize=2)
    muscle_ax.set_ylim([-0.2, 0.2])
    learnMSA_ax.bar(xtick_range, agg_data_df['learnMSA_diff'])
    learnMSA_ax.set_xticks(xtick_range, custom_xtick_labels, rotation=90, ha='center', fontsize=2)
    learnMSA_ax.set_ylim([-0.2, 0.2])

    #save figures to file
    
    muscle_fig.savefig('muscle_hist.png', dpi=1200)
    learnMSA_fig.savefig('learnMSA_hist.png', dpi=1200)

def main():

    set_script_path()

    outputs_path = '../plots/'
    create_dir_safely(outputs_path)
    jobs_path = '../results_ref/'
    alns_path = '../jobs_to_run/'
    agg_data_full_path = os.path.join(outputs_path, 'jacc_aggregated_scores_ref.tsv')

    argParser = argparse.ArgumentParser()
    argParser.add_argument("-ps", "--plotstem", default=0, action="count", help="Creates stem plots from aggregated data.", required=False)
    argParser.add_argument("-ph", "--plothist", default=0, action="count", help="Creates histograms from aggregated data.", required=False)
    args = argParser.parse_args()

    if args.plotstem>0:
        plot_aggregated_data_stem(agg_data_full_path)

    if args.plothist>0:
        plot_hist(agg_data_full_path)
    
    job_name_list = get_job_name_list(jobs_path)

    red_index_path_dict = get_red_index_path_dict(jobs_path, job_name_list)

    aln_len_dict = get_aln_lengths(alns_path, job_name_list)

    sample_num_dict = get_sample_num_dict(red_index_path_dict)

    score_file_path_dict = get_score_file_path_dict(jobs_path, job_name_list)

    total_score_df = write_aggregated_score_file(score_file_path_dict, sample_num_dict, outputs_path, aln_len_dict)


    
if __name__ == "__main__":
    main()

