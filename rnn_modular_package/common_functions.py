"""
common functions used by the scripts.

"""

import datetime
import os
import pickle
from collections import OrderedDict

import numpy as np
import pandas as pd
from dateutil import relativedelta

np.warnings.filterwarnings('ignore')


def date_fmt():
    """
    function to set the time format
    """
    return '%Y-%m-%d %H:%M:%S'


def get_time():
    """
    function to get time in a set format.
    :return: time in str format.
    
    """
    return datetime.datetime.now().strftime(date_fmt())


def add_args(given_func, *args, **keywords):
    """
    create a new function with additional arguments and key-arguments
    
    :param given_func: original function
    :param args: arguments supplied
    :param keywords: key-arguments supplied
    :return: new function
    """
    given_func.func = []
    given_func.args = args
    given_func.keywords = keywords
    
    return given_func


def write_now(wr_file, statement, rt, msg=''):
    """
    function to write files.
    """
    # print(statement, rt, msg)
    wr_file.write("%s %s %s\n" % (statement, rt, msg))
    wr_file.flush()


def cal_runtime(st, ed):
    """
    function to get time difference.
    :param st: start time
    :param ed: end time
    :return: time difference
    """
    d1 = datetime.datetime.strptime(st, date_fmt())
    d2 = datetime.datetime.strptime(ed, date_fmt())
    
    td = relativedelta.relativedelta(d2, d1)
    # print(td)
    return "%s minutes" % str(
            np.round((td.days * 24 + td.hours * 60 + td.minutes + float(td.seconds) / 60), decimals=2))


def write_common_log(list_file_names=(), clog_filename='combined_rnn_log.log'):
    """

    :param list_file_names: list of individual log files
    :param clog_filename: name of combined log file
    """
    common_log = open(clog_filename, 'w')
    list_file_names = sorted(list_file_names)
    arr_time = []
    for fn in list_file_names:
        ind_log = open(fn, 'r').readlines()
        write_now(common_log, "Individual log File", fn)
        for line in ind_log:
            line = line.rstrip('\n')
            write_now(common_log, "", line)
            if "Total run time " in line:
                arr_line = [line for line in line.split(" ") if line != ""]
                arr_time.append(float(arr_line[-2]))
    write_now(common_log, "\nAverage iteration run time\t", np.round(np.mean(arr_time), decimals=2), 'minutes')


def output_in_text(input_array=(), work_dir=''):
    """
    module read pre-saved *.pkl files and convert them in text files
    :param input_array: list of input pkl files to be converted to txt files
    :param work_dir: directory from which calculations are launched
    :return: tuple containing txt file names. These file contains predicted sequences and RNN score arranged in
    ascending manner.
    """
    
    arr_files = [x.replace('.pkl', '.txt') for x in input_array]
    
    for po_file in input_array:
        i_filename = os.path.join(work_dir, po_file)
        o_filename = os.path.join(work_dir, po_file.replace('../', '').replace('pkl', 'txt'))
        if not os.path.exists(o_filename):
            d = OrderedDict(pickle.load(open(i_filename, 'rb')))
            
            # sort the values in descending order
            d_sorted_by_value = OrderedDict(sorted(d.items(), key=lambda x: (x[1], x[0]), reverse=True))
            
            sorted_dict = OrderedDict()
            for k, v in d_sorted_by_value.items():
                sorted_dict.update({k: v})
            
            o_file_h = open(o_filename, 'w')
            for b in sorted_dict:
                o_file_h.write('%s\t%.5f\n' % (b, sorted_dict[b][0][0]))
            o_file_h.close()
    
    return arr_files


def ranking_bootstrapping(list_random_files, ref=''):
    """
    ranking_wrt_reference module: read multiple input files and return the ranks for each of the sequence.
    
    :param list_random_files: list containing name of output text files (randomized set)
    :param ref: reference to reference file, among the list_random_files.
    :return: dictionary consisting sequences and their ranks in each of the file in list_random_files.
    """
    
    rank_dic = OrderedDict()
    ref_file_data = open(ref, 'r').readlines()
    
    for line in ref_file_data:
        seq_i = line.rstrip('\n').split('\t')
        rank_dic.update({seq_i[0]: []})
    
    for j in range(0, len(ref_file_data)):
        seq_o = ref_file_data[j].rstrip('\n').split('\t')
        rank_dic[seq_o[0]].append([j + 1, seq_o[1]])
    
    for fn in list_random_files[1:]:
        comp_file_data = open(fn, 'r').readlines()
        
        for k in range(0, len(comp_file_data)):
            seq_m = comp_file_data[k].rstrip('\n').split('\t')
            rank_dic[seq_m[0]].append([k + 1, seq_m[1]])
    
    return rank_dic


def average_ranking_and_rnn_score(odict, o_filename='rank_across_trials.csv', txt_file_names=(), header=(),
                                  file_format='csv'):
    """
    function to calculate average ranks and RNN score
    :param odict:
    :param o_filename:
    :param txt_file_names:
    :param header:
    :param file_format:
    :return:
    """
    
    stack_ranks = []
    for seq in list(odict.keys()):
        loc_arr = [seq, np.hstack(odict[seq])]
        stack_ranks.append(np.hstack(loc_arr))
    
    if not header:
        header = ['seq']
        for i in range(0, len(txt_file_names)):
            hn = txt_file_names[i].replace('all_predict_', '').replace('.txt', '')
            header.append("%s_%s" % ("rank", hn))
            header.append("%s_%s" % ("predicted_score", hn))
    
    df = pd.DataFrame(stack_ranks, columns=header)
    
    average_score = []
    average_rank = []
    
    if len(df.columns) > 3:
        for col in df.columns:
            if 'predicted_score_' in col:
                if len(average_score) == 0:
                    average_score = [float(x) for x in df[col].values]
                else:
                    average_score = np.vstack((average_score, [float(x) for x in df[col].values]))
            elif 'rank_' in col:
                if len(average_rank) == 0:
                    average_rank = [float(x) for x in df[col].values]
                else:
                    average_rank = np.vstack((average_rank, [float(x) for x in df[col].values]))
        
        df['Av_RNN_score'] = np.round(np.mean(average_score, axis=0), decimals=5)
        df['Av_RNN_rank'] = np.round(np.mean(average_rank, axis=0), decimals=5)
    
    elif len(df.columns) == 3:
        for col in df.columns:
            if 'predicted_score_' in col:
                df['Av_RNN_score'] = df[col]
            if 'rank_' in col:
                df['Av_RNN_rank'] = df[col]
    
    df = df.sort_values(by=['Av_RNN_score', 'Av_RNN_rank'], ascending=[False, True])
    df = df.set_index(['seq'])
    
    # write file in csv format. Append ".csv" to file name, if not provided in o_filename.
    if file_format not in o_filename[-3:]:
        df.to_csv("%s.%s" % (o_filename, "csv"), sep='\t', float_format='%.5f')
    else:
        df.to_csv(o_filename, sep='\t', float_format='%.5f')
    return df
