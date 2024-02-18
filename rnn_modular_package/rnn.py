"""
core script to create, train, and test the model.

Input: Read output dataframe file generated using input_preparation.py. File name is defined as output_df in
path_parameters*.py file. See sample path file (path_parameters_sample.py). parameters for constructing RNN are
described in rnn_parameters.py.

This is a machine learning script developed to create a model to predict "Tg Response".
    workflow:
        1. Read output dataframe file generated using input_preparation file.
        2. Train and test the RNN model using 10:1 ratio of the data given.
        3. Prediction of new sequences and their "Tg Response" values based on the model.
        4. Save the trained model in checkpoints directory defined using -c flag in command line.
        5. Checkpoints are read again and sequences are filtered (a.a. that appears more than twice in the data).
        6. Output is stored in filename defined by prediction_out_file_name variable in path_parameters*.py file.

AUTHOR

    Ravi Kumar Verma

VERSION/DATE

    31st July, 2018

"""

import getopt
import os
import sys
from collections import Counter
from itertools import product

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

import importlib

working_dir = os.getcwd()
sys.path.append(working_dir)

rnn_parameters = importlib.import_module('rnn_parameters', package=None)
common_functions = importlib.import_module('common_functions', package=None)
initialize_graph = importlib.import_module('initialize_graph', package=None)

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '4'
np.warnings.filterwarnings('ignore')


def get_input_from_df(input_df=None, input_aa2num2=None, sequence_key=None, response_key=None):
    """
    :param input_df:
    :param input_aa2num2:
    :return: [l_whole_x, l_whole_y] where
            l_whole_x = sequences in one-hot encoding
            l_whole_y = their corresponding response values in excel sheet
    
    :param sequence_key: header for sequence data in excel sheet
    :param response_key: header for response data in excel sheet
    """
    
    df = pd.read_pickle(input_df)
    
    # print(df.columns, df[sequence_key])
    # our data-set
    
    l_whole_x = np.array([[input_aa2num2[s] for s in seq] for seq in df[sequence_key]])
    if response_key in df.columns:
        print(response_key)
        l_whole_y = np.array(df[response_key].values)
    else:
        l_whole_y = []
    
    return [l_whole_x, l_whole_y, df]


def convert_data_to_hot_encoding(input_aas=None):
    """
    :type input_aas: amino-acid sequence array defined in parameters file

    :return:
        l_aa2num2 : dictionary containing amino acids as keys and the one-hot encoding as values.
        l_num2aa: reverse dictionary, containing one-hot encoding as keys and amino-acid identities as values.
    """
    
    l_aa2num2 = {k: v for v, k in enumerate(input_aas)}
    l_num2aa = {v: k for k, v in l_aa2num2.items()}
    
    # print(l_aa2num2)
    # print(l_num2aa)
    return [l_aa2num2, l_num2aa]


def get_all_possibilities(input_aa2num2=None, input_df=None):
    # get the a.a. that appears more than twice in the data
    
    temp = [[] for _ in range(7)]
    for i in range(7):
        al = [input_aa2num2[x[i]] for x in input_df['tail']]
        counter = Counter(al)
        temp[i] = [key for key in counter if counter[key] > 2]
    
    all_comb = product(*temp)
    return all_comb


def print_parameters_used(b, c, d, e, f, g, h, i):
    print("\n\tUSED PARAMETERS:")
    print("\n\trnn :", b)
    print("\tjobtype:", c)
    print("\tepochs:", d)
    print("\tcheckpoint dir:", e)
    print("\tcheckpoint meta file:", f)
    print("\tprediction output file:", g)
    print("\taa2num2:", h)
    print("\tnum2aa:", i, "\n\n")


def show_help():
    print("\n\n\t(-:  RNN SCRIPT  :-)\n\n")
    print("Written by:RAVI KUMAR VERMA\n")
    print("\nUSAGE:")
    print("\tpython3 rnn.py -p path_parameters.py -r rnn_type -j jobtype -e epochs -c checkpoint_dir")
    print("\n\nOPTIONS:\n")
    print("\t-h or --help:\tprint help")
    print("\t-p or --pfile:\tfile with paths to input files. See path_parameters_sample.py")
    print("\t-r or --rnn:\tuse either GRU or LSTM\n\t\t\tDEFAULT: LSTM")
    print("\t-e or --epochs:\tno of optimization iterations\n\t\t\tDEFAULT: 3000")
    print("\t-c or --ckdir:\tname of checkpoint directory where checkpoints to be saved\n\t\t\tDEFAULT: checkpoints")
    print("\t-j or --jobtype: use prediction\n\t\t\tDEFAULT: prediction")
    
    print(
            "\t\t\tif jobtype== validation:\tcode performs cross validation")
    print("\t\t\tif jobtype== prediction:\ttrain on whole dataset and perform analysis")
    print("\t-g or --gpu: id for the gpu:")
    print("\t\tuse any integer between 0 and 7")
    
    print("\nIMPORTANT CONSIDERATIONS:")
    print("\tUse different checkpoint directory names while doing multiple runs.")
    print("\nWORKFLOW:")
    print("\trnn.py calls initialize_graph.py script contain functions to generate and analyze RNN model.\n")
    print("\tCode performs following operations:")
    print(
            "\t1. Read dataframe file generated by input_preparation file as input. The dataframe filename is read "
            "from "
            "rnn_parameters.py file.")
    print("\t2. Train and test the RNN model using 10:1 ratio of the input data given.")
    print("\t3. Prediction of new sequences and their 'Tg Response' values based on the model.")
    print("\t4. Save the trained model in checkpoints directory defined as using -c or --ckdir flags.")
    print(
            "\t5. Checkpoints are read again and sequences are filtered (a.a. that appears more than twice in the "
            "data).")
    print(
            "\t6. Output is stored in filename defined by prediction_out_file_name variable. File name is derived from"
            "-c or --ckdir flags.\n")


def main(argv):
    if argv:
        try:
            opts, argument = getopt.getopt(argv, "hp:r:j:e:c:l:g:t:a::",
                                           ["pfile=", "rnn=", "jobtype=", "epochs=", "ckdir=", "log=", "gpu=",
                                            "training=", "artificial_sequences="])
        except getopt.GetoptError:
            print('\n', 'ERROR:', '\n\t', 'Check your arguments', '\n')
            show_help()
            sys.exit(2)
        
        rt = "LSTM"
        jt = "prediction"
        ep = 3000
        ck = "checkpoints"
        lg = "log.log"
        g = ''
        pfile = ''
        tr = ''
        afs = ''
        for opt, arg in opts:
            if opt == "-h":
                show_help()
                sys.exit()
            elif opt in ("-p", "--pfile"):
                pfile = arg
            elif opt in ("-r", "--rnn"):
                rt = arg
            elif opt in ("-j", "--jobtype"):
                jt = arg
            elif opt in ("-e", "--epochs"):
                ep = arg
            elif opt in ("-c", "--ckdir"):
                ck = arg
            elif opt in ("-l", "--log"):
                lg = arg
            elif opt in ("-g", "--gpu"):
                g = arg
            elif opt in ("-t", "--training"):
                tr = arg
            elif opt in ("-a", "--artificial_sequences"):
                afs = arg
        return [pfile, rt, jt, ep, ck, lg, g, tr, afs]


if __name__ == "__main__":
    args = main(sys.argv[1:])
    print(args)
    parms = __import__(args[0].replace('.py', ''))
    rnn_type = args[1]
    jobtype = args[2]
    epochs = int(args[3])
    checkpoints_dir = args[4].rstrip('/')
    log = open(args[5], 'w')
    gpu_id = args[6]
    training = args[7]
    artificial_sequences = args[8]
    
    if not os.path.exists(parms.checkpoint_store_dir):
        #print("\n### ERROR : checkpoint_store_dir MENTIONED IN %s DOES NOT EXIST. CREATE IT MANUALLY.###\n" % args[0])
        #sys.exit(1)
        print ('\n No Checkpoint Saving Directory exists, create one')
        os.system('mkdir -p '+parms.checkpoint_store_dir)
    checkpoint = "%s_%s_%s" % (checkpoints_dir, rnn_type, epochs)
    
    aa2num2, num2aa = convert_data_to_hot_encoding(input_aas=rnn_parameters.aas)
    
    whole_x, whole_y, data_frame = get_input_from_df(input_df=os.path.join(parms.meta_dir, parms.output_df),
                                                     input_aa2num2=aa2num2,
                                                     sequence_key=parms.sequence_key,
                                                     response_key=parms.response_key)
    
    common_functions.write_now(log, "X length", str(len(whole_x)), '')
    common_functions.write_now(log, "Y length", str(len(whole_y)), '')
    common_functions.write_now(log, "batch_size", str(rnn_parameters.batch_size), '')
    
    all_combinations = ()
    
    # check if training and prediction to be done on separate data-sets
    if parms.output_df != parms.prediction_df:
        whole_x_prediction, whole_y_prediction, data_frame_prediction = get_input_from_df(
                input_df=os.path.join(parms.meta_dir, parms.prediction_df),
                input_aa2num2=aa2num2,
                sequence_key=parms.sequence_key,
                response_key=parms.response_key)
        
        # check if yes is given for -a or --artificial_sequences flags
        # if "yes" the code will generate all possible combination of sequences using >2 occurrence criteria per
        # position and will calculate RNN scores for them.
        if artificial_sequences == "yes":
            all_combinations = get_all_possibilities(input_aa2num2=aa2num2,
                                                     input_df=data_frame_prediction)
        
        # if "no" the code will use sequences given in the prediction data-frame and will calculate the scores for them.
        elif artificial_sequences == "no":
            all_combinations = whole_x_prediction
        
        # print(len(list(all_combinations)))
        common_functions.write_now(log, "X length prediction", str(len(whole_x_prediction)), '')
        common_functions.write_now(log, "Y length prediction", str(len(whole_y_prediction)), '')
        common_functions.write_now(log, "batch_size prediction", str(rnn_parameters.batch_size), '')
    
    # check if training and prediction data-frames are exactly the same the code generate all possible combination of
    #  sequences using >2 occurrence criteria per position and will calculate RNN scores for them.
    else:
        whole_x_prediction, whole_y_prediction, data_frame_prediction = get_input_from_df(
                input_df=os.path.join(parms.meta_dir, parms.prediction_df),
                input_aa2num2=aa2num2,
                sequence_key=parms.sequence_key,
                response_key=parms.response_key)
        # check if yes is given for -a or --artificial_sequences flags
        # if "yes" the code will generate all possible combination of sequences using >2 occurrence criteria per
        # position and will calculate RNN scores for them.
        if artificial_sequences == "yes":
            all_combinations = get_all_possibilities(input_aa2num2=aa2num2,
                                                     input_df=data_frame_prediction)
        
        # if "no" the code will use sequences given in the prediction data-frame and will calculate the scores for them.
        elif artificial_sequences == "no":
            all_combinations = whole_x_prediction
        
        common_functions.write_now(log, "prediction X length ", str(len(whole_x_prediction)), '')
        common_functions.write_now(log, "prediction Y length", str(len(whole_y_prediction)), '')
        common_functions.write_now(log, "prediction batch_size", str(rnn_parameters.batch_size), '')
    
    if jobtype == "validation":
        # k-fold cross validation
        """
            n_splits: refers to number of folds
            number of training and test data-sets.
    
        """
        
        prediction_out_file_name = "%s_%s%s" % ('all_predict', checkpoint, '.pkl')
        print_parameters_used(rnn_type, jobtype, epochs, os.path.join(parms.checkpoint_store_dir, checkpoints_dir),
                              checkpoint, prediction_out_file_name, aa2num2, num2aa)
        
        kf = KFold(n_splits=10, shuffle=True, random_state=33)
        
        iteration = 0
        log.write("Cross Validation:\n")
        # k-fold cross validation
        arr_validation_runtime = []
        for train_index, test_index in kf.split(whole_x):
            common_functions.write_now(log, "Training Iteration", iteration)
            s_time = common_functions.get_time()
            common_functions.write_now(log, "Iteration Start time", s_time)
            
            mse = initialize_graph.graph_initialization(lstm_size=rnn_parameters.lstm_size,
                                                        lstm_layers=rnn_parameters.lstm_layers,
                                                        batch_size=rnn_parameters.batch_size,
                                                        learning_rate=rnn_parameters.learning_rate,
                                                        epochs=epochs,
                                                        whole_x=whole_x,
                                                        whole_y=whole_y,
                                                        train_index=train_index,
                                                        test_index=test_index,
                                                        ckpt_dir=os.path.join(parms.checkpoint_store_dir,
                                                                              checkpoints_dir),
                                                        ckpt=checkpoint,
                                                        caltype='model',
                                                        cell=rnn_type,
                                                        fixed_seed=rnn_parameters.random_seed,
                                                        gpu_id=gpu_id,
                                                        cv=True)
            # calculate run time
            common_functions.write_now(log, "\nMSE:", mse, "\n")
            e_time = common_functions.get_time()
            common_functions.write_now(log, "Iteration end time", e_time)
            run_validation_time = common_functions.cal_runtime(s_time, e_time)
            common_functions.write_now(log, "run time", run_validation_time)
            arr_validation_runtime.append(run_validation_time)
            
            iteration += 1
        
        # calculate average run time
        common_functions.write_now(log, "average run time",
                                   np.mean(
                                           np.array(
                                                   [float(x.replace(' minutes', '')) for x in arr_validation_runtime])),
                                   'minutes')
    
    if jobtype == "prediction":
        log.write("Prediction:\n")
        # arr_training_runtime = []
        # arr_prediction_runtime = []
        
        checkpoint_dir_prediction = checkpoints_dir
        checkpoint_prediction = checkpoint
        
        prediction_out_file_name = "%s_%s%s" % ('all_predict', checkpoint_prediction, '.pkl')
        print_parameters_used(rnn_type, jobtype, epochs, checkpoint_dir_prediction, checkpoint_prediction,
                              prediction_out_file_name, aa2num2, num2aa)
        
        t_st_time = common_functions.get_time()
        if training == "yes":
            common_functions.write_now(log, "Training start time", t_st_time)
            
            # train on whole data-set
            mse = initialize_graph.graph_initialization(lstm_size=rnn_parameters.lstm_size,
                                                        lstm_layers=rnn_parameters.lstm_layers,
                                                        batch_size=rnn_parameters.batch_size,
                                                        learning_rate=rnn_parameters.learning_rate,
                                                        epochs=epochs,
                                                        whole_x=whole_x,
                                                        whole_y=whole_y,
                                                        train_index=range(len(whole_x)),
                                                        test_index=range(len(whole_x)),
                                                        ckpt_dir=os.path.join(parms.checkpoint_store_dir,
                                                                              checkpoints_dir),
                                                        ckpt=checkpoint_prediction,
                                                        caltype='model',
                                                        cell=rnn_type,
                                                        fixed_seed=rnn_parameters.random_seed,
                                                        gpu_id=gpu_id,
                                                        cv=True)# should be True for saving checkpoint file
            # calculate run time for training
            common_functions.write_now(log, "MSE:", mse)
            t_ed_time = common_functions.get_time()
            common_functions.write_now(log, "Training end time", t_ed_time)
            run_training_time = common_functions.cal_runtime(t_st_time, t_ed_time)
            common_functions.write_now(log, "Training run time", run_training_time)
            # arr_training_runtime.append(run_training_time)
        elif training == "no":
            pass
        
        p_st_time = common_functions.get_time()
        common_functions.write_now(log, "Prediction start time", p_st_time)
        
        initialize_graph.graph_initialization(lstm_size=rnn_parameters.lstm_size,
                                              lstm_layers=rnn_parameters.lstm_layers,
                                              batch_size=rnn_parameters.batch_size,
                                              learning_rate=rnn_parameters.learning_rate,
                                              epochs=epochs,
                                              whole_x=whole_x_prediction,
                                              whole_y=whole_y_prediction,
                                              ckpt_dir=os.path.join(parms.checkpoint_store_dir, checkpoints_dir),
                                              ckpt=checkpoint_prediction,
                                              caltype='data',
                                              cell=rnn_type,
                                              sequences_for_prediction=all_combinations,
                                              num2aa=num2aa,
                                              out_file_name=prediction_out_file_name,
                                              gpu_id=gpu_id,
                                              fixed_seed=rnn_parameters.random_seed)
        print(prediction_out_file_name)
        # calculate run time for prediction
        p_ed_time = common_functions.get_time()
        common_functions.write_now(log, "Prediction end time", p_ed_time)
        run_prediction_time = common_functions.cal_runtime(p_st_time, p_ed_time)
        common_functions.write_now(log, "Prediction run time", run_prediction_time)
        # arr_prediction_runtime.append(run_training_time)
        
        # calculate average run time for training
        common_functions.write_now(log, "\nTotal run time", common_functions.cal_runtime(t_st_time, p_ed_time))

