"""
parameters file consists of all the parameters need to build the RNN.

specify type of RNN needed in rnn_type.
    currently two types of RNN can be constructed.
        1. LSTM: refers to tf.contrib.rnn.BasicLSTMCell
        2. GRU: refers to tf.contrib.rnn.GRUCell(lstm_size)
 

"""
import os
import sys

home_dir = os.environ['HOME']
code_dir = 'PycharmProjects/kdel/rnn_modular_package'
sys.path.append(os.path.join(home_dir, code_dir))

# parameters for input_preparation script:
meta_dir = home_dir + "/PycharmProjects/kdel/raw_data"  # directory where input raw data file is stored.
checkpoint_store_dir = home_dir + "/PycharmProjects/kdel/stored_checkpoints/batch1/" # store checkpoints in this directory
input_excel = "SERCaMP_Library_Gene_List_with_mass_spec.xlsx"
sequence_key = "tail"  # header for sequence data in excel sheet
response_key = "Tg Response - 8hr"  # header for response data in excel sheet
output_df = "old_tail_response.df"  # file name for the processed output dataframe, to be used in training and testing
prediction_df = "new_SERCaMP_Response_Data.df"  # file name for the processed output dataframe, to be used for prediction
