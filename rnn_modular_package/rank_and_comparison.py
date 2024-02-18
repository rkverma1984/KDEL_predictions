"""
script to rank output predictions across various trials

"""
import collections
import getopt
import os
import sys

user = os.path.expanduser('~')

working_dir = os.getcwd()
sys.path.append(working_dir)

try:
    from common_functions import output_in_text, average_ranking_and_rnn_score
except ImportError:
    from .common_functions import output_in_text, average_ranking_and_rnn_score

home_dir = os.environ['HOME']


def ranking_wrt_reference(data, ref=''):
    """
    ranking_wrt_reference module: read multiple input files and compare original data with the new trial data.
    """
    
    rank_dic = collections.OrderedDict()
    ref_file_data = open(ref, 'r').readlines()
    
    for line in ref_file_data:
        line = line.rstrip('\n')
        seq_i = line.split('\t')
        rank_dic.update({seq_i[0]: []})
    
    for j in range(0, len(ref_file_data)):
        seq_o = ref_file_data[j].rstrip('\n').split('\t')
        rank_dic[seq_o[0]].append([j + 1, seq_o[1]])
    
    for fnc in data:
        if ref not in fnc:
            comp_file_data = open(fnc, 'r').readlines()
            
            for k in range(0, len(comp_file_data)):
                seq_m = comp_file_data[k].rstrip('\n').split('\t')
                rank_dic[seq_m[0]].append([k + 1, seq_m[1]])
    
    return rank_dic


# following section calls functions defined above and do analysis ######


def show_help():
    print("\n\n\t(-:  PICKLE TO TXT CONVERSION and RANKING SCRIPT  :-)\n\n")
    print("Written by:RAVI KUMAR VERMA\n")
    print("\nUSAGE:")
    print("\tpython3 rank_and_comparison.py -o outputfilename -i inputfile1 -i inputfile2 -i inputfile3")
    print("\n\nOPTIONS:")
    print("\t-h or --help:\tprint help")
    print("\t-o or --ofile:\tOutput csv filename")
    print("\t-i or --ifile:\tInput file names")
    print("\t\t\tIf using multiple input files provide each of their name using \n\t\t\tflag -i separately. see "
          "COMMAND LINE.\n")
    print("\nDESCRIPTION:")
    print("\tScript read the outputs from the rnn.py script saved as pickle (.pkl) and convert them to txt file.")
    print("\tFile names are derived by replacing '.pkl' with '.txt'.")
    print("\tIf a single file name is given, the prediction sequences are shorted and printed to the text file based "
          "on their score.")
    print("\tIf more than one pickle file is given, the script will compare the sequence ranks in each of them and "
          "will provide an average ranks and scores.\n")


def main(argv):
    if argv:
        try:
            opts, argument = getopt.getopt(argv, "ho:i::", ["ofile=", "ifile="])
            # print opts, argument
        except getopt.GetoptError:
            print('\n', 'ERROR:', '\n\t', 'Check your arguments', '\n')
            show_help()
            sys.exit(2)
        
        ipf = []
        ofn = ''
        for opt, arg in opts:
            if opt == '-h':
                show_help()
                sys.exit()
            elif opt in ("-i", "--ifile"):
                ipf.append(arg)
            elif opt in ("-o", "--ofile"):
                ofn = arg
        
        print('\n', 'input files :', ipf)
        print('output file name:', ofn)
        
        return [ipf, ofn]


if __name__ == "__main__":
    args = main(sys.argv[1:])
    pkl_files = args[0]
    out_filename = args[1]
    
    # convert .pkl files to human readable form
    # arr_prediction_txt_file_names = output_in_text(input_array=pkl_files, work_dir=working_dir)
    
    arr_prediction_txt_file_names = [x.replace('.pkl', '.txt') for x in pkl_files]
    ranking_out = ranking_wrt_reference(arr_prediction_txt_file_names, ref=arr_prediction_txt_file_names[0])
    average_ranking_and_rnn_score(ranking_out, o_filename=out_filename, txt_file_names=arr_prediction_txt_file_names)
