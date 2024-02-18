"""
Script process the excel sheet provided and stores the data in a DataFrame.

Input files:
    Input file is an excel sheet.
        "SERCaMP_Library_Gene_List_with_mass_spec.xlsx" in present case

Output:
    processed output is stored in file defined in test_rnn_parameters.py by output_df variable.

"""

import os
import sys
import getopt
import pandas as pd
import numpy as np
from Bio import SeqIO
import requests
from io import StringIO
from IPython.display import display

user = os.path.expanduser('~')

working_dir = os.getcwd()
sys.path.append(working_dir)


def get_fasta(code):
    """
    module to pull sequence from the UNIPROT database

    """
    if len(code.split(',')) == 2:
        code = code.split(',')[0]
    params = {"query": code, "format": "fasta"}
    response = requests.get("http://www.uniprot.org/uniprot/", params)
    
    for r in SeqIO.parse(StringIO(response.text), "fasta"):
        return str(r.seq)


def read_excel_and_save_df(meta_dir, input_excel, output_df, work_dir):
    """
    read the excel file provided and store the information in a DataFrame.
   
    """
    
    df = pd.read_excel(os.path.join(meta_dir, input_excel), header=1, parse_cols='B:Q')
    
    df = df.dropna(axis=0, how='all')
    
    s = pd.Series([''.join(df.iloc[i, 4:11]) for i in range(df.shape[0])])
    df['tail'] = s.values
    
    n_pos = np.where(df.index.values == 'SERCaMP (+)')[0][0]
    label = [0] * n_pos + [1] * (df.shape[0] - n_pos)
    df['Label'] = label
    
    df = df.drop(['Number', 'Amino Acids'] + ['Unnamed: %d' % (i + 5) for i in range(7)], axis=1)
    
    df = df.set_index('UniProt')
    display(df)
    
    # the sequence entry in the 34th line is wrong and is modified accordingly
    # df['tail'][34] = 'PQPRDEL'
    
    # get full sequence from UNIPROT
    full_sequence = [get_fasta(code) for code in df.index]
    
    # verify all the sequences are correct. last 7 residues from the UNIPROT sequence are matched with excel entry.
    for i, code in enumerate(df.index):
        # print out UNIPROT id if the last 7 residues does not match
        if full_sequence[i][-7:] != df['tail'][i]:
            print(i, code)
    
    """
    update the full sequence in the DataFrame
    """
    df['fullSeq'] = full_sequence
    
    # save processed DataFrame in pickle format
    df.to_pickle(os.path.join(work_dir, "%s_%s" % ("complete", output_df)))
    df.to_csv(os.path.join(work_dir, "%s_%s" % ("complete", output_df.replace('.df', '.txt'))), sep='\t',
              header=df.columns, index=False)
    
    # only write the relevant two columns (tail sequence and Tg response at 8hrs)
    my_df = df[['tail', 'Tg Response - 8hr']]
    my_df = my_df.set_index([np.arange(len(my_df))])
    decimals = pd.Series([5], index=['Tg Response - 8hr'])
    my_df = my_df.round(decimals)
    my_df.to_pickle(os.path.join(work_dir, "%s" % output_df))
    my_df.to_csv(os.path.join(work_dir, output_df.replace('.df', '.txt')), sep='\t',
                 header=['tail', 'Tg Response - 8hr'], index=False, float_format='%.5f')


def show_help():
    print("Usage:\n")
    print("ipython input_preparation.py -p parameters.py")


def main(argv):
    if argv:
        try:
            opts, argument = getopt.getopt(argv, "hp:j::", ["pfile="])
        except getopt.GetoptError:
            print('\nERROR:\n    Check your arguments', '\n')
            show_help()
            sys.exit(2)
        
        pf = ''
        jn = ''
        for opt, arg in opts:
            if opt == '-h':
                show_help()
                sys.exit()
            elif opt in ("-p", "--pfile"):
                pf = arg
        
        print('input parameters file :', pf)
        
        return [pf, jn]


if __name__ == "__main__":
    args = main(sys.argv[1:])
    parms = __import__(args[0].replace('.py', ''))
    
    read_excel_and_save_df(parms.meta_dir, parms.input_excel, parms.output_df, working_dir)
