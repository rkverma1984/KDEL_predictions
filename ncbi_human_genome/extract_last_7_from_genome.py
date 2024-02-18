"""
script to extract last seven residues from the human proteome
Input: human proteome in fasta format
Output: dataframe with 'fasta_header' and 'tail' as columns

"""
import glob
import pandas as pd

genome_files_in_dir = glob.glob('*.faa')
print("Available genome files:")
for fn in genome_files_in_dir:
    print("\t", fn)
genome_fasta_file = input("\nEnter the genome file name: ")
genome_f_data = open(genome_fasta_file, 'r').read()

out_prefix = 'last_7_residues'
out_file_name = '%s_%s' % (out_prefix, genome_fasta_file)


def separate_each_fasta_file(data):
    return [x for x in data.split(">") if x != '']


individual_fasta_sequences = separate_each_fasta_file(genome_f_data)


def get_sequences_dictionary(data):
    seq_dict = {}
    temp_data = [x for x in data.split("\n") if x != '']
    full_seq = "".join(temp_data[1:])
    last_7_residues = full_seq[-7:]
    
    seq_dict.update({temp_data[0]: last_7_residues})
    
    return seq_dict


temp_array = []
for block in individual_fasta_sequences:
    for key, value in dict.items(get_sequences_dictionary(block)):
        temp_array.append([key.replace('>', ''), value])

df = pd.DataFrame(temp_array, columns=['fasta_header', 'tail'])
df = df.set_index('fasta_header')
df.to_csv(out_file_name, sep='\t')
df.to_pickle(out_file_name.replace('.faa', '.df'))
print(df)
