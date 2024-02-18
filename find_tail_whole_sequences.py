import os
import pandas as pd
workdir='/home/bxie/PycharmProjects/benchmark/building_new_model/199_training_database'
sequences_txt =os.path.join(workdir, 'raw_data','first_set_tail_response.txt')
excel_file =os.path.join(workdir,'rnn_modular_package/199data_training_batch1_prediction_P13_human_proteins','all_top_100_unique.xlsx')

# read input tail sequences file #
def input_tail(tail_txt):
  f = open(tail_txt,'r')
  lines=f.read().split('\n')
  f.close()
  seqs = []
  for line in lines:
    if not line=='' and not 'tail' in line:
      seqs.append(line.split()[0])
  return seqs

# search tails from excel file to get uniprot ID # 
def search_xlsx(xlsx_file, seqs):
  data =pd.read_excel(xlsx_file)
  IDs = []
  for s in seqs:
    for i in range(len(data['tail'])):
      if s ==data['tail'][i]:
        IDs.append(data['uniprot_id'][i])
  IDs = list(set(IDs)) 
  return IDs

# Uniprot ID was not included in the first batch, search from downloaded database #
def search_in_downloaded_uniprot_database(seqs):
  f = open('/home/bxie/PycharmProjects/benchmark/uniprot-human-review.fasta','r')
  lines=f.read().split('\n')
  f.close()
  mark = []
  for i in range(len(lines)):
    if '>' in lines[i]:
      mark.append(i) ### starting index
  dictionary={} #value is uniprot id, key is sequence
  for i in range(len(mark)-1):
    dictionary[(''.join(x for x in lines[mark[i]:mark[i+1]] if not '>' in x))[-7:]]=lines[mark[i]].split('|')[1]
  dictionary[(''.join(x for x in lines[mark[-1]:] if not '>' in x and x!=''))[-7:]] = lines[mark[-1]].split('|')[1]
  full_sequences_dictionary={}
  for i in range(len(mark)-1):
    full_sequences_dictionary[lines[mark[i]].split('|')[1]] = ''.join(x for x in lines[mark[i]:mark[i+1]] if not '>' in x)
  full_sequences_dictionary[lines[mark[-1]].split('|')[1]] = ''.join(x for x in lines[mark[-1]:] if not '>' in x)
  IDs = []
  for key in dictionary.keys():
    if key in seqs:
      IDs.append(dictionary[key])
  f = open('relative_sequences.txt','w')
  full_sequences =[]
  for i in IDs:
    f.write('>'+i+'\n')
    f.write(full_sequences_dictionary[i]+'\n')
    full_sequences.append(full_sequences_dictionary[i])
  f.close()
  return IDs


seqs = input_tail(sequences_txt)
uniprot_ids = search_in_downloaded_uniprot_database(seqs)
#uniprot_ids = search_xlsx(excel_file, seqs)

