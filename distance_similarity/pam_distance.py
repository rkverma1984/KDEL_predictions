#!/usr/bin/env python
"""
Calculate sequences similarity through PAM250 scoring matrix

We compare each sequence to the reference sequence and get PAM250 score. LIST the sequences that have top 3 lowest PAM250 scores.
Input file: reference sequence, a string format
            a txt file that contained both sequences and corresponding experimental tg values .

To run this code: python pam_distance.py --ref_sequence 'KKKKDEL' --sequences_txt 'cutoff2_full_dataset_199_exp.txt'
Author: Bing Xie
"""

import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--ref_sequence', default = 'KKKKDEL', help ='The sequence/sequences that can be used as the reference')
parser.add_argument('--sequences_txt', default = 'cutoff2_full_dataset_199_exp.txt', help ='A txt file contains sequences that are needed to be compare')
args = parser.parse_args()

reference_seq = args.ref_sequence
test_seq_file =args.sequences_txt 


from Bio.SubsMat.MatrixInfo import *
def score_match(pair, matrix):
  if pair in matrix.keys():
    return matrix[pair]
  else:
    return matrix[(tuple(reversed(pair)))]

def score_pairwise(seq1, seq2, matrix): # no gap penalty concerned
  score =0
  assert len(seq1)==len(seq2)
  for i in range(len(seq1)):
    pair = (seq1[i],seq2[i])
    score+=score_match(pair, matrix)
  return score

def read_sequences_from_txt(seq_file):
  f = open(seq_file,'r')
  lines=f.read().split('\n')
  f.close()
  sequences={}
  for line in lines[1:]:
    if not line=='' and line !='>':
      if float(line.split()[1])>2: ## Use threshold 2 as cutoff, any sequence tg value above 2, considered as positive ones
        sequences[line.split()[0]]= float(line.split()[1])
  return sequences


human_seqs=read_sequences_from_txt(test_seq_file)
pam_options = [250] ######## Using PAM250 matrix to calculate sequences similarity
for num in pam_options:
  dictionary ={}
  for h in human_seqs.keys():
    score = score_pairwise(reference_seq,h,locals()['pam'+str(num)])
    dictionary[reference_seq+'_'+h] = score
  dictionary = sorted(dictionary.items(), key=lambda item:item[1])
  values = sorted(list(set([d[1] for d in dictionary])))
  for score_rank in range(3): #### Select top 3 lowest score
    lowest_score=values[score_rank]
    length=0
    for top_index in range(len(dictionary)):
      if dictionary[top_index][1]==lowest_score:
        print ('\n****** The Furthest Distance between  experimental sequences and human predictions is ******')
        print ('****** '+dictionary[top_index][0].split('_')[0]+' and '+dictionary[top_index][0].split('_')[1]+'******')
        print ('Tg score: '+str(human_seqs[dictionary[top_index][0].split('_')[1]]), dictionary[top_index][0].split('_')[1])
        print ('****** The Pam'+str(num)+' score is '+str(dictionary[top_index][1])+' ******')
        length+=1
    print ('total '+str(length)+' sequences have the '+str(score_rank+1)+'th lowest score\n')
