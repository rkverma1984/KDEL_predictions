#!/usr/bin/env python
"""
Calculate the frequency of all amino acids shown up in each position.

We count each position amino acid frequency, then sum up the total of 7 positions frequency. List top 10 highest frequent sequences
Input file: A .txt format file that only contained sequences.

To run this code:  python calculate_occurrence_frequency.py --sequences_txt 'only_seq_cutoff2_first_set.txt'

Author: Bing Xie
"""
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--sequences_txt', default = 'only_seq_cutoff2_first_set.txt', help ='A txt file only contains sequences that are needed to be compare')
args = parser.parse_args()
input_file = args.sequences_txt

def read_sequences_from_txt(seq_file):
  f = open(seq_file,'r')
  lines=f.read().split('\n')
  f.close()
  sequences=[]
  for line in lines:
    if not line=='':
      sequences.append(line)
  return sequences
def extract_sequences(data_file):
  f = open(data_file,'r')
  lines=f.read().split('\n')
  f.close()
  sequences =[]
  for line in lines:
    if line !='' and line!='>':
      sequences.append(line)
  return sequences

def _count_frequency(sequences,position):
  keys = [s[position] for s in sequences]
  uniq_keys = list(set(keys))
  dictionary=dict.fromkeys(uniq_keys, 0)#initialize dictionary
  for s in sequences:
    dictionary[s[position]]+=1
  dictionary = sorted(dictionary.items(), key=lambda item:item[1])[::-1]## sort frequence from high to low
  return dictionary  ## return in a particular position all possibilities of AA

def score_last_7_amino_acids(sequences):
  position_0=_count_frequency(sequences,0)
  position_1=_count_frequency(sequences,1)
  position_2=_count_frequency(sequences,2)
  position_3=_count_frequency(sequences,3)
  position_4=_count_frequency(sequences,4)
  position_5=_count_frequency(sequences,5)
  position_6=_count_frequency(sequences,6)
  min_length = min(len(position_0), len(position_1),len(position_2),len(position_3),len(position_4),len(position_5),len(position_6))
  print ('***** TOP Amino Acid in each position ****')
  for i in range(min_length):
    print (position_0[i], position_1[i],position_2[i],position_3[i],position_4[i],position_5[i],position_6[i])

  # filter as top 5,  or can try using cutoff >10 to reduce calculation
  all_possibilies={}
  for a in position_0[:5]:
    for b in position_1[:5]:
      for c in position_2[:5]:
        for d in position_3[:5]:
          for e in position_4[:5]:
            for f in position_5[:5]:
              for g in position_6[:5]:
                all_possibilies[a[0]+b[0]+c[0]+d[0]+e[0]+f[0]+g[0]] = a[1]+b[1]+c[1]+d[1]+e[1]+f[1]+g[1]

  sorted_sequences = sorted(all_possibilies.items(), key=lambda item:item[1])[::-1]## sort frequence from high to low
  del all_possibilies
  output=[]
  for s in sorted_sequences:
    #if s[1]>900:
      #print s
      output.append(s)
  return output

sequences  = extract_sequences(input_file)
output = score_last_7_amino_acids(sequences)
print ('****** all combination sequences and total count number  ******')
print (output[:10])

