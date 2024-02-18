import pandas as pd
import numpy as np
import os,sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as stats
#import argparse
#parser=argparse.ArgumentParser()
#parser.add_argument('--input_experimental_txt', default='second_set_SERCaMP_Response_Data.txt')
#parser.add_argument('--run_type', default='ranking')
#args=parser.parse_args()

#input_experimental_txt = args.input_experimental_txt
#run_type=args.run_type
#home_dir ='/home/bxie/PycharmProjects/KDEL/benchmark/building_new_model/95_training_database_predict_89batch'
#experimental_data = os.path.join(home_dir,'raw_data', input_experimental_txt)

def read_excel(xslx_file): #extract sequences from positive and negtive xlsx files
  human_seq=[]
  arti_seq = []
  d = pd.read_excel(xslx_file)
  if 'HIT' in d.keys():
    for x in d['HIT']:
      human_seq.append(x)
  if 'QUERY' in d.keys():
    for y in d['QUERY']:
      arti_seq.append(y)
  else:
    print ('no human or aritificial equence found in this excel')
  return human_seq, arti_seq

#if os.path.isfile(home_dir+'/negatives_39_sequences.xlsx') and os.path.isfile(home_dir+'/positives_29_sequences.xlsx'):
#  neg_sequences_h, neg_sequences_a = read_excel(home_dir+'/negatives_39_sequences.xlsx')
#  print ('#--- read from negative_39_sequences.xlsx, total%4d sequences'%len(neg_sequences_h))
#  pos_sequences_h, pos_sequences_a = read_excel(home_dir+'/positives_29_sequences.xlsx')
#  print ('#--- read from poitives_29_sequences.xlsx, total%4d sequences'%len(pos_sequences_h))
#else:
#  print ('CANNOT FIND REFERENCE EXCEL FILES! TERMINATE PROCESS')
#  sys.exit()
def _IO(input_file):
  f = open(input_file,'r')
  line=f.read().split('\n')
  f.close()
  return line

def read_experimental_txt(experimental_file):
  experimental_dictionary={}
  line=_IO(experimental_file)
  for line in line[1:]:
    if line!='':
      experimental_dictionary[line.split()[0]] = float(line.split()[1])
  experimental_dictionary = sorted(experimental_dictionary.items(), key=lambda x: x[1])[::-1]
  return experimental_dictionary

def read_average_prediction_reults(prediction_file,comparision_type):
  d = pd.read_cv(prediction_file,delimiter='\t')
  sequences = d['seq']
  if comparision_type =='score':
    RNN = d['Av_RNN_core']
  elif comparision_type =='ranking':
    RNN = d['Av_RNN_rank']
  else:
   print ('Unrecognize compariion type')
   sys.exit()
  prediction ={}
  for i in range(len(sequences)):
    prediction[sequences[i]] = RNN[i]
  prediction =sorted(prediction.items(), key=lambda x: x[1])
  return prediction

def read_seperate_prediction_results(prediction_file_comman_name, model_inx, comparision_type,  main_dir):
  mintype='LSTM'
  epoch=3000
  single_prediction_file=os.path.join(main_dir, 'rnn_modular_package', 'first32',prediction_file_comman_name+'_'+str(model_inx)+'_'+mintype+'_'+str(epoch)+'.txt')
  if os.path.isfile(single_prediction_file):
    lines=_IO(single_prediction_file)
    sequences = []
    RNN = []
    for i in range(len(lines)):
      if lines[i]!='':
        sequences.append(lines[i].split()[0])
        if comparision_type =='score':
          RNN.append(float(lines[i].split()[1]))
        elif comparision_type =='ranking':
          RNN.append(i)
        else:
          print ('Unrecognize compariion type')
          sys.exit()
    prediction = {}
    for i in range(len(sequences)):
      prediction[sequences[i]] = RNN[i]
    prediction =sorted(prediction.items(), key=lambda x: x[1])
  else:
    print ('Check individual prediction file name. Format hould be all_predict_checkpoints_$INDEX_LSTM_3000.txt')
    print ('Current input file name i '+single_prediction_file)
    sys.exit()
  return prediction

def matching_value(experimental_data, predictionary_data, comparision_type):
  #core type: experimental_data =experimental_tg, predictionary_data= predictionary_score
  #ranking type: experimental_data =experimental_tg, predictionary_data= predictionary_rank 
  #output: prediction_value is a list type, experimental_values is a list type, sequences is a list type, ordered by experimental_data
  prediction_values=[]
  experimental_values=[]
  sequences=[]
  if comparision_type =='score':
    for i in range(len(experimental_data)):
      for j in range(len(predictionary_data)):
        if experimental_data[i][0] == predictionary_data[j][0]:
          prediction_values.append(predictionary_data[j][1])
          experimental_values.append(experimental_data[i][1])
          sequences.append(experimental_data[i][0])
  elif comparision_type =='ranking':
    for i in range(len(experimental_data)):
      for j in range(len(predictionary_data)):
        if experimental_data[i][0] == predictionary_data[j][0]:
          prediction_values.append(j)
          experimental_values.append(i)
          sequences.append(experimental_data[i][0])
  return prediction_values, experimental_values, sequences

def plot_value(pred_rnn, tg_value, sequences, comparision_type, std_pred, std_exp, neg_sequences_h, neg_sequences_a,pos_sequences_h, pos_sequences_a ):
  #tg_value hould be a list only with tg value
  #pred_rnn hould be a list only with rnn score
  human_neg=0 #count number
  human_pos=0 #counts number
  artificial_pos =0 #counts number 
  artificial_neg = 0#count number
  other =0    #count number
  nh_index = None
  na_index = None
  ph_index = None
  pa_index = None
  o_index = None
  plt.figure() 
  for i in range(len(tg_value)):
    if sequences[i] in neg_sequences_h:
      plt.errorbar(tg_value[i], pred_rnn[i], xerr=std_exp[i],yerr=std_pred[i],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='red')
      human_neg+=1
      nh_index = i
    elif sequences[i] in neg_sequences_a:
      plt.errorbar(tg_value[i], pred_rnn[i],xerr=std_exp[i],yerr=std_pred[i],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='pink')
      artificial_neg+=1
      na_index= i
    elif sequences[i] in pos_sequences_h:
      plt.errorbar(tg_value[i],pred_rnn[i],xerr=std_exp[i],yerr=std_pred[i],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='blue')
      human_pos+=1
      ph_index = i
    elif sequences[i] in pos_sequences_a:
      plt.errorbar(tg_value[i],pred_rnn[i],xerr=std_exp[i],yerr=std_pred[i],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='cyan')
      artificial_pos+=1
      pa_index = i
    else:
      plt.errorbar(tg_value[i],pred_rnn[i],xerr=std_exp[i],yerr=std_pred[i],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='black')
      other+=1
      o_index =i
  if pa_index:
    plt.errorbar(tg_value[pa_index], pred_rnn[pa_index],xerr=std_exp[pa_index],yerr=std_pred[pa_index],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5,color='cyan',label='positive artificial')
  if ph_index:
    plt.errorbar(tg_value[ph_index], pred_rnn[ph_index],xerr=std_exp[ph_index],yerr=std_pred[ph_index],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='blue', label='positive human')
  if nh_index:
    plt.errorbar(tg_value[nh_index], pred_rnn[nh_index],xerr=std_exp[nh_index],yerr=std_pred[nh_index],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='red', label='negtative human')
  if o_index:
    plt.errorbar(tg_value[o_index], pred_rnn[o_index],xerr=std_exp[o_index],yerr=std_pred[o_index],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='black',label='other')
  if na_index:
    plt.errorbar(tg_value[na_index], pred_rnn[na_index],xerr=std_exp[na_index],yerr=std_pred[na_index],fmt='o',markersize= 6, capthick=2,capsize=2,alpha=0.5, color='pink', label='negtative artificial')
  print ('#--- total number of human negative i ', human_neg,'-------*')
  print ('#--- total number of human poitive is ', human_pos,'-------*')
  print ('#--- total number of artificial negative i ', artificial_neg,'---*')
  print ('#--- total number of artificial poitive is ', artificial_pos,'---*')
  print ('#--- total number of other is ', other,'---------------*')
  pearson, p_value=stats.pearsonr(tg_value,pred_rnn )
  coeffs=np.polyfit(tg_value, pred_rnn, 1)
  corre_coeff = np.corrcoef([tg_value, pred_rnn])[0][-1]
  intercept = coeffs[-1]
  slope=coeffs[-2]
  yl = slope*np.array([min(tg_value)-1, max(tg_value)+1])+intercept
  line_text = "y = %0.2fx"%slope
  if intercept >=0:
    line_text+="+"
  else:
    line_text+="-"
  line_text += "%0.2f"%np.abs(intercept)
  print ('********** Pearson R is '+str(pearson)+' ******************************')
  print ('********** correlation i '+line_text+' *******************************')
  plt.legend(loc='lower right')
  plt.plot([min(pred_rnn)-2,max(pred_rnn)+2], yl,'--',lw=1,color='black')
  if comparision_type =='ranking':
    plt.xlabel('Experimental Ranking')
    plt.ylabel('Prediction Ranking')
    plt.xticks([0,9,19,29,39,49,59,69,79,88],['1','10','20','30','40','50','60','70','80','89'])
    plt.yticks([0,9,19,29,39,49,59,69,79,88],['1','10','20','30','40','50','60','70','80','89'])
  if comparision_type =='score':
    plt.xlabel('Experimental Tg Reponse')
    plt.ylabel('RNN Average Score')
    plt.xlim(0.,7.)
    plt.ylim(0.5,8.0)
  plt.legend(loc='lower right')
  plt.tight_layout()
  plt.savefig('combined_'+comparision_type+'_comparison_based_RNN.png')
  plt.show()
  plt.clf()
  return
#------------------- running section -------------------------------#
#run_type='score' or 'ranking', define in front of this code
#calc_all =False
#experimental_tg = read_experimental_txt(experimental_data)
#n_models = 32

#if calc_all:
#  prediction_data = os.path.join(home_dir, 'rnn_modular_package', 'combine_32_model.csv')
#  predictionary = read_average_prediction_reults(prediction_data, run_type)
#  prediction_values, experimental_values, sequences = matching_value(experimental_tg, predictionary,run_type)
#  plot_value(prediction_values,experimental_values, sequences, run_type,[0.0 for i in range(len(prediction_values))],[0.0 for i in range(len(experimental_values))])
#else:
#  pred_value =[]
#  exp_value =[]

#  for i in range(n_models):#32 model
#    predictionary = read_seperate_prediction_results('all_predict_checkpoints', i, run_type)
#    prediction_values, experimental_values, sequences = matching_value(experimental_tg, predictionary,run_type)
#    pred_value.append(prediction_values)
#    exp_value.append(experimental_values)
#  ave_pred_value = []
#  std_pred_value =[]
#  ave_exp_value=[]
#  std_exp_value=[]
#  for j in range(len(pred_value[0])):
#    ave_pred_value.append(np.mean([pred_value[x][j] for x in range(n_models)]))
#    std_pred_value.append(np.std([pred_value[x][j] for x in range(n_models)]))
#    ave_exp_value.append(np.mean([exp_value[x][j] for x in range(n_models)]))
#    std_exp_value.append(np.std([exp_value[x][j] for x in range(n_models)]))
#  sem_exp_value=[s/np.sqrt(n_models) for s in std_exp_value]
#  sem_pred_value=[s/np.sqrt(n_models) for s in std_pred_value]
  #plot_value(ave_pred_value,ave_exp_value, sequences, run_type,std_pred_value,std_exp_value)
#  plot_value(ave_pred_value,ave_exp_value, sequences, run_type,sem_pred_value,sem_exp_value)  #SEM error bar 
  #f = open('32_models_ave_experiment_value.txt','w')
  #for a in range(len(ave_exp_value)):
  #  f.write("%6.3f"%(ave_exp_value[a])+'  '+"%6.3f"%(std_exp_value[a])+'\n')
#---------------------------------------------------------------------------------------------#
  #f.close()
 
#---------------------------------------------------------------------------------------------#
