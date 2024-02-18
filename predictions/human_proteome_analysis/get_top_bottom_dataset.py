"""
script extract top and bottom sequences from output file generated from rank_and_comparison.py and compare it with
experimental data and all_processed_last7_human_proteome.csv.

output:
    all_top.xlsx, containing top scored sequences from n iterations, with all the relevant data.
    all_bottom.xlsx, containing top scored sequences from n iterations, with all the relevant data.
    
    ** rows containing sequences that are in experimental data-set will be highlighted in pink.

"""
import os
import pandas as pd
import numpy as np

top_dataset_cutoff = 4.0
bottom_dataset_cutoff = 1.5

refseq_dir = "../../raw_data/ncbi_human_genome/gene2refseq/"
refseq_file_name = "gene2refseq_08_10_2018_prot_acc_human.txt"
# output file names
out_excel_name_top = 'top_dataset.xlsx'
out_excel_name_top_unique = 'top_dataset_unique.xlsx'
out_excel_name_bottom = 'bottom_dataset.xlsx'
out_excel_name_bottom_unique = 'top_dataset_unique.xlsx'

# input exp and predicted data
df_exp = pd.read_pickle('../../raw_data/full_dataset_199_exp.df')
# df_exp = pd.read_csv('../../raw_data/full_dataset_200_exp.txt', sep='\t')
df_av_rank_score = pd.read_csv('rank_comparisons_0_31.csv', sep='\t',
                               index_col='seq')  # rank_and_comparison.py output
df_processed_data = pd.read_csv('../../raw_data/ncbi_human_genome/proteome/all_processed_last7_human_proteome.txt',
                                sep='\t')

df_gene2refseq = pd.read_csv(os.path.join(refseq_dir, refseq_file_name), sep='\t')
df_mass_spec = pd.read_excel("../../raw_data/Copy of Harvey Lab Mass Spec(2).xlsx")


def write_xml(data_frame=None, count_gene_id_exp_data=None, out_xml=None):
    # print(count_gene_id_exp_data)
    xml = ['<table>']
    arr_unique_gene_id = []
    for a in data_frame.index:
        row = data_frame.iloc[a]
        # current_gene_id = row['count_gene_id']
        
        # try:
        #     next_gene_id = data_frame.iloc[a + 1]['count_gene_id']
        # except:
        #     next_gene_id = ''
        
        # if current_gene_id not in arr_unique_gene_id:
        #     if current_gene_id in count_gene_id_exp_data:
        #         # print(current_gene_id)
        #         xml.append('  <gene_count' + str(current_gene_id) + '_EXP>')
        #     elif current_gene_id not in count_gene_id_exp_data:
        #         xml.append('  <gene_count' + str(current_gene_id) + '_PRED>')
        #     arr_unique_gene_id.append(current_gene_id)
        
        xml.append('    <row' + str(a + 1) + '>')
        for i, col_name in enumerate(row.index):
            # xml.append('      <field name="{0}">{1}</field>'.format(col_name, row.iloc[i]))
            xml.append('      <{0}>{1}</{0}>'.format(col_name, row.iloc[i]))
        xml.append('    </row' + str(a + 1) + '>')
        
        # if next_gene_id != '':
        #     if current_gene_id != next_gene_id:
        #         if current_gene_id in count_gene_id_exp_data:
        #             xml.append('  </gene_count' + str(current_gene_id) + '_EXP>')
        #         if current_gene_id not in count_gene_id_exp_data:
        #             xml.append('    </gene_count' + str(current_gene_id) + '_PRED>')
        # else:
        #     if current_gene_id in count_gene_id_exp_data:
        #         xml.append('  </gene_count' + str(current_gene_id) + '_EXP>')
        #     if current_gene_id not in count_gene_id_exp_data:
        #         xml.append('  </gene_count' + str(current_gene_id) + '_PRED>')
    
    xml.append('</table>')
    
    with open(out_xml, 'w') as f:
        f.write('\n'.join(xml))


def write_excel(data_frame=None, out_excel=None, worksheet_name='Worksheet1', count_gene_id_exp_data=None):
    """
    module to write output into excel sheets. rows with sequences given in data_frame_exp will be colored.
    :param count_gene_id_exp_data:
    :param data_frame: input DataFrame
    :param worksheet_name: output worksheet name
    :param out_excel: output excel sheet name
    :param count_gene_id_exp_data: array containing list of count_gene_id positions in data_frame to be colored
    :return:
    """
    
    count_gene_id_exp_data = ['"%s"' % x for x in count_gene_id_exp_data]
    # print(count_gene_id_exp_data)
    writer = pd.ExcelWriter(out_excel)
    data_frame.to_excel(writer, worksheet_name, index=False)
    
    workbook = writer.book
    worksheet = writer.sheets[worksheet_name]
    number_rows = len(data_frame.index)
    
    format_all = workbook.add_format({'num_format': '0'})
    worksheet.set_column('A:Z', None, format_all)
    
    color_range1 = "A1:Z{}".format(number_rows + 1)
    format1 = workbook.add_format({
        'bg_color'  : '#FFC7CE',  # pink
        'font_color': '#000000',
        # 'num_format': '0.0000',
    })
    
    for entry in count_gene_id_exp_data:
        worksheet.conditional_format(color_range1, {
            'type'    : 'formula',
            'criteria': '=LEFT($A1, 20)=' + entry,
            'format'  : format1
        })
    
    format_five_dec = workbook.add_format({'num_format': '0.00000'})
    worksheet.set_column('E:G', None, format_five_dec)
    
    writer.save()


def mapping(data_frame=None, data_frame_av_rank_score=None, data_frame_processed_data=None, data_frame_exp=None,
            col=None, col_exp_response=None, above=None, below=None):
    arr_seq_top = []
    arr_seq_bottom = []
    if above is not None:
        arr_seq_top = [data_frame.index[x] for x in range(0, len(data_frame.index)) if
                       data_frame_av_rank_score['Av_RNN_score'][x] >= above]
        print('\nunique top sequences', len(arr_seq_top))
    if below is not None:
        arr_seq_bottom = [data_frame.index[x] for x in range(0, len(data_frame.index)) if
                          data_frame_av_rank_score['Av_RNN_score'][x] < below]
        print('unique bottom sequences', len(arr_seq_bottom))
    arr_seq = arr_seq_top + arr_seq_bottom
    
    temp_arr_all = []
    temp_score_all = []
    temp_exp_score = []
    
    for seq in arr_seq:
        seq_loc_average_score = list(np.where(data_frame_av_rank_score.index == seq))[0]
        average_score = round(data_frame_av_rank_score['Av_RNN_score'].iloc[seq_loc_average_score[0]], 5)
        seq_loc_processed = list(np.where(data_frame_processed_data[col] == seq)[0])
        for loc in seq_loc_processed:
            temp_score_all.append(average_score)
            temp_arr_all.append(list(data_frame_processed_data.iloc[loc])[1:])
            
            seq_loc_exp_score = list(np.where(data_frame_exp[col] == seq))[0]
            if len(seq_loc_exp_score) == 0:
                temp_exp_score.append('')
            else:
                temp_exp_score.append(round(data_frame_exp[col_exp_response].iloc[seq_loc_exp_score[0]], 5))
    
    header = list(data_frame_processed_data.columns[1:])
    temp_df_all = pd.DataFrame(temp_arr_all, columns=header, index=None)
    
    # Insert average RNN score and exp_score to the data-frame
    temp_df_all.insert(2, 'av_RNN_score', temp_score_all)
    temp_df_all.insert(2, 'exp_score', temp_exp_score)
    
    # copy pubmed id column
    df_pub_header = temp_df_all['pubmed_id']
    
    # delete pubmed id column
    temp_df_all.drop(labels=['pubmed_id'], axis=1, inplace=True)
    
    # count pubmed ids
    arr_pubmed_counts = []
    
    for pc in df_pub_header:
        if pd.isnull(pc):
            arr_pubmed_counts.append(0)
        else:
            if "," in pc:
                arr_pubmed_counts.append(len(pc.split(",")))
            elif pc.isdigit():
                arr_pubmed_counts.append(1)
    # insert pumbed id counts to the data-frame
    temp_df_all.insert(2, 'pubmed_id_counts', arr_pubmed_counts)
    
    arr_mass_spec = []
    for protein_gi in list(temp_df_all['protein_gi']):
        gi_loc_mass_spec = list(np.where(df_mass_spec['Accession'] == protein_gi)[0])
        if gi_loc_mass_spec:
            response = list(set(df_mass_spec['Tg Fold Induction'].iloc[gi_loc_mass_spec]))[0]
            arr_mass_spec.append(response)
        else:
            arr_mass_spec.append('')
    
    unique_tails = []
    count = 0
    arr_count = []
    exp_count = []
    for tail in list(temp_df_all['tail']):
        if tail in unique_tails:
            arr_count.append(count)
            if tail in list(data_frame_exp[col]):
                exp_count.append(count)
        else:
            count += 1
            arr_count.append(count)
            unique_tails.append(tail)
            if tail in list(data_frame_exp[col]):
                exp_count.append(count)
    
    temp_df_all.insert(0, 'count_gene_id', arr_count)
    temp_df_all.insert(6, 'mass_spec_response', arr_mass_spec)
    return [temp_df_all, list(set(exp_count))]


def unique(data_frame=None):
    temp_data_frame = data_frame.copy()
    for col in list(temp_data_frame.columns)[1:6]:
        # print(col)
        unique_entries = []
        for entry in temp_data_frame[col]:
            if entry not in unique_entries:
                unique_entries.append(entry)
            else:
                unique_entries.append('')
        temp_data_frame[col] = unique_entries
    return temp_data_frame


# map top and bottom average rnn scored sequences to data extracted from ncbi

df_mapped_data_all_top_all, list_exp_data_row_top_all = mapping(data_frame=df_av_rank_score,
                                                                data_frame_av_rank_score=df_av_rank_score,
                                                                data_frame_processed_data=df_processed_data,
                                                                data_frame_exp=df_exp,
                                                                col='tail',
                                                                col_exp_response='Tg Response - 8hr',
                                                                above=top_dataset_cutoff,
                                                                below=None)

mapped_data_all_top_all_unique = unique(data_frame=df_mapped_data_all_top_all)

df_mapped_data_all_bottom_all, list_exp_data_row_bottom_all = mapping(data_frame=df_av_rank_score,
                                                                      data_frame_av_rank_score=df_av_rank_score,
                                                                      data_frame_processed_data=df_processed_data,
                                                                      data_frame_exp=df_exp,
                                                                      col='tail',
                                                                      col_exp_response='Tg Response - 8hr',
                                                                      above=None,
                                                                      below=bottom_dataset_cutoff)

mapped_data_all_bottom_all_unique = unique(data_frame=df_mapped_data_all_bottom_all)

# write top and bottom sequences to excel sheet, color rows with experimental values.
write_excel(data_frame=df_mapped_data_all_top_all,
            out_excel=out_excel_name_top,
            worksheet_name='top',
            count_gene_id_exp_data=list_exp_data_row_top_all
            )

write_excel(data_frame=mapped_data_all_top_all_unique,
            out_excel=out_excel_name_top_unique,
            worksheet_name='top',
            count_gene_id_exp_data=list_exp_data_row_top_all
            )

write_excel(data_frame=df_mapped_data_all_bottom_all,
            out_excel=out_excel_name_bottom,
            worksheet_name='bottom',
            count_gene_id_exp_data=list_exp_data_row_bottom_all
            )

write_excel(data_frame=mapped_data_all_bottom_all_unique,
            out_excel=out_excel_name_bottom_unique,
            worksheet_name='bottom',
            count_gene_id_exp_data=list_exp_data_row_bottom_all
            )

# write xml files
write_xml(data_frame=mapped_data_all_top_all_unique,
          count_gene_id_exp_data=list_exp_data_row_top_all,
          out_xml=out_excel_name_top_unique.replace('.xlsx', '.xml'))

write_xml(data_frame=mapped_data_all_bottom_all_unique,
          count_gene_id_exp_data=list_exp_data_row_bottom_all,
          out_xml=out_excel_name_bottom_unique.replace('.xlsx', '.xml'))


# stats block #

def stats(data_frame=None, col1=None, count_gene_id_exp_data=(), sets='yes'):
    temp_all = []
    temp_pred = []
    tails_all = []
    tails_pred = []
    data_frame = data_frame.fillna('')
    for i in range(0, len(data_frame['count_gene_id'])):
        count_gene_id = list(data_frame['count_gene_id'])[i]
        if list(data_frame[col1])[i] != '':
            tails_all.append(count_gene_id)
            temp_all.append(list(data_frame[col1])[i])
            if count_gene_id not in count_gene_id_exp_data:
                tails_pred.append(count_gene_id)
                temp_pred.append(list(data_frame[col1])[i])
    
    # print(temp_all)
    if sets != 'no':
        print('counts', str(len(list(set(temp_all)))) + ' from ' + str(len(list(set(tails_all)))) + ' tails, ',
              str(len(list(set(temp_pred)))) + ' from ' + str(len(list(set(tails_pred)))) + ' tails')
    if sets == 'no':
        print('counts', str(len(temp_all)) + ' from ' + str(len(list(set(tails_all)))) + ' tails, ',
              str(len(temp_pred)) + ' from ' + str(len(list(set(tails_pred)))) + ' tails')


print('\nTotal tails:top')
stats(data_frame=mapped_data_all_top_all_unique, col1='count_gene_id',
      count_gene_id_exp_data=list_exp_data_row_top_all)

print('\nTotal refseq ids:top')
stats(data_frame=mapped_data_all_top_all_unique, col1='refseq_id', count_gene_id_exp_data=list_exp_data_row_top_all)
print('\nTotal Gene ids:top')
stats(data_frame=mapped_data_all_top_all_unique, col1='gene_id', count_gene_id_exp_data=list_exp_data_row_top_all)

print('\nSequences with Exp top')
stats(data_frame=mapped_data_all_top_all_unique, col1='exp_score', count_gene_id_exp_data=list_exp_data_row_top_all)
print('\nSequences with MassSpec top')
stats(data_frame=mapped_data_all_top_all_unique, col1='mass_spec_response',
      count_gene_id_exp_data=list_exp_data_row_top_all)

print('\nSequences signal peptides:top')
stats(data_frame=mapped_data_all_top_all_unique, col1='signal_peptide',
      count_gene_id_exp_data=list_exp_data_row_top_all, sets='no')

print('\nSequences tm proteins:top')
stats(data_frame=mapped_data_all_top_all_unique, col1='tm_region_info',
      count_gene_id_exp_data=list_exp_data_row_top_all, sets='no')

# bottom stats
print('\nTotal tails:bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='count_gene_id',
      count_gene_id_exp_data=list_exp_data_row_bottom_all)

print('\nTotal refseq ids:bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='refseq_id',
      count_gene_id_exp_data=list_exp_data_row_bottom_all)
print('\nTotal Gene ids:bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='gene_id',
      count_gene_id_exp_data=list_exp_data_row_bottom_all)

print('\nSequences with Exp :bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='exp_score',
      count_gene_id_exp_data=list_exp_data_row_bottom_all)
print('\nSequences with MassSpec :bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='mass_spec_response',
      count_gene_id_exp_data=list_exp_data_row_bottom_all)

print('\nSequences signal peptides:bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='signal_peptide',
      count_gene_id_exp_data=list_exp_data_row_bottom_all, sets='no')

print('\nSequences tm proteins:bottom')
stats(data_frame=mapped_data_all_bottom_all_unique, col1='tm_region_info',
      count_gene_id_exp_data=list_exp_data_row_bottom_all, sets='no')


def stats_signal_tm(data_frame=None, col1=None, count_gene_id_exp_data=()):
    temp_all = []
    temp_pred = []
    tails_all = []
    tails_pred = []
    data_frame = data_frame.fillna('')
    for i in range(0, len(data_frame['count_gene_id'])):
        count_gene_id = list(data_frame['count_gene_id'])[i]
        if list(data_frame[col1])[i] != '':
            tails_all.append(count_gene_id)
            temp_all.append(list(data_frame[col1])[i])
            if count_gene_id not in count_gene_id_exp_data:
                tails_pred.append(count_gene_id)
                temp_pred.append(list(data_frame[col1])[i])
    
    pub_count = 0
    uni_count = 0
    for a in temp_pred:
        if 'PubMed' in a:
            pub_count += 1
        elif 'UniProtKB' in a:
            uni_count += 1
    print('pubmed ref', pub_count, 'uniprot ref', uni_count)


print('\n**********\nSequences signal peptides ref:top')
stats_signal_tm(data_frame=mapped_data_all_top_all_unique, col1='signal_peptide',
                count_gene_id_exp_data=list_exp_data_row_top_all)

print('\nSequences tm proteins ref :top')
stats_signal_tm(data_frame=mapped_data_all_top_all_unique, col1='tm_region_info',
                count_gene_id_exp_data=list_exp_data_row_top_all)

print('\nSequences signal peptides ref:bottom')
stats_signal_tm(data_frame=mapped_data_all_bottom_all_unique, col1='signal_peptide',
                count_gene_id_exp_data=list_exp_data_row_bottom_all)

print('\nSequences tm proteins ref:bottom ')
stats_signal_tm(data_frame=mapped_data_all_bottom_all_unique, col1='tm_region_info',
                count_gene_id_exp_data=list_exp_data_row_bottom_all)
