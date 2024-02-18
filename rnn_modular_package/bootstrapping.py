import os
import sys
import numpy as np
import pandas as pd

sys.path.append(os.getcwd())

try:
    from inputs_outputs import *
except ImportError:
    from .inputs_outputs import *

# read data #
proteome_f_data = open(os.path.join(proteome_dir, proteome_fasta_file_name), 'r').read()
df_gene2refseq = pd.read_csv(os.path.join(refseq_dir, refseq_file_name), sep='\t')
pubmed_file_name = pd.read_csv(os.path.join(gene2pubmed_dir, pubmed_file_name), sep='\t')
df_gene2go = pd.read_csv(os.path.join(gene2go_dir, go_file_name), sep='\t')
df_geneinfo = pd.read_csv(os.path.join(geneinfo_dir, gene_info_file_name), sep='\t')


# functions for data processing #

def download_updated(out_f_name=None, seq_to_download=None):
    # batch entrez url link
    joined_entries = ','.join(seq_to_download)
    print(joined_entries)
    batch_entrez_url_link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta" \
                            "&retmode=text&id="
    # os.system('wget -O %s "%s%s"' % (os.path.join(proteome_dir, out_f_name), batch_entrez_url_link, joined_entries))
    os.system('curl -o {} "{}{}"'.format(os.path.join(proteome_dir, out_f_name), batch_entrez_url_link, joined_entries))


def get_embl_table(refseq_id=None):
    embl_link_prefix = "https://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query="
    embl_link_suffix = "&sort=score&columns=id%2Creviewed%2Cprotein%20names%2Corganism%2Clength%2Cfeature(SIGNAL)" \
                       "%2Csequence%2Ccitation%2Cfeature(" \
                       "TRANSMEMBRANE)%2Ccomment(SUBCELLULAR%20LOCATION)%2Cfeature(" \
                       "TOPOLOGICAL%20DOMAIN)&format=tab&force=yes"
    
    if os.path.exists(uniprot_dir):
        pass
    else:
        os.makedirs(uniprot_dir)
    
    # os.system('wget -O %s/%s "%s%s%s"' % (uniprot_dir, refseq_id + '.txt', embl_link_prefix, refseq_id,
    #                                       embl_link_suffix))
    os.system('curl -s -o %s/%s "%s%s%s"' % (uniprot_dir, refseq_id + '.txt', embl_link_prefix, refseq_id,
                                             embl_link_suffix))
    try:
        embl_df = pd.read_csv(os.path.join(uniprot_dir, refseq_id + '.txt'), sep='\t')
    except pd.io.common.EmptyDataError:
        embl_df = pd.DataFrame()
        # os.system('rm -fr '+ refseq_id + '.txt')
    return embl_df


def separate_each_fasta_file(data):
    return [x for x in data.split(">") if x != '']


def get_sequences_dictionary_all(data):
    seq_dict = {}
    temp_data = [x for x in data.split("\n") if x != '']
    full_seq = "".join(temp_data[1:])
    last_7_residues = full_seq[-7:]
    
    seq_dict.update({temp_data[0]: last_7_residues})
    return seq_dict


def get_fasta(seq_block=None):
    """
    seq_block = sequence block
    """
    return "".join([x for x in seq_block.split("\n")[1:] if x != ''])


def process_fasta_file(f_handle=None, refseq_data_frame=None, pubmed_data_frame=None, go_data_frame=None,
                       gene_info_data_frame=None):
    individual_fasta_sequences = separate_each_fasta_file(f_handle)
    
    temp_array = []
    temp_updated = []
    temp_removed = []
    
    i = 0
    for block in individual_fasta_sequences:
        i += 1
        fasta_seq = get_fasta(seq_block=block)
        for key, value in dict.items(get_sequences_dictionary_all(block)):
            refseq_id = key.split(' ')[0]
            print(i, refseq_id)
            if refseq_id in list(refseq_data_frame['protein_accession.version']):
                gene_id_loc = (np.where(refseq_data_frame['protein_accession.version'] == refseq_id))[0]
                gene_id = list(set(refseq_data_frame['GeneID'].iloc[gene_id_loc]))[0]
                protein_gi = list(set(refseq_data_frame['protein_gi'].iloc[gene_id_loc]))[0]
                
                gene_status = []
                for gbl in gene_id_loc:
                    gene_status.append(refseq_data_frame['status'].iloc[gbl])
                
                go_loc = (np.where(go_data_frame['GeneID'] == gene_id))[0]
                arr_localization = []
                arr_function = []
                for gl in go_loc:
                    category = go_data_frame['Category'].iloc[gl]
                    go_term = go_data_frame['GO_term'].iloc[gl]
                    if 'Component' in category:
                        arr_localization.append(go_term)
                    if 'Function' in category:
                        arr_function.append(go_term)
                
                type_loc = (np.where(gene_info_data_frame['GeneID'] == gene_id))[0]
                gene_type = []
                
                for tl in type_loc:
                    gene_type.append(gene_info_data_frame['type_of_gene'].iloc[tl])
                
                pubmed_loc = (np.where(pubmed_data_frame['GeneID'] == gene_id))[0]
                pubmed_ids = []
                
                for pid in pubmed_loc:
                    pubmed_ids.append(pubmed_data_frame['PubMed_ID'].iloc[pid])
                
                # print(i, refseq_id, gene_id, list(set(gene_status)), list(set(gene_type)))
                
                epub_count = 0
                uniprot_id = ''
                tm_protein_type = ''
                signal_peptide = ''
                tm_region_info = ''
                try:
                    embl_data_frame = get_embl_table(refseq_id)
                    if embl_data_frame.empty:
                        refseq_id_new = (refseq_id.split('.'))[0]
                        embl_data_frame = get_embl_table(refseq_id_new)
                        print(i, refseq_id, refseq_id_new)
                except:
                    pass
                
                if not embl_data_frame.empty:
                    for epub in range(0, len(list(embl_data_frame['PubMed ID']))):
                        if not pd.isnull(embl_data_frame['PubMed ID'][epub]):
                            line = str(list(embl_data_frame['PubMed ID'])[epub]).rstrip(';')
                            if ';' in line:
                                embl_pubmed_ids = [int(line) for line in line.split(';') if line != ' ']
                            else:
                                embl_pubmed_ids = [line]
                            intersection_list = list(set(embl_pubmed_ids).intersection(set(pubmed_ids)))
                            if len(intersection_list) >= epub_count and fasta_seq == str(
                                    embl_data_frame['Sequence'][epub]):
                                uniprot_id = list(embl_data_frame['Entry'])[epub]
                                epub_count = len(intersection_list)
                                
                                signal_peptide = embl_data_frame['Signal peptide'][epub]
                                
                                tm_region_info = embl_data_frame['Transmembrane'][epub]
                                # print(uniprot_id, tm_region_info)
                                
                                if pd.isnull(tm_region_info):
                                    tm_protein_type = ''
                                elif not pd.isnull(tm_region_info):
                                    if tm_region_info.count("TRANSMEM") == 1:
                                        tm_protein_type = 'Single-pass'
                                    elif tm_region_info.count("TRANSMEM") > 1:
                                        tm_protein_type = 'Multi-pass'
                                    else:
                                        tm_protein_type = ''
                            else:
                                pass
                else:
                    # print("embl_data_frame is empty")
                    pass
                temp_array.append(
                        [value,
                         gene_id,
                         refseq_id,
                         protein_gi,
                         uniprot_id,
                         signal_peptide,
                         tm_protein_type,
                         tm_region_info,
                         " ".join(key.split(' ')[1:]),
                         ", ".join(list(set(gene_status))),
                         ", ".join(list(set(gene_type))),
                         ", ".join(list(set(arr_localization))),
                         ", ".join(list(set(arr_function))),
                         ", ".join(str(x) for x in sorted(list(set(pubmed_ids)))),
                         ])
            
            else:
                new_gb_id = '.'.join((refseq_id.split('.')[0], str(int(refseq_id.split('.')[1]) + 1)))
                if new_gb_id in list(refseq_data_frame['protein_accession.version']):
                    temp_updated.append(new_gb_id)
                else:
                    temp_removed.append(refseq_id)
    df_temp = pd.DataFrame(temp_array,
                           columns=['tail',
                                    'gene_id',
                                    'refseq_id',
                                    'protein_gi',
                                    'uniprot_id',
                                    'signal_peptide',
                                    'tm_protein_type',
                                    'tm_region_info',
                                    'fasta_header',
                                    'status',
                                    'type_of_gene',
                                    'component',
                                    'function',
                                    'pubmed_id'])
    df_temp_updated = pd.DataFrame(temp_updated, columns=['tail'])
    df_temp_removed = pd.DataFrame(temp_removed, columns=['tail'])
    
    return [df_temp, df_temp_updated, df_temp_removed]


# data processing #

df_processed, df_updated, df_removed = process_fasta_file(f_handle=proteome_f_data,
                                                          refseq_data_frame=df_gene2refseq,
                                                          go_data_frame=df_gene2go,
                                                          pubmed_data_frame=pubmed_file_name,
                                                          gene_info_data_frame=df_geneinfo)

# save processed data to txt
df_processed.to_csv(os.path.join(proteome_dir, processed_data_f_name), sep='\t')

# write  list of updated refseq ids to txt file
df_updated.to_csv(os.path.join(proteome_dir, out_updated_ids_name), sep='\t', header=False, index=False)

# download fasta files for the updated entries
download_updated(out_f_name=updated_seq_fasta_f_name, seq_to_download=df_updated['tail'].tolist())

# write removed refseq ids to txt file
df_removed.to_csv(os.path.join(proteome_dir, out_removed_ids_name), sep='\t', header=False, index=False)
