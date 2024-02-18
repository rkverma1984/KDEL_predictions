import os,glob
import mdtraj as md
import numpy as np
import _heatmap
from _heatmap import *
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument('--ref_pdb', type=str,default='human_6ddf.DAMGO.xtal_pose_rep1.1.pdb')
#parser.add_argument('--dcd_file',type=str,default='human_6ddf.DAMGO.xtal_pose_rep1.1_e0-120.protein.wrapped.dcd')
#args = parser.parse_args()

zone_threshold = 0.5 # within this zone, the contact residues will be considered.
frequency_threshold = 0.005 # only above this percent, the contact residues will be considered.

#reference_pdb = args.ref_pdb
#traj_dcd = args.dcd_file

def load_traj_by_heavy(top_pdb, traj_dcd):
    structure=md.load(top_pdb)
    atom_slices = structure.topology.select_atom_indices(selection = 'heavy') # consider heavy atoms
    heavy_structure = structure.atom_slice(atom_slices)
    trajectory = md.load(traj_dcd,top = structure,atom_indices = atom_slices,stride=10 )
    trajectory.superpose(trajectory, frame = 0)
    
    chain_list = []
    for chain in trajectory.topology.chains:
        chain_list.append(chain)
    return chain_list,trajectory

def load_traj_by_all(top_pdb, traj_dcd):
    structure=md.load(top_pdb)
    trajectory = md.load(traj_dcd,top = structure,stride=10 )
    trajectory.superpose(trajectory, frame = 0)

    chain_list = []
    for chain in trajectory.topology.chains:
        chain_list.append(chain)
    return chain_list,trajectory



def load_traj_by_bb(top_pdb, traj_dcd):
    structure=md.load(top_pdb)
    atom_slices = structure.topology.select('backbone')
    trajectory = md.load(traj_dcd,top = structure,atom_indices = atom_slices,stride=10 )
    trajectory.superpose(trajectory, frame = 0)

    chain_list = []
    for chain in trajectory.topology.chains:
        chain_list.append(chain)
    return chain_list,trajectory

def _index_info(chain_list, chain_id_pair):
    """
    Default order: The first chain is not the receptor!
    e.g. chain_id_pair = [0, -1] for human_6ddf.DAMGO.xtal
    """
    rec_id = chain_id_pair[1]
    o_id = chain_id_pair[0]

    receptor = chain_list[rec_id]
    other_chain = chain_list[o_id]

    receptor_starting_index = 0
    receptor_ending_index = 0
    real_receptor_index = int(str(chain_list[rec_id].residue(0))[3:]) # the first residue index shows up in the pdb

    other_starting_index = 0
    other_ending_index = 0
    real_other_index = int(str(chain_list[o_id].residue(0))[3:]) # the first residue index shows up in the pdb

    for c in chain_list[:rec_id]:
        receptor_starting_index+= c.n_residues # the first mdtraj residues, accumulated by multiple chains (0-based)
    for c in chain_list[:rec_id+1]:
        receptor_ending_index+= c.n_residues   # the first mdtraj residues, accumulated by multiple chains (0-based)
        
    if o_id ==0:
        other_starting_index = 0
    elif o_id>0:
        for c in chain_list[:o_id]:
            other_starting_index+= c.n_residues # accumulated mdtraj residues (0-based)
    else:
        print ('wrong other chain index')
    other_ending_index = 0
    for c in chain_list[:o_id+1]:
        other_ending_index+= c.n_residues    # accumulated mdtraj residues (0-based)
    return real_receptor_index,receptor_starting_index,receptor_ending_index,real_other_index,other_starting_index,other_ending_index


def convert_Indexes2Labels( chain_list,chain_id_pair, receptor_indexes, real_receptor_index,receptor_starting_index, other_chains_indexes, real_other_index, other_starting_index, gap):
    """
    receptor_indexes:        the list that only contain contacted residues, the mdtraj index(0-based)
    real_receptor_index:     the first residue index shows up in the pdb
    receptor_starting_index: the first mdtraj residues, accumulated by multiple chains (0-based)

    other_chains_indexes:    the list that only contain contacted residues, the mdtraj index(0-based)
    real_other_index:        the first residue index shows up in the pdb
    other_starting_index:    the first mdtraj residues, accumulated by multiple chains (0-based)

    gap: here is a bool type, particularly for Galpha. If any other chains have gap, this part needs to be modified.
    """

    rec_id = chain_id_pair[1]
    o_id = chain_id_pair[0]

    receptor = chain_list[rec_id]
    other_chain = chain_list[o_id]

    receptor_residues_labels = []
    receptor_TMs_labels = []
    other_chain_residues_labels = []

    for r in receptor_indexes:
        resn_r = str(receptor.residue(r-receptor_starting_index))[:3]
        index = r + real_receptor_index - receptor_starting_index
        receptor_residues_labels.append(index)
        #if (index+2) in human_BW_dic.keys():
        #    tm_index = human_BW_dic[index+2]
        #    receptor_TMs_labels.append(str(tm_index)+' '+resn_r)
        #else:
        receptor_TMs_labels.append(str(index)+' '+resn_r) # for kdel
    if gap:
        for o in other_chains_indexes:
            resn_o = str(other_chain.residue(o-other_starting_index))[:3]
            o_i = o+real_other_index-other_starting_index
            if o_i < 56:
                other_chain_residues_labels.append(str(o_i)+' '+resn_o)
            if o_i > 55 and o_i < 108:
                other_chain_residues_labels.append(str(o_i+126)+' '+resn_o)
            if o_i > 107:
                other_chain_residues_labels.append(str(o_i+133)+' '+resn_o)
    else:
        for o in other_chains_indexes:
            resn_o = str(other_chain.residue(o-other_starting_index))[:3]
            index = o + real_other_index - other_starting_index
            other_chain_residues_labels.append(str(index)+' '+resn_o)
            
    return receptor_residues_labels,receptor_TMs_labels,other_chain_residues_labels


def calc_frequency_matrix(chain_list, trajectory, chain_id_pair, zone_threshold, frequency_threshold, gap):
    connect_dic={}
    r_rec_index,rec_start_i,rec_end_i,r_oth_index,oth_start_index,oth_end_index = _index_info(chain_list, chain_id_pair)

    for i in tqdm(range(rec_start_i, rec_end_i)): # loop receptor index
        for j in range(oth_start_index, oth_end_index): # loop other protein index
            counts = 0
            dist = md.compute_contacts(trajectory,[[i,j]])
            data = dist[0]
            for d in data:
                if d[0] < zone_threshold:
                    counts+=1
            contact_fre = float(counts)/float(len(data))
            #print (len(data))
            if contact_fre>frequency_threshold:
                connect_dic['R:'+str(i)+'_O:'+str(j)] = contact_fre #key is receptor resid + other_chain resid, value is frequency

    receptor_indexes = []
    other_chains_indexes = []
    for key in connect_dic.keys():
        receptor_indexes.append(int(key.split('_')[0][2:]))
        other_chains_indexes.append(int(key.split('_')[1][2:]))
    receptor_indexes = sorted(list(set(receptor_indexes)))
    other_chains_indexes = sorted(list(set(other_chains_indexes)))


    frequency_matrix = np.zeros((len(receptor_indexes), len(other_chains_indexes)))
    for i in range(len(receptor_indexes)):
        for j in range(len(other_chains_indexes)):
            key = 'R:'+str(receptor_indexes[i])+'_O:'+str(other_chains_indexes[j])
            if key in connect_dic.keys():
                frequency_matrix[i][j] = connect_dic[key]
    return frequency_matrix, receptor_indexes, other_chains_indexes



def calc_freq_matrix_fixed_indexes(chain_list,trajectory, chain_id_pair, receptor_indexes, other_chains_indexes, zone_threshold, frequency_threshold, gap):
    rec_id = chain_id_pair[1]
    o_id = chain_id_pair[0]
    receptor = chain_list[rec_id]
    other_chain = chain_list[o_id]

    connect_dic = {}
    
    for i in tqdm(range(len(receptor_indexes))): # loop receptor index
        for j in range(len(other_chains_indexes)): # loop other protein index
            counts = 0
            dist = md.compute_contacts(trajectory,[[i,j]])
            data = dist[0]
            for d in data:
                if d[0] < zone_threshold:
                    counts+=1
            contact_fre = float(counts)/float(len(data))
            #print (len(data))
            if contact_fre>frequency_threshold:
                connect_dic['R:'+str(i)+'_O:'+str(j)] = contact_fre #key is receptor resid + other_chain resid, value is frequency


    frequency_matrix = np.zeros((len(receptor_indexes), len(other_chains_indexes)))
    for i in range(len(receptor_indexes)):
        for j in range(len(other_chains_indexes)):
            key = 'R:'+str(receptor_indexes[i])+'_O:'+str(other_chains_indexes[j])
            if key in connect_dic.keys():
                frequency_matrix[i][j] = connect_dic[key]

    return frequency_matrix

def cal_freq_binding_site(ligand_top,top_pdb, traj_dcd, fixed_list, zone_threshold, frequency_threshold):
    f = open(top_pdb,'r')
    lines=f.read().split('\n')
    f.close()
    start_index = None
    end_index = None
    for line in lines:
        if line.startswith('ATOM') and ' N ' in line:
            start_index = int(line[23:27].strip())
            break
    print('Starting residue index is',start_index)    
    structure=md.load(top_pdb)
    trajectory = md.load(traj_dcd,top = structure,stride=10 )
    trajectory.superpose(trajectory, frame=0)
    print ('There are',structure.n_chains, 'chains in this system.')
    atom_slices = structure.topology.select_atom_indices(selection = 'heavy') # consider heavy atoms
    heavy_structure = structure.atom_slice(atom_slices)

    rec_range = []
    if ligand_top:
        ligand_resid = 0
        rec_range = [1,heavy_structure.n_residues+1]
    else:
        ligand_resid =structure.n_residues-1
        rec_range = [0,heavy_structure.n_residues+0]
    lig_slices = structure.topology.select('resid '+str(ligand_resid))
    lig_traj = md.load(traj_dcd,top = structure, atom_indices=lig_slices)  
    connect_dic={}#key is residue index, value is frequency

    if fixed_list is not None:
        for i in [(f+rec_range[0]) for f in fixed_list]:
            counts = 0
            dist = md.compute_contacts(trajectory,[[i,ligand_resid]])
            data = dist[0]
            for d in data:
                if d[0] < zone_threshold:
                    counts+=1
            contact_fre = float(counts)/float(len(data))
            
            resn_r = str(heavy_structure.topology.residue(i-rec_range[0]))[:3]
            connect_dic[resn_r+' '+str(i+start_index-rec_range[0])] = contact_fre
        print ('use  fixed_list')
    else:
        for i in range(rec_range[0],rec_range[1]):
            counts = 0
            dist = md.compute_contacts(trajectory,[[i,ligand_resid]])
            data = dist[0]
            for d in data:
                if d[0] < zone_threshold:
                    counts+=1
            contact_fre = float(counts)/float(len(data))
            resn_r = str(heavy_structure.topology.residue(i-rec_range[0]))[:3]

            if contact_fre>frequency_threshold:
                connect_dic[resn_r+' '+str(i+start_index-rec_range[0])] = contact_fre

    return connect_dic


