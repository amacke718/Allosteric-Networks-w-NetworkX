import Bio.PDB
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import time

# Load a PDB of your Protein Structure
pdb_code = "6ugd"
pdb_filename = "./6UGD_APO.pdb"
MAP_FOLD = 'STARTING'

PROT_STATE ='APO'

HEX_CHAIN = ['A', 'B', 'C', 'D', 'E', 'F']
NO_RES = 317*6
STARTING_RES = 156
NO_RES_PRO = 317

# Establish the differences between the chains
chainA_cij_index=range(156, 473)
chainB_cij_index=range(473, 790)
chainC_cij_index=range(790, 1107)
chainD_cij_index=range(1107, 1424)
chainE_cij_index=range(1424, 1741)
chainF_cij_index=range(1741, 2058)
  
####################################################################################################
# Upload Contact Data
####################################################################################################
# Contact Map computed with VMD using "measure" function
# CHAIN1 RES1 CHAIN2 RES2 CONTACTS
# Commented out the Empty dataframes
ChainA_ChainA = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_A_A.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainA_ChainB = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_A_B.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainA_ChainC = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_A_C.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainA_ChainD = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_A_D.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainA_ChainE = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_A_E.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainA_ChainF = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_A_F.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)

ChainB_ChainB = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_B_B.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainB_ChainC = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_B_C.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainB_ChainD = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_B_D.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainB_ChainE = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_B_E.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainB_ChainF = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_B_F.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)

ChainC_ChainC = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_C_C.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainC_ChainD = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_C_D.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainC_ChainE = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_C_E.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainC_ChainF = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_C_F.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)

ChainD_ChainD = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_D_D.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainD_ChainE = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_D_E.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
#ChainD_ChainF = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_D_F.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)

ChainE_ChainE = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_E_E.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)
ChainE_ChainF = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_E_F.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)

ChainF_ChainF = pd.read_csv('./{0}_CONTACTS/atom_pair_contacts_F_F.dat'.format(MAP_FOLD), delimiter=',', comment='#', header=None)

# Adjust the list based on which chains were found to be in contact w/ e/o
Non_Empty_for_Concat=[ChainA_ChainA, ChainA_ChainB, ChainA_ChainE, ChainA_ChainF, ChainB_ChainB, ChainB_ChainC, ChainB_ChainF, ChainC_ChainC, ChainC_ChainD, ChainD_ChainD, ChainD_ChainE, ChainE_ChainE, ChainE_ChainF, ChainF_ChainF]
        
All_Contacts = pd.concat(Non_Empty_for_Concat, axis=0, ignore_index=True)

# Identify the col for easier looping
CHAIN1=All_Contacts[0]
RES1=All_Contacts[1]
CHAIN2=All_Contacts[2]
RES2=All_Contacts[3]

RES_CONTACTS=All_Contacts[4]

####################################################################################################
# This section uploads the PDB by which a needed dic is made
####################################################################################################

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# Upload the PDB
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
all_resname_by1 = []
all_resnum = []
for residue in structure.get_residues():
    name = residue.get_resname()
    all_resname_by1.append(d3to1[name])
    
    num = residue.id[1]
    all_resnum.append(num)
    
#print(all_resnum)

keys = all_resnum
values = all_resname_by1

# Make dictionary for easy resid/resname calling
RESNUM_RESNAME_DIC = {}
for x, y in zip(keys, values):
    RESNUM_RESNAME_DIC[x] = y

####################################################################################################
# This section defines necessary functions for handling the DCCM data
####################################################################################################
# DCCM pairs for Hexamer
def PARSE_Cij_PAIRS (RESIDUE_INDEX):
    for residue in range(0, len(I_NUM)):
        if I_NUM[residue] == RESIDUE_INDEX:
            res1=I_NUM[residue]+STARTING_RES   # Index converted to a Res ID
            res2=J_NUM[residue]+STARTING_RES
            # Assign Chain 1 Identity
            if res1 in chainA_cij_index:
                chain1='A'
                res1 = int(res1)
            elif res1 in chainB_cij_index:
                chain1='B'
                res1 = int(res1 - NO_RES_PRO)
            elif res1 in chainC_cij_index:
                chain1='C'
                res1 = int(res1 - NO_RES_PRO - NO_RES_PRO)
            elif res1 in chainD_cij_index:
                chain1='D'
                res1 = int(res1 - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO)
            elif res1 in chainE_cij_index:
                chain1='E'
                res1 = int(res1 - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO)
            elif res1 in chainF_cij_index:
                chain1='F'
                res1 = int(res1 -NO_RES_PRO - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO)
            # Assign Chain 2 Identity & Write it to Output
            if res2 in chainA_cij_index:
                chain2='A'
                res2 = int(res2)
                dummy_list=[res1, res2, Cij[residue]]
                OUTPUT.write('{0}{1},{2}{3},{4} \n'.format(chain1, dummy_list[0], chain2, dummy_list[1], dummy_list[2]))
            elif res2 in chainB_cij_index:
                chain2='B'
                res2 = int(res2 - NO_RES_PRO)
                dummy_list=[res1, res2, Cij[residue]]
                OUTPUT.write('{0}{1},{2}{3},{4} \n'.format(chain1, dummy_list[0], chain2, dummy_list[1], dummy_list[2]))
            elif res2 in chainC_cij_index:
                chain2='C'
                res2 = int(res2 - NO_RES_PRO - NO_RES_PRO)
                dummy_list=[res1, res2, Cij[residue]]
                OUTPUT.write('{0}{1},{2}{3},{4} \n'.format(chain1, dummy_list[0], chain2, dummy_list[1], dummy_list[2]))
            elif res2 in chainD_cij_index:
                chain2='D'
                res2 = int(res2 - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO)
                dummy_list=[res1, res2, Cij[residue]]
                OUTPUT.write('{0}{1},{2}{3},{4} \n'.format(chain1, dummy_list[0], chain2, dummy_list[1], dummy_list[2]))
            elif res2 in chainE_cij_index:
                chain2='E'
                res2 = int(res2 - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO)
                dummy_list=[res1, res2, Cij[residue]]
                OUTPUT.write('{0}{1},{2}{3},{4} \n'.format(chain1, dummy_list[0], chain2, dummy_list[1], dummy_list[2]))
            elif res2 in chainF_cij_index:
                chain2='F'
                res2 = int(res2 - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO - NO_RES_PRO)
                dummy_list=[res1, res2, Cij[residue]]
                OUTPUT.write('{0}{1},{2}{3},{4} \n'.format(chain1, dummy_list[0], chain2, dummy_list[1], dummy_list[2]))

####################################################################################################
# This section handles the DCCM data
####################################################################################################

print('~ Loading Correlation Data ~ \n')
# Cross Correlation Matrix Parsing
Cij_DATA=np.genfromtxt('DCCM_{}_list.txt'.format(PROT_STATE), comments='#')
I_NUM=Cij_DATA[:, 0]
J_NUM=Cij_DATA[:, 1]
Cij=Cij_DATA[:, 2]

# Cross Correlation Pairs

if os.path.exists('{0}_Cij_pairs.txt'.format(PROT_STATE)):
    os.remove('{0}_Cij_pairs.txt'.format(PROT_STATE))
print('~ About to Parse through the Cij to correct the indices, this can take a while! ~ \n')
OUTPUT = open(r'{0}_Cij_pairs.txt'.format(PROT_STATE), 'a')
OUTPUT.write('### Cij PAIRS for {0} ###\n'.format(PROT_STATE))
count=0
for i in range(0, NO_RES-1):
    PARSE_Cij_PAIRS(i)
    print('### {}% Progress ###'.format(round((count/len(range(0, NO_RES-1)))*100), 2))
    count=count+1
OUTPUT.close()
print('~ Cij Pairs Renumbered & New List Exported ~')

####################################################################################################
# This section evaluates each residue pair based on the contact map
####################################################################################################

# For Parsing
print('~ Loading Corrected Correlation Data ~ \n')
# All Chains
# Upload the Corrected DCCM Output File
CORRECTED_DCCM = np.genfromtxt('{0}_Cij_pairs.txt'.format(PROT_STATE), delimiter=',', comments='#', dtype=str)
CORRECTED_I_NUM=(CORRECTED_DCCM[:, 0])
CORRECTED_J_NUM=(CORRECTED_DCCM[:, 1])
Cij=(CORRECTED_DCCM[:, 2]).astype(float)

# Make File to hold Node IDs & Cij
if os.path.exists('START_RES_PAIRS.txt'):
    os.remove('START_RES_PAIRS.txt')

RSCB_RES_PAIRS = open(r'START_RES_PAIRS.txt', 'a')
RSCB_RES_PAIRS.write('### Residue Pairs based on CA Contacts ###\n')
RSCB_RES_PAIRS.write('### Contacts computed by VMD measure contacts with a cutoff of 10A \n')
RSCB_RES_PAIRS.write('### NODE1, NODE2, Contacts\n')

for pair in range(0,len(RES_CONTACTS)):   # len(RES_CONTACTS)
    #print('~~~~~~~~~~~~~~~~~~~~~')
    check_res=RES1[pair]
    next_covalent_res_id=check_res+1
    #print(check_res, next_covalent_res_id)
    #print(CHAIN1[pair], CHAIN2[pair])
    if CHAIN1[pair] == CHAIN2[pair]:
        if RES1[pair] == RES2[pair]:
            # skip pairs between same residues 
            continue        
        elif RES2[pair] == next_covalent_res_id:
            # skip pairs between residues that are covalently bonded
            continue
    if RES1[pair] < RES2[pair]:        
	    res1 = RES1[pair]
	    res2 = RES2[pair]
	    #print(res1, res2)
	    
	    resname1 = RESNUM_RESNAME_DIC[res1]
	    resname2 = RESNUM_RESNAME_DIC[res2]
	    #print(resname1, resname2)
	    
	    contact = RES_CONTACTS[pair]
	    
	    RSCB_RES_PAIRS.write('{0}{1},{2}{3},{4}\n'.format(CHAIN1[pair], res1, CHAIN2[pair], res2, contact))
    else:
        # skip the repetative Residue Pairs
        continue
RSCB_RES_PAIRS.close()

####################################################################################################
# This section matches the residues pairs found in contact with their Cij from MD
####################################################################################################

# Upload the Res Pairs
print('~ Loading Res Pairs from CA Contacts ~ \n')
NODES_TO_PLOT = np.genfromtxt('START_RES_PAIRS.txt', delimiter=',', comments='#', dtype=str)
NODE_RES1 = (NODES_TO_PLOT[:, 0])
NODE_RES2 = (NODES_TO_PLOT[:, 1])
NODE_Iij = NODES_TO_PLOT[:, 2]   # Don't need the strengths for the Graph


print('~ Creating File Combining the Residue Pairs from Contacts with Cij ~ \n')
if os.path.exists('START_{0}_Cij.txt'.format(PROT_STATE)):
    os.remove('START_{0}_Cij.txt'.format(PROT_STATE))

NETWORK_INFO = open(r'START_{0}_Cij.txt'.format(PROT_STATE), 'a')
NETWORK_INFO.write('### All Graph Network Information for {0} ###\n'.format(PROT_STATE))
NETWORK_INFO.write('### Node1, Node2, Correlation\n')

count=0
print('~ Match Nodes with Correlation Data ~ \n')
print('Length of Cij: ', len(Cij))
for corr_pair in range(0, len(Cij)):   # len(Cij)
    #print(CORRECTED_I_NUM[corr_pair], CORRECTED_J_NUM[corr_pair])
    for node_pair in range(0, len(NODE_Iij)):
        #print(NODE_RES1[node_pair], NODE_RES2[node_pair])
        if NODE_RES1[node_pair] == CORRECTED_I_NUM[corr_pair] and NODE_RES2[node_pair] == CORRECTED_J_NUM[corr_pair]:
            #print('~ Nest 1 ~')
            my_cij = Cij[corr_pair]
            NETWORK_INFO.write('{0},{1},{2}\n'.format(NODE_RES1[node_pair], NODE_RES2[node_pair], my_cij))
            continue
    print('### {}% Progress ###'.format(round((count/len(Cij))*100, 1)))
    count=count+1
NETWORK_INFO.close()           
print('~ Graph Network Information for {0} Extracted ~ \n'.format(PROT_STATE))

print('~ Done ~')
