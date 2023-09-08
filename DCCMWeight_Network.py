import Bio.PDB
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import time

PROT_STATE ='APO'
HEX_CHAIN = ['A', 'B', 'C', 'D', 'E', 'F']
NO_RES = 317*6
STARTING_RES = 156
NO_RES_PRO = 317
# Set the Cij cutoff - It is necessary to determine an appropriate cutoff for your system of interest
Cij_CUTOFF = 0.5
GRAPH_COLOR = 'honeydew'

####################################################################################################
# This section Filters the left over pairs with a specified Cij cutoff for the Graph Network
####################################################################################################

Cij_TO_FILTER = np.genfromtxt('START_{}_Cij.txt'.format(PROT_STATE), delimiter=',', comments='#', dtype=str)
CORRECTED_I_NUM=(Cij_TO_FILTER[:, 0])
CORRECTED_J_NUM=(Cij_TO_FILTER[:, 1])
Cij=(Cij_TO_FILTER[:, 2]).astype(float)

# Create lists to hold Cij with a given cutoff
cutoff_I = []
cutoff_J = []
cutoff_Cij = []
for cij_pair_I, cij_pair_J, corr_value in zip(CORRECTED_I_NUM, CORRECTED_J_NUM, Cij):
    if np.abs(corr_value) > Cij_CUTOFF:
        cutoff_I.append(cij_pair_I)
        cutoff_J.append(cij_pair_J)
        cutoff_Cij.append(corr_value)
# Weight is determined by w=-log(|c_ij|)
# Create File to hold the Weight dataframe with appropriate cutoff
if os.path.exists('{0}_Cij_{1}cutoff_WEIGHTS.txt'.format(PROT_STATE, Cij_CUTOFF)):
    os.remove('{0}_Cij_{1}cutoff_WEIGHTS.txt'.format(PROT_STATE, Cij_CUTOFF))

WEIGHT_OUTPUT = open(r'{0}_Cij_{1}cutoff_WEIGHTS.txt'.format(PROT_STATE, Cij_CUTOFF), 'a')
WEIGHT_OUTPUT.write('### Graph Network Pairs & Weights for {0} from calculated DCCM of MD ###\n'.format(PROT_STATE))
WEIGHT_OUTPUT.write('# Cij cutoff = {0}\n'.format(Cij_CUTOFF))
WEIGHT_OUTPUT.write('# Residue 1, Residue 2, -log(|Cij|)\n')
time.sleep(1)

# Need to carry out some calculations for the weights
print('~ Calculating the Weights based on the Cij ~ \n')
count=0
for cij_pair_I, cij_pair_J, corr_value in zip(cutoff_I, cutoff_J, cutoff_Cij):
    my_weight=np.log((np.abs(corr_value)))*-1
    # the avg dccm will compute X/0 as 'nan' which should be a 0 so it is filtered out here along with the infinities that result in the log(0)
    if my_weight != float('inf') and np.isnan(my_weight) == False:
        WEIGHT_OUTPUT.write('{0},{1},{2}\n'.format(cij_pair_I, cij_pair_J, my_weight))

    count=count+1

WEIGHT_OUTPUT.close()

####################################################################################################
# This section uploads the Graph Network data & Evaluates the Network
####################################################################################################

df = pd.read_csv('{0}_Cij_{1}cutoff_WEIGHTS.txt'.format(PROT_STATE, Cij_CUTOFF), sep=",", comment='#')
df.columns = ['RESNUM1', 'RESNUM2', 'WEIGHTS']
check_res1 = df['RESNUM1'].astype(str)
check_res2 = df['RESNUM2'].astype(str)
NET_WEIGHTS = df['WEIGHTS']
df['PAIRS'] = df[['RESNUM1', 'RESNUM2']].apply(tuple, axis=1)
NET_EDGES = list(df['PAIRS'])

####################################################################################################
# Build the Graph Network
# Undirected Graph
####################################################################################################

GRAPH = nx.Graph()

####################################################################################################

print('~~~~~~~~~~ Hexamer {0} Cij cutoff used: {1} ~~~~~~~~~~'.format(PROT_STATE, Cij_CUTOFF))
node_id=[]
count=0
for i in range(0,len(NET_WEIGHTS)):
    if NET_EDGES[i][0] not in node_id:   # and NET_EDGES[i][1] not in node_id
        node_id.append(NET_EDGES[i][0])
    if NET_EDGES[i][1] not in node_id:
        node_id.append(NET_EDGES[i][1])
    count=count+1

# Dictionary with Strings first for Tree Parsing
NODES_DICT={}
for y,z in zip(node_id,range(0,len(node_id))):
   NODES_DICT[y]=z  
   
# Dictionary with Int first for Graph Labeling  
REV_NODES_DICT=dict(list(enumerate(node_id)))
NODES_USED=list(range(0,len(node_id)))

GRAPH.add_nodes_from(NODES_USED)

for i in range(0,len(NET_EDGES)):
    KEY1=NODES_DICT.get(NET_EDGES[i][0])
    KEY2=NODES_DICT.get(NET_EDGES[i][1])
    GRAPH.add_weighted_edges_from([(KEY1, KEY2, NET_WEIGHTS[i])])    

for n in GRAPH.nodes():
    GRAPH.nodes[n]['label'] = REV_NODES_DICT[n]

FIG = plt.figure(figsize=(30,30))
nx.draw(GRAPH, with_labels=True, font_size=8, font_weight='bold', labels=REV_NODES_DICT, node_color=GRAPH_COLOR, node_size=700, edge_color='black', width=2.0, pos=nx.spring_layout(GRAPH,seed=71896))
plt.savefig('{0}_Cij_{1}cutoff_PSG.png'.format(PROT_STATE, Cij_CUTOFF), dpi=100)  
#plt.show()

print('No. of Network Nodes: ', GRAPH.number_of_nodes())
print('No. of Network Edges: ', GRAPH.number_of_edges())
####################################################################################################
# Evaluate some graph network characteristics
####################################################################################################
# Betweeness is exported as a dictionary - each node has a value
if os.path.exists('{0}_Cij_{1}cutoff_PSG_Betweeness.txt'.format(PROT_STATE, Cij_CUTOFF)):
    os.remove('{0}_Cij_{1}cutoff_PSG_Betweeness.txt'.format(PROT_STATE, Cij_CUTOFF))

betweeness = nx.betweenness_centrality(GRAPH, normalized=True, endpoints=True)
with open ('{0}_Cij_{1}cutoff_Betweeness.txt'.format(PROT_STATE, Cij_CUTOFF), 'w') as f:
    f.write('# PSG Betweeness Centrality \n# Node, Betweeness\n')
    for key, value in betweeness.items():
        #f.write('%s:%s\n' % (key, value))
        KEY_NAMES=REV_NODES_DICT.get(key)
        f.write('%s,%s\n' % (KEY_NAMES, value))
with open ('KEYS+RESID.txt', 'w') as f:
    for key, value in betweeness.items():
        KEY_NAMES=REV_NODES_DICT.get(key)
        f.write('%s,%s\n' % (key, KEY_NAMES))
        
# Neighbors 
if os.path.exists('{0}_Cij_{1}cutoff_PSG_Neighbors.txt'.format(PROT_STATE, Cij_CUTOFF)):
    os.remove('{0}_Cij_{1}cutoff_PSG_Neighbors.txt'.format(PROT_STATE, Cij_CUTOFF))

neighbors = nx.average_neighbor_degree(GRAPH)
with open ('{0}_Cij_{1}cutoff_PSG_Neighbors.txt'.format(PROT_STATE, Cij_CUTOFF), 'w') as f:
    f.write('# PSG Nodes with K Links\n# Manually Added to the File by running the percentages individually \n# Node, No. Links\n')
    for n in GRAPH.nodes():
        f.write('{0},{1}\n'.format(REV_NODES_DICT[n],len(GRAPH.edges(n))))

# Print out the Assortativity for logging
assortativity = nx.degree_pearson_correlation_coefficient(GRAPH)
print('Assortativity for Entire Graph: ', assortativity)

# Print out the Largest Cluster for logging
largest_cluster = max(nx.connected_components(GRAPH), key=len)
print('Largest Cluster: ', (len(largest_cluster)))
print('Largest Cluster/N: ', (len(largest_cluster))/NO_RES)

####################################################################################################
# Path Analysis:
# Shortest path by weight determined with Dijkstra Path
# Suboptimal paths are determined with what is assumed to be Yen's Algorithm
####################################################################################################
# Network X version of Yen?
# https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.simple_paths.shortest_simple_paths.html
from itertools import islice

def k_shortest_paths(G, source, target, k):
    return list(
        islice(nx.shortest_simple_paths(G, source, target, weight='weight'), k)
    )

# Groups to be Mapped
# all output will be saved to the following location
GROUP_FOLD = 'atp_to_pl1'

chain1 = 'A'
chain2 = 'F'

# 2 Positions in ATP pocket 
group_mon1 = ['{}235'.format(chain1)]

# 2 Positions in PL2
group_mon2 = ['{}312'.format(chain2)]

# If you want to take paths between a larger group, just include a longer list for looping
#group_mon2 = ['{}739'.format(chain2), '{}740'.format(chain2), '{}741'.format(chain2), '{}742'.format(chain2),'{}743'.format(chain2), '{}744'.format(chain2), '{}745'.format(chain2), '{}746'.format(chain2), '{}747'.format(chain2), '{}748'.format(chain2), '{}749'.format(chain2), '{}750'.format(chain2), '{}751'.format(chain2), '{}752'.format(chain2), '{}753'.format(chain2), '{}754'.format(chain2), '{}755'.format(chain2), '{}756'.format(chain2)]

# Changable values for evaluating 
no_paths = 20000
my_bins = 50
####################################################################################################
# Make File to hold the paths to evaluate the node degeneracy
if os.path.exists('./{0}/{1}_Cij_{2}cutoff_{0}_dijkstra_paths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF)):
    os.remove('./{0}/{1}_Cij_{2}cutoff_{0}_dijkstra_paths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF))

SHORT_PATHS = open(r'./{0}/{1}_Cij_{2}cutoff_{0}_dijkstra_paths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF), 'a')
SHORT_PATHS.write('### Identified dijsktra paths for {} to CT Hlx ###\n'.format(GROUP_FOLD))
    
# Loop through all of the pairs
for i in range(0,len(group_mon1)):
    mon1_to_check=group_mon1[i]
    for j in range(0, len(group_mon2)):
        mon2_to_check=group_mon2[j]

        # The Network ID is extracted for referencing 
        FROM_RES = NODES_DICT['{}'.format(mon1_to_check)]
        TO_RES = NODES_DICT['{}'.format(mon2_to_check)]
####################################################################################################
        print()
        print('~~~~~~~~~~ {0} to {1} ~~~~~~~~~~'.format(mon1_to_check, mon2_to_check))

        # Dijkstra Path Check
        print('\n~ Dijkstra Path Check ~')
        DP=(nx. dijkstra_path(GRAPH, source=FROM_RES, target=TO_RES, weight='weight'))
        dijkstra_shortest_path = []
        for i in DP:
            dijkstra_shortest_path.append(REV_NODES_DICT[i])
        print(dijkstra_shortest_path)
        dp_path_length=(nx.path_weight(GRAPH, path=DP, weight='weight'))
        print(dp_path_length)
        SHORT_PATHS.write('{0}-{1}\n'.format(round(dp_path_length,4), dijkstra_shortest_path))

####################################################################################################
        # Supoptimal Paths with Yen
  
        # Make File to hold the paths to evaluate the node degeneracy
        if os.path.exists('./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimal_Paths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths)):
            os.remove('./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimal_Paths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths))

        ALL_PATHS = open(r'./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimal_Paths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths), 'a')
        ALL_PATHS.write('### Identified paths between {0} and {1} ###\n'.format(mon1_to_check, mon2_to_check))
        ALL_PATHS.write('### No. of Paths Identified: {} ###\n'.format(no_paths))

        # Make File to hold the path lengths in case you want histograms later
        if os.path.exists('./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimals_Lengths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths)):
            os.remove('./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimals_Lengths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths))

        LENGTHS = open(r'./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimals_Lengths.txt'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths), 'a')
        LENGTHS.write('### Path Length for Each Path Identified from CiJ with cutoff PSG ###\n')
        LENGTHS.write('### No. of Paths Identified: {} ###\n'.format(no_paths))

        yen_paths = []
        yen_path_lengths = []
        count=0
        print('~~~~~~~ Finding {0} Paths ~~~~~~~'.format(no_paths))
        for my_path in k_shortest_paths(GRAPH, FROM_RES, TO_RES, no_paths):
            #print('## Progress {}% ###'.format(round((count/no_paths)*100, 2)))
            # Identified Path
            yen_paths.append(my_path)
            # Extract ResID for saving for each path
            resid_path = []
            for i in my_path:
                resid_path.append(REV_NODES_DICT[i])
                #ALL_PATHS.write("'{}',".format(REV_NODES_DICT[i]))
            # Save the Extracted Path    
            ALL_PATHS.write('{}\n'.format(resid_path))
            #ALL_PATHS.write('\n')
    
            # Get Path Length & Save
            my_path_length=(nx.path_weight(GRAPH, path=my_path, weight='weight'))
            LENGTHS.write('{}\n'.format(my_path_length))
            yen_path_lengths.append(my_path_length)
    
            count=count+1
    
        print('\n~ Yen Path Check ~')
        YP=[]
        for i in yen_paths[0]:
            YP.append(REV_NODES_DICT[i])
        #print(YP)
        #print(yen_path_lengths[0])
        # PLOT THE INDIVIDUAL HISTOGRAM
        plt.figure(figsize=(10,10))   
        # density=True, stacked=True
        plt.hist(yen_path_lengths, bins=my_bins, weights=np.ones(len(yen_path_lengths))/len(yen_path_lengths), align='mid', facecolor='white', edgecolor='black') 
        plt.xlabel('Path Weights', fontsize=28, fontweight='bold')
        plt.xticks(fontsize=16, fontweight='bold')
        plt.yticks(fontsize=16, fontweight='bold')
        plt.savefig('./{0}/{1}_Cij_{2}cutoff_{3}to{4}_{5}paths_YenSuboptimals.png'.format(GROUP_FOLD, PROT_STATE, Cij_CUTOFF, mon1_to_check, mon2_to_check, no_paths))
        plt.close()
        #plt.show()
SHORT_PATHS.close()

print('~ Done ~')
