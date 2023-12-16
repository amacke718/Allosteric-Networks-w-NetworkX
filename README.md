# Allosteric-Networks-w-NetworkX
How I mapped protein allostery with NetworkX
This procedure was developed for the analysis carried out in the the description of Spastin Allostery (https://doi.org/10.1063/5.0139273)

## Step 1: Calculate a Contact Map

We first generate a contact map based on the Ca connectivity (10A) with the VMD measure function (in a .tcl). This connectivity will determine the edge connectivity of the network. 

## Step 2: Calculate the DCCM from MD simulation

We calculated the average DCCM across multiple trajectories with MD Traj

## Step 3: Convert the DCCM Indices to Res IDs

The DCCM converts the resid to indices starting at 0. While managable as is, it is much more convenient to convert the indices to Res ID (especially since we work with quaternary protein structures with multiple chains). That procedure has been included in the "ContactMap+DCCM.py"

This will define the weights of the edges in the network. 

It is important to note that in our procedure we do not include covalently bound edges.

## Step 4: Build the Network & Evaluate an appropriate cutoff

Weights that fall below a determined c_ij cutoff are excluded from the network as they have insignificant contribution to the porpagation of allosteric signals but contribute to extended computational expense. It is therefore advised to evaluate (based on the ratio of the largest cluster to no. of residues) an appropriate c_ij cutoff for your system. As a part of the provided code, the network topology is described by providing the no. of nodes, no. of edges, the assortativity (described in the following publication: https://doi.org/10.1021/acs.jpcb.1c04792 ) and the Largest Cluster/N ratio.

While evaluating it is recommended to comment out the path analysis portion of the code. 

## Step 5: Evaluate the Network and Compute and Path Analysis of interest

As a part of the provided code, the number of neighbors for each residue (Cij_0.6cutoff_PSG_Neighbors.txt) and the betweenness centrality by weight for each residue (Cij_0.6cutoff_Betweeness.txt) is provided as output. For the dynamic network, we identified those residues that are in the top 10% betweenness centrality which are the most important nodes for the network. 

If you are interested in collecting paths, you can uncomment the path analysis portion. In this section, the shortest path by weight between the ndicated residues (at the top of the code) are computed. The shortest path describes the optimal path for the signal to propogate via Dijkstra's algorithm (DOI: 10.1109/ICSMC.2000.886462). However, this path is not always the best indication of network shifts when comparing a ligand bound v. unbound state of your protein. It may then be necessary to extract suboptimal paths (described in the following publication: https://doi.org/10.1021/ct4008603 ). This code will utilze Yen's Algorithm to extract an indicated number of suboptimal paths (we use 20,000 per residue combination). Keep in mind that the larger your system is the more paths you should take to be appropriately descriptive. This does have an amount of computational cost associated with it. We also utilize a compbination of residues to descibe paths rather than a single residue pair as residues that are covalently bonded will feel the region similarly, especially when looking at protein binding regions. 
