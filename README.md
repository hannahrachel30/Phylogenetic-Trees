
UPGMA and Neighbor-Joining Phylogenetic Tree Construction: 

Overview
This repository contains implementations of two phylogenetic tree construction algorithms—UPGMA (Unweighted Pair Group Method with Arithmetic Mean) and Neighbor Joining—based on 
distance matrix inputs. Both implementations output phylogenetic trees in Newick format. The repository also includes sample input files (DM-p127.txt,DM-p139.txt) to demonstrate how
the programs work.

Key Features
UPGMA: Centroid-linkage clustering using average inter-cluster distances.
Neighbor Joining: Systematic calculation of the “transition” distance matrix, recomputing cluster distances after each merge.

Implementation Details
In this project, the UPGMA implementation starts by treating each taxon as a separate cluster and then repeatedly merges the two closest clusters based on their average inter-cluster 
distances, recalculating the relevant distances at each step, until only one cluster remains. Throughout, the program logs which clusters are merged along with their average distance.
The Neighbor Joining algorithm, on the other hand, calculates the average distance for each taxon (or cluster), constructs a transition matrix  that helps identify the pair of clusters
most likely to be neighbors, then merges them and updates distances to form a new cluster. Both algorithms ultimately produce a single phylogenetic tree, which is represented in Newick
format; UPGMA outputs centroid-linkage trees, while Neighbor Joining includes branch lengths inferred from the updated transition matrices.

Troubleshooting and Notes
Efficiency: Ensure the cluster distance updates are computed without re-checking all original pairwise distances unnecessarily.
Floating-Point Precision: You may see slight differences in decimal outputs compared to example scripts if you’re using different rounding or data types.
Console Output Consistency: The instructions require specific step-by-step logs. Verify your logs match the expected format.
