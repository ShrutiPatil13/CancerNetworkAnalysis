#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import louvain
from scipy.stats import fisher_exact
from collections import Counter

max_tissue = df_gene_seq_comb1_wildtype_uni['max_tissue'].values.tolist()
donors = df_gene_seq_comb1_wildtype_uni['sample id'].values.tolist()

# cluster information
def get_tissue_gene_clusters(clust_list,max_tissue, gene_wildtype,donors):
    clust_tissue_list=[]
    clust_gene_list=[]
    clust_donor_list=[]
    clust_gene_list_len=[]
    clust_donor_list_len=[]
    clust_gene_set=[]
    for clist in clust_list:
        clist_tissue=[]
        clist_gene=[]
        clist_donor=[]
        clist_gene_len=[]
        for i in clist:
            clist_tissue.append(max_tissue[i])
            clist_gene.extend(gene_wildtype[i])
            clist_donor.append(donors[i])
            clist_gene_len.append(len(gene_wildtype[i]))
        clust_tissue_list.append(clist_tissue)
        clust_gene_list.append(clist_gene)
        clust_gene_set.append(list(set(clist_gene)))
        clust_donor_list.append(clist_donor)
        clust_gene_list_len.append(clist_gene_len)
        clust_donor_list_len.append(len(clist_donor))
        
    from collections import Counter
    # get count of tissues in cluster list
    clust_tissue_list_counter=[]
    for clist in clust_tissue_list:
        sorted_counter = Counter(clist).most_common()
        clust_tissue_list_counter.append(dict(sorted_counter))
     
    # get count of genes in cluster list
    clust_gene_list_counter=[]
    for clist in clust_gene_list:
        sorted_counter = Counter(clist).most_common()
        clust_gene_list_counter.append(dict(sorted_counter))
        
    # get count of num_mut in cluster list
    clust_num_mut_counter=[]
    for clist in clust_gene_list_len:
        sorted_counter = Counter(clist).most_common()
        clust_num_mut_counter.append(dict(sorted_counter))
        
    mut_load=[]
    for i in range(len(clust_donor_list)):
        mut_load.append(len(clust_gene_list[i])/len(clust_donor_list[i]))
                
    #return clust_tissue_list_counter, clust_gene_list_counter, clust_gene_set, mut_load, clust_gene_list_len,clust_num_mut_counter, clust_donor_list_len
    return clust_tissue_list, clust_gene_list, clust_gene_set, mut_load, clust_gene_list_len,clust_num_mut_counter, clust_donor_list_len

# louvain clustering - resolution specified
def get_clusters(ssn,tissue_wildtype,gene_wildtype,donors,resolution):
    # louvain clustering
    partition = louvain.find_partition(ssn672, louvain.RBConfigurationVertexPartition, resolution_parameter = resolution)
    ssn672.vs['membership'] = partition.membership
    plot(partition, 'ssn672_louvain_rbc_resolution_{resolution}.png'.format(resolution=resolution),mark_groups=True, **style)    
    clust_list_louvain = list(partition)
    # clusters info
    clust_list_tissue_ml, clust_list_gene_ml, gene_set_ml, mut_load_ml, mut_num_ml, mut_num_ml_counter, num_donors_ml = get_tissue_gene_clusters(clust_list_louvain,tissue_wildtype, gene_wildtype, donors)
    louvain_clusters = pd.DataFrame()
    louvain_clusters['cancer_types'] = clust_list_tissue_ml
    louvain_clusters['genes'] = clust_list_gene_ml
    louvain_clusters['gene_set'] = gene_set_ml
    louvain_clusters['mutational_load'] = mut_load_ml
    louvain_clusters['mut_num_counter'] = mut_num_ml_counter
    louvain_clusters['num_samples'] = num_donors_ml
    #louvain_clusters.to_csv('ssn672_louvain_rbc_clusters_resolution_{resolution}_info.csv'.format(resolution=resolution), index = None, sep=',')
    cols_keep = ['cancer_types','genes','mutational_load','num_samples']
    louvain_clusters = louvain_clusters[cols_keep]
    return louvain_clusters

# Annotation enrichment check - using Fischers exact test
louvain_clusters1 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 1)
louvain_clusters09 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 0.9)
louvain_clusters08 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 0.8)
louvain_clusters05 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 0.5)
louvain_clusters02 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 0.2)
louvain_clusters01 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 0.1)
louvain_clusters2 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 2)
louvain_clusters108 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 1.8)
louvain_clusters105 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 1.5)
louvain_clusters102 = get_clusters(ssn672, tissue_wildtype, gene_wildtype, donors, 1.2)


def gen_dict(cluster_resolution_df):
    # create clusters_resolutionX dict
    clusters_resolution1_dict = {}
    # Iterate over each row (cluster) in the DataFrame
    for idx, row in cluster_resolution_df.iterrows():
        cluster_label = idx  # Assuming the index is the cluster label
        cancer_types = row['cancer_types']
        genes = row['genes']
        mutational_load = row['mutational_load']
        num_samples = row['num_samples']
        
        # Create a dictionary for the cluster
        cluster_info = {
            'cluster_label': cluster_label,
            'cancer_types': cancer_types,
            'genes': genes,
            'mutational_load': mutational_load,
            'num_samples': num_samples
        }  
        # Add the cluster to the clusters_resolution1_dict
        clusters_resolution1_dict[cluster_label] = cluster_info
    return clusters_resolution1_dict

clusters_resolution1_dict = gen_dict(louvain_clusters1)
clusters_resolution09_dict = gen_dict(louvain_clusters09)
clusters_resolution08_dict = gen_dict(louvain_clusters08)
clusters_resolution05_dict = gen_dict(louvain_clusters05)
clusters_resolution02_dict = gen_dict(louvain_clusters02)
clusters_resolution01_dict = gen_dict(louvain_clusters01)
clusters_resolution2_dict = gen_dict(louvain_clusters2)
clusters_resolution108_dict = gen_dict(louvain_clusters108)
clusters_resolution105_dict = gen_dict(louvain_clusters105)
clusters_resolution102_dict = gen_dict(louvain_clusters102)


# Check for enriched annotations using Fischer exact test

resolutions = ['resolution1', 'resolution09', 'resolution08', 'resolution05', 'resolution02', 'resolution01', 'resolution2', 'resolution108', 'resolution105', 'resolution102']  # Add more resolutions as needed
# Dictionary to store p-values for each annotation and each resolution
p_values = {}
# Calculate total number of clusters for resolution 1
total_clusters_resolution1 = len(clusters_resolution1_dict)

# Loop through each annotation
for annotation in ['cancer_types', 'genes']:
    # Initialize dictionary to store counts for each annotation
    annotation_counts = {r: {'in_cluster': Counter(), 'not_in_cluster': Counter()} for r in resolutions}
    
    # Loop through each resolution
    p_values[annotation] = {r: {} for r in resolutions}
    for r in resolutions:
        for cluster_label, cluster_info in globals()[f'clusters_{r}_dict'].items():
            annotations = cluster_info[annotation]
            
            for a in annotations:
                annotation_counts[r]['in_cluster'][a] += 1
                
                for other_cluster_label, other_cluster_info in globals()[f'clusters_{r}_dict'].items():
                    if other_cluster_label != cluster_label:
                        other_annotations = other_cluster_info[annotation]
                        for oa in other_annotations:
                            annotation_counts[r]['not_in_cluster'][oa] += 1
        
        # Perform the Fisher exact test for each annotation
        for a in annotation_counts[r]['in_cluster']:
            in_cluster_count = annotation_counts[r]['in_cluster'][a]
            not_in_cluster_count = annotation_counts[r]['not_in_cluster'][a]
            # Calculate the counts for the contingency table
            in_cluster_remaining = total_clusters_resolution1 - in_cluster_count
            not_in_cluster_remaining = total_clusters_resolution1 - not_in_cluster_count
            # Ensure the remaining counts are non-negative
            in_cluster_remaining = max(in_cluster_remaining, 0)
            not_in_cluster_remaining = max(not_in_cluster_remaining, 0)
            contingency_table = [
                [in_cluster_count, in_cluster_remaining],
                [not_in_cluster_count, not_in_cluster_remaining]
            ]
            odds_ratio, p_value = fisher_exact(contingency_table)
            p_values[annotation][r][a] = p_value

# Initialize a dictionary to store the annotations that are enriched in each resolution compared to resolution 1
enriched_annotations = {r: [] for r in resolutions}

# Compare p-values from other resolutions to resolution 1
for annotation in p_values:
    for r in resolutions:
        if r != 'resolution1':  # Skip resolution 1
            for a in p_values[annotation][r]:
                if a is not None:
                    if a in p_values[annotation]['resolution1']:
                        if p_values[annotation][r][a] < p_values[annotation]['resolution1'][a]:
                            enriched_annotations[r].append((annotation, a))
                        else:
                            print(f"Key {a} not found in resolution 1 for annotation {annotation}")
                
# Print the p-values for each annotation and each resolution
for annotation in p_values:
    print(f"Annotation: {annotation}")
    for r in p_values[annotation]:
        print(f"Resolution {r}: p-values = {p_values[annotation][r]}")

# Print the enriched annotations for each resolution
for r in enriched_annotations:
    if enriched_annotations[r]:
        print(f"Resolution {r}:")
        for annotation, a in enriched_annotations[r]:
            print(f"- {annotation}: {a}")
    else:
        print(f"Resolution {r}: No annotations enriched more than resolution 1")



# Plot % of terms enriched
enriched_counts_cancer_type = {}
enriched_counts_genes = {}

# Calculate the count of enriched terms for each annotation for each resolution
for res, annotations in enriched_annotations.items():
    for annotation, term in annotations:
        if annotation == 'cancer_types':
            if res not in enriched_counts_cancer_type:
                enriched_counts_cancer_type[res] = 0
            enriched_counts_cancer_type[res] += 1
        elif annotation == 'genes':
            if res not in enriched_counts_genes:
                enriched_counts_genes[res] = 0
            enriched_counts_genes[res] += 1
            
# resolution profile
optimiser = louvain.Optimiser()
profile_rbc = optimiser.resolution_profile(ssn672, louvain.RBConfigurationVertexPartition,resolution_range=(0,2))

# Calculate the percentages
total_cancer_types = 16
total_genes = 1264

percentage_enriched_cancer_type = {res: (count / total_cancer_types) * 100 for res, count in enriched_counts_cancer_type.items()}
percentage_enriched_genes = {res: (count / total_genes) * 100 for res, count in enriched_counts_genes.items()}

# Sorting resolutions in ascending order
sorted_resolutions = sorted(percentage_enriched_cancer_type.keys())
custom_resolutions = ['resolution01','resolution02','resolution05','resolution08','resolution09','resolution1','resolution102','resolution105','resolution108','resolution2']
# Plotting the percentages
plt.figure(figsize=(10, 8))

# Cancer Type
cancer_type_data = [percentage_enriched_cancer_type.get(res, 0) for res in custom_resolutions]
plt.plot(custom_resolutions, cancer_type_data, marker='o', label='Cancer Types', color='b')

# Genes
genes_data = [percentage_enriched_genes.get(res, 0) for res in custom_resolutions]
plt.plot(custom_resolutions, genes_data, marker='o', label='Genes', color='g')

plt.xlabel('Resolution parameter', fontsize=20)
plt.ylabel('Percentage of Enriched Terms', fontsize=20)
plt.title('Percentage of Enriched Terms across Resolutions', fontsize=24, y = 1.05, x = 0.45)
custom_xticks = ['0.1','0.2','0.5','0.8','0.9','1','1.2','1.5','1.8','2']
plt.xticks(custom_resolutions, custom_xticks, fontsize=18)  # Rotate x-axis ticks by 90 degrees
plt.yticks(fontsize=18)
plt.legend(fontsize=18)
plt.grid(False)  # Remove grid lines
plt.show()


