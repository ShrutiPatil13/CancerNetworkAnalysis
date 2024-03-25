#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from igraph import Graph, mean
from igraph import VertexSeq
from igraph import EdgeSeq
from igraph import summary
from igraph import plot
from igraph import GraphBase
import networkx as nx

icgc_data = pd.read_csv('/Users/shrutipatil/Downloads/final_consensus_passonly.snv_mnv_indel.icgc.public.maf', sep='\t')
cols_keep = ['Hugo_Symbol','Variant_Classification', 'Genome_Change','Project_Code','Donor_ID']
icgc_data = icgc_data[cols_keep]

census_info = pd.read_csv("Census_all_gene_info.csv",header=0)
census_genes = census_info['Gene Symbol'].values.tolist()

# filter variants 
vars_exclude = ['IGR','Intron',"5'Flank"]
icgc_data_var = icgc_data[~icgc_data['Variant_Classification'].isin(vars_exclude)]
icgc_data_var = icgc_data_var.reset_index(drop = True)

# filter for cgc genes
icgc_data_var = icgc_data_var[~icgc_data_var['Hugo_Symbol'].isin(census_genes)]
icgc_data_var = icgc_data_var.reset_index(drop = True)
icgc_data_var['Project_Code'] = icgc_data_var['Project_Code'].str.split('-').str[0]

# check for recurrent mutations
gene_mut = icgc_data_var['Genome_Change']
gene_mut = gene_mut.values.tolist()
hotspots = []
for change in gene_mut:
    change_count = gene_mut.count(change)
    if change_count > 1:
        hotspots.append(change)        
# icgc data var hotspots
icgc_data_hotspots = icgc_data_var[icgc_data_var['Genome_Change'].isin(hotspots)]
icgc_data_hotspots = icgc_data_hotspots.reset_index(drop=True)

icgc_data_hotspot_donor = icgc_data_hotspots.groupby('Donor_ID').aggregate({'Hugo_Symbol': list, 'Genome_Change': list, 'Variant_Classification': list, 'Project_Code': set})
icgc_data_hotspot_donor = icgc_data_hotspot_donor.reset_index()
donor_tissues923 = icgc_data_hotspot_donor['Tissue'].values.tolist()
donor_genes923 = icgc_data_hotspot_donor['Hugo_Symbol'].values.tolist()
donor_ids923 = icgc_data_hotspot_donor['Donor_ID'].values.tolist()
# by gene
icgc_data_hotspot_gene = icgc_data_hotspots.groupby('Hugo_Symbol').aggregate({'Donor_ID': list, 'Genome_Change': list, 'Variant_Classification': list, 'Project_Code': list})
icgc_data_hotspot_gene = icgc_data_hotspots.reset_index()

eval_genes = icgc_data_hotspot_gene.iloc[:,0]
eval_genes = eval_genes.values.tolist()
cds_mut_all = icgc_data_hotspot_gene.iloc[:,2]
cds_mut_all = cds_mut_all.values.tolist()

eval_sites = []
sites_len = []
for mut in cds_mut_all:
    mut_uni = list(set(mut))
    eval_sites.append(mut_uni)
    sites_len.append(len(mut_uni))
    
sites = pd.DataFrame()
sites['gene'] = eval_genes
sites['sites'] = eval_sites
sites['number_sites'] = sites_len

# create mut_list for cds_mut
def create_dict(cds_mut_uni):
    mut_list=[]
    mutated_seq=[]
    for c,mut1 in enumerate(cds_mut_uni):
        split1 = mut1.split(':')
        mut = split1[1]
        #new_seq=[]
        if '+' in mut or '-' in mut or '*' in mut or '(' in mut or ')' in mut or '?' in mut:
            continue
        elif 'del' in mut:
            x = mut.split('del')
            if len(x[1]) > 5:
                continue
            elif "_" in x[0]:
                y = x[0].split('_')  
                i = int(y[0])
                j = int(y[1])
                #print(i,",",j)
                mut_list.append({'position':i,'org_n':x[1],'new_n':"", 'cds_mut':mut})
            else:
                i = int(x[0])
                mut_list.append({'position':i,'org_n':x[1],'new_n':"", 'cds_mut':mut})
        elif 'ins' in mut:
            x = mut.split('ins')
            aa_ins = x[1]
            y = x[0].split('_')  
            i = int(y[0])
            j = int(y[1])
            mut_list.append({'position':i,'org_n':"",'new_n':aa_ins, 'cds_mut':mut})
        elif '>' in mut:
            x = mut.split('>')
            aa2 = x[1] 
            if '_' in x[0]:
                y = x[0].split('_')  
                i = int(y[0])
                match = re.match(r"([0-9]+)([a-z]+)", y[1], re.I)
                if match:
                    items = match.groups()
                j = int(items[0])
                aa1 = items[1]
                #j = re.sub("[^0-9]","", y[1])
                #j = int(j)-1
                mut_list.append({'position':i,'org_n':aa1,'new_n':aa2, 'cds_mut':mut})
            else:
                ind = x[0]
                aa1 = ind[-1]
                ind = ind[:len(ind)-1]
                ind = int(ind)
                mut_list.append({'position':ind,'org_n':aa1,'new_n':aa2, 'cds_mut':mut})    
    org_seq = ""
    org_len=0
    for d in mut_list:
        org_seq = org_seq + d["org_n"]
    return mut_list, org_seq

mut_list_all=[]
org_seq_all=[]
for i in range(sites.shape[0]):
    gene = sites.iloc[i,0]
    cds_mut_uni = sites.iloc[i,1]
    mut_list, org_seq = create_dict(cds_mut_uni)
    mut_list_all.append(mut_list)
    org_seq_all.append(org_seq)

icgc_hotspot_donor_gene = icgc_data_hotspots.astype(str).groupby(['Donor_ID','Hugo_Symbol']).aggregate({'Genome_Change': list, 'Project_Code': set})
icgc_hotspot_donor_gene = icgc_hotspot_donor_gene.reset_index()

# add wildtype info for genes
donor_id_uni = list(set(icgc_hotspot_donor_gene['Donor_ID']))
for uni_id in donor_id_uni:
    donor_df = pd.DataFrame()
    donor_df = icgc_hotspot_donor_gene[icgc_hotspot_donor_gene['Donor_ID'] == uni_id]
    #print(donor_df)
    genes_donor = donor_df.iloc[:,1]
    genes_donor = genes_donor.values.tolist()
    #print(genes_donor)
    tissue = donor_df.iloc[0,3]
    for gene in eval_genes:
        if gene in genes_donor:
            print(gene)
        else:
            donor_df = donor_df.append({'Donor_ID': uni_id, 'Hugo_Symbol': gene, 'Genome_Change': ['none'], 'Project_Code': tissue}, ignore_index=True)
            donor_df = donor_df.sort_values(by=['Hugo_Symbol'])
    icgc_hotspot_donor_gene = icgc_hotspot_donor_gene.append(donor_df, ignore_index = True)
    icgc_hotspot_donor_gene = icgc_hotspot_donor_gene.drop_duplicates(subset=['Donor_ID', 'Hugo_Symbol'])

icgc_hotspot_donor_gene = icgc_hotspot_donor_gene.sort_values(by = ['Donor_ID', 'Hugo_Symbol'])

# create mut list and art sequence
gene_seq=[]
sample_id_seq = []
gene_name_seq = []
sample_tissue = []
for i in range(len(icgc_hotspot_donor_gene.iloc[:,:])):
    gene = icgc_hotspot_donor_gene.iloc[i,1]
    #print('Gene: ', gene)
    sampleid = icgc_hotspot_donor_gene.iloc[i,0]
    #print('ID: ', sampleid)
    tissue = icgc_hotspot_donor_gene.iloc[i,3]
    cds_mut = icgc_hotspot_donor_gene.iloc[i,2]
    index = eval_genes.index(gene)
    org_seq = org_seq_all[index]
    mut_list_gene = mut_list_all[index]
    cds_mut_gene_uni = list(set(cds_mut))
    if len(cds_mut_gene_uni) == 1:
        #print('len = 1')
        for mut1 in cds_mut_gene_uni:
            if mut1 == 'none':
                #print("len org: ", len(org_seq))
                gene_seq.append(org_seq)
                sample_id_seq.append(sampleid)
                gene_name_seq.append(gene)
                sample_tissue.append(tissue)
            else:
                split1 = mut1.split(':')
                mut = split1[1]
                i = next((index for (index, d) in enumerate(mut_list_gene) if d["cds_mut"] == mut), None)
                if i == None:
                    continue
                else:
                    seq3 = ""
                    trial1 = mut_list_gene[:i]
                    trial2 = mut_list_gene[i+1:]
                    for d in trial1:
                        seq3 = seq3 + d["org_n"]
                    seq3 = seq3 + mut_list_gene[i]['new_n']
                    for d in trial2:
                        seq3 = seq3 + d["org_n"]
                    gene_seq.append(seq3)
                    sample_id_seq.append(sampleid)
                    gene_name_seq.append(gene)
                    sample_tissue.append(tissue)
    else:
        index_list = []
        for mut1 in cds_mut_gene_uni:
            split1 = mut1.split(':')
            mut = split1[1]
            i = next((index for (index, d) in enumerate(mut_list_gene) if d["cds_mut"] == mut), None)
            if i == None:
                continue
            index_list.append(i)
            #print("index_list", index_list)
            index_list = sorted(index_list)
            #print("sorted index_list", index_list)
            seq3 = ""
            for ind in index_list:
                if seq3=="":
                    trial1 = mut_list_gene[:ind]
                    for d in trial1:
                        seq3 = seq3 + d["org_n"]
                    seq3 = seq3 + mut_list_gene[ind]['new_n']
                    update_ind = ind
                else:
                    trial1 = mut_list_gene[update_ind+1:ind]
                    for d in trial1:
                        seq3 = seq3 + d["org_n"]
                    seq3 = seq3 + mut_list_gene[ind]['new_n']
                    update_ind = ind
            trial1 = mut_list_gene[update_ind+1:]
            for d in trial1:
                seq3 = seq3 + d["org_n"]
        gene_seq.append(seq3)
        sample_id_seq.append(sampleid)
        gene_name_seq.append(gene)
        sample_tissue.append(tissue)
        
df_gene_seq1 = pd.DataFrame()
df_gene_seq1['sample id'] = sample_id_seq
df_gene_seq1['gene'] = gene_name_seq
df_gene_seq1['seq'] = gene_seq
df_gene_seq1['tissue'] = sample_tissue

seq_len_gene=[]
for seq in df_gene_seq1['seq']:
    seq_len_gene.append(len(seq))
df_gene_seq1['seq_len'] = seq_len_gene

df_gene_seq_comb1 = df_gene_seq1.groupby('sample id')['seq'].apply(lambda x: ''.join(x)).reset_index()
seq_len1=[]
for seq in df_gene_seq_comb1['seq']:
    seq_len1.append(len(seq))
df_gene_seq_comb1['seq_len'] = seq_len1

# wildtype sequence to be added to network
wildtype_seq=''
for seq in org_seq_all:
    wildtype_seq += seq
# make df with wildtype
wildtype_donor = ['D0']
wildtype_seq = [wildtype_seq]
wildtype_seq_len = list(len(wildtype_seq))

df_gene_seq_comb1_wildtype = pd.DataFrame()
df_gene_seq_comb1_wildtype['sample id'] = wildtype_donor
df_gene_seq_comb1_wildtype['seq'] = wildtype_seq
df_gene_seq_comb1_wildtype['seq_len'] = wildtype_seq_len

df_gene_seq_comb1_wildtype = df_gene_seq_comb1_wildtype.append(df_gene_seq_comb1)
df_gene_seq_comb1_wildtype_uni = df_gene_seq_comb1_wildtype.drop_duplicates(subset=['seq'])

#df_gene_seq_comb1_wildtype_uni.to_csv('df_gene_seq_comb1_wildtype_uni_293g.csv', index = None, sep=',')

# Bipartite network -------------------------------------------
tissue923_list = []
for t in donor_tissues923:
    t_list=[]
    t_list.append(t)
    tissue923_list.append(t_list)
    
from itertools import product
comb_list1=[]
for i in range(len(donor_tissues923)):
    l1 , l2 = tissue923_list[i] , donor_genes923[i]
    comb = list(product(l1,l2))
    comb_list1.append(comb)

edge_list1 = []
for lst in comb_list1:
    for item in lst:
        edge_list1.append(item)

tissue923_set = list(set(donor_tissues923))

import networkx as nx
from networkx.algorithms import bipartite

B = nx.Graph()
B.add_nodes_from(eval_genes, bipartite=1)
B.add_nodes_from(tissue923_set, bipartite=0)
B.add_edges_from(edge_list1)
edges = list(B.edges())
print(edges)
u = [n for n in B.nodes if B.nodes[n]['bipartite'] == 0]
v = [n for n in B.nodes if B.nodes[n]['bipartite'] == 1]

pos = dict()
pos.update( (n, (1, 3*i)) for i, n in enumerate(u) ) # put nodes from X at x=1
pos.update( (n, (2, 2*i)) for i, n in enumerate(v) ) # put nodes from Y at x=2
plt.figure(3,figsize=(40,40))
nx.draw_networkx(B, pos=pos, node_size = 100, width = 0.5, font_size = 8)
plt.show()

# Tissue projection
Bt = bipartite.weighted_projected_graph(B, u)
edgest = Bt.edges()
weights = [Bt[u][v]['weight'] for u,v in edgest]
weight = [w/10 for w in weights]
weight=[]
for w in weights:
    if w >= 30:
        weight.append(3.5)
    elif w >= 25:
        weight.append(2.5)
    elif w >= 15:
        weight.append(1.0)
    elif w >= 12 and w < 15:
        weight.append(0.8)
    elif w >= 9 and w < 12:
        weight.append(0.5)
    elif w >= 6 and w < 9:
        weight.append(0.1)
    elif w >= 3 and w < 6:
        weight.append(0)
    else:
        weight.append(0)
#post = nx.circular_layout(Bt)
post = nx.circular_layout(Bt)
plt.figure(3,figsize=(12,12))
nx.draw_networkx(Bt, pos = post, node_size = 1500, font_size = 20, width = weight)

# Gene projection
Bg = bipartite.weighted_projected_graph(B, v)
edgesg = Bg.edges()
weights_g = [Bg[u][v]['weight'] for u,v in edgesg]                                                                                                                                                            
weight=[]
for w in weights_g:
    if w >= 8:
        weight.append(2.5)
    elif w == 7:
        weight.append(1.0)
    elif w >= 5 and w<= 6:
        weight.append(0.3)
    elif w >= 3 and w <= 4:
        weight.append(0.1)
    elif w >= 2 and w < 3:
        weight.append(0.05)
    else:
        weight.append(0)
#weight = [w/10 for w in weights_g]
posg = nx.spring_layout(Bg, k=0.25)
plt.figure(3,figsize=(20,18))
nx.draw_networkx(Bg, pos = posg, node_size = 500, font_size = 12, width = weight)


