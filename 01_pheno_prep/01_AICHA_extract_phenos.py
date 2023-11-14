import pandas as pd
import numpy as np
import glob
import os
import time
import argparse

tic = time.perf_counter()


"""
#Commented section for running parallellizing

# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--total_parts', type=int,
                    help='Number of partitions of the data')
parser.add_argument('--no_part', type=int,
                    help='Which partition is this')
args = parser.parse_args()
args.no_part = args.no_part-1


#load QC'ed subs

#chunk into parts
chunk_size = int(np.ceil(len(subj_list)/args.total_parts))
chunked_list = [subj_list[j:j+chunk_size] for j in range(0, len(subj_list), chunk_size)]
subj_list = chunked_list[args.no_part]

print("Total number of subjects :", subs_no)
print("Chunk {0} out of {1}.".format(args.no_part+1, args.total_parts))
print("Number of subjects in chunk: ", len(subj_list))

"""

print('Loading stuff')
base_path = "/data/clusterfs/lag/users/jitame/SENT_CORE/"

#load QC'ed subs
s_list = open("/data/clusterfs/lag/users/jitame/SENT_CORE/subj_sent_N30652_exome_final_pass_sex.txt")
subj_list = sorted(list(s_list.read().split('\n')))[1:]

#get total no of subs
subs_no = len(subj_list)

## SET UP INDICES
aicha_len = 384

#SENT_core indices:
#specify regions of interest based on this: https://link.springer.com/article/10.1007/s00429-018-1810-2/tables/2: SENT_CORE
sent_core_l_ind = [2, 30, 32, 40, 56, 98, 102, 146, 148, 166, 168, 170, 172, 174, 182, 184, 222, 224]
sent_core_r_ind = [3, 31, 33, 41, 57, 99, 103, 147, 149, 167, 169, 171, 173, 175, 183, 185, 223, 225]
sent_core_bi_ind = sorted(sent_core_l_ind + sent_core_r_ind)

#GET LABELS
aicha_test = "/data/clusterfs/lag/projects/lg-ukbiobank/working_data/imaging_data/AICHA/1000099/AICHA_timeseries_NO_GSR_cormat.txt"
aicha2 = pd.read_csv(aicha_test, sep=";", index_col=0)
sent_node_names_l = list(aicha2.columns[sent_core_l_ind])
sent_node_names_r = list(aicha2.columns[sent_core_r_ind])
sent_node_names_bi = list(aicha2.columns[sent_core_bi_ind])
print(sent_node_names_bi)
sent_node_names_asym = [i[:-2] for i in sent_node_names_l]

def get_edge_names(list_names):
    df_names=pd.DataFrame(index=list_names, columns=list_names)
    for i in range(len(list_names)):
        for j in range(len(list_names)):
            df_names.iloc[i, j] = list_names[i]+"*"+list_names[j]
    return list(df_names.to_numpy()[np.triu_indices(len(list_names), k=1)].flatten())

print("Specifying empty dataframes")
#SPECIFY EMPTY DATAFRAMES
sent_edge_names = get_edge_names(sent_node_names_bi)
sent_edge_names_asym = get_edge_names(sent_node_names_asym)

#SENT CORE
#edges
df_sent_edges = pd.DataFrame(index=subj_list, columns=sent_edge_names)
df_sent_edges_asym = pd.DataFrame(index=subj_list, columns=sent_edge_names_asym)

if not os.path.exists(os.path.join(base_path, "pheno", "asym_mats_sent")): os.makedirs(os.path.join(base_path, "pheno", "asym_mats_sent"))

print("Start collecting nodes, edges and asymmetries")
i=0
for sub in subj_list: 
    if i%300 == 0:
        print("{0} % done".format((i/len(subj_list))*100))
    ## PER SUBJECT
    #load z-transformed correlation matrix
    aicha = np.loadtxt("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/multilateral/FuncNet_AICHA/mats_AICHA/fcmatz/{0}_fcmatz.txt".format(sub))

    #set negative values and diagonals to zero
    aicha[aicha < 0] = 0
    np.fill_diagonal(aicha, 0)

    ## SPECIFY NETWORKS
        #SENT CORE
    l_sent = aicha[np.ix_(sent_core_l_ind, sent_core_l_ind)]
    r_sent = aicha[np.ix_(sent_core_r_ind, sent_core_r_ind)]
    inter_sent = aicha[np.ix_(sent_core_l_ind, sent_core_r_ind)]
    bi_sent = aicha[np.ix_(sent_core_bi_ind, sent_core_bi_ind)]

    ## EDGE ASYMMETRIES
    #for edges we use L-R (hemispheric difference) -non-normalized as correlations are semi-normalized and normalizing may bias metrics
        #SENT CORE
    #calculate HD
    sent_edges_asym = l_sent - r_sent

    ## STORE DATA IN DF
    #save asymmetry matrices to folder
    np.savetxt(os.path.join(base_path, "pheno", "asym_mats_sent", "{0}_sent_hdmat.txt".format(sub)), sent_edges_asym)

        #SENT CORE
        #edges
    df_sent_edges.loc[sub] = bi_sent[np.triu_indices(len(sent_core_bi_ind), k=1)] #k=1 for skipping diagonal
    df_sent_edges_asym.loc[sub] = sent_edges_asym[np.triu_indices(len(sent_core_l_ind), k=1)]

    i=i+1

print("Saving all dataframes")
#SAVE ALL DFS
#edges
df_sent_edges.to_csv(os.path.join(base_path, "pheno", "sent_edges_N{0}.csv".format(subs_no)))
df_sent_edges_asym.to_csv(os.path.join(base_path, "pheno", "sent_edges_asym_N{0}.csv".format(subs_no)))
                          
toc = time.perf_counter()

print(f"Done in {(toc - tic)/60:0.4f} minutes")