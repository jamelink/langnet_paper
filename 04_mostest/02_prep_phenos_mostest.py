import os
import pandas as pd

"""
Preps pheno file for MOSTEST:
- keep only heritable columns
- put subs in same order as FAM-file
- also separately for pseudo-autosomal region

Written by: J. S. Amelink
Date: March 8, 2023
(based on earlier code)
"""

def load_column_names(fn):
    return [str(x) for x in open( fn ).read().split('\n')[0:-1] ]

def prep_data_frame(fname, fname_col_in, her_cols, idx, out, xy=False):
    print("File: ", fname)
    if xy:
        print("Pseudoautosomal regions")
        
    #set column names
    cols = ['subid', 'famid'] + list(pd.read_csv(fname_col_in).columns[1:])
    heritable_cols_to_keep = load_column_names(her_cols)
    print("Heritable columns: ", len(heritable_cols_to_keep))
    
    #reorganize dataframe
    df = pd.read_csv(fname+".txt",sep="\t", names=cols)
    
    #filter heritable columns because more genetic signal
    df = df[['subid', 'famid'] + heritable_cols_to_keep]
    
    #reindex   
    df.set_index(df.columns[0], inplace=True)
    df = df.iloc[:, 1:]
    df = df.reindex(idx)

    print(df.shape)
    if xy:
        df.to_csv(out+"_mostest_XY.txt", sep="\t", index=False, header=True)
    else:
        df.to_csv(out+"_mostest.txt", sep="\t", index=False, header=True)

        
## RUN STUFF
base_path = "/data/clusterfs/lag/users/jitame/SENT_CORE/"

subs_no_in=29682
subs_no_out=29681

fnames = [os.path.join(base_path, "pheno", "sent_edges_N{0}_resid".format(subs_no_out)),
         os.path.join(base_path, "pheno", "sent_edges_asym_N{0}_resid".format(subs_no_out))]

fnames_col_in = [os.path.join(base_path, "pheno", "sent_edges_N{0}.csv".format(subs_no_in)),
         os.path.join(base_path, "pheno", "sent_edges_asym_N{0}.csv".format(subs_no_in))]

heritable_cols = [os.path.join(base_path, "heritable_edges.txt"),
                 os.path.join(base_path, "heritable_edges_asym.txt")]

#load order subjects
fam = pd.read_csv(os.path.join(base_path, "geno", "mostest", "in", "mostest_geno_in_c14.fam"), sep="\t", header=None)
fam_idx = list(fam.iloc[:, 0])

fam_XY = pd.read_csv(os.path.join(base_path, "geno", "mostest", "in", "mostest_geno_in_cXY.fam"), sep="\t", header=None)
fam_idx_XY = list(fam_XY.iloc[:, 0])


#for fname, fname_col_in, her_cols  in zip(fnames, fnames_col_in, heritable_cols):
#    prep_data_frame(fname, fname_col_in, her_cols, fam_idx_XY, xy=True)
#    prep_data_frame(fname, fname_col_in, her_cols, fam_idx)
    
sub_cats = ["heritable_edges_L", "heritable_edges_R", "heritable_edges_inter", "heritable_edges_intra"]

for sub_cat  in sub_cats:
    her_cols = os.path.join(base_path, sub_cat+".txt")
    prep_data_frame(fnames[0], fnames_col_in[0], her_cols, fam_idx_XY, fnames[0]+"_"+sub_cat, xy=True)
    prep_data_frame(fnames[0], fnames_col_in[0], her_cols, fam_idx, fnames[0]+"_"+sub_cat)