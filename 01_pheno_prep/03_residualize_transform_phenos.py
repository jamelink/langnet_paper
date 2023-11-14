import statsmodels.api as sm
from sklearn.preprocessing import quantile_transform
import os
import pandas as pd
import numpy as np
import time

tic = time.perf_counter()

base_path = "/data/clusterfs/lag/users/jitame/SENT_CORE/"
subs_no = str(29682)

"""
Residualize all files from (2) with covars from (3)

python_path="/home/jitame/bin/envs/std_work_env/bin/python"
file_name="/home/jitame/bin/code/AICHA/01_pheno_prep/03_residualize_transform_phenos.py"
$python_path $file_name

"""
print("Load covariates")
#set up covariates and impute data
covs=pd.read_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_pc10_for_correction_N30660_gwas_batch.txt", sep="\t")
covs=covs.set_index('subject_id')
covs=covs.sort_index()
covs=covs.dropna()

print("Load subjects")
exome_sub_list = os.path.join(base_path, "subj_sent_N30652_exome_final_pass_sex.txt")
exome_subs = [int(x) for x in open( exome_sub_list ).read().split('\n')[:-1] ]

def residualize(file_name, covs, exome_subs, drop_all=True):
    input_file = file_name + '.csv'
    data = pd.read_csv(input_file)
    data = data.set_index(data.columns[0]) 
    data = data.loc[exome_subs]

    if drop_all:
        data = data.dropna()
        data = data.sort_index()
        #data[data.isnull()] = 0
        missing_sub = list((set(list(covs.index.values)).difference(list(data.index.values))))
        print("No. subjects missing from subject file: ", len(missing_sub))
        covs=covs.loc[data.index.values]

        #define new dataframe
        data_new=pd.DataFrame(columns=data.columns, index=data.index.values)

        #residualize
        for dep_var in data.columns: 
            model = sm.OLS(data[dep_var], exog=covs)
            results = model.fit()
            data_new[dep_var] = results.resid
    
            
    #quantile transformation
    X = data_new.to_numpy()
    data_new2 = pd.DataFrame(data=quantile_transform(X, n_quantiles=1000, output_distribution='normal', random_state=0, copy=True), columns=data_new.columns, index=data_new.index.values)

    #reorder and save
    initial_cols = data_new.columns
    data_new['Subject_ID'] = data.index.values.astype(int)
    data_new['Family_ID'] = data.index.values.astype(int)
    data_new = data_new[['Subject_ID', 'Family_ID', *initial_cols]]
    data_new.to_csv(file_name[:-5] + '{0}_resid.txt'.format(len(data_new)), na_rep="NA", sep="\t", index=False, header=False)
    data_new2 = pd.concat([data_new[['Family_ID', 'Subject_ID']], data_new2], axis=1)
    data_new2.to_csv(file_name[:-5] + '{0}_resid_norm.txt'.format(len(data_new2)), na_rep="NA", sep="\t", index=False, header=False)
    list_w = list(data_new2.columns)
    with open(file_name[:-5]+"_resid_col_names.txt", "w") as file:
        for row in list_w:
            file.write(str(row)+'\n')

#files:
f_names=[os.path.join(base_path, "pheno", "sent_edges_N{0}".format(subs_no)),
         os.path.join(base_path, "pheno", "sent_edges_asym_N{0}".format(subs_no))]


#for-loop where the magic happens
run_all = True
if run_all:
    for fn in f_names:
        print("Input: ", fn)
        try:
            residualize(fn, covs, exome_subs, drop_all=True)
        except ValueError:
            print("Value Error, skipping this one")
            continue

toc = time.perf_counter()

print(f"Done in {(toc - tic)/60:0.4f} minutes")