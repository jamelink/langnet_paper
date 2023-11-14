import pandas as pd
import os
import numpy as np
import time

tic = time.perf_counter()

"""
This file transforms a tab file from UKB to covars used for this study.

Covars included:
- Sex
- Age
- Age^2
- Sex*age
- Site
- Genotype array
- Exome array
- MR variables:
    - head motion rs
    - inv_tSNR_rs
    - X-Y-Z and table position
- 10 Genetic PCs

For running:

python_path="/home/jitame/bin/envs/std_work_env/bin/python"
file_name="/home/jitame/bin/code/AICHA/01_pheno_prep/02_generate_covars.py"
$python_path $file_name

"""

with open("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars.txt") as f:
    firstline = f.readline().rstrip()
header_covs = firstline.split()

#create dictionary with fieldnames. Also handy as reference/overview
field_names = dict([('f.eid', 'subject_id'),
                    ('f.31.0.0', 'sex'),
                    ('f.34.0.0', 'year-of-birth'),
                    ('f.52.0.0', 'month-of-birth'),
                    ('f.53.2.0', 'date-site' ),
                    ('f.54.2.0', 'site'),
                    ('f.25741.2.0', 'MR_head_motion_rs'),
                    ('f.25743.2.0', 'MR_inv_tSNR_rs'),
                    ('f.25756.2.0', 'MR_X_brain_pos'),
                    ('f.25757.2.0', 'MR_Y_brain_pos'),
                    ('f.25758.2.0', 'MR_Z_brain_pos'),
                    ('f.25759.2.0', 'MR_table_pos'),
                    ('f.22008.0.0', 'geno_well'),
                    ('f.22000.0.0', 'geno_batch')])
field_names_2=dict([('f.22009.0.{0}'.format(i), 'geno_PC_{0}'.format(i)) for i in range(1, 41)])
field_names = field_names | field_names_2

print("Load data")

#add to dataframe
covs = pd.read_csv('/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_rs_N30660.txt', sep=" ", header=None)
covs = covs.drop(covs.columns[[0]], axis=1)
covs.columns = header_covs
covs.rename(columns=field_names, inplace=True)
covs = covs.set_index('subject_id')

print("Convert age and create genotype array and site dummies")
#compute interactions
def get_age(born_year, born_month, visit):
    visit_date=str(visit).split('-')
    return float(visit_date[0])-float(born_year) - ((float(born_month) - float(visit_date[1]))/12)

covs['age'] = [get_age(covs['year-of-birth'][x], covs['month-of-birth'][x], covs['date-site'][x]) for x in covs.index.values]
covs['age_sq'] = np.square(covs['age'])
covs['age_sex'] = covs['age']*covs['sex']

#make dummies
covs['geno_array_dummy'] = 1
covs.loc[covs['geno_batch'] < 0, 'geno_array_dummy'] = 0

site_dummies = pd.get_dummies(covs['site'])
site_dummies.columns = ["site_dummy_{0}".format(x) for x in site_dummies.columns]

exome_batch = pd.read_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/exome_batch_N30660.txt",
                          sep=" ",
                         skiprows=1,
                         names=["exome_batch", "subject_id"])
exome_batch.sort_values(by="subject_id", inplace=True)
exome_batch.set_index("subject_id", inplace=True)

print("Save complete file")

#add together and save
covs = pd.concat([covs, site_dummies], axis=1)
covs = pd.concat([covs, exome_batch], axis=1)

covs.to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_all_N{0}.txt".format(len(covs)), sep="\t")

print("Keep only relevant files")
covs_qual_for_gen = covs[['sex', 'geno_array_dummy', 'site']]
covs_qual_for_gen['Family_ID'] = covs_qual_for_gen.index.values.astype(int)
covs_qual_for_gen['Subject_ID'] = covs_qual_for_gen.index.values.astype(int)
covs_qual_for_gen = covs_qual_for_gen[['Family_ID', 'Subject_ID', 'sex', 'geno_array_dummy', 'site']]
covs_qual_for_gen.to_csv('/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_trait_gen_N{0}.txt'.format(len(covs_qual_for_gen)), sep="\t", index=False, header=False)

#get rid of vars that can't be used for residualization
lose_vars=['year-of-birth','month-of-birth','date-site','site','geno_batch', 'geno_well']
lose_vars2=["geno_PC_{0}".format(i) for i in range(11, 41)]
lose_vars3=["exome_batch"]
lose_vars4=['geno_array_dummy']

print("Save files")
covs.drop(lose_vars, axis=1, inplace=True)
covs.to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_pc40_for_correction_N{0}.txt".format(len(covs)), sep="\t")
covs.drop(lose_vars2, axis=1, inplace=True)
covs.to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_pc10_for_correction_N{0}_exome_gwas_batch.txt".format(len(covs)), sep="\t")
covs.drop(lose_vars3, axis=1).to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_pc10_for_correction_N{0}_gwas_batch.txt".format(len(covs)), sep="\t")
covs.drop(lose_vars4, axis=1).to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_pc10_for_correction_N{0}_exome_batch.txt".format(len(covs)), sep="\t")

toc = time.perf_counter()

print(f"Done in {(toc - tic)/60:0.4f} minutes")