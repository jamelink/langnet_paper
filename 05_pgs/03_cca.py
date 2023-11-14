import os
import re
import glob
from tqdm import tqdm

import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.cross_decomposition import CCA
from scipy.stats import pearsonr
from sklearn.preprocessing import quantile_transform
from numpy.random import default_rng

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from nilearn import plotting, image
from matplotlib import cm


## SETTINGS ##
n_perm = 10000
base_path = "/data/clusterfs/lag/users/jitame/SENT_CORE/"
results_path = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/rest-multimodal/results/"
plot_path = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/rest-multimodal/results"

fig_args = {"figsize":(18,10),"dpi":300}


## Helper functions

def read_pgs_chr(fn, chr_no):
    data = pd.read_csv(fn, sep="\s+", engine="python", na_values='-9')
    return data.rename(columns={"SCORESUM":"score_c{}".format(chr_no)})

def get_pgs(fn, pheno_name):
    file_list = [fn.format(pheno_name, chr_no) for chr_no in range(1,23)]
    data = pd.concat(map(read_pgs_chr, file_list, range(1,23)), join="inner", axis=1).T.drop_duplicates().T
    chrom_scores = ["score_c{0}".format(i) for i in range(1,23)]
    data["score_{0}".format(pheno_name)] = np.sum(data[chrom_scores].to_numpy(),axis=1)
    return data[["FID", "score_{0}".format(pheno_name)]]

def get_all_pgs(fn, pheno_list):
    return pd.concat(map(get_pgs, [fn]*len(pheno_list), pheno_list), join="inner", axis=1).T.drop_duplicates().T    

def load_column_names(fn):
    return [str(x) for x in open( fn ).read().split('\n')[:-1] ]

def residualize(data, covs, drop_all=True):
    #input_file = file_name + '.csv'
    #data = pd.read_csv(input_file)
    #data = data.set_index(data.columns[0]) 
    #data = data.loc[exome_subs]

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
    return data_new2

def add_row(df, new_data, column_names, index_name):
    df_new = pd.DataFrame(data=new_data.reshape((1, len(column_names))), columns=column_names, index=[index_name])
    return pd.concat([df, df_new], axis=0)

## CCA ##

def CCA_core(X,Y):
    #fit data
    ca = CCA(n_components=1)
    ca.fit(X, Y)
    X_c, Y_c = ca.transform(X, Y)
    
    #add to dataframes
    r, p = pearsonr(X_c[:, 0], Y_c[:, 0])
    loadings = ca.coef_
    return r, loadings

def greater(null_distribution, observed, n_resamples):
    """
    from scipy.stats.permutation_test
    https://github.com/scipy/scipy/blob/v1.11.1/scipy/stats/_resampling.py#L1234-L1707
    """
    cmps = null_distribution >= observed
    pvalues = (cmps.sum(axis=0) + 1) / (n_resamples + 1)  # see [1]
    return pvalues

def CCA_wrap(all_data, df_loadings, df_corrs, df_null_dist, pheno_prs, brain_phenos, n_perm, seed=2023):
    """
    """
    #get right data
    X = all_data[pheno_prs].to_numpy().reshape(-1, 1)
    Y = all_data[brain_phenos].to_numpy()
    
    rng = default_rng(seed)
    
    r_distribution = np.zeros(n_perm)
    
    for x in tqdm(range(n_perm)):
        X_perm = rng.permutation(X)
        r, _ = CCA_core(X_perm, Y)
        r_distribution[x] = r
        
    #final loadings
    r, loadings = CCA_core(X, Y)
    p = greater(r_distribution, r, n_perm)
    
    df_null_dist = add_row(df_null_dist, r_distribution, range(n_perm), pheno_prs)
    df_corrs = add_row(df_corrs, np.array([r, p]), "R P".split(), pheno_prs)
    df_loadings = add_row(df_loadings, np.array(loadings), brain_phenos, pheno_prs)
    
    return df_corrs, df_loadings, df_null_dist

## PLOTTING ##

def load_column_names(fn):
    return [str(x) for x in open( fn ).read().split('\n')[:-1] ]  

def get_edge_names(list_names):
    """
    Returns edge combinations for edges
    """
    df_names=pd.DataFrame(index=list_names, columns=list_names)
    for i in range(len(list_names)):
        for j in range(len(list_names)):
            df_names.iloc[i, j] = list_names[i]+"*"+list_names[j]
            
    #heritable = load_column_names("/data/clusterfs/lag/users/jitame/SENT_CORE/heritable_edges.txt") + load_column_names("/data/clusterfs/lag/users/jitame/SENT_CORE/heritable_edges_asym.txt")
    #heritable = ["sent_edges_"+x for x in heritable]
    
    return [x for x in list(df_names.to_numpy()[np.triu_indices(len(list_names), k=1)].flatten()) if x in heritable]


def get_idx(cat):
    """
    Returns indices for category
    """
    sent_core_l_ind = [2, 30, 32, 40, 56, 98, 102, 146, 148, 166, 168, 170, 172, 174, 182, 184, 222, 224]
    sent_core_r_ind = [3, 31, 33, 41, 57, 99, 103, 147, 149, 167, 169, 171, 173, 175, 183, 185, 223, 225]
    sent_core_bi_ind = sorted(sent_core_l_ind + sent_core_r_ind)
    
    if cat == "edges":
        return sent_core_bi_ind
    elif cat == "edges_HD":
        return sent_core_l_ind
    
    
def set_up_coordinates(cat):
    """
    get coordinates for plot
    """
    #set up atlas + coordinates
    aicha_atlas = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/multilateral/FuncNet_AICHA/segs/AICHA.nii"
    aicha_img = image.load_img(aicha_atlas)
    aicha_coords = plotting.find_parcellation_cut_coords(aicha_img)
    
    #indices
    idx = get_idx(cat)

    #get coordinates and names
    return aicha_coords[idx, :]

def get_names_aicha(cat):
    """
    Return node names
    """
    aicha_test = "/data/clusterfs/lag/projects/lg-ukbiobank/working_data/imaging_data/AICHA/1000099/AICHA_timeseries_NO_GSR_cormat.txt"
    aicha2 = pd.read_csv(aicha_test, sep=";", index_col=0)
    
    idx = get_idx(cat)
    
    if cat == "edges":
        return list(aicha2.columns[idx]) #[x for x in list(aicha2.columns[idx]) if x in heritable]
    elif cat == "edges_HD":
        return list(x[:-2] for x in aicha2.columns[idx]) #[x for x in list(aicha2.columns[idx]) if x in heritable]


def make_nice_mat(input_data, cat, var_name):
    """
    put data from column into right spot in matrix
    """
    #get indices
    ind_names = get_names_aicha(cat)
        
    #set up dataframe
    df_out = pd.DataFrame(index=ind_names, columns=ind_names)
    df_out.iloc[np.diag_indices(len(ind_names)), np.diag_indices(len(ind_names))] = 0
    
    #put data in right place
    for x in list(input_data.columns):
        nodes = x.split("*", 1)
        df_out.loc[nodes[0], nodes[1]] = input_data.loc[var_name, x]
        df_out.loc[nodes[1], nodes[0]] = input_data.loc[var_name, x]
    
    return df_out
        
    
def plot_results_brain_cca(data, cat, brain_pheno, beh_pheno, ax, out=None):
    """
    Main function:
    - loads data
    - significance testing
    - plot using nilearn
    """
    #get correlation values
    load_mat = make_nice_mat(data, cat, beh_pheno)  
    
    #make colormap
    colors = ["mediumblue", "cornflowerblue", "lightgrey", "whitesmoke", "lightgrey", "lightcoral", "indianred"]
    cmap1 = LinearSegmentedColormap.from_list("mycmap", colors)
    
    #plot
    if np.abs(data.sum().sum()) < 0.0001:
        plotting.plot_markers([1]*len(load_mat),
                              set_up_coordinates(cat), 
                              node_cmap="binary",
                              node_size=30,
                              alpha=0.8, 
                              display_mode='lyrz',
                              axes=ax, 
                              title=beh_pheno[6:], 
                              colorbar=False)
        
    else:    
        plotting.plot_connectome(load_mat.to_numpy(dtype=float),
                             set_up_coordinates(cat),
                             title=beh_pheno[6:],
                             #edge_cmap="coolwarm", 
                             edge_cmap=cmap1,
                             node_color="dimgrey",
                             node_size=30,
                             display_mode="lyrz",
                             colorbar=True,
                             alpha=0.8,
                             axes=ax)
                             #output_file=out)    

def load_column_names(fn):
    return [str(x) for x in open( fn ).read().split('\n')[:-1] ]  

def get_edge_names(list_names):
    """
    Returns edge combinations for edges
    """
    df_names=pd.DataFrame(index=list_names, columns=list_names)
    for i in range(len(list_names)):
        for j in range(len(list_names)):
            df_names.iloc[i, j] = list_names[i]+"*"+list_names[j]
            
    heritable = load_column_names("/data/clusterfs/lag/users/jitame/SENT_CORE/heritable_edges.txt") + load_column_names("/data/clusterfs/lag/users/jitame/SENT_CORE/heritable_edges_asym.txt")
    
    return [x for x in list(df_names.to_numpy()[np.triu_indices(len(list_names), k=1)].flatten()) if x in heritable]

def get_names_aicha(cat):
    """
    Return node names
    """
    aicha_test = "/data/clusterfs/lag/projects/lg-ukbiobank/working_data/imaging_data/AICHA/1000099/AICHA_timeseries_NO_GSR_cormat.txt"
    aicha2 = pd.read_csv(aicha_test, sep=";", index_col=0)
    
    idx = get_idx(cat)
    
    if cat == "edges":
        return list(aicha2.columns[idx]) #[x for x in list(aicha2.columns[idx]) if x in heritable]
    elif cat == "edges_HD":
        return [x[:-2] for x in list(aicha2.columns[idx])] 
    
def get_hd_name(edge_name):
    print(edge_name)
    edge_split = edge_name.split("*")
    return edge_split[0][:-2]+"*"+edge_split[1][:-2]


## PROCESS PGS ## 
print("Process PGS")
#get pgs
pgs = get_all_pgs(fn="/data/clusterfs/lag/users/jitame/SENT_CORE/geno/polygenic-scores/prs_out/{0}/{0}_prs_chr{1}.profile",
                  pheno_list=["read", "dyslexia", "hand"])
pgs["FID"] = pgs["FID"].astype("int")
pgs.set_index("FID", inplace=True)
pgs.to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/geno/polygenic-scores/all_scores_uncor.csv")

# get covariates
covs=pd.read_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/covars/covars_pc10_gwas.tsv", sep="\t")
covs["FID"] = covs["FID"].astype("int")
covs.set_index("FID", inplace=True)
covs.drop(["IID"], axis=1)
covs=covs.dropna()

# normalize PGS
normalized_pgs = residualize(pgs, covs)
normalized_pgs.to_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/geno/polygenic-scores/all_scores_norm.csv")

print("Run CCA Language Network")
#load phenotypes
brain_phenos = pd.read_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/pheno/sent_edges_N29681_resid_norm.txt",
                     header=None,
                     sep="\t",
                     names=load_column_names("/data/clusterfs/lag/users/jitame/SENT_CORE/pheno/sent_edges_N_resid_col_names.txt")
                    )
brain_phenos.set_index("Family_ID", inplace=True)
brain_phenos.drop("Subject_ID",axis=1, inplace=True)

her_names = load_column_names(os.path.join("/data/clusterfs/lag/users/jitame/SENT_CORE", "heritable_edges.txt"))

print(brain_phenos.shape)
brain_phenos = brain_phenos[her_names]
print(brain_phenos.shape)

#set up stuff for CCA
all_data = brain_phenos.join(normalized_pgs)

brain_pheno_names = list(brain_phenos.columns)
prs_pheno_names = list(normalized_pgs.columns)   

df_null_dist = pd.DataFrame(columns=range(n_perm))
df_corrs = pd.DataFrame(columns="R P".split())
df_loadings = pd.DataFrame(columns=brain_pheno_names)

#run CCA
for prs_pheno in prs_pheno_names:
    df_corrs, df_loadings, df_null_dist = CCA_wrap(all_data,
                                     df_loadings,
                                     df_corrs,
                                     df_null_dist,
                                     prs_pheno,
                                     brain_pheno_names,
                                     n_perm)


print("Run CCA asymmetries")
#load asymmetries
brain_phenos_asym = pd.read_csv("/data/clusterfs/lag/users/jitame/SENT_CORE/pheno/sent_edges_asym_N29681_resid_norm.txt",
                     header=None,
                     sep="\t",
                     names=load_column_names("/data/clusterfs/lag/users/jitame/SENT_CORE/pheno/sent_edges_asym_N_resid_col_names.txt")
                    )
brain_phenos_asym.set_index("Family_ID", inplace=True)
brain_phenos_asym.drop("Subject_ID",axis=1, inplace=True)

her_names_asym = load_column_names(os.path.join("/data/clusterfs/lag/users/jitame/SENT_CORE", "heritable_edges_asym.txt"))

print(brain_phenos_asym.shape)
brain_phenos_asym = brain_phenos_asym[her_names_asym]
print(brain_phenos_asym.shape)

#set up stuff
all_data = brain_phenos_asym.join(normalized_pgs)

brain_pheno_names = list(brain_phenos_asym.columns)
prs_pheno_names = list(normalized_pgs.columns)   

df_null_dist_asym = pd.DataFrame(columns=range(n_perm))
df_corrs_asym = pd.DataFrame(columns="R P".split())
df_loadings_asym = pd.DataFrame(columns=brain_pheno_names)

#run magic
for prs_pheno in prs_pheno_names:
    df_corrs_asym, df_loadings_asym, df_null_dist_asym = CCA_wrap(all_data,
             df_loadings_asym,
             df_corrs_asym,
             df_null_dist_asym,
             prs_pheno,
             brain_pheno_names,
             n_perm)

print("Save files")
df_corrs_asym.to_csv(os.path.join(base_path, "geno", "CCA_corr_asym.csv"))
df_loadings_asym.to_csv(os.path.join(base_path, "geno", "CCA_loadings_asym.csv"))
df_null_dist_asym.to_csv(os.path.join(base_path, "geno", "CCA_null_dist_asym.csv"))
df_corrs.to_csv(os.path.join(base_path, "geno", "CCA_corr.csv"))
df_loadings.to_csv(os.path.join(base_path, "geno", "CCA_loadings.csv"))
df_null_dist.to_csv(os.path.join(base_path, "geno", "CCA_null_dist.csv"))
    

print("Making density plots null distributions")
print(df_corrs)
sns.kdeplot(df_null_dist.T, bw_method=0.25)
print(df_corrs_asym)
sns.kdeplot(df_null_dist_asym.T, bw_method=0.25)
plt.savefig(fname=os.path.join(plot_path, "CCA_null_distributions.png"), bbox_inches="tight")


print("Make main plot")

fig, ax = plt.subplots(len(prs_pheno_names), 2, gridspec_kw={'height_ratios': [1]*len(prs_pheno_names)}, **fig_args)

for x, pheno_name in enumerate(prs_pheno_names):
    plot_results_brain_cca(data=df_loadings,
                           cat="edges",
                           brain_pheno="edge",
                           beh_pheno=pheno_name,
                           ax=ax[x, 0],
                           out=None)
    ax[x, 0].text(0.2, 0.94, "R={:.3f}, P={:.3e}".format(df_corrs.loc[pheno_name, "R"], df_corrs.loc[pheno_name, "P"]), fontsize=12)

for x, pheno_name in enumerate(prs_pheno_names):
    plot_results_brain_cca(data=df_loadings_asym,
                           cat="edges_HD",
                           brain_pheno="edge hemispheric difference",
                           beh_pheno=pheno_name,
                           ax=ax[x, 1],
                           out=None)
    
    ax[x, 1].text(0.2, 0.94, "R={:.3f}, P={:.3e}".format(df_corrs_asym.loc[pheno_name, "R"], df_corrs_asym.loc[pheno_name, "P"]), fontsize=12)

fig.suptitle("Multivariate associations of polygenic scores with brain phenotypes", y=0.98, fontsize=24)
ax[0, 0].text(0, 1.13, "A. language network", fontsize=16)
ax[0, 1].text(0, 1.13, "B. hemispheric differences", fontsize=16)

plt.savefig(fname=os.path.join(plot_path, "CCA_results.png"), bbox_inches="tight")