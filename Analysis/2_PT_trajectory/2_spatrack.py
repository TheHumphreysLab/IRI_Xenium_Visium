import spaTrack as spt
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import warnings
import os
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
plt.rcParams['figure.dpi'] = 100 

adata = sc.read_h5ad("pt_ad_BBKNN.h5ad")

### subsample cells 
def obs_key_wise_subsampling(adata, obs_key, N):

    counts = adata.obs[obs_key].value_counts()
    # subsample indices per group defined by obs_key
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=False) for group in counts.index]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()

adata = obs_key_wise_subsampling(adata, "ident", 5000)
adata.obs = adata.obs[['x_centroid', 'y_centroid','celltype_bbknn','ident']]
adata.obs['cell_id'] = adata.obs.index
adata.obs.columns = ['x','y','annotation', 'time','cell_id']

adata1_sub=adata[adata.obs['time']=='ShamR']
adata2_sub=adata[adata.obs['time']=='Hour4R']
adata3_sub=adata[adata.obs['time']=='Hour12R']
adata4_sub=adata[adata.obs['time']=='Day2R']
adata5_sub=adata[adata.obs['time']=='Day14R']
adata6_sub=adata[adata.obs['time']=='Week6R']

import pandas as pd
import itertools

adata_list = [adata1_sub, adata2_sub, adata3_sub, adata4_sub, adata5_sub, adata6_sub]
sample_ids = ['ShamR', 'Hour4R', 'Hour12R', 'Day2R', 'Day14R', 'Week6R']
epsilon_list = [0.01, 0.01, 0.01, 0.01, 0.002]
alpha_list = [0.3, 0.5, 0.5, 0.1, 0.1]

results = []
pi_list = []
for i in range(len(adata_list) - 1):
    
    adata_start = adata_list[i]
    adata_end = adata_list[i + 1]
    
    sample_start_id = sample_ids[i]
    sample_end_id = sample_ids[i + 1]
    print(f'calculating transition between {sample_start_id} and {sample_end_id}')

    pi = spt.transfer_matrix(
        adata_start, 
        adata_end, 
        epsilon=epsilon_list[i], 
        alpha=alpha_list[i],
    )
    pi_list.append(pi)
    # Generate coordinate information for Sankey input
    pi_matrix_coord = spt.generate_animate_input(
        pi_list=[pi], 
        adata_list=[adata_start, adata_end], 
        time='time', 
        annotation='annotation'
    )
    
    # Create the transition matrix
    transition_matrix = pd.crosstab(
        pi_matrix_coord['slice1_annotation'], 
        pi_matrix_coord['slice2_annotation']
    )
    
    # Store results with corresponding epsilon, alpha, and the sample pair IDs
    results.append({
        'sample_start_id': sample_start_id,
        'sample_end_id': sample_end_id,
        'transition_matrix': transition_matrix
    })
    df = transition_matrix.div(transition_matrix.sum(axis=1), axis=0)
    file_name = os.path.join('spatrack_result',f'Xenium_5000_{sample_ids[i]}_{sample_ids[i+1]}_TM.csv')
    print(file_name)
    transition_matrix.to_csv(file_name)
    

## Get the transition between day2 and day14 data
pi_matrix_coord = spt.generate_animate_input(
        pi_list=[pi_list[3]], 
        adata_list=[adata4_sub, adata5_sub], 
        time='time', 
        annotation='annotation'
    )

## Define Day2 PT state based on destination in day14 
inj_FRPT = pi_matrix_coord[pi_matrix_coord.slice1_annotation.isin(['Inj_S1','Inj_S2','Inj_S3']) & (pi_matrix_coord.slice2_annotation=="Failed_repair")].slice1.values
inj_healthy = pi_matrix_coord[pi_matrix_coord.slice1_annotation.isin(['Inj_S1','Inj_S2','Inj_S3']) & (pi_matrix_coord.slice2_annotation.isin(['Healthy_S1','Healthy_S2','Healthy_S3']))].slice1.values
healthy = pi_matrix_coord[pi_matrix_coord.slice1_annotation.isin(['Healthy_S1','Healthy_S2','Healthy_S3']) & (pi_matrix_coord.slice2_annotation.isin(['Healthy_S1','Healthy_S2','Healthy_S3']))].slice1.values

adata4_sub.obs['Inj_des'] = 'Other'
adata4_sub.obs.loc[adata4_sub.obs.index.isin(inj_FRPT), 'Inj_des'] = 'Maladaptive injured PT'
adata4_sub.obs.loc[adata4_sub.obs.index.isin(inj_healthy), 'Inj_des'] = 'Recovering injured PT'
adata4_sub.obs.loc[adata4_sub.obs.index.isin(healthy), 'Inj_des'] = 'Stable healthy PT'
adata4_sub = adata4_sub[adata4_sub.obs['Inj_des'] != 'Other']

adata4_sub.obs['Inj_des'] = adata4_sub.obs['Inj_des'].cat.reorder_categories(['Stable healthy PT','Recovering injured PT','Maladaptive injured PT']) 
custom_colors = ["#b8d4a2", "#fbad75", "#CDB4DB"]  
adata4_sub.uns["Inj_des_colors"] = custom_colors  

## Figure 3D
with rc_context({"figure.figsize": (5.5, 5), "figure.dpi":300}):
    ax = sc.pl.scatter(
        adata4_sub, basis="pca", color=["Inj_des"], components=[1,2], size = 50, frameon=False, show=False,
    )
    plt.savefig("Day2PT_fate_PC.png", bbox_inches = "tight", transparent=True)
    plt.show()
 
## Figure 3E
with rc_context({"figure.figsize": (3,2), 'figure.dpi':300}):
    ax = sc.pl.violin(adata4_sub, ['Krt20','Il34','Vcam1','Serpine1','Cd44','Klf5','Haao','Hmgcs2','Cxcl12'], stripplot=False, show=False,
                 groupby="Inj_des", use_raw=False, jitter=False)
    for a in ax:
        a.set_xticklabels([])
        a.set_xlabel('')
        for sp in ['top', 'right', 'bottom']:
            a.spines[sp].set_visible(False)
        
        a.set_yticks([0, 1, 2])
        a.set_yticklabels([0, 1, 2])
    
    plt.savefig("iPT_day2_compare.png", bbox_inches = "tight", transparent=True)
    plt.show()
    
## Supplementary Figure 5   
import matplotlib.colors as mcolors

cmap = mcolors.LinearSegmentedColormap.from_list('WhRd',['#eeeeee', "#fffacd", "red", "darkred"], N=256) 
with rc_context({"figure.figsize": (5, 4)}):
    ax = sc.pl.embedding(
            adata4_sub, basis='pca', color=["Inj_des",'Havcr1','Vcam1','Cxcl12','Haao','Hmgcs2','Kynu','Ppara',
                                           'Serpine1','Cd44','Klf5','Gprc5a','Runx1'], frameon=False, size=50, 
                legend_loc=False, color_map=cmap, show = False
        )
    plt.savefig("marker_umap.png", bbox_inches = "tight", transparent=True)
    plt.show()    