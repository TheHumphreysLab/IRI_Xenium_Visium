import scanpy as sc
import pandas as pd
import numpy as np

from scipy.sparse import csr_matrix, coo_matrix
from sklearn.neighbors import radius_neighbors_graph
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
    
def permute(adata, groupby = "celltype_raw", 
            n_permutations = 100, 
            niche_radius = 15.0,
            permute_radius = 100.0, 
            spatial_key = "spatial",
            seed = 123):
    '''
    Permutation test: Randomizing the actual spatial locations of the cells, then computed the proximity 
    between all possible pairs of cell types under this randomized spatial arrangement.
    '''
    categories_str_cat = list(adata.obs[groupby].cat.categories)
    categories_num_cat = range(len(categories_str_cat))
    map_dict = dict(zip(categories_num_cat, categories_str_cat))

    categories_str = adata.obs[groupby]
    categories_num = adata.obs[groupby].replace(categories_str_cat, categories_num_cat)
    labels = categories_num.to_numpy()
    print(len(labels))
    
    max_index = len(categories_num_cat)
    pair_counts = np.zeros((max_index, max_index))
    
    ## calculate the true interactions
    con = radius_neighbors_graph(adata.obsm[spatial_key], niche_radius,  mode='connectivity', include_self=False)
    con_coo = coo_matrix(con)
    print(f"Calculating observed interactions...")
    for i, j, val in zip(con_coo.row, con_coo.col, con_coo.data):
        if i >= j:  
            continue

        type1 = labels[i]
        type2 = labels[j]

        if val:  
            if type1 != type2:
                pair_counts[type1, type2] += 1
                pair_counts[type2, type1] += 1
            else:
                pair_counts[type1, type2] += 1
                
    #pair_counts = pair_counts/pair_counts.sum()  ## normalize by total counts                        
    ## calculate the null hypothesis
    coords = adata.obsm[spatial_key]

    pair_counts_null = np.zeros((max_index, max_index, n_permutations))

    for perm in range(n_permutations):
        np.random.seed(seed=None if seed is None else seed + perm)
        
        permuted_coords = coords + np.random.uniform(-permute_radius, permute_radius, size=coords.shape)
        permuted_con = radius_neighbors_graph(permuted_coords, niche_radius, mode='connectivity', include_self=False)
        permuted_con_coo = coo_matrix(permuted_con)
        
        if perm%100==1:
            print(f"Performing permutation {perm}...")

        pair_counts_permuted = np.zeros((max_index, max_index))

        for i, j, val in zip(permuted_con_coo.row, permuted_con_coo.col, permuted_con_coo.data):
            if i >= j:  
                continue

            type1 = labels[i]
            type2 = labels[j]

            if val:  
                if type1 != type2:
                    pair_counts_permuted[type1, type2] += 1
                    pair_counts_permuted[type2, type1] += 1
                else:
                    pair_counts_permuted[type1, type2] += 1
        
        #pair_counts_permuted = pair_counts_permuted/pair_counts_permuted.sum()  ## normalize by total counts 
        pair_counts_null[:, :, perm] = pair_counts_permuted
                
    return categories_str_cat, pair_counts, pair_counts_null


adata = sc.read_h5ad("../Xenium_all.h5ad")

for ident in adata.obs['ident'].unique():
    ad = adata[adata.obs["ident"] == ident].copy()
    cell_cnt = pd.DataFrame(ad.obs["celltype"].value_counts()).reset_index()
    cell_cnt.columns = ['CellType', "Count"]

    celltypes, pair_counts, pair_counts_null = permute(ad, "celltype", n_permutations = 1000, niche_radius = 27.5, permute_radius = 50)
    pair_counts_permuted_means = np.mean(pair_counts_null, axis=2)
    pair_counts_permuted_stds = np.std(pair_counts_null, axis=2)

    z_scores = (pair_counts - pair_counts_permuted_means) / pair_counts_permuted_stds
    z_scores[z_scores<0]=0
    z_score_df = pd.DataFrame(z_scores,index = celltypes, columns = celltypes)
    z_score_df.to_csv(f"{ident}_zscores_d55_perm.csv")  
