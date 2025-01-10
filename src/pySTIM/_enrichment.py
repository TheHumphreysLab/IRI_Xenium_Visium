import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from sklearn.neighbors import radius_neighbors_graph
from anndata import AnnData
import requests
import json
import time
import seaborn as sns
sns.set_style("ticks")
from matplotlib.ticker import MaxNLocator

def calculate_enrichment(
    adata: "AnnData", 
    groupby: str = "celltype", 
    n_permutations: int = 100, 
    niche_radius: float = 15.0,
    permute_radius: float = 50.0, 
    spatial_key: str = "spatial",
    seed: int = 123
) -> pd.DataFrame:
    """
    Permutation test to calculate proximity enrichment frequency between cell type pairs.

    Parameters
    ----------
    adata : AnnData
        An anndata object containing the data.
    groupby : str, optional
        Column name to group by, by default "celltype".
    n_permutations : int, optional
        Number of permutations for the test, by default 100.
    niche_radius : float, optional
        Radius to define the niche, by default 15.0.
    permute_radius : float, optional
        Radius for permutation, by default 50.0.
    spatial_key : str, optional
        Key for accessing spatial coordinates in `adata.obsm`, by default "spatial".
    seed : int, optional
        Random seed for reproducibility, by default 123.

    Returns
    -------
    pd.DataFrame
        DataFrame containing z-scores for cell type enrichment.
    """
    categories_str_cat = list(adata.obs[groupby].cat.categories)
    categories_num_cat = range(len(categories_str_cat))
    map_dict = dict(zip(categories_num_cat, categories_str_cat))

    categories_str = adata.obs[groupby]
    categories_num = adata.obs[groupby].replace(categories_str_cat, categories_num_cat)
    labels = categories_num.to_numpy()
    print(f"Total number of {len(labels)} cells")
    
    max_index = len(categories_num_cat)
    pair_counts = np.zeros((max_index, max_index))
    
    # Calculate the true interactions
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
                
    coords = adata.obsm[spatial_key]

    pair_counts_null = np.zeros((max_index, max_index, n_permutations))

    for perm in range(n_permutations):
        np.random.seed(seed=None if seed is None else seed + perm)
        
        permuted_coords = coords + np.random.uniform(-permute_radius, permute_radius, size=coords.shape)
        permuted_con = radius_neighbors_graph(permuted_coords, niche_radius, mode='connectivity', include_self=False)
        permuted_con_coo = coo_matrix(permuted_con)
        
        if (perm + 1) % 100 == 1:
            print(f"Permutation iterations: {perm}...")

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
        
        pair_counts_null[:, :, perm] = pair_counts_permuted

    pair_counts_permuted_means = np.mean(pair_counts_null, axis=2)
    pair_counts_permuted_stds = np.std(pair_counts_null, axis=2)
    
    z_scores = (pair_counts - pair_counts_permuted_means) / pair_counts_permuted_stds
    z_scores[z_scores < 0] = 0
    z_score_df = pd.DataFrame(z_scores, index=categories_str_cat, columns=categories_str_cat)
    print("Finished!")
    
    return z_score_df


def plot_connectivity(df, cmap, dpi = 100, figsize = (6,5)):

    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    labels = df.columns
    np.fill_diagonal(df.values, 0)
    
    fig, ax = plt.subplots(figsize=figsize, dpi = dpi)

    # Use the mask in the heatmap
    sns.heatmap(df, mask=mask, cbar=False, cmap=cmap, ax=ax, annot_kws={"fontsize":6}, linewidths=1, linecolor='#D3D3D3')

    ax.set(xticks=np.arange(df.shape[1])+0.5,
           yticks=np.arange(df.shape[0])+0.5,
           xticklabels=labels, yticklabels=labels)

    ax.tick_params(axis='both', which='both', length=0, labelsize=12)

    plt.setp(ax.get_xticklabels(), rotation=40, ha="right", rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor")

    norm = mcolors.Normalize(vmin=df.min().min(), vmax=df.max().max())
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, location="right")

    triangle_coords = [(0, 0), (len(df), 0), (len(df), len(df))]
    triangle = plt.Polygon(triangle_coords, facecolor='white', edgecolor='white')
    ax.add_patch(triangle)

    triangle_coords = [(0, 0), (0, len(df)), (len(df), len(df))]
    triangle = plt.Polygon(triangle_coords, facecolor='white', edgecolor='#D3D3D3', fill=False)
    ax.add_patch(triangle)

    ax.grid(False)
    plt.show()
    
def gesa(gene_list, enrichr_library, show_number = 10, bar_color = '#fccde5'):
    """
    This function is adapted from the Enrichr Appyter notebooks.
    For more information, visit: https://maayanlab.cloud/Enrichr/
    """
    terms = []
    pvalues = []
    adjusted_pvalues = []
    
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(gene_list)
    description = ''
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise APIFailure
        
    data = json.loads(response.text)
    time.sleep(0.5)
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    short_id = data["shortId"]
    gene_set_library = enrichr_library
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise APIFailure
        
    data = json.loads(response.text)
    
    if len(data[enrichr_library]) == 0:
        raise NoResults
        
    results_df = pd.DataFrame(data[enrichr_library][0:show_number])
    terms = list(results_df[1])
    terms = [i.split(' (')[0] for i in terms]
    pvalues = list(results_df[2])
    adjusted_pvalues = list(results_df[6])

    ### bar plot
    bar_color_not_sig = 'white'
    edgecolor=None
    linewidth=0
    
    plt.figure(figsize=(4,3), dpi = 150)
    
    bar_colors = [bar_color if (x < 0.05) else bar_color_not_sig for x in pvalues]
    fig = sns.barplot(x=np.log10(pvalues)*-1, y=terms, palette=bar_colors, edgecolor=edgecolor, linewidth=linewidth)
    fig.set_xlabel('−log₁₀p‐value', fontsize=12, labelpad=10)
    fig.set_ylabel(f'{enrichr_library}', fontsize=12, labelpad=10)
    fig.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.tick_params(axis='x', which='major', labelsize=10)
    fig.set_yticks([])
    fig.tick_params(axis='y', which='both', length=0)
    fig.yaxis.label.set_visible(True)
    
    if max(np.log10(pvalues)*-1)<1:
        fig.xaxis.set_ticks(np.arange(0, max(np.log10(pvalues[i])*-1), 0.1))
        
    for ii, annot in enumerate(terms):
        if adjusted_pvalues[ii] < 0.05:
            annot = '  *'.join([annot]) 
        else:
            annot = '  '.join([annot])
            
        title_start= max(fig.axes.get_xlim())/200
        fig.text(title_start, ii, annot, ha='left', va='center', wrap = True, fontsize = 10)
        
    fig.spines['right'].set_visible(False)
    fig.spines['top'].set_visible(False)
    
    plt.show() 