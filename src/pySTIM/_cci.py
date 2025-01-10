import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.stats.multitest
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
from descartes import PolygonPatch
from matplotlib.patches import Patch, Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.colors import LinearSegmentedColormap
from typing import Dict, List, Optional, Tuple
from anndata import AnnData

from ._dataset import lr_db
from ._sp_plot import plot_scatter
from ._utils import to_hex, crop, get_color_map, split_field

def compute_cci(
    adata: "AnnData",
    group: str,
    sender: str,
    receiver: str,
    contact_radius: float = 30,
    p_value_threshold: float = 0.05,
    fold_change_threshold: float = 1,
    exp_fraction_threshold: float = 0.1,
    expression_threshold: float = 0,
    spatial_key: str = 'spatial'
) -> Dict[str, pd.DataFrame]:
    """
    Infer statistically significant ligand-receptor (LR) pairs by analyzing their co-expression levels within proximal sender-receiver cell pairs.

    Parameters
    ----------
    adata : AnnData
        An anndata object containing the data.
    group : str
        Column name to group cells by (e.g., cell type).
    sender : str
        Group name of the sender cells.
    receiver : str
        Group name of the receiver cells.
    contact_radius : float, optional
        Radius to define contact between cells, by default 30.
    p_value_threshold : float, optional
        P-value threshold for filtering significant interactions, by default 0.05.
    fold_change_threshold : float, optional
        Fold change threshold for filtering significant interactions, by default 1.
    exp_fraction_threshold : float, optional
        Expression fraction threshold for filtering significant interactions, by default 0.1.
    expression_threshold : float, optional
        Expression level threshold to consider a gene expressed, by default 0.
    spatial_key : str, optional
        Key for accessing spatial coordinates in `adata.obsm`, by default 'spatial'.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary with two DataFrames: "cell_pair" for cell pairs and "lr_pair" for ligand-receptor pairs.
    """

    # Load ligand-receptor pairs
    lr_network = lr_db.load_dataframe()
    lr_network['lr_pair'] = lr_network['from'].str.cat(lr_network['to'], sep='-')

    ligand = lr_network["from"].unique()
    expressed_ligand = list(set(ligand) & set(adata.var_names))
    print("Number of expressed ligand genes: ", len(expressed_ligand))
    lr_network = lr_network[lr_network["from"].isin(expressed_ligand)]
    
    receptor = lr_network["to"].unique()
    expressed_receptor = list(set(receptor) & set(adata.var_names))
    print("Number of expressed receptor genes: ", len(expressed_receptor))
    lr_network = lr_network[lr_network["to"].isin(expressed_receptor)]
    ## drop duplicates
    lr_network = lr_network[["from","to","lr_pair"]]
    lr_network = lr_network.drop_duplicates()

    # Prepare sender and receiver datasets
    sender_cells = adata[adata.obs[group] == sender]
    receiver_cells = adata[adata.obs[group] == receiver]

    # Calculate contact cell pairs
    contact_pairs = get_close_contact_cell_pairs(sender_cells.obs, receiver_cells.obs, contact_radius)
    noncontact_pairs = random_select_non_contact_cell_pairs(sender_cells.obs, receiver_cells.obs, contact_pairs)

    results = []

    for lr_pair in tqdm(lr_network.itertuples(), total=len(lr_network), desc='Processing LR Pairs'):
        
        ligand, receptor = lr_pair._1, lr_pair.to

        contact_lr_exp, co_exp_count = get_lr_pair_expressions_for_cell_pairs(
            sender_cells, receiver_cells, contact_pairs, ligand, receptor, expression_threshold)
        
        noncontact_lr_exp, _ = get_lr_pair_expressions_for_cell_pairs(
            sender_cells, receiver_cells, noncontact_pairs, ligand, receptor, expression_threshold)

        contact_mean_exp = np.mean(contact_lr_exp)
        noncontact_mean_exp = np.mean(noncontact_lr_exp)

        fold_change = contact_mean_exp / noncontact_mean_exp
        exp_fraction = np.sum(contact_lr_exp > 0) / len(contact_lr_exp)
        LR_score = np.log2(contact_mean_exp + 1)
        
        _, p_value = scipy.stats.ttest_ind(contact_lr_exp, noncontact_lr_exp, equal_var=False, alternative='greater')
        
        if p_value < p_value_threshold and fold_change > fold_change_threshold and exp_fraction > exp_fraction_threshold:
            results.append([sender, receiver, ligand, receptor, lr_pair.lr_pair, LR_score, contact_mean_exp, exp_fraction, co_exp_count, fold_change, p_value])
        
    results_df = pd.DataFrame(results, columns=[
        'Sender', 'Receiver', 'Ligand', 'Receptor', 'LR_pair', 'LR_Score', 'mean_contact_lr_exp', 'exp_fraction', 'co_exp_count', 'fold_change', 'pval'])
    
    results_df = results_df[~np.isnan(results_df['pval'])]
    
    results_df['pval_adjusted'] = statsmodels.stats.multitest.multipletests(results_df['pval'], method='fdr_bh')[1]
    
    results_df_filtered = results_df[(results_df['pval_adjusted'] < p_value_threshold) &
                                     (results_df['fold_change'] > fold_change_threshold) &
                                     (results_df['exp_fraction'] > exp_fraction_threshold)]
    
    cell_pair_df = pd.DataFrame(contact_pairs, columns = ["cell_sender", "cell_receiver"])
    res = {"cell_pair": cell_pair_df, "lr_pair": results_df_filtered}
    
    return res

def get_close_contact_cell_pairs(query_df, ref_df, contact_radius, xy_cols=['x_centroid', 'y_centroid']):
    query_coords = np.array(query_df[xy_cols])
    ref_coords = np.array(ref_df[xy_cols])

    contact_cell_pairs = []
    point_tree = scipy.spatial.cKDTree(ref_coords)
    neighbor_ids = point_tree.query_ball_point(query_coords, contact_radius)

    for i in range(len(neighbor_ids)):
        for j in neighbor_ids[i]:
            contact_cell_pairs.append((query_df.index[i], ref_df.index[j]))
    return np.array(contact_cell_pairs)


def random_select_non_contact_cell_pairs(query_df, ref_df, contact_cell_pairs):
    non_contact_query_cells = query_df.index[~query_df.index.isin(contact_cell_pairs[:, 0])]
    non_contact_ref_cells = ref_df.index[~ref_df.index.isin(contact_cell_pairs[:, 1])]

    n_pairs = len(contact_cell_pairs)
    return np.stack((np.random.choice(non_contact_query_cells, n_pairs), np.random.choice(non_contact_ref_cells, n_pairs)), axis=1)


def get_lr_pair_expressions_for_cell_pairs(adata_ligand_cells, adata_receptor_cells, cell_pairs, ligand_gene, receptor_gene, expression_threshold):
    ligand_map = pd.DataFrame({'cell_id': adata_ligand_cells.obs.index, 'index': range(adata_ligand_cells.shape[0])}).set_index('cell_id')
    receptor_map = pd.DataFrame({'cell_id': adata_receptor_cells.obs.index, 'index': range(adata_receptor_cells.shape[0])}).set_index('cell_id')

    lig_cell_ids = np.array(ligand_map['index'].loc[cell_pairs[:, 0]])
    rec_cell_ids = np.array(receptor_map['index'].loc[cell_pairs[:, 1]])

    lig_gene_id = np.where(adata_ligand_cells.var.index == ligand_gene)[0][0]
    rec_gene_id = np.where(adata_receptor_cells.var.index == receptor_gene)[0][0]

    lig_exp = adata_ligand_cells.X[lig_cell_ids, lig_gene_id]
    rec_exp = adata_receptor_cells.X[rec_cell_ids, rec_gene_id]

    co_exp_count = np.sum((lig_exp > expression_threshold) & (rec_exp > expression_threshold))

    return lig_exp * rec_exp, co_exp_count

# Example usage
# result = compute_cci(adata, group, 'Sender', 'Receiver')

def vis_cci(adata, res, cmap_dict):
    
    celltype_sender = np.unique(res["lr_pair"]["Sender"])[0]
    celltype_Receiver = np.unique(res["lr_pair"]["Receiver"])[0]
    
    key_added = f"{celltype_sender}-{celltype_Receiver}_ccipair"
    ta.obs[key_added] = "Rest"
    adata.obs.loc[adata.obs.index.isin(res["cell_pair"]["cell_sender"].tolist()), key_added] = celltype_sender
    adata.obs.loc[adata.obs.index.isin(res["cell_pair"]["cell_receiver"].tolist()), key_added] = celltype_Receiver
    
    if cmap_dict:
        plot_scatter(adata, color_by = key_added, cmap = cmap_dict, ptsize=1, dpi=150, legend_col=1)
    else:
        cmap_dict = {
            celltype_sender: '#2596be',
            celltype_Receiver: "#e8847d",
            "Rest": '#eeeeee'
        }
        plot_scatter(adata, color_by = key_added, cmap = cmap_dict, ptsize=1, dpi=150, legend_col=1)


def vis_lr_cci(adata, group, ligand, receptor, sender, receiver, bg_color, cmap1, cmap2, 
                   ptsize, xlims = None, ylims = None, ticks = True, poly=False):

    adata = adata.copy()
    adata.obs["cell_idx"] = adata.obs.index.astype(str)
    subset_idx, new_coord = crop(adata, xlims, ylims)
    colors_background = [bg_color] * len(subset_idx)

    cell_sender_set = set(adata[adata.obs[group] == sender].obs.cell_idx)
    subset_idx_set = set(subset_idx)
    sender_idx = [idx for (idx, cell) in enumerate(adata.obs["cell_idx"]) if cell in cell_sender_set and idx in subset_idx_set]
    sender_coord = adata.obs.iloc[sender_idx]

    print(f"Number of sender cells: {len(sender_idx)}")
    
    counts = np.array(sc.get.obs_df(adata, keys=ligand, use_raw=False).to_list())
    c_max = np.quantile(counts, 0.995)
    norm1 = matplotlib.colors.Normalize(vmin=min(counts), vmax=c_max)
    all_colors = cmap1(counts / c_max)
    all_colors = np.clip(all_colors, 0, 1)
    sender_colors = [all_colors[i] for i in sender_idx]
    sender_color_hex = [to_hex(i) for i in sender_colors]

    cell_receiver_set = set(adata[adata.obs[group] == receiver].obs.cell_idx)
    subset_idx_set = set(subset_idx)
    receiver_idx = [idx for (idx, cell) in enumerate(adata.obs["cell_idx"]) if cell in cell_receiver_set and idx in subset_idx_set]
    receiver_coord = adata.obs.iloc[receiver_idx]

    print(f"Number of receiver cells: {len(receiver_idx)}")
    
    counts2 = np.array(sc.get.obs_df(adata, keys=receptor, use_raw=False).to_list())
    c_max = np.quantile(counts2, 0.995)
    norm2 = matplotlib.colors.Normalize(vmin=min(counts2), vmax=c_max)
    all_colors = cmap2(counts2 / c_max)
    all_colors = np.clip(all_colors, 0, 1)
    receiver_colors = [all_colors[i] for i in receiver_idx]
    receiver_color_hex = [to_hex(i) for i in receiver_colors]

    fig, ax = plt.subplots(dpi=150, figsize=(5,5))
    plt.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)

    if poly:
        poly_dict = {idx + 1: adata.uns['poly'][key] for idx, key in enumerate(adata.uns['poly'].keys())}
        new_lib = [*poly_dict.values()]
        [ax.add_patch(PolygonPatch(new_lib[i], fc=colors_background[j], ec="#bfbfbf", alpha=1, zorder=0.1, linewidth=0.1)) for j, i
                            in enumerate(subset_idx)]
        [ax.add_patch(PolygonPatch(new_lib[i], fc=sender_color_hex[j], ec="#bfbfbf", alpha=1, zorder=0.1, linewidth=0.1)) for j, i
                            in enumerate(sender_idx)]
        [ax.add_patch(PolygonPatch(new_lib[i], fc=receiver_color_hex[j], ec="#bfbfbf", alpha=1, zorder=0.1, linewidth=0.1)) for j, i
                            in enumerate(receiver_idx)]
    else:
        ax.scatter(np.array(new_coord.x_centroid), np.array(new_coord.y_centroid), marker='o', linewidth=0,
                   alpha=1, color='#f2f2f2', s=ptsize)
        ax.scatter(np.array(sender_coord.x_centroid), np.array(sender_coord.y_centroid), marker='o', linewidth=0,
                   alpha=1, color=sender_color_hex, s=ptsize)
        ax.scatter(np.array(receiver_coord.x_centroid), np.array(receiver_coord.y_centroid), marker='o', linewidth=0,
                   alpha=1, color=receiver_color_hex, s=ptsize)

    if xlims and ylims:
        cbar_ax1 = fig.add_axes([1.0, 0.1, 0.01, 0.8]) 
        cbar_ax2 = fig.add_axes([1.08, 0.1, 0.01, 0.8]) 
    else:
        cbar_ax1 = fig.add_axes([0.8, 0.2, 0.01, 0.4])
        cbar_ax2 = fig.add_axes([.88, 0.2, 0.01, 0.4]) 
        
    cbar1 = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap=cmap1), cax=cbar_ax1, orientation='vertical')
    cbar1.set_ticks([min(counts), np.quantile(counts, 0.995)])
    cbar1.set_ticklabels(['', ''], fontsize=10)
    cbar1.ax.set_title(ligand, fontsize=10)
    cbar1.ax.tick_params(length=0)

    cbar2 = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm2, cmap=cmap2), cax=cbar_ax2, orientation='vertical')
    cbar2.set_ticks([min(counts2), np.quantile(counts2, 0.995)])
    cbar2.set_ticklabels(['', ''], fontsize=10)
    cbar2.ax.set_title(receptor, fontsize=10)
    cbar2.ax.tick_params(length=0)
      

    scalebar = AnchoredSizeBar(ax.transData,
               20, '', 'lower left', 
               bbox_to_anchor=(0.1, 0.05),
               borderpad=0, 
               bbox_transform=ax.transAxes, 
               sep=0,
               color='black',
               frameon=False,
               size_vertical=0.5,
               fontproperties=fm.FontProperties(size=12))

    ax.add_artist(scalebar)

    ax.axis('scaled')
    ax.axis("off")
        
    ax.grid(False)
    ax.invert_yaxis()
    plt.show()