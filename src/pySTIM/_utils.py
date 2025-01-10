import random
import tifffile as tiff
from PIL import Image
#import cv2
import scanpy as sc
from tqdm import tqdm
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pandas as pd
from typing import Any, List, Tuple, Union, Dict, Optional
from scipy.sparse import csc_matrix, csr_matrix, coo_matrix
from shapely.geometry import Polygon
from sklearn.neighbors import radius_neighbors_graph, kneighbors_graph
from sklearn.cluster import MiniBatchKMeans
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import seaborn as sns
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def normalize_data(adata, scale_factor=100):
    
    data_matrix = adata.X.toarray()
    
    total_counts = np.sum(data_matrix, axis=1).reshape(-1, 1)
    norm_mtx = (data_matrix / total_counts) * scale_factor
    log_mtx = np.log1p(norm_mtx)

    adata.X = csc_matrix(log_mtx.astype(np.float32))
    
    return adata
    
def to_hex(rgb_tuple: Tuple[float, float, float]) -> str:
    """
    Convert an RGB tuple to its hex color.
    
    Parameters:
    - rgb_tuple (Tuple[float, float, float]): A tuple representing RGB values (range 0 to 1).
    
    Returns:
    - str: Hexadecimal string representing the color.
    """
    r = int(rgb_tuple[0] * 255)
    g = int(rgb_tuple[1] * 255)
    b = int(rgb_tuple[2] * 255)
    hex_code = '#{0:02x}{1:02x}{2:02x}'.format(r, g, b)
    return hex_code


def crop(adata: Any, 
         xlims: Optional[Tuple[float, float]], 
         ylims: Optional[Tuple[float, float]]) -> Tuple[List[int], pd.DataFrame]:
    """
    Crop adata based on x and y limits.

    Parameters:
    - adata: An anndata object.
    - xlims: Tuple containing the x-axis limits.
    - ylims: Tuple containing the y-axis limits.

    Returns:
    - Tuple[List[int], pd.DataFrame]: A tuple containing the subset indices and the subsetted cell coordinates DataFrame.
    """
    
    cell_coord = adata.obs
    x_new = adata.obsm["spatial"][:,0]
    y_new = adata.obsm["spatial"][:,1]

    expand = 1.05
    if xlims is None:
        minx = x_new.min()
        maxx = x_new.max()
        xlims = [(minx + maxx) / 2.0 - (maxx - minx) / 2.0 * expand,
                 (minx + maxx) / 2.0 + (maxx - minx) / 2.0 * expand]

    if ylims is None:
        miny = y_new.min()
        maxy = y_new.max()
        ylims = [(miny + maxy) / 2.0 - (maxy - miny) / 2.0 * expand,
                 (miny + maxy) / 2.0 + (maxy - miny) / 2.0 * expand]

    x_lim_idx = [idx for idx, value in enumerate(x_new) if value > xlims[0] and value < xlims[1]]
    y_lim_idx = [idx for idx, value in enumerate(y_new) if value > ylims[0] and value < ylims[1]]

    subset_idx = list(set(x_lim_idx) & set(y_lim_idx))
    cell_coord2 = cell_coord.iloc[subset_idx]

    return subset_idx, cell_coord2


def get_color_map(
    adata: Any,
    seed: int = 42,
    color_by: Optional[str] = None,
    cmap: Optional[Union[str, Dict[str, str]]] = None,
    genes: Optional[Union[str, List[str]]] = None,
    subset_idx: Optional[List[int]] = None,
):
    """
    Get a color map for the given adata.

    Parameters:
    - adata: An anndata object.
    - color_by (Optional[str]): Category from adata.obs by which to color the output.
    - cmap (Optional[Union[str, Dict[str, str]]]): Color map or a dictionary to map values to colors.
    - seed (int): Seed for random color generation.
    - genes (Optional[Union[str, List[str]]]): Whether the map is based on genes or not.
    - subset_idx (Optional[List[int]]): Indices of the subset of data to be considered.

    Returns:
    - Union[str, Dict[str, str]]: A color map or a dictionary mapping values to colors.
    """

    def generate_unique_random_colors(num_colors: int):
        random.seed(seed)
        colors = set()
        while len(colors) < num_colors:
            color = "#" + "".join(random.choices("0123456789ABCDEF", k=6))
            colors.add(color)
        return list(colors)

    def generate_cell_type_map():
        if color_by not in adata.obs:
            raise ValueError(f"{color_by} not found in adata.obs")
        adata.obs[color_by] = adata.obs[color_by].astype('category')
        cell_types = np.array(adata.obs[color_by].cat.categories)
        key = f"{color_by}_colors"
        if key in adata.uns.keys():
            cell_colors = adata.uns[key]
        else:
            cell_colors = generate_unique_random_colors(len(cell_types))
        return dict(zip(cell_types, cell_colors))

    if genes is None:
        genes = False

    if subset_idx is None:
        subset_idx = list(range(adata.shape[0]))

    if not genes:
        if color_by is None:
            map_dict = ["#6699cc"] * len(subset_idx)
        elif cmap is None:
            map_dict = generate_cell_type_map()
        else:
            map_dict = cmap
    else:
        if cmap is None:
            colors = ["#f2f2f2", "#ffdbdb", "#fc0303"]
            map_dict = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
        else:
            map_dict = cmap

    return map_dict



def split_field(df: pd.DataFrame, 
                n_fields_x: int, 
                n_fields_y: int, 
                x_col: str = 'x_centroid', 
                y_col: str = 'y_centroid') -> Tuple[List[List[Tuple[float, float]]], List[Tuple[float, float]]]:
    """
    Splits a spatial region into a grid of FOV based on x and y centroids in obs.
    
    Parameters:
    - df: AnnData.obs containing the x and y centroids.
    - n_fields_x: Number of FOV horizontally.
    - n_fields_y: Number of FOV vertically.
    - x_col: Column name for x centroids in df. Default is 'x_centroid'.
    - y_col: Column name for y centroids in df. Default is 'y_centroid'.
    
    Returns:
    - rectangles: each rectangle is represented by a list of its vertices.
    - centroids: centroids of these rectangles.
    """
    
    # Expand the bounding box slightly for better visualization
    expand = 1.1
    minx = df[x_col].min()
    maxx = df[x_col].max()
    x_range = expand * (maxx - minx)
    
    miny = df[y_col].min()
    maxy = df[y_col].max()
    y_range = expand * (maxy - miny)
    
    # Calculate the dimensions of each rectangle
    x_step = x_range / n_fields_x
    y_step = y_range / n_fields_y
    
    # Calculate the starting coordinates for the grid
    x_start = (minx + maxx) / 2.0 - (maxx - minx) / 2.0 * expand
    y_start = (miny + maxy) / 2.0 - (maxy - miny) / 2.0 * expand
    
    # Loop to create the grid of rectangles and calculate their centroids
    rectangles = []
    centroids = []
    for i in range(n_fields_x):
        for j in range(n_fields_y):
            top_left = (x_start + x_step * i, y_start + y_step * j)
            bottom_right = (top_left[0] + x_step, top_left[1] + y_step)
            centroid = (top_left[0] + x_step / 2, top_left[1] + y_step / 2)
    
            rectangle = [top_left, (bottom_right[0], top_left[1]), bottom_right, (top_left[0], bottom_right[1])]
            rectangles.append(rectangle)
            centroids.append(centroid)
    
    return rectangles, centroids

def filter_polygon(adata, poly_key = "poly", threshold = 10, view_summary = True):
    '''
    Filters polygons based on a size threshold
    
    Parameters:
    - adata: AnnData object.
    - threshold: Minimum area threshold for polygons to be kept.
    - view_summary: If True, displays a histogram of polygon areas.
    - inplace: If True, modifies `adata` in place, otherwise returns modified data.
    '''

    poly_dict = adata.uns.get(poly_key, {})
    polygon_areas = np.array([Polygon(poly['coordinates'][0]).area for poly in poly_dict.values()])
    
    if view_summary:
        print("Polygon area distribution")
        plt.figure(figsize=(5, 3))
        plt.hist(polygon_areas, bins=100, align='left', color='skyblue')
        plt.xlim(-10, 250)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.tight_layout()
        plt.show()
        
    print(f"Original polygon number: {len(adata.uns['poly'])}")
    filtered_poly_dict = {pid: pd for pid, pd in poly_dict.items() if Polygon(pd['coordinates'][0]).area >= threshold}
    print(f"Keep polygon number: {len(filtered_poly_dict)}")
    
    print("Saving result...")
    keep_poly_number = [int(i) for i in filtered_poly_dict.keys()]
    adata.uns["seg"] = adata.uns["seg"][adata.uns["seg"].polygon_number.isin(keep_poly_number)]
    adata.uns[poly_key] = filtered_poly_dict
    adata = adata[adata.obs.polygon_number.isin(keep_poly_number)]
    return adata

def smooth_polygon(adata, window_size=5):
    """
    Smoothens the coordinates of polygons using a moving average filter.

    Parameters:
    - adata: AnnData object containing the dataset.
    - window_size: The size of the moving window for smoothing.
    """
    def smooth_coordinates(polygon, window_size):
        n = len(polygon)
        polygon = np.array(polygon)
        smoothed = polygon.copy()
        for i in range(n):
            indices = [(i + j - window_size // 2) % n for j in range(window_size)]
            smoothed[i] = np.mean(polygon[indices], axis=0)
        return smoothed.tolist()

    seg = adata.uns["seg"]
    all_polygons = sorted(set(seg["polygon_number"]))

    new_polygons = []
    for polygon_number in tqdm(all_polygons, desc="Smoothing polygons"):
        group = seg[seg.polygon_number == polygon_number]
        raw_coords = group.iloc[:, 3:5].values.tolist()
        smoothed_coords = smooth_coordinates(raw_coords, window_size = window_size)
        new_polygons.append({
            'coordinates': [smoothed_coords],
            'type': 'Polygon'
        })

    adata.uns['poly'] = dict(zip(map(str, all_polygons), new_polygons))
    
def align_adata(transformation, adata_xe, adata_vis):
    
    ### Convert Visium spot coordinates to Xenium DAPI coordinates and scale
    points_vis = np.array(adata_vis.obs[["X_coords","Y_coords"]])
    points_vis2 = np.hstack((points_vis, np.ones((points_vis.shape[0], 1))))
    points_xe = np.dot(points_vis2, H.T)
    points_xe2 = points_xe[:, :2] / points_xe[:, 2, np.newaxis]
    points_xe2 = points_xe2*0.2125

    ### Check the alignment
    fig, axs = plt.subplots(1,2,dpi=100)
    axs[0].scatter(points_xe2[:,0], points_xe2[:,1],  marker='o', linewidth=0, alpha=1, color="black", s=15)
    axs[0].axis('scaled')
    axs[0].set_title("Visium")
    axs[1].scatter(adata_xe.obs["x_centroid"], adata_xe.obs["y_centroid"],  marker='o', linewidth=0, alpha=1, color="black", s=1)
    axs[1].axis('scaled')
    axs[1].set_title("Xenium")
    axs[0].invert_yaxis()
    axs[1].invert_yaxis()
    plt.show()

    return points_xe2

def add_img(adata, img_path, dapi_path, transformation=None, scale_factor=0.2125, library_key="spatial"):
    # Read the histology image
    if img_path.lower().endswith(('.tif', '.tiff')):
        img = tiff.imread(img_path)
    else:
        img = np.array(Image.open(img_path))
    
    # Read the DAPI image
    if dapi_path.lower().endswith(('.tif', '.tiff')):
        dapi_image = tiff.imread(dapi_path)
    
    # Apply transformation if provided
    if transformation is not None:
        height = dapi_image.shape[-1] * scale_factor
        width = dapi_image.shape[-2] * scale_factor

        if isinstance(transformation, (np.ndarray, list)):
            img = cv2.warpPerspective(img, transformation, output_shape=(int(height), int(width)))
        
        elif isinstance(transformation, str):

            lm_df = pd.read_csv(transformation)
            lm_df["X"] = lm_df["fixedX"] * scale_factor
            lm_df["Y"] = lm_df["fixedY"] * scale_factor
            fixed_points = lm_df[['X', 'Y']].values
            alignment_points = lm_df[['alignmentX', 'alignmentY']].values
            fixed_points = np.float32(fixed_points)
            alignment_points = np.float32(alignment_points)

            # Compute the transformation matrix
            transformation_matrix, status = cv2.findHomography(alignment_points, fixed_points)
            img = cv2.warpPerspective(img, transformation_matrix, output_shape=(int(height), int(width)))
    
    # Add the image to the anndata object
    if library_key not in adata.uns:
        adata.uns[library_key] = {}

    adata.uns[library_key]["image"] = img

    return adata

   
def register_xenium_visium(adata_xe, adata_vis, spot_diameter=100, agg='celltype'):
    from scipy.spatial import cKDTree

    # Create k-d tree for Visium spots
    visium_coords = adata_vis.obsm["spatial"]
    visium_tree = cKDTree(visium_coords)
    
    cell_to_spot_mapping = {}

    # Iterate over Xenium cells and find the corresponding Visium spots within the specified diameter
    for idx, row in adata_xe.obs.iterrows():
        x_centroid, y_centroid = row['x_centroid'], row['y_centroid']
        cell_id = row.name

        indices = visium_tree.query_ball_point([x_centroid, y_centroid], spot_diameter / 2)
        for spot_idx in indices:
            spot_id = adata_vis.obs.index[spot_idx]
            if spot_id not in cell_to_spot_mapping:
                cell_to_spot_mapping[spot_id] = []
            cell_to_spot_mapping[spot_id].append(cell_id)

    print(f"Original Xenium cell #: {adata_xe.shape[0]}, mapped Xenium cell #: {sum(len(cells) for cells in cell_to_spot_mapping.values())}")

    counts = [len(cells) for cells in cell_to_spot_mapping.values()]
    plt.figure(figsize=(5, 3))
    plt.hist(counts, bins=15, edgecolor='k', alpha=0.7)
    plt.title('Distribution of Xenium Cells Mapped to Visium Spots')
    plt.xlabel('Number of Xenium Cells per Visium Spot')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.show()
    
    # Aggregate gene expression for cells within each spot
    print("Calculating pseudo-Visium data.....")
    gene_overlap = [g for g in adata_xe.var_names if g in adata_vis.var_names]
    count_df = sc.get.obs_df(adata_xe, keys=gene_overlap, use_raw=True)
    agg_gene = []

    for spot_id in adata_vis.obs_names:
        cells_in_spot = cell_to_spot_mapping.get(spot_id, [])
        counts = count_df.loc[cells_in_spot,]
        sum_counts = counts.sum(axis=0).to_numpy()
        agg_gene.append(sum_counts)

    agg_gene_df = pd.DataFrame(agg_gene, columns=gene_overlap, index=adata_vis.obs_names)
    agg_gene_df.fillna(0, inplace=True)

    # Calculating correlation
    agg_gene_df = agg_gene_df.loc[adata_vis.obs.index,]
    Xenium_pseudo_gene_exp = agg_gene_df.loc[:, gene_overlap]
    Visium_gene_df = sc.get.obs_df(adata_vis, keys=gene_overlap, use_raw=True)

    corr = pd.DataFrame(Xenium_pseudo_gene_exp.corrwith(Visium_gene_df, method="pearson"), columns=["pearson_r"])
    corr[corr["pearson_r"] <= 0] = 0
    corr = corr.fillna(0)
    print("Top 5 correlated spatial gene expression .....")
    print(corr.sort_values(by="pearson_r", ascending=False).head(10))

    return cell_to_spot_mapping, Xenium_pseudo_gene_exp


def calculate_cn(adata, mode='radius', group = "ident", cluster = "celltype", radius=22.5, k=20, n_neighborhoods=20, plot = True, figsize = (5,5), dpi=100):
    
    samples = list(adata.obs[group].cat.categories)
    neighbors_dict = {}
    cell_idx = []

    for sample in samples:
        tmp = adata[adata.obs[group] == sample]
        cell_idx += list(tmp.obs_names)
        x_y_coordinates = tmp.obsm["spatial"]
        cts = np.array(tmp.obs[cluster])

        if mode == 'radius':
            connectivity_matrix = radius_neighbors_graph(x_y_coordinates, radius, mode='connectivity', include_self=False)
        elif mode == 'knn':
            connectivity_matrix = kneighbors_graph(x_y_coordinates, k, mode='connectivity', include_self=False)
        else:
            raise ValueError("Mode should be either 'radius' or 'knn'.")

        one_hot_matrix = pd.get_dummies(cts).values
        cell_by_neighborhood_matrix = connectivity_matrix.dot(one_hot_matrix)

        cell_type_columns = pd.get_dummies(cts).columns
        cell_by_neighborhood_df = pd.DataFrame(cell_by_neighborhood_matrix, columns=cell_type_columns)

        row_sums = cell_by_neighborhood_df.sum(axis=1)
        cell_by_neighborhood_df = cell_by_neighborhood_df.div(row_sums, axis=0)
        cell_by_neighborhood_df.fillna(0, inplace=True)
        cell_by_neighborhood_df["cell_type"] = cts
        cell_by_neighborhood_df[["x", "y"]] = x_y_coordinates
        cell_by_neighborhood_df[group] = sample
        neighbors_dict[sample] = cell_by_neighborhood_df

    all_neighbors_df = pd.concat(neighbors_dict.values())
    
    values = np.array(all_neighbors_df.iloc[:, :-4])
    km = MiniBatchKMeans(n_clusters=n_neighborhoods, random_state=0)
    labelskm = km.fit_predict(values)
    k_centroids = km.cluster_centers_
    
    cell_types = all_neighbors_df.columns[:-4]
    tissue_avgs = values.mean(axis = 0)
    fc = np.log2(((k_centroids+tissue_avgs)/(k_centroids+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
    fc = pd.DataFrame(fc, columns = cell_types, index = ["CN"+ str(i) for i in range(n_neighborhoods)])
    
    if plot:
        ### Visualize the cell type enrichment in each CN
        sns.set(font_scale=1, rc={"figure.dpi":dpi})
        g = sns.clustermap(fc, vmin=-2, vmax=2, cmap="RdYlBu_r", cbar_pos=None, 
                       row_cluster=False, col_cluster=True, figsize=figsize)
    
        #g.ax_heatmap.yaxis.set_label_position("right")
        #g.ax_heatmap.yaxis.tick_right()
        #plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        
        cbar_ax =  g.fig.add_axes([1.0, 0.1, 0.2, 0.02])  # Adjust these values as needed for your layout
        norm = Normalize(vmin=-2, vmax=2)
        cbar = ColorbarBase(cbar_ax, cmap="RdYlBu_r", norm=norm, orientation='horizontal')
        cbar.set_label('Enrichment Score', labelpad=-40, horizontalalignment='center')
        
        plt.tight_layout()
        plt.show()
        
    adata.obs["CN"] = pd.Categorical(pd.Series(labelskm))
    
def add_metadata(adata_vis, adata_xe, anchor_dict, col = "CN"):

    from collections import Counter

    metadata = pd.DataFrame(index=adata_vis.obs.index, columns=[col])
    
    for spot_id, cell_ids in anchor_dict.items():
        values = np.array(adata_xe.obs.loc[cell_ids, col])
        if values.size > 0:
            mapped_value = Counter(values).most_common(1)[0][0]
            metadata.loc[spot_id, col] = mapped_value
            
    metadata[col] = pd.Categorical(metadata[col], categories=adata_xe.obs[col].cat.categories)
    adata_vis.obs[col] = metadata[col]

    valid_spots = metadata.dropna().index
    n_filtered = len(adata_vis.obs.index) - len(valid_spots)
    print(f"Filtering {n_filtered} Visium spots due to no cells mapping")
    adata_vis = adata_vis[adata_vis.obs.index.isin(valid_spots)]
    
    return adata_vis
    