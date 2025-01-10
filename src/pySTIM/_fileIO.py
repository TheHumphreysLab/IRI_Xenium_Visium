import os
import json
import anndata
import tarfile
import scanpy as sc
import pandas as pd
import numpy as np
from typing import Any
from tqdm import tqdm
from scipy.sparse import csr_matrix
from PIL import Image
from typing import Any, Dict
import logging

def read_xenium(xenium_path: str) -> sc.AnnData:
    """
    Reads and processes Xenium data from the specified path.

    Parameters
    ----------
    xenium_path : str
        Path to the directory containing Xenium data.

    Returns
    -------
    sc.AnnData
        Annotated data matrix containing the processed Xenium data.
    """

    def read_data(file_name, file_type='csv'):
        file_path = os.path.join(xenium_path, file_name)
        tqdm.write(f"Reading {file_name}...")
        if file_type == 'csv':
            return pd.read_csv(file_path, compression='gzip')
        elif file_type == 'parquet':
            return pd.read_parquet(file_path)
    
    try:
        tqdm.write("Starting to read and process Xenium data.")

        # Read cell feature matrix using Scanpy
        tqdm.write("Reading cell feature matrix...")
        xenium_adata = sc.read_10x_h5(os.path.join(xenium_path, 'cell_feature_matrix.h5'))

        # Read spatial and molecular data
        cell_coord = read_data('cells.parquet', 'parquet')
        molecule_data = read_data('transcripts.parquet', 'parquet')
        cell_seg = read_data('cell_boundaries.parquet', 'parquet')

        # Remove unassigned transcripts
        molecule_data = molecule_data[~(molecule_data['cell_id'] == "UNASSIGNED")]

        # Ensure the 'analysis' directory is present
        analysis_dir = os.path.join(xenium_path, 'analysis')
        if not os.path.exists(analysis_dir):
            tqdm.write("Extracting analysis.tar.gz...")
            tar_path = os.path.join(xenium_path, 'analysis.tar.gz')
            with tarfile.open(tar_path, 'r:gz') as tar:
                tar.extractall(path=xenium_path)

        # Read UMAP and clustering results
        tqdm.write("Reading UMAP data...")
        umap_path = os.path.join(analysis_dir, 'umap/gene_expression_2_components/projection.csv')
        clusters_path = os.path.join(analysis_dir, 'clustering/gene_expression_graphclust/clusters.csv')
        xenium_umap = pd.read_csv(umap_path)
        xenium_clusters = pd.read_csv(clusters_path)

        # Update AnnData
        xenium_adata.obs = cell_coord.copy()

        # Filter cells based on UMAP results
        tqdm.write("Filtering cells based on UMAP results...")
        cell_id_keep = xenium_umap['Barcode'].tolist()
        xenium_adata = xenium_adata[xenium_adata.obs['cell_id'].isin(cell_id_keep)]
        cell_seg = cell_seg[cell_seg['cell_id'].isin(cell_id_keep)]

        # Update spatial information
        xenium_adata.obsm['spatial'] = xenium_adata.obs[['x_centroid', 'y_centroid']].copy().to_numpy()

        # Update UMAP data
        xenium_umap.index = xenium_adata.obs_names
        xenium_adata.obsm['X_umap'] = xenium_umap.iloc[:, 1:3].copy().to_numpy()

        # Update cluster information
        clusters = xenium_clusters['Cluster'].to_list()
        xenium_adata.obs['cluster'] = list(map(str, clusters))

        # Filter and update molecular data
        tqdm.write("Updating molecular data...")
        gene_names = xenium_adata.var_names.astype(str).tolist()
        molecule_data = molecule_data[molecule_data['feature_name'].isin(gene_names)]
        molecule_data = molecule_data[molecule_data['cell_id'].isin(cell_id_keep)]

        xenium_adata.uns['seg'] = cell_seg
        xenium_adata.uns['transcript'] = molecule_data
        xenium_adata.uns['datatype'] = "Xenium"

        tqdm.write("Data loading completed.")
        return xenium_adata
    
    except Exception as e:
        raise RuntimeError(f"Failed to read and process Xenium data: {str(e)}")


def read_layers(hd_dir: str, bin_size: int = 2) -> sc.AnnData:
	"""
	Reads VisiumHD layers, optionally including UMAP and cluster data if available.

	Parameters:
	- hd_dir (str): The VisiumHD directory.
	- bin_size (int): The size of the bin to determine spatial resolution.

	Returns:
	- sc.AnnData: An annotated data object.
	"""
	
	def compute_corner_points(df: pd.DataFrame, px_width: float, cell_col='barcode',
								x_col='pxl_row_in_fullres', y_col='pxl_col_in_fullres') -> pd.DataFrame:
		corner_coordinates = {
			'new_x': df[x_col] - 0.5 * px_width,
			'new_y': df[y_col] - 0.5 * px_width,
			'barcode': df[cell_col]
		}
		return pd.DataFrame(corner_coordinates)
	
	layer_dir = os.path.join(hd_dir, f"binned_outputs/square_{bin_size:03}um")
	h5_file = os.path.join(layer_dir, "filtered_feature_bc_matrix.h5")
	pos_file = os.path.join(layer_dir, "spatial/tissue_positions.parquet")
	json_file = os.path.join(layer_dir, "spatial/scalefactors_json.json")
	
	layer_data = sc.read_10x_h5(h5_file)
	layer_data.var_names_make_unique()
	
	pos = pd.read_parquet(pos_file)
	pos.index = pos['barcode']

	all_cells = set(layer_data.obs_names)
	
	pos = pos[(pos['pxl_row_in_fullres'] > 0) & (pos['pxl_col_in_fullres'] > 0)]
	cells_keep = pos.index.intersection(all_cells)
	layer_data = layer_data[cells_keep]
    
	pos = pos.loc[layer_data.obs.index, ]
	layer_data.obs = layer_data.obs.join(pos)
	
	with open(json_file) as f:
		json_data = json.load(f)
		
	# Adjust for bin_size in spatial resolution
	px_width = bin_size / json_data["microns_per_pixel"]
	corner_coordinates = compute_corner_points(pos, px_width)
	layer_data.obsm["spatial"] = corner_coordinates[['new_y', 'new_x']].values
	
	# Load UMAP and clusters if available
	if not os.path.exists(layer_dir + "/analysis"):
		print(f"No analysis directory for bin size {bin_size}.")
		
	umap_path = os.path.join(layer_dir, 'analysis/umap/gene_expression_2_components/projection.csv')
	clusters_path = os.path.join(layer_dir, 'analysis/clustering/gene_expression_graphclust/clusters.csv')
	
	if os.path.exists(umap_path) and os.path.exists(clusters_path):
		umap = pd.read_csv(umap_path)
		clusters = pd.read_csv(clusters_path)
		umap.index = umap['Barcode']
		clusters.index = clusters['Barcode']
		
		# Filter AnnData based on cells present in UMAP
		cells_keep = umap.index.intersection(layer_data.obs_names)
		layer_data = layer_data[cells_keep]
		layer_data.obsm["X_umap"] = umap.loc[cells_keep, ["UMAP-1", "UMAP-2"]].copy().to_numpy()
		layer_data.obs["cluster"] = clusters.loc[cells_keep, "Cluster"].copy().to_numpy()
		layer_data.obs["cluster"] = layer_data.obs["cluster"].astype('category')
		
		
	layer_data.uns['spatial'] = {
			'scalefactors': json_data,
		}	
	layer_data.uns["spatial"]["images"] = {
		res: np.asarray(Image.open(os.path.join(layer_dir, f"spatial/tissue_{res}_image.png"))) for res in ["hires", "lowres"]
	}
	
	return layer_data

def read_visiumHD(hd_dir: str, bins: Any = "all") -> Dict[str, sc.AnnData]:
    """
    Loads and processes Visium HD data for specified bin sizes.

    Parameters
    ----------
    hd_dir : str
            Directory containing VisiumHD data.
    bins : Any
            Specific bin sizes to load, can be 'all', a single integer, or a list of integers.

    Returns
    -------
    Dict[str, sc.AnnData]
            A dictionary of annotated data objects indexed by bin size.
    """
    if isinstance(bins, str) and bins.lower() == "all":
        bin_sizes = [2, 8, 16]
    elif isinstance(bins, int):
        bin_sizes = [bins]
    elif isinstance(bins, list) and all(isinstance(bin_size, int) for bin_size in bins):
        bin_sizes = bins
    else:
        raise ValueError(
            "Invalid 'bins' parameter. It must be 'all', a single integer, or a list of integers."
        )

    adata_dict = {}

    for bin_size in bin_sizes:
        logging.info(f"Loading {bin_size}um binned data...")
        adata_dict[f"bin_{bin_size}um"] = read_layers(hd_dir, bin_size=bin_size)

    return adata_dict


def load_xenium_baysor(data_dir: str) -> sc.AnnData:
    """
    Load Baysor segmented Xenium data.

    Parameters:
    - data_dir (str): Path to the Baysor output directory.

    Returns:
    - two Anndata object, one saves cells, the other saves polygons.
    """

    print("Loading baysor segmentation output...")
    counts = pd.read_csv(os.path.join(data_dir, "segmentation_counts.tsv"), sep="\t")

    cells_meta = pd.read_csv(os.path.join(data_dir, "segmentation_cell_stats.csv"))
    cells_meta.index = cells_meta.cell

    molecule_data = pd.read_csv(os.path.join(data_dir, "segmentation.csv"))
    molecule_data = molecule_data[molecule_data["cell"] != 0]
    molecule_data["x_location"] = molecule_data["x"] * 0.2125
    molecule_data["y_location"] = molecule_data["y"] * 0.2125
    
    features = list(counts["gene"])
    counts = csr_matrix(counts.iloc[:, 1:].values, dtype=np.int32)
    adata = anndata.AnnData(counts.T, obs=cells_meta)
    adata.var_names = features
    adata.uns["datatype"] = "xenium"

    adata.obs["x_centroid"] = adata.obs["x"] * 0.2125 ## convert pixel to micron
    adata.obs["y_centroid"] = adata.obs["y"] * 0.2125
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
    
    print("Loading polygons...")
    polygon_df = pd.read_csv(os.path.join(data_dir, "polygon.csv"))
    polygon_meta = pd.read_csv(os.path.join(data_dir, "polygon_meta.csv"))
    polygon_count = pd.read_csv(os.path.join(data_dir, "polygon_count.csv"))
    polygon_df["x_centroid"] = polygon_df["x_location"] * 0.2125
    polygon_df["y_centroid"] = polygon_df["y_location"] * 0.2125

    poly_counts = csr_matrix(polygon_count.iloc[:, 1:].values, dtype=np.int32)
    polygon_centroid = polygon_df.groupby('polygon_number').mean().reset_index()
    polygon_meta = polygon_meta.merge(polygon_centroid, on="polygon_number", how="inner")
    polygon_meta.index = polygon_meta.polygon_number.astype(str)
    poly_adata = anndata.AnnData(poly_counts, obs=polygon_meta)
    poly_adata.var_names = list(polygon_count.columns[1:])

    poly_adata.uns["seg"] = polygon_df
    poly_adata.uns["transcript"] = molecule_data

    print("Creating polygons...")
    seg = poly_adata.uns['seg']
    grouped = seg.groupby("polygon_number")
    all_polygons = sorted(list(set(seg["polygon_number"])))
    new_lib = [{
        'coordinates': [grouped.get_group(cell).iloc[:, 3:5].to_numpy().tolist()],
        'type': 'Polygon'
    } for cell in tqdm(all_polygons)]

    new_poly = dict(zip(map(str, all_polygons), new_lib))
    poly_adata.uns['poly'] = new_poly

    return adata, poly_adata

def load_polygon(data_dir):

    print("Loading baysor segmentation output...")
    molecule_data = pd.read_csv(os.path.join(data_dir, "segmentation.csv"))
    molecule_data = molecule_data[molecule_data["cell"] != 0]
    molecule_data["x_location"] = molecule_data["x"] * 0.2125
    molecule_data["y_location"] = molecule_data["y"] * 0.2125
    print("Loading polygons...")
    polygon_df = pd.read_csv(os.path.join(data_dir, "polygon.csv"))
    polygon_meta = pd.read_csv(os.path.join(data_dir, "polygon_meta.csv"))
    polygon_count = pd.read_csv(os.path.join(data_dir, "polygon_count.csv"))
    polygon_df["x_centroid"] = polygon_df["x_location"] * 0.2125
    polygon_df["y_centroid"] = polygon_df["y_location"] * 0.2125

    poly_counts = csr_matrix(polygon_count.iloc[:, 1:].values, dtype=np.int32)
    polygon_centroid = polygon_df.groupby('polygon_number').mean().reset_index()
    polygon_meta = polygon_meta.merge(polygon_centroid, on="polygon_number", how="inner")
    polygon_meta.index = polygon_meta.polygon_number.astype(str)
    poly_adata = anndata.AnnData(poly_counts, obs=polygon_meta)
    poly_adata.var_names = list(polygon_count.columns[1:])

    poly_adata.uns["seg"] = polygon_df
    poly_adata.uns["transcript"] = molecule_data

    print("Creating polygons...")
    seg = poly_adata.uns['seg']
    grouped = seg.groupby("polygon_number")
    all_polygons = sorted(list(set(seg["polygon_number"])))
    new_lib = [{
        'coordinates': [grouped.get_group(cell).iloc[:, 3:5].to_numpy().tolist()],
        'type': 'Polygon'
    } for cell in tqdm(all_polygons)]

    new_poly = dict(zip(map(str, all_polygons), new_lib))
    poly_adata.uns['poly'] = new_poly

    return  poly_adata

    
