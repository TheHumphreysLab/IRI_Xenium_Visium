import os
import sys
import math
import numpy as np
import scanpy as sc
import pandas as pd

import pySTIM as pst

### Loading all dataset and integration
sample_list = ['ShamL', 'ShamR', 'Hour4L', 'Hour4R', 'Hour12L', 'Hour12R', 
				'Day2L', 'Day2R', 'Day14L', 'Day14R', 'Week6L', 'Week6R']

spadata_dict = {}
folder_path = 'run_baysor/'
for sample in sample_list:
	data_dir = os.path.join(folder_path, sample)
	adata, _ = pcs.xn.load_xenium_baysor(data_dir)
	spadata_dict[sample] = adata 
	
	
for key, adata in spadata_dict.items():
	print('='*50)
	print(key)
	print(adata.shape)
	sc.pp.filter_cells(adata, min_genes=3)
	sc.pp.filter_cells(adata, min_counts=10)

	
for key, adata in spadata_dict.items():
	print(key)
	adata.obs_names = [f"{key}_{name}" for name in adata.obs_names]
	adata.obs["ident"] = key
	
	
kidney_all = anndata.concat(list(spadata_dict.values()), join = "inner", fill_value = 0)
kidney_all.raw = kidney_all

from scipy.sparse import csc_matrix
def normalize_data(adata, scale_factor=100):
	
	data_matrix = adata.X.toarray()
	
	total_counts = np.sum(data_matrix, axis=1).reshape(-1, 1)
	norm_mtx = (data_matrix / total_counts) * scale_factor
	log_mtx = np.log1p(norm_mtx)
	
	adata.X = csc_matrix(log_mtx.astype(np.float32))
	
	return adata

kidney_all = normalize_data(kidney_all)

sc.pp.pca(kidney_all)
kidney_all.obs['ident'] = kidney_all.obs['ident'].astype('category')
kidney_all.obs['ident'] = kidney_all.obs['ident'].cat.reorder_categories(order)

## Harmony integration
import scanpy.external as sce
sce.pp.harmony_integrate(kidney_all, "ident", max_iter_harmony = 15)

### Save anndata object and read in R (For label transfer using Seurat integration)
from scipy.io import mmwrite

mmwrite("IRI_all_count.mtx", adata.raw.X)
pd.DataFrame(adata.obs_names, columns=["cells"]).to_csv("py2r/IRI_all_cells.csv")
pd.DataFrame(adata.var_names, columns=["genes"]).to_csv("py2r/IRI_all_genes.csv")
pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1","UMAP2"]).to_csv("py2r/IRI_all_umap.csv")
pd.DataFrame(adata.obsm["X_pca"]).to_csv("py2r/IRI_all_pca.csv")
pd.DataFrame(adata.obsm["X_pca_harmony"]).to_csv("py2r/IRI_all_pca_harmony.csv")
adata.obs[["x_centroid", "y_centroid"]].to_csv("py2r/IRI_all_spatial.csv")
adata.obs.to_csv("py2r/IRI_all_meta.csv")