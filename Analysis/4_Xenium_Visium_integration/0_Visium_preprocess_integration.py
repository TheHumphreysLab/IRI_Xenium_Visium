import numpy as np
import scanpy as sc
import pandas as pd
import anndata

sham = sc.read_visium(path = 'spaceranger_outs/mouse_sham/outs/')
hour4 = sc.read_visium(path = 'spaceranger_outs/mouse_hour4/outs/')
hour12 = sc.read_visium(path = 'spaceranger_outs/mouse_hour12/outs/')
day2 = sc.read_visium(path = 'spaceranger_outs/mouse_day2/outs/')
day14 = sc.read_visium(path = 'spaceranger_outs/mouse_day14/outs/')
week6 = sc.read_visium(path = 'spaceranger_outs/mouse_week6/outs/')

ad_list = [sham,hour4,hour12,day2,day14,week6]

for item in ad_list:
    item.var_names_make_unique()
    sc.pp.filter_cells(item, min_counts = 500)
    item.obsm["spatial"] = item.obsm["spatial"].astype(int)
    
for item in ad_list:
    sc.pp.calculate_qc_metrics(item, inplace=True)
    item.var['mt'] = [gene.startswith('MT-') for gene in item.var_names]
    item.obs['mt_frac'] = item[:, item.var['mt']].X.sum(1).A.squeeze()/item.obs['total_counts']
    
    
for item in ad_list:
    item.obs["X_coords"] = item.obsm["spatial"][:,0]
    item.obs["Y_coords"] = item.obsm["spatial"][:,1]
    
adata = anndata.concat(ad_list, join = "inner", fill_value = 0, uns_merge = "unique")
adata.raw = adata

sc.pp.normalize_total(adata, inplace = True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=3000, inplace=True)

sc.pp.pca(adata)
import scanpy.external as sce
sce.pp.harmony_integrate(adata, "ident", max_iter_harmony = 15)
adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]

sc.pp.neighbors(adata, n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=2, key_added="res2")
