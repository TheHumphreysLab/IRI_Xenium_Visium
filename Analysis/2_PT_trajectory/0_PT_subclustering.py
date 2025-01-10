import numpy as np
import scanpy as sc
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt


adata = sc.read_h5ad("../Xenium_all.h5ad")
pt_ad = adata[adata.obs.celltype_raw.isin(["PTS1", "PTS2", "PTS3", 'Inj_PT', 'FR_PT'])]

pt_ad_BBKNN = pt_ad.raw.to_adata()
pt_ad_BBKNN = normalize_data(pt_ad_BBKNN)

sc.pp.pca(pt_ad_BBKNN)

import scanpy.external as sce
sce.pp.bbknn(pt_ad_BBKNN, batch_key="ident", n_pcs=30)
sc.tl.leiden(pt_ad_BBKNN, resolution=2, key_added="bbknn_res2")

