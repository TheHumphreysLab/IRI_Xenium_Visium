import os
import sys
import anndata
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns

adata = sc.read_h5ad('Visium_IRI_R.h5ad')

import decoupler as dc
net = dc.get_collectri(organism='mouse', split_complexes=False)
dc.run_ulm(
	mat=adata,
	net=net,
	source='source',
	target='target',
	weight='weight',
	verbose=True
)
acts = dc.get_acts(adata, obsm_key='ulm_estimate')
df = dc.rank_sources_groups(acts, groupby='CN_rename', reference='rest', method='wilcoxon')
df = df[(df['pvals_adj']<0.05)&(df['meanchange']>0)]

sc.tl.rank_genes_groups(adata, groupby='CN_rename', reference='rest', method='t-test_overestim_var', use_raw=False, pts=True)
dedf = sc.get.rank_genes_groups_df(adata, group=None)
dedf = dedf[(dedf['pvals_adj']<0.001)&(dedf['logfoldchanges']>0.5)]
dedf = dedf[dedf['pct_nz_group']>0.4]

filtered_tf = {}

for CN in adata.obs['CN_rename'].unique():
	tmp = df[(df['group']==CN) & (df['names'].isin(dedf[dedf['group']==CN].names))]
	filtered_tf[CN] = tmp
	
source_markers = ['Hivep2','Nfat5','Nfatc2','Spi1','Nfkb1', 'Irf1','Irf5','Irf8',
				'Runx1','Runx2','Ets1','Stat2']

## Supplementary Fig 9D
with plt.rc_context({"figure.dpi": (300)}):
	dp = sc.pl.dotplot(acts, source_markers, 'assigned_cn_label', cmap = 'RdBu_r', figsize=(4.5,4), show=False,
					standard_scale='var', dendrogram=False, return_fig=True, swap_axes=True)
	dp.legend(colorbar_title='Z-scaled scores')
	
	dp.savefig("Visium_CN7_TF.pdf", transparent=True, bbox_inches="tight", dpi=300)
	plt.show()
	
fc = dedf[(dedf['group']=='CN7: Fibro-inflammatory Niche') & (dedf['names'].isin(source_markers))][['names','logfoldchanges']]
fc.index = fc['names']
fc = fc.iloc[:,1:].T
fc.index = ['CN7']

aggregated = sc.get.aggregate(adata, by='CN_rename', func=["mean",'sum'])
aggregated = pd.DataFrame(aggregated['CN7: Fibro-inflammatory Niche',:].layers['mean'], 
             columns = aggregated.var_names, index = ['CN7'])

subset_net = net[net.target.isin(dedf[(dedf['group']=='CN7: Fibro-inflammatory Niche') & (dedf['logfoldchanges']>1)]['names'])]

dc.plot_network(
	net=subset_net,
	obs = aggregated,
	act = fc,
	n_sources=['Hivep2','Nfat5','Nfatc2','Spi1','Nfkb1', 'Irf1','Irf5','Irf8',
					'Runx1','Runx2','Ets1','Stat2'],
	n_targets=12,
	node_size=80,
	c_pos_w='grey',
	c_neg_w='red',
	s_cmap='Reds',
	t_cmap='Greens',
	figsize=(8,8),
	dpi=300,
	save = 'GRN.pdf'
)