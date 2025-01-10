import anndata
import numpy as np
import scanpy as sc
import pandas as pd

### Imputation using SpaGE (https://github.com/tabdelaal/SpaGE)
from SpaGE.main import SpaGE

adata_xe = sc.read_h5ad("Xenium_all.h5ad")

def Norm_rna(x):
	return np.log(((x/np.sum(x))*1000000)+1)

def Norm_spatial(x):
	return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

db_path = 'CellChat_db.csv'
lr_network = pd.read_csv(db_path, index_col=0)
lr_network.columns = ['from', 'to', 'pathway_name', 'annotation']
lr_genes = list(set(list(lr_network['from'].unique()) + list(lr_network['to'].unique())))


meta = pd.read_table("sc_counts/meta_data.csv.gz")

### Example code of imputation on Sham data
chunks = pd.read_table("sc_counts/Control_count.dge.txt.gz",index_col=0,chunksize=1000000)
rna_data = pd.concat(chunks)
meta = pd.read_table("sc_counts/meta_data.csv.gz")
meta = meta[(meta['Group'] == "Control") & (meta['Replicates'].isin(["1_1",'1_2']))]
rna_data = rna_data[meta.index.tolist()]
gene_count = np.sum(rna_data>0, axis=1)
rna_data = rna_data.loc[gene_count >=10, :]
rna_data = rna_data.apply(Norm_rna,axis=0)

ad = adata[adata.obs.ident == "ShamR"]

spatial_data = pd.DataFrame.sparse.from_spmatrix(ad.raw.X.transpose())
spatial_data.columns = ad.obs_names.tolist()
genes = ad.var_names.tolist()
spatial_data.index = genes
cell_count = np.sum(spatial_data,axis=0)
cell_count = cell_count[cell_count>0]
spatial_data = spatial_data[cell_count.index.tolist()]
spatial_data = spatial_data.apply(Norm_spatial,axis=0)
spatial_data = spatial_data.sparse.to_dense()

selected_genes = list(set(ad.var_names.tolist() + list(lr_genes)))
gene_to_impute = rna_data.index.intersection(selected_genes)
print(len(gene_to_impute))

Imp_Genes = SpaGE(spatial_data.T, 
			rna_data.T, 
			n_pv=30, 
			genes_to_predict = gene_to_impute)

Imp_Genes.index = spatial_data.columns
ct_mtx = scipy.sparse.csr_matrix(Imp_Genes.values)
ad_impute = anndata.AnnData(ct_mtx)
ad_impute.var_names = Imp_Genes.columns
ad_impute.obs_names = spatial_data.columns

ad_impute.obs = ad.obs

ad_impute.write_h5ad("impuation_h5ad/sham_imputed.h5ad")

### Load all imputation data
ad_path = [
'hour4_imputed.h5ad',
'hour12_imputed.h5ad',
'day2_imputed.h5ad',
'week6_imputed.h5ad',
'sham_imputed.h5ad',
'day14_imputed.h5ad']

ad_list = dict()
for file in ad_path:
	ident = file.split('_')[0]
	print(ident)
	adata = sc.read_h5ad(os.path.join("impuation_h5ad", file))
	ad_list[ident] = adata

imputed_all = anndata.concat(ad_list,join = "outer")

## Supplementary Figure 10A
genes = ['Nphs2', 'Ehd3', 'Slc5a2', 'Slc22a6', 'Slc7a13', 'Havcr1', 
			'Vcam1', 'Bst1', 'Slc12a1', 'Slc12a3', 'Scnn1g', 'Aqp2', "Clnk",
			'Slc4a9','Krt19', "Akap12",'Emcn', 'Mrc2', 'Acta2', 'Cd74']
		
mat = np.zeros((len(genes), len(genes)))

for i, gene1 in enumerate(genes):
	for j, gene2 in enumerate(genes):
		expr1 = imputed_all[:, gene1].X.toarray().flatten()
		expr2 = adata_xe[:, gene2].X.toarray().flatten()
		
		corr, _ = scipy.stats.pearsonr(expr1, expr2)
		mat[i, j] = corr
		
df = pd.DataFrame(mat, index=genes, columns=genes)

plt.figure(figsize=(8,7), dpi=200)
plt.pcolor(df, cmap='RdYlBu_r')
plt.yticks(np.arange(0.5, len(df.columns), 1), df.columns)
plt.xticks(np.arange(0.5, len(df.index), 1), df.index, rotation=90)
plt.colorbar(label='Pearson Correlation', shrink=0.5)
plt.show()
