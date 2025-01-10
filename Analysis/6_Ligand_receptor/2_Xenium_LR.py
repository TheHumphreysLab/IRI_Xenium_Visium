import pySTIM as pst
import scanpy as sc
import numpy as np
import pandas as pd
import concurrent.futures

celltyps = ['FR_PT', 'Fib', 'Immune','Inj_PT']

df = pd.DataFrame({
    "celltype_sender": np.repeat(celltyps, len(celltyps)),
    "celltype_receiver": list(celltyps)*len(celltyps),
})
df = df[df['celltype_sender'] != df['celltype_receiver']]
df["celltype_pair"] = df["celltype_sender"].str.cat(
    df["celltype_receiver"], sep="-")

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
    
    
def compute_cci_for_pair(pair, adata, lr_network):
    
    s, r = pair.split(sep='-')
    print(s, '->', r)
    result = pst.compute_cci(
        adata=adata,
        group='celltype',
        lr_network=lr_network,
        sender=s,
        receiver=r,
        contact_radius=30,
        p_value_threshold=0.05,
        spatial_key='spatial'
    )
    return pair, result

all_result = {}

for i in ['hour4', 'hour12', 'day2', 'week6', 'sham', 'day14']:
    adata = ad_list[i]
    
    pairs = df['celltype_pair'].tolist()
    
    res = {}
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_pair = {executor.submit(compute_cci_for_pair, pair, adata, lr_network): pair for pair in pairs}
        
        for future in concurrent.futures.as_completed(future_to_pair):
            pair = future_to_pair[future]
            try:
                pair_result = future.result()
                res[pair_result[0]] = pair_result[1]
            except Exception as exc:
                print(f'{pair} generated an exception: {exc}')
                
    all_result[i] = res
    

df_list = []

for i in ['hour4', 'hour12', 'day2', 'week6', 'sham', 'day14']:
    res = all_result[i]
    result = pd.DataFrame(columns=all_result['week6']['FR_PT-Immune'].columns)
    for l in df.index:
        if res[df['celltype_pair'][l]] is not None and not res[df['celltype_pair'][l]].empty:
            res[df['celltype_pair'][l]] = res[df['celltype_pair'][l]].sort_values('LR_Score', ascending=False)
            result = pd.concat([result, res[df['celltype_pair'][l]]], axis=0, join='outer')
    result['time'] = i
    df_list.append(result)
    
final_df = pd.concat(df_list, axis=0, join='outer')
final_df = final_df[(final_df.LR_Score>0.1)&(final_df.co_exp_count>500)&(final_df.fold_change>1.2)]
final_df = final_df.reset_index(drop=True)