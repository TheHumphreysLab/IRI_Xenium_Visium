import pandas as pd
import numpy as np
from tqdm import tqdm
from typing import Any

def create_polygons(adata: Any) -> Any:
    """
    Create polygons for each unique cell ID based on segmentation info.
    
    Parameters:
    - adata (AnnData): An AnnData object containing segmentation data under the key 'seg' in its 'uns' slot.

    Returns:
    - AnnData: The input AnnData object, updated with polygons for each cell under the key 'poly' in its 'uns' slot.
    """

    tenx_seg = adata.uns['seg']
    grouped = tenx_seg.groupby('cell_id')

    all_cells = sorted(list(set(tenx_seg.cell_id)))
 
    new_lib = [{
        'coordinates': [grouped.get_group(cell).iloc[:, 1:].to_numpy().tolist()],
        'type': 'Polygon'
    } for cell in tqdm(all_cells)]
 
    new_poly = dict(zip(map(str, all_cells), new_lib))
    adata.uns['poly'] = new_poly
    
    return adata
