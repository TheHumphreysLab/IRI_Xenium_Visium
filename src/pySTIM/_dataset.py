import os
import pandas as pd
import urllib.request
import anndata as ad
import json

class AMetadata:
    def __init__(self, name, doc_header, url, file_path):
        self.name = name
        self.doc_header = doc_header
        self.url = url
        self.file_path = file_path
        
    def download(self):
        os.makedirs(os.path.dirname(self.file_path), exist_ok=True)
        if not os.path.exists(self.file_path):
            print(f"Downloading {self.name} dataset...")
            urllib.request.urlretrieve(self.url, self.file_path)
            print(f"Downloaded {self.name} dataset to {self.file_path}")
        else:
            print(f"{self.name} dataset already exists at {self.file_path}")
            
    def load_anndata(self):
        self.download()
        adata = ad.read_h5ad(self.file_path)
        return adata
    
    def load_dataframe(self):
        self.download()
        df = pd.read_csv(self.file_path)
        return df
    
    def load_json(self):
        self.download()
        with open(self.file_path, 'r') as json_file:
            poly_data = json.load(json_file)
        return poly_data
    

Xenium_iri_poly = AMetadata(
    name="Xenium-polygon",
    doc_header="Smoothed Baysor cell segmentation polygon data of Xenium day14 sample",
    url="https://figshare.com/ndownloader/files/47106217",
    file_path="data/Day14R_poly_data_smoothed.json"  
)

Visium_sham = AMetadata(
	name="Visium-Sham",
	doc_header="Pre-processed mouse Visium sham (rep1) kidney sample.",
	url="https://figshare.com/ndownloader/files/47095234",
	file_path="data/Visium_Sham.h5ad"  
)

Xenium_sham = AMetadata(
	name="Xenium-Sham",
	doc_header="Pre-processed mouse Xenium sham (rep1) kidney sample.",
	url="https://figshare.com/ndownloader/files/47095570",
	file_path="data/Xenium_Sham.h5ad"  
)

Xenium_iri = AMetadata(
	name="Xenium-iri",
	doc_header="Pre-processed mouse Xenium Day14 (rep1) kidney sample.",
	url="https://figshare.com/ndownloader/files/47108383",
	file_path="data/Xenium_IRI.h5ad"  
)

lr_db = AMetadata(
    name="CellChat ligand receptor database (Secreted Signaling only)",
    doc_header="Ligand-Receptor pairs database.",
    url="https://figshare.com/ndownloader/files/47095957",
    file_path="data/CellChatDB_lr_network_mouse.csv"  
)