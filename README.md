# IRI_Xenium_Visium

Multimodal spatial transcriptomic characterization of mouse kidney injury and repair

## Reproduce
|  Figure |Description   |  Code |
| :------------: | :------------: | :------------: |
|  Fig 2 | Xenium data preprocessing and visualization |  [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/tree/main/Analysis/1_Xenium_cell_mapping) |
| Fig 3|Xenium PT spatial trajectory analysis| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/tree/main/Analysis/2_PT_trajectory) |
| Fig 4|Xenium cell neighborhood analysis| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/tree/main/Analysis/3_Xenium_cell_neighborhood) |
| Fig 5 A,B|Xenium Visium registration and alignment| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/tree/main/Analysis/4_Xenium_Visium_integration) |
| Fig 5 E,F|Visium cell neighborhood functional characterization| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/5_Functional_cell_neighborhood/Figure5E.ipynb) |
| Fig 5 G|TF activity inference| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/5_Functional_cell_neighborhood/DecoupleR.py) |
| Fig 6 A|Xenium LR analysis| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/6_Ligand_receptor/2_Xenium_LR.py) |
| Fig 6 C-H|Visium cell-cell communication| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/6_Ligand_receptor/Figure6C-H.ipynb) |
|Supple Fig 1|Xenium and Visium QC | [link]() |
|Supple Fig 3|Xenium snRNA-seq integration | [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/1_Xenium_cell_mapping/1_Xenium_snRNA_integration.R) |
|Supple Fig 4|Xenium cell composition analysis| [link]() |
|Supple Fig 8|Pseudo-Visium analysis| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/4_Xenium_Visium_integration/Supplementary_Fig8.ipynb) |
|Supple Fig 10|Xenium imputation | [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/Analysis/6_Ligand_receptor/1_Xenium_imputation.py) |
|Supple Fig 13|Comparative analysis| [link]() |
|Supple Fig 17|Benchmark segmentation algorithm on Xenium| [link]() |
|Supple Fig 18|Benchmark spot deconvolution algorithm on Visium| [link](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/tree/main/Analysis/7_Benchmark) |


## Package

The [**pySTIM**](https://github.com/qiaoxy0/STIM) provides a comprehensive computational pipeline for the integration and analysis of spatial transcriptomics data. This tool is designed to integrate tissue histology with both high-definition spatial transcriptomics (e.g., Xenium) and transcriptome-wide spatial transcriptomics (e.g., Visium). 

### Installation

Installation using pip:\
`pip install pySTIM` 

### Import
`import pySTIM as pst`

### Usage
- **[Cross-Modality Integration of Single-Cell Resolution Xenium In Situ Data and Whole Transcriptome Visium Data](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/notebooks/Integration.ipynb)**
- **[Cell Neighborhoods Mapping and Functional Characterization](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/notebooks/CN_analysis.ipynb)**
- **[Inference of Ligand-Receptor Interactions](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/notebooks/LR_analysis.ipynb)**
- **[Core Visualization Functions](https://github.com/TheHumphreysLab/IRI_Xenium_Visium/blob/main/notebooks/Visualizations.ipynb)**
