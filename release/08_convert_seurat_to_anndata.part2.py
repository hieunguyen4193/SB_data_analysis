import scanpy as sc
import numpy as np
import os
import anndata2ri
import pathlib
from scipy import io
import anndata#
import pandas as pd
from tqdm import tqdm
import argparse
import sys

# Activate the anndata2ri conversion between SingleCellExperiment and AnnData
# anndata2ri.activate()

#Loading the rpy2 extension enables cell magic to be used
#This runs R code in jupyter notebook cells
# %load_ext rpy2.ipython

sc.settings.verbosity = 3
# sc.logging.print_versions()

import warnings
warnings.filterwarnings("ignore")

outdir = "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20250102"

path_to_main_src = "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/release"
samplesheet = pd.read_csv(os.path.join(path_to_main_src, "SampleSheet_all_seurat_objects.csv"))
samplesheet["dataset_name"] = samplesheet.apply(
    lambda row: f"{row['dataset_name']}_reIntegration" if row["reIntegration"] == "yes" else row["dataset_name"],
    axis=1
)

for index, row in samplesheet.iterrows():
    PROJECT = row['PROJECT']
    dataset_name = row['dataset_name']
    path_to_s_obj = row['path']
    re_integration = row['reIntegration']
    
    if re_integration in ["yes", ""]:
        reduction_name = "cca_UMAP"
    else:
        reduction_name = "SCT_UMAP"
    
    print(f"Working on PROJECT {PROJECT}, dataset name {dataset_name}")
    
    path_to_s_obj = path_to_s_obj.replace(".rds", ".addedInfo.rds")
    
    path_to_main_output = os.path.join(outdir, PROJECT, "data_analysis")
    path_to_08_output = os.path.join(path_to_main_output, "08_output", dataset_name)
    os.makedirs(path_to_08_output, exist_ok=True)
    
    path_to_seurat2anndata = os.path.join(path_to_08_output, "seurat2anndata")
    os.makedirs(path_to_seurat2anndata, exist_ok=True)

    object_name = f"{PROJECT}_{dataset_name}"
    ##### MAIN CONVERSION FUNCTIONS
    print(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(object_name)))
    X = io.mmread(os.path.join(path_to_seurat2anndata, "counts_{}.mtx".format(object_name)))

    # create anndata object
    adata = anndata.AnnData(X=X.transpose().tocsr())

    # load cell metadata:
    cell_meta = pd.read_csv(os.path.join(path_to_seurat2anndata, "metadata_{}.csv".format(object_name)))

    # load gene names:
    with open(os.path.join(path_to_seurat2anndata, "gene_names_{}.csv".format(object_name)), 'r') as f:
        gene_names = f.read().splitlines()

    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    adata.var.index = gene_names

    # load dimensional reduction:
    pca = pd.read_csv(os.path.join(path_to_seurat2anndata, "pca_{}.csv".format(object_name)))
    pca.index = adata.obs.index

    # set pca and umap
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

    # save dataset as anndata format
    adata.write(os.path.join(path_to_seurat2anndata, '{}.h5ad'.format(object_name)))
