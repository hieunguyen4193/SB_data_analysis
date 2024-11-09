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

outdir = "/media/hieunguyen/HD01/outdir/CRC1382/SBharadwaj_20240318"

path_to_main_src = "/home/hieunguyen/CRC1382/src_2023/SBharadwaj/deep_seq_batch"
samplesheet = pd.read_csv(os.path.join(path_to_main_src, "SampleSheet_for_DGE_and_CellChat.csv"))
samplesheet["full_dataset_name"] = samplesheet[["PROJECT", "dataset_name"]].apply(
    lambda x: f"{x[0]}_{x[1]}", axis = 1
)
samplesheet = samplesheet[~samplesheet["full_dataset_name"].duplicated()]

for full_name in samplesheet["full_dataset_name"].unique():
    print(f"Working on dataset {full_name}")
    project = samplesheet[samplesheet["full_dataset_name"] == full_name]["PROJECT"].values[0]
    dataset_name = samplesheet[samplesheet["full_dataset_name"] == full_name]["dataset_name"].values[0]
    path_to_input_s_obj = samplesheet[samplesheet["full_dataset_name"] == full_name]["path"].values[0]
    path_to_main_output = os.path.join(outdir, project, "data_analysis")
    
    path_to_seurat2anndata = os.path.join(path_to_main_output, "09_output", "seurat2anndata")

    object_name = f"{full_name}"
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
