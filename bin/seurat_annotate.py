import os
import sys
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import anndata as ad
import scipy.io
from pathlib import Path
import numpy as np

if len(sys.argv) != 3:
    print("Usage: seurat_annotate.py <Path_to_humansolo_Solo.out> <Path_to_Output_CSV>")
    sys.exit(1)

# Input paths
sample_dir = sys.argv[1]
output_csv = sys.argv[2]
filtered_path = os.path.join(sample_dir, "Gene", "filtered")

# Load 10X data
mtx_file = Path(filtered_path) / "matrix.mtx"
feature_file = Path(filtered_path) / "features.tsv"
barcode_file = Path(filtered_path) / "barcodes.tsv"

# Read expression matrix and annotations
X = scipy.io.mmread(mtx_file).T.tocsr()
var_names = pd.read_csv(feature_file, header=None, sep="\t")[1].astype(str).values
obs_names = pd.read_csv(barcode_file, header=None)[0].astype(str).values

adata = ad.AnnData(X=X)
adata.var_names = var_names
adata.obs_names = obs_names
adata.var_names_make_unique()

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Preserve raw normalized data for annotation
adata.raw = adata.copy()

# Optional scaling and PCA
if adata.n_vars >= 2:
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=5, n_pcs=min(adata.n_vars, 10))
    sc.tl.leiden(adata, resolution=0.5)

# CellTypist annotation
annotation_adata = adata.raw.to_adata()
if annotation_adata.n_vars < 5:
    print("Too few genes for reliable annotation; assigning 'Unknown'.")
    adata.obs['cell_type'] = 'Unknown'
else:
    models.download_models(model='Immune_All_Low.pkl')
    predictions = celltypist.annotate(annotation_adata, model='Immune_All_Low.pkl')
    # Flatten predicted labels array to 1D
    labels = predictions.predicted_labels.values
    labels = labels.flatten()
    adata.obs['cell_type'] = labels

# Save barcode -> cell type mapping
cell_info = pd.DataFrame({
    'Barcode': adata.obs_names,
    'CellType': adata.obs['cell_type'].values
})
cell_info.to_csv(output_csv, index=False)

print("✅ CellTypist annotation complete.")
