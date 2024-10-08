import scanpy as sc
import numpy as np
from .utils import (
    add_atac_data_to_sc, 
    add_cell_types_to_sc, 
    filter_doublets, 
    add_regions
)
from .config import (
    d_sample_to_condition_vz
)


def preprocess_data_sc(filename, sample, use_annotations=True, filename_annotation=None):
    """
    Preprocess single-cell data by adding ATAC data and filtering.

    Parameters:
    - filename (str): Path to the 10X H5 file.
    - sample (dict): Information about the sample, including ATAC filename.
    - use_annotations (bool): Whether to use cell annotations.
    - filename_annotation (str): Path to the annotation file.

    Returns:
    - AnnData: Preprocessed AnnData object.
    """
    adata = sc.read_10x_h5(filename)
    adata = add_atac_data_to_sc(adata, sample["atac_filename"])

    if use_annotations:
        adata = add_cell_types_to_sc(
            adata, sample["library_name"], sample["sample_name"], filename_annotation
            )
    else:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        adata = filter_doublets(adata)
   
    duplicates = adata.var_names[adata.var_names.duplicated()]
    if len(duplicates) > 0:
        print("Duplicates found.")
        print("Removing the data duplicated...")
        adata = adata[:, ~adata.var_names.duplicated()].copy()
        print("Size AFTER dropping duplicated:")
        adata.var_names_make_unique()

    adata.obs.index = adata.obs.index + "." + sample["library_name"]
    adata.obs["Sample"] = sample["library_name"]
    return adata


def preprocess_data_vz(folder_name, idx_sample, filename_morph):
    """
    Preprocess spatial data and add region metadata.
    
    Parameters:
    - folder_name (str): Path to the 10X mtx folder.
    - idx_sample (str): Identifier for the sample.
    - filename_morph (str): Path to the morphological file for regions.
    
    Returns:
    - AnnData: Preprocessed spatial AnnData object.
    """
    adata = sc.read_10x_mtx(folder_name)
    
    sc.pp.filter_cells(adata, min_genes=100)

    adata.var['mt'] = adata.var_names.str.startswith('Mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, .02)
    upper_lim_pct = np.quantile(adata.obs.pct_counts_mt.values, .98)
    lower_lim_pct = np.quantile(adata.obs.pct_counts_mt.values, .02)

    print(f"genes_by_counts: {lower_lim} to {upper_lim}, pct_mt: {lower_lim_pct} to {upper_lim_pct}")
    print(idx_sample + ": Original size: " + str(len(adata)))
    adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]
    print(idx_sample + ": After filtering per n_genes_by_counts: " + str(len(adata)))
    adata = adata[adata.obs.pct_counts_mt < upper_lim_pct, :]
    print(idx_sample + ": After filtering per pct_counts_mt: " + str(len(adata)))

    adata.obs['Sample'] = idx_sample
    adata.obs['Condition'] = d_sample_to_condition_vz[idx_sample]
    adata = add_regions(adata, filename_morph)
    return adata

