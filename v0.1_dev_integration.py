import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
sys.path.append('./')  # uncomment for local import

import tangram as tg
tg.__version__

import scvelo as scv

import warnings

warnings.filterwarnings('ignore')
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

import scvi

import doubletdetection
from scipy.stats import median_abs_deviation as mad
from typing import List


clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=10
)

from filenames import (
    ROOT_FOLDER_DIR,
    FILE_ANNOTATIONS_MA,
    FILE_ANNOTATIONS_FA, 
    FILE_ANNOTATIONS_MP,  
    FILE_ANNOTATIONS_FP,  
    FILE_ANNOTATIONS_VZ_R, 

    FILE_ATAC_ANNOTATIONS_MC1, 
    FILE_ATAC_ANNOTATIONS_MC2, 
    FILE_ATAC_ANNOTATIONS_ME1, 
    FILE_ATAC_ANNOTATIONS_ME2, 

    FILE_ATAC_ANNOTATIONS_FC1, 
    FILE_ATAC_ANNOTATIONS_FE1, 

    FILE_ATAC_ANNOTATIONS_PMC1, 
    FILE_ATAC_ANNOTATIONS_PMC2, 
    FILE_ATAC_ANNOTATIONS_PME1, 
    FILE_ATAC_ANNOTATIONS_PME2, 

    FILE_ATAC_ANNOTATIONS_PFC, 
    FILE_ATAC_ANNOTATIONS_PFE,
)


d_sample_to_dirname_sc = {
    "3886": "P7_filtered/P7.SC_3886_",
    "3887": "P7_filtered/P7.SC_3887_",
    "3888": "P7_filtered/P7.SC_3888_",
    "3889": "P7_filtered/P7.SC_3889_",
    "4203": "P7_filtered/P7.SC_4203_",
    "4204": "P7_filtered/P7.SC_4204_",
    "4205": "P7_filtered/P7.SC_4205_",
    "4206": "P7_filtered/P7.SC_4206_",
    "4207": "P7_filtered/P7.SC_4207_",
    "4208": "P7_filtered/P7.SC_4208_",
    "4209": "P7_filtered/P7.SC_4209_",
    "4210": "P7_filtered/P7.SC_4210_"
}

d_sample_to_categ_sc = {
    "3886": "CM",
    "3887": "CM",
    "3888": "EcigM",
    "3889": "EcigM",
    "4203": "CM",
    "4204": "CM",
    "4205": "EcigM",
    "4206": "EcigM",
    "4207": "CF",
    "4208": "EcigF",
    "4209": "CF",
    "4210": "EcigF"
}

sample_to_dataset = {
    '3886': 'MA',
    '3887': 'MA',
    '3888': 'MA',
    '3889': 'MA',
    '4203': 'MP',
    '4204': 'MP',
    '4205': 'MP',
    '4206': 'MP',
    '4207': 'FA',
    '4208': 'FA',
    '4209': 'FP',
    '4210': 'FP',
}

d_sample_to_dirname_vz = {
    "P94A1":"p94_a1_run2",
    "P95A1":"p95_a1_run2",
    "P100A1":"p100_a1_run2",
    "P81A2":"p81_a2_run2",
    "P81A1":"p81_a1_run2",
    "P83A1":"p83_a1_run2",
    "P112A1":"p112_a1_run2",
    "P113A1":"p113_a1_run2",
    "P80A1":"p80_a1_run2",
    "P96A1":"p96_a1_run2",
    "P82A1":"p82_a1_run2",
    "P96A1":"p96_a1_run2",
    "P80P1":"p80_p1_run2",
    "P94P1":"p94_p1_run2",
    "P96P1":"p96_p1_run2",
    "P82P1":"p82_p1_run2",
}

d_sample_to_condition_vz = {
    "P94A1":"CM",
    "P95A1":"CM",
    "P100A1":"CM",
    "P81A2":"EcigM",
    "P81A1":"EcigM",
    "P83A1":"EcigM",
    "P112A1":"EcigM",
    "P113A1":"EcigM",
    "P80A1":"EcigM",
    "P96A1":"CF",
    "P94P1":"CM",
    "P80P1":"EcigM",
    "P82A1":"EcigF",
    "P96P1":"CF",
    "P82P1":"EcigF",
}


def add_atac_data_to_sc(adata, filename_annotation):
    """
    Adds atac information per cell to an AnnData object by merging metadata 
    from result folder.
    """
    print(f"Loading metadata: {filename_annotation}")
    adata_atac = sc.read_h5ad(filename_annotation)
    print(adata_atac.obs) 
    #adata_atac = adata_atac.obs.drop(columns=["total_counts","n_counts","leiden"])

    # Join with the AnnData object's observations (adata.obs)
    print(f"Size before merging: {adata.shape}")
    # Step 1: Get the intersection of indices
    common_index = adata.obs.index.intersection(adata_atac.obs.index)

    # Step 2: Subset both AnnData objects to only keep common indices
    adata_subset = adata[common_index, :].copy()
    adata_atac_subset = adata_atac[common_index, :].copy()

    # Step 3: Merge the obs dataframes from both AnnData objects on the common index
    print("Merging ATAC with scRNA data...")
    merged_obs = pd.concat([adata_subset.obs, adata_atac_subset.obs], axis=1)

    # Step 4: Create a new AnnData object with the merged obs
    # We use the X matrix from adata1_subset, but you can modify this to use adata2's or combine as needed
    adata_merged = sc.AnnData(X=adata_subset.X, obs=merged_obs, var=adata_subset.var)
    adata_merged.obs = adata_merged.obs.drop(columns=["total_counts", "n_counts", "leiden"])
    
    print(adata_merged.obs)

    return adata_merged


def add_cell_types_to_sc(adata, idx_name, sample_name, filename_annotation):
    """
    Adds cell types to an AnnData object by merging metadata annotations from a CSV file.
    
    Parameters:
    - adata: AnnData object to be annotated.
    - idx_name: Name of the index column in the AnnData object.
    - sample_name: Sample identifier to filter the annotation file.
    - filename_annotation: Path to the CSV file with annotations.
    
    Returns:
    - Updated AnnData object with cell types added.
    """
    posterior = filename_annotation == FILE_ANNOTATIONS_MP or filename_annotation == FILE_ANNOTATIONS_FP 
    print(f"Loading metadata: {filename_annotation}")
    
    # Load CSV with proper header handling
    annotations = pd.read_csv(filename_annotation)
    #annotations = annotations.tail(annotations.shape[0] -1) 
    # Filter the annotations by sample_name
    annotations = annotations[annotations['Sample'] == sample_name]
    print(annotations)
    # Properly extract and clean CellID
    prefix_index_size = 5 if filename_annotation == FILE_ANNOTATIONS_MP else 4
    column_index_name = "Unnamed: 0"
    if filename_annotation == FILE_ANNOTATIONS_MA or filename_annotation == FILE_ANNOTATIONS_FA:
        column_index_name = "X"
    
    annotations["CellID"] = annotations[column_index_name].str[prefix_index_size:]
    annotations["CellID"] = annotations["CellID"].str.strip()
    annotations["CellID"] = annotations["CellID"].astype(str)
    # Rename 'cell.type' column to 'cell_type' for consistency
    annotations.rename(columns={'cell.type': 'cell_type'}, inplace=True)
    
    # Remove unnecessary columns before merging
    annotations.drop(columns=[column_index_name, "Sample", "Group"], inplace=True)
    
    # Set 'CellID' as the index
    annotations.set_index("CellID", inplace=True)
    
    # Join with the AnnData object's observations (adata.obs)
    print(f"Size before merging: {adata.shape}")
    adata.obs.index = adata.obs.index.str.strip()
    adata.obs = adata.obs.join(annotations, how="left")

    # Check how many cells were successfully annotated
    annotated_cells = adata.obs['cell_type'].notna().sum()

    # Update annotations:
    if not posterior:
        adata.obs['cell_type'] = adata.obs['cell_type'].replace('exLayer4','exLayer5/6-tmp')
        adata.obs['cell_type'] = adata.obs['cell_type'].replace('exLayer5/6','exLayer4')
        adata.obs['cell_type'] = adata.obs['cell_type'].replace('exLayer5/6-tmp','exLayer5/6')

    print(f"Number of annotated cells: {annotated_cells}")
    
    # Remove NaN and Print the updated AnnData object
    adata = adata[adata.obs['cell_type'].notna()].copy()
    print(f"Size after merging: {adata.shape}")
    
    return adata


def add_regions(intg_adata_slice, filename):
# Read the CSV file
    csv_file = pd.read_csv(filename)
# Extract the first column (index column)
    index_column = csv_file.iloc[:, 0].astype(str)
    filtered_adata = intg_adata_slice[intg_adata_slice.obs_names.isin(index_column)]
    filtered_adata_df = filtered_adata.to_df()
    csv_file = csv_file.set_index(csv_file.columns[0])
    #csv_file.index = csv_file.index.astype(str)
    result_df = csv_file.join(filtered_adata_df, how='inner')
    result_df_col6 = result_df[["region"]]
# Ensure both have the same index (common key)
    intg_adata_slice.obs_names = intg_adata_slice.obs_names.astype(str)
    result_df_col6.index = result_df_col6.index.astype(str)
# Subset the AnnData object to include only rows that match the metadata DataFrame's index
    intg_adata_res = intg_adata_slice[intg_adata_slice.obs_names.isin(result_df_col6.index)]
# Add the metadata columns to the AnnData object's obs

    print(intg_adata_res)
    for col in result_df_col6.columns:
        intg_adata_res.obs[col] = result_df_col6.loc[intg_adata_res.obs_names, col]

    return intg_adata_res


def preprocess_data_sc(filename, sample, use_annotations = True, filename_annotation = FILE_ANNOTATIONS_MA):
    
    adata = sc.read_10x_h5(filename) 
    adata = add_atac_data_to_sc(adata, sample["atac_filename"])

    if use_annotations:
        adata = add_cell_types_to_sc(adata, sample["library_name"], sample["sample_name"], filename_annotation)
    else:
        sc.pp.filter_cells(adata, min_genes=200)
        adata.var['mt'] = adata.var_names.str.startswith('Mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
        lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, .02)
        upper_lim_pct = np.quantile(adata.obs.pct_counts_mt.values, .98)
        lower_lim_pct = np.quantile(adata.obs.pct_counts_mt.values, .02)
        print(f"genes_by_counts: {lower_lim} to {upper_lim}, pct_mt: {lower_lim_pct} to {upper_lim_pct}")
        print(sample["library_name"]+ ": Original size: " + str(len(adata)))
        adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]
        print(sample["library_name"] + ": After filtering per n_genes_by_counts: " + str(len(adata)))
        adata = adata[adata.obs.pct_counts_mt < upper_lim_pct, :]
        print(sample["library_name"] + ": After filtering per pct_counts_mt: " + str(len(adata)))

        # doublet detection:
        doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.50)
        doublet_score = clf.doublet_score()
        adata.obs["doublet"] = doublets
        adata.obs["doublet_score"] = doublet_score
        adata.uns['doublets_removed'] = adata.obs.doublet.sum()
        print("Doublets count: " + str((adata.uns['doublets_removed'])) + " final adata size: ")
        print(str(len(adata)))
        adata = adata[adata.obs.doublet == 0]
        print("Size of AnnData:")
        print(str(len(adata)))

    print("single-cell data loaded and annotated...")
    print(adata.obs)
    duplicates = adata.var_names[adata.var_names.duplicated()]
    if len(duplicates) > 0:
        print("Duplicates found.")
        print(duplicates)
        print("Removing the data duplicated...")
        adata = adata[:, ~adata.var_names.duplicated()].copy()
        print("Size AFTER dropping duplicated.")
        print(str(len(adata)))

    # index update:
    adata.obs.index = adata.obs.index + "." + sample["library_name"]
    # define sample attribute
    adata.obs["Sample"] = sample["library_name"]

    return adata


def preprocess_data_vz(folder_name, idx_sample, filename_morph):
    adata = sc.read_10x_mtx(folder_name)
    sc.pp.filter_cells(adata, min_genes=100)
    #sc.pp.filter_genes(adata, min_cells=3)
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



EPOCHS=600
TOPGENESCOUNT=5000
VERSION = "v5.cell.g" + str(TOPGENESCOUNT) + ".e" + str(EPOCHS)
FILE_SC_ADATA = VERSION + ".h5ad"


def add_coords(adata_lcl, sample_dirname, reverse_coords, angle = 0):
    filename_positions = ROOT_FOLDER_DIR +  "SpaceRangerTissuePositions/" + sample_dirname + "tissue_positions.csv"
    print(f"Loading tissue positions: {filename_positions}")
    positions = pd.read_csv(filename_positions)
    positions = positions.tail(positions.shape[0] -1)
    positions.index = positions.index.rename("CellID")
    positions["CellID"] = positions["barcode"]
    positions.drop("barcode", inplace=True, axis=1)
    positions.index = positions["CellID"]
    
    positions['pxl_col_in_fullres'] = positions['pxl_col_in_fullres'].astype('float32')
    positions['pxl_row_in_fullres'] = positions['pxl_row_in_fullres'].astype('float32')

    adata_lcl.obs = adata_lcl.obs.join(positions, how="left")
    adata_lcl.obs['x'] = adata_lcl.obs['pxl_row_in_fullres']
    adata_lcl.obs['y'] = adata_lcl.obs['pxl_col_in_fullres']
    
    if reverse_coords:
        x_coords = adata_lcl.obs['x']
        y_coords = adata_lcl.obs['y']
        x, y = generate_symmetric_coordinates(x_coords, y_coords, False)
        adata_lcl.obs['x'] = x
        adata_lcl.obs['y'] = y
    if angle > 0:
        x_coords = adata_lcl.obs['x']
        y_coords = adata_lcl.obs['y']
        x, y = rotate_coordinates(x_coords, y_coords, angle)
        adata_lcl.obs['x'] = x
        adata_lcl.obs['y'] = y    

    adata_lcl.obs.drop(
        columns=['CellID', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
        inplace=True,
    )
    return adata_lcl


def rotate_coordinates(x, y, angle_degrees = 180):
    # Ensure x and y are numpy arrays
    x = np.array(x)
    y = np.array(y)

    # Calculate the midpoint of the dataset on the x-axis and y-axis
    midpoint_x = np.mean(x)
    midpoint_y = np.mean(y)

    # Convert the angle from degrees to radians
    angle_radians = np.radians(angle_degrees)

    # Create the rotation matrix
    rotation_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians), np.cos(angle_radians)]
        ])

    # Center the points around the midpoint
    centered_points = np.vstack((x - midpoint_x, y - midpoint_y))

    # Rotate the points
    rotated_points = rotation_matrix.dot(centered_points)

    # Shift the points back to the original midpoint
    rotated_x = rotated_points[0, :] + midpoint_x
    rotated_y = rotated_points[1, :] + midpoint_y

    # Convert the rotated coordinates back to lists
    return rotated_x.tolist(), rotated_y.tolist()


def generate_symmetric_coordinates(x, y, is_x_axis=True):
    # Ensure x and y are numpy arrays
    x = np.array(x)
    y = np.array(y)

    if is_x_axis:
        # Calculate the midpoint of the dataset on the x-axis
        midpoint_x = np.mean(x)

        # Generate the symmetric x coordinates
        symmetric_x = 2 * midpoint_x - x

        # The y coordinates remain the same for symmetry around the x-axis midpoint
        symmetric_y = y
    else:
        # Calculate the midpoint of the dataset on the x-axis
        midpoint_y = np.mean(y)

        # Generate the symmetric x coordinates
        symmetric_y = 2 * midpoint_y - y

        # The y coordinates remain the same for symmetry around the x-axis midpoint
        symmetric_x = x

    return symmetric_x.tolist(), symmetric_y.tolist()


def run_integration(ds_name: str, vz_sample, sc_sample_ids, file_annotations = None):

    adata_vz_list = []
    file_mapped_sp_sc = "mapped.sc.vz." + ds_name + "." + vz_sample["name"] + "." + VERSION + ".tangram.h5ad"
    file_sc = "sc." + ds_name + "." + FILE_SC_ADATA
    file_sp = "VZ." + ds_name + "." + vz_sample["name"] + "." + FILE_SC_ADATA

    # Creating ANNDATA for VZ data:
    print("[VZ] Adding " + vz_sample["name"])
    filename = ROOT_FOLDER_DIR + d_sample_to_dirname_vz[vz_sample["name"]] + "/outs/filtered_feature_bc_matrix"
    # Loading VZ dataset
    adata_vz_lcl = preprocess_data_vz(filename, vz_sample["name"], vz_sample["file_morph_mapping"])
    adata_vz_list.append(
        add_coords(
            adata_vz_lcl, d_sample_to_dirname_vz[vz_sample["name"]], vz_sample["flipped"], vz_sample["rotation_angle"] 
        )
    )
        
    adata_sp= sc.concat(adata_vz_list)
    print("Merging VZ data done.")
    print(adata_sp.obs)
    #adata_vz.write_h5ad("VZ." + FILE_SC_ADATA)
    use_annotations = file_annotations is not None
    cluster_type = "cell_type" if use_annotations else "leiden"
    if not os.path.exists(file_mapped_sp_sc):
        print("Saved mapped spatial vs single-cell data not found... It will be created from scratch.")
        if not os.path.exists(file_sc):
            print("Saved single-cell data not found... It will be created from scratch.")
            # Read sc data
            adata_sc_list = []
            for sample_code in sc_sample_ids:
                print("[SC] Adding " + sample_code["library_name"])
                filename = ROOT_FOLDER_DIR + d_sample_to_dirname_sc[sample_code["library_name"]] + "filtered_feature_bc_matrix.h5"
                adata_sc_list.append(preprocess_data_sc(filename, sample_code, use_annotations, file_annotations))
            
            adata_sc = sc.concat(adata_sc_list)
            adata_sc.obs["Condition"] = adata_sc.obs.Sample.map(d_sample_to_categ_sc)

            sc.pp.filter_genes(adata_sc, min_cells=50)
            print("Size after removing genes in low number of cells: " + str(adata_sc))
            print(f"Observation attributes: {adata_sc.obs_names}")
            adata_sc.obs.groupby("Sample").count()
            adata_sc.raw_unnormalized = adata_sc.copy()

            sc.pp.normalize_total(adata_sc)
            sc.pp.log1p(adata_sc)
            adata_sc.raw = adata_sc.copy()
            sc.pp.highly_variable_genes(adata_sc, n_top_genes=TOPGENESCOUNT)

            # Run MODEL SCVI & leiden
            if False:
                cluster_type = "leiden"
                sc.settings.n_jobs = 100
                scvi.model.SCVI.setup_anndata(
                    adata_sc,
                    categorical_covariate_keys = ["Sample"],
                    continuous_covariate_keys = ["pct_counts_mt"]
                )
                model = scvi.model.SCVI(adata_sc)
                model.train(max_epochs = 50)
                model.get_latent_representation().shape
                adata_sc.obsm["X_scVI"] = model.get_latent_representation() 
                sc.pp.neighbors(adata_sc, use_rep = "X_scVI")

                # Clustering
                print("Clustering sc data...")
                sc.tl.leiden(adata_sc, resolution= 2.5)
            else:
                cluster_type = "cell_type"

            print("Saving sc-data h5ad")
            adata_sc.write_h5ad(file_sc)
        else:
            print("Loading the saved single-cell h5ad data ...")
            adata_sc = sc.read_h5ad(file_sc)

        hvg_genes = adata_sc.var.index[adata_sc.var['highly_variable']].tolist()
        #adata
        tg.pp_adatas(adata_sc, adata_sp, genes=hvg_genes)

        print(f"#training genes in sc: {adata_sc.uns['training_genes'][:100]} ...")
        print(f"#training genes in sp: {adata_sp.uns['training_genes'][:100]} ...")

        adata_sc.write_h5ad(file_sc)
        adata_sp.write_h5ad(file_sp)

        print("Running mapping...")
        ad_map = tg.map_cells_to_space(
            adata_sc=adata_sc,
            adata_sp=adata_sp,
            device='cpu',
            num_epochs=EPOCHS,
            #mode='clusters',
            #cluster_label=cluster_type,
            #device='cuda:0',
        )
        
        ad_map.write_h5ad(file_mapped_sp_sc)
        adata_sc.write_h5ad(file_sc)
        adata_sp.write_h5ad(file_sp)
        print("Saved ad_map, ad_sc, and ad_sp.")
    else:
        print("Found saved ad_map, ad_sc, and ad_sp. Loading now...")
        adata_sc = sc.read_h5ad(file_sc)
        ad_map = sc.read_h5ad(file_mapped_sp_sc)
        adata_sp = sc.read_h5ad(file_sp)

    try:
        print("Running spatial cell projections...")
        tg.project_cell_annotations(ad_map, adata_sp, annotation=cluster_type)
        annotation_list = sorted(list(pd.unique(adata_sc.obs[cluster_type])))
        tg.plot_cell_annotation_sc(
            adata_sp, 
            annotation_list,
            x='x', 
            y='y',
            spot_size= 200, 
            scale_factor=0.1, 
            perc=0.001, 
            filename=".tg." + ds_name + "." + VERSION + "." + cluster_type + ".scvi.pdf"
        )

        ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
        print(f"Gene expression on overlap data: {ad_ge}")

        # Saving files:
        #adata_sc.write_h5ad(file_sc)
        #adata_sp.write_h5ad(file_sp)

        genes_astrocyte = ['gfap', 'slc1a2', 'aqp4']
        genes_OPC = ['pdgfra', 'itpr2', 'tcf7l2']
        genes_oligodendrocyte = ['olig2', 'bmp4']
    
        print("Printing figures for some regions/cell types...")
        fig = tg.plot_genes_sc(
            genes_astrocyte, 
            adata_measured=adata_sp, 
            adata_predicted=ad_ge, 
            spot_size= 200, 
            scale_factor=0.1, 
            perc = 0.001, 
            return_figure=True
        )
        fig.savefig("_" + ds_name + "." + VERSION + "_genes.astr.pdf")
        fig = tg.plot_genes_sc(
            genes_OPC, 
            adata_measured=adata_sp, 
            adata_predicted=ad_ge, 
            spot_size= 200, 
            scale_factor=0.1, 
            perc = 0.001, 
            return_figure=True
        )
        fig.savefig("_" + ds_name + "." + VERSION + "_genes.opc.pdf")
        fig = tg.plot_genes_sc(
            genes_oligodendrocyte, 
            adata_measured=adata_sp, 
            adata_predicted=ad_ge, 
            spot_size= 200, 
            scale_factor=0.1, 
            perc = 0.001, 
            return_figure=True
        )
        fig.savefig("_" + ds_name + "." + VERSION + "_genes.olig.pdf")

        df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_sp, adata_sc)
        print(f"Similarity scores on all genes: {df_all_genes}")

        print("Plotting AUC graph...")
        tg.plot_auc(df_all_genes, "_" + ds_name + "." + VERSION + ".auc.pdf") 

    except Exception as e:
        print(e)
    finally:
        print(f"Exception occured. But continuing to next task.")
    


datasets =[
    {
        "name": "MAC_SC_VZ_ATAC", 
        "params_vz_sample_ids": {"name":"P94A1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/MA94.cor.csv",  "rotation_angle": 0, "flipped": True},
        "params_sc_sample_ids": [
            {"library_name":"3886", "sample_name":"MC1", "atac_filename": FILE_ATAC_ANNOTATIONS_MC1}, 
            {"library_name":"3887", "sample_name":"MC2", "atac_filename": FILE_ATAC_ANNOTATIONS_MC2}
        ],
        "sc_annotation": FILE_ANNOTATIONS_MA,
    },
    {
        "name":"MAE_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P80A1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/MA80.cor.csv",  "rotation_angle": 0, "flipped": True}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"3888", "sample_name":"ME1", "atac_filename": FILE_ATAC_ANNOTATIONS_ME1},
            {"library_name":"3889", "sample_name":"ME2", "atac_filename": FILE_ATAC_ANNOTATIONS_ME2}
        ],
        "sc_annotation": FILE_ANNOTATIONS_MA,
    },
    {
        "name": "FAC_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P96A1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/FA96.cor.csv",  "rotation_angle": 0, "flipped": True}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"4207", "sample_name":"FC1", "atac_filename": FILE_ATAC_ANNOTATIONS_FC1},
        ],
        "sc_annotation": FILE_ANNOTATIONS_FA,
    },
    {
        "name": "FAE_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P82A1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/FA82.cor.csv",  "rotation_angle": 0, "flipped": True}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"4208", "sample_name":"FE1", "atac_filename": FILE_ATAC_ANNOTATIONS_FE1},
        ],
        "sc_annotation": FILE_ANNOTATIONS_FA,
    }
    ,
    {
        "name": "MPC_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P94P1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/MP94.cor.csv",  "rotation_angle": 0, "flipped": True}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"4203", "sample_name":"MC1", "atac_filename": FILE_ATAC_ANNOTATIONS_PMC1},
            {"library_name":"4204", "sample_name":"MC2", "atac_filename": FILE_ATAC_ANNOTATIONS_PMC2},
        ],
        "sc_annotation": FILE_ANNOTATIONS_MP
    },
    {
        "name": "MPE_SC_VZ_ATAC",        
        "params_vz_sample_ids": {"name":"P80P1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/MP80.cor.csv",  "rotation_angle": 0, "flipped": True}, 
        "params_sc_sample_ids": [
            {"library_name":"4205", "sample_name":"ME1", "atac_filename": FILE_ATAC_ANNOTATIONS_PME1},
            {"library_name":"4206", "sample_name":"ME2", "atac_filename": FILE_ATAC_ANNOTATIONS_PME2},
        ],
        "sc_annotation": FILE_ANNOTATIONS_MP
    },
    {
        "name": "FPC_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P96P1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/FP96.cor.csv", "rotation_angle": 270, "flipped": False}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"4209", "sample_name":"PFC", "atac_filename": FILE_ATAC_ANNOTATIONS_PFC},
        ],
        "sc_annotation": FILE_ANNOTATIONS_FP
    },
    {
        "name": "FPE_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P82P1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/FP82.cor.csv", "rotation_angle": 0, "flipped": True}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"4210", "sample_name":"PME", "atac_filename": FILE_ATAC_ANNOTATIONS_PFE},
        ],
        "sc_annotation": FILE_ANNOTATIONS_FP
    }
]


for dataset in datasets:
    run_integration(
        dataset["name"], 
        dataset["params_vz_sample_ids"], 
        dataset["params_sc_sample_ids"], 
        dataset["sc_annotation"]
    )


