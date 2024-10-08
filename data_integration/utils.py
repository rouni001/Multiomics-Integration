import pandas as pd
import numpy as np
import scanpy as sc

from .config import (    
    ROOT_FOLDER_DIR,
    FILE_ANNOTATIONS_MP, 
    FILE_ANNOTATIONS_FP, 
    FILE_ANNOTATIONS_FA, 
    FILE_ANNOTATIONS_MA
)

def add_atac_data_to_sc(adata, filename_annotation):
    """
    Add ATAC-seq data to an AnnData object containing the scRNA data.
    
    Parameters:
    - adata (AnnData): AnnData object containing the scRNA data
    - filename_annotation (str): Path to ATAC metadata file.
    
    Returns:
    - AnnData: Updated AnnData object.
    """
    
    adata_atac = sc.read_h5ad(filename_annotation)
    common_index = adata.obs.index.intersection(adata_atac.obs.index)
    adata_subset = adata[common_index, :].copy()
    adata_atac_subset = adata_atac[common_index, :].copy()
    merged_obs = pd.concat([adata_subset.obs, adata_atac_subset.obs], axis=1)
    adata_merged = sc.AnnData(X=adata_subset.X, obs=merged_obs, var=adata_subset.var)
    adata_merged.obs = adata_merged.obs.drop(columns=["total_counts", "n_counts", "leiden"])
    return adata_merged


def add_cell_types_to_sc(adata, idx_name, sample_name, filename_annotation):
    """
    Adds cell types to an AnnData object by merging metadata annotations from a CSV file.

    Parameters:
    - adata (AnnData): AnnData object to be annotated.
    - idx_name (str): Name of the index column in the AnnData object.
    - sample_name (str): Sample identifier to filter the annotation file.
    - filename_annotation (str): Path to the CSV file with annotations.

    Returns:
    - Updated AnnData object with cell types added.
    """
    posterior = filename_annotation == "FILE_ANNOTATIONS_MP" or filename_annotation == "FILE_ANNOTATIONS_FP"
    
    # Load CSV file and filter by sample name
    annotations = pd.read_csv(filename_annotation)
    annotations = annotations[annotations['Sample'] == sample_name]

    # Determine prefix size based on the annotation type
    prefix_index_size = 5 if filename_annotation == "FILE_ANNOTATIONS_MP" else 4
    column_index_name = "Unnamed: 0"  # Default column name
    
    # Adjust for other file types
    if filename_annotation == "FILE_ANNOTATIONS_MA" or filename_annotation == "FILE_ANNOTATIONS_FA":
        column_index_name = "X"
    
    # Extract and clean CellID
    annotations["CellID"] = annotations[column_index_name].str[prefix_index_size:]
    annotations["CellID"] = annotations["CellID"].str.strip()
    annotations["CellID"] = annotations["CellID"].astype(str)
    
    # Rename 'cell.type' to 'cell_type' for consistency
    annotations.rename(columns={'cell.type': 'cell_type'}, inplace=True)
    
    # Drop unnecessary columns
    annotations.drop(columns=[column_index_name, "Sample", "Group"], inplace=True)
    
    # Set 'CellID' as the index and join with AnnData object
    annotations.set_index("CellID", inplace=True)
    adata.obs.index = adata.obs.index.str.strip()
    adata.obs = adata.obs.join(annotations, how="left")
    
    # Check the number of successfully annotated cells
    annotated_cells = adata.obs['cell_type'].notna().sum()
    
    # If not posterior, replace certain cell type labels
    if not posterior:
        adata.obs['cell_type'] = adata.obs['cell_type'].replace('exLayer4', 'exLayer5/6-tmp')
        adata.obs['cell_type'] = adata.obs['cell_type'].replace('exLayer5/6', 'exLayer4')
        adata.obs['cell_type'] = adata.obs['cell_type'].replace('exLayer5/6-tmp', 'exLayer5/6')
    
    # Filter out rows where 'cell_type' is NaN
    adata = adata[adata.obs['cell_type'].notna()].copy()
    
    return adata



def add_cell_types_to_sc(adata, idx_name, sample_name, filename_annotation):
    """
    Add cell types to scRNA AnnData object from a CSV file.

    Parameters:
    - adata (AnnData): AnnData object to be annotated.
    - idx_name (str): Name of the index column in the AnnData object.
    - sample_name (str): Sample identifier to filter the annotation file.
    - filename_annotation (str): Path to the CSV file with annotations.

    Returns:
    - AnnData: Updated AnnData object with cell types added.
    """
    posterior = filename_annotation == FILE_ANNOTATIONS_MP or filename_annotation == FILE_ANNOTATIONS_FP
    print(f"Loading metadata: {filename_annotation}")

    # Load CSV with proper header handling
    annotations = pd.read_csv(filename_annotation)
    #annotations = annotations.tail(annotations.shape[0] -1)
    # Filter the annotations by sample_name
    annotations = annotations[annotations['Sample'] == sample_name]
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


def filter_doublets(adata):
    """
    Filter doublets from an AnnData object.
    
    Parameters:
    - adata (AnnData): AnnData object to filter.
    
    Returns:
    - AnnData: Filtered AnnData object.
    """
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs.pct_counts_mt < np.quantile(adata.obs.pct_counts_mt, 0.98))]
    return adata


def add_regions(intg_adata_slice, filename):
    """
    Add metadata about the morphological regions from a CSV file to an AnnData object
    representation the spatial visium data of a slice tissue.

    Parameters:
    - intg_adata_slice (AnnData): AnnData object to which region metadata will be added.
    - filename (str): Path to the CSV file containing region information.

    Returns:
    - AnnData: Updated AnnData object with region information added.
    """
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
    intg_adata_res = intg_adata_slice[intg_adata_slice.obs_names.isin(result_df_col6.index)].copy()
    # Add the metadata columns to the AnnData object's obs

    for col in result_df_col6.columns:
        intg_adata_res.obs[col] = result_df_col6.loc[intg_adata_res.obs_names, col]

    return intg_adata_res


def add_coords(adata, sample_dirname, reverse_coords, angle = 0):
    filename_positions = ROOT_FOLDER_DIR + f"SpaceRangerTissuePositions/{sample_dirname}tissue_positions.csv"

    positions = pd.read_csv(filename_positions, header=None)
    positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    positions = positions[positions["in_tissue"] == 1]  # Filter to include only positions that are in tissue

    # Set index by barcode and join with adata
    positions.set_index("barcode", inplace=True)

    # Join positions with AnnData obs
    adata.obs = adata.obs.join(positions, how="left")

    # Extract 'x' and 'y' coordinates
    adata.obs['x'] = adata.obs['pxl_row_in_fullres']
    adata.obs['y'] = adata.obs['pxl_col_in_fullres']

    if reverse_coords:
        adata.obs['x'] = -adata.obs['x']
        adata.obs['y'] = -adata.obs['y']

    if angle > 0:
        adata.obs['x'], adata.obs['y'] = rotate_coordinates(adata.obs['x'], adata.obs['y'], angle)

    adata.obs.drop(columns=['pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)

    return adata


def rotate_coordinates(x, y, angle_degrees):
    """
    Rotate coordinates by a given angle.

    Parameters:
    - x (Series): The x coordinates to rotate.
    - y (Series): The y coordinates to rotate.
    - angle_degrees (float): The angle in degrees to rotate the coordinates.

    Returns:
    - tuple: Rotated x and y coordinates.
    """
    angle_radians = np.radians(angle_degrees)
    rotation_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians)],
        [np.sin(angle_radians), np.cos(angle_radians)]
    ])

    coords = np.vstack([x, y])
    rotated_coords = rotation_matrix @ coords
    return rotated_coords[0, :], rotated_coords[1, :]


def generate_symmetric_coordinates(x, y, is_x_axis=True):
    """
    Generate symmetric coordinates around the midpoint.

    Parameters:
    - x (array-like): X coordinates.
    - y (array-like): Y coordinates.
    - is_x_axis (bool): Whether to generate symmetry around the X axis. If False, it generates around the Y axis.

    Returns:
    - tuple: Symmetric X and Y coordinates.
    """
    x = np.array(x)
    y = np.array(y)

    if is_x_axis:
        midpoint_x = np.mean(x)
        symmetric_x = 2 * midpoint_x - x
        symmetric_y = y
    else:
        midpoint_y = np.mean(y)
        symmetric_y = 2 * midpoint_y - y
        symmetric_x = x

    return symmetric_x.tolist(), symmetric_y.tolist()

