from typing import Optional

import scanpy as sc
import tangram as tg
import os

from .preprocessing import preprocess_data_sc, preprocess_data_vz
from .utils import add_coords
from .config import ROOT_FOLDER_DIR, d_sample_to_dirname_vz, d_sample_to_dirname_sc, d_sample_to_categ_sc
import pandas as pd


def run_integration(ds_name: str, vz_sample, sc_sample_ids, file_annotations: Optional[str], epochs: int, top_genes_count: int):
    """
    Run the integration of spatial visium data and single-cell data, 
    mapping cells to space using Tangram, and save the results.

    Parameters:
    - ds_name (str): Dataset name.
    - vz_sample (dict): Metadata for the Visium (VZ) sample.
    - sc_sample_ids (list of dict): List of metadata for single-cell samples.
    - file_annotations (str, optional): File path for cell type annotations (default: None).

    Returns:
        None
    """
    version = "v5.cell.g" + str(top_genes_count) + ".e" + str(epochs)
    file_sc_adata = f"{version}.h5ad"
    adata_vz_list = []
    file_mapped_sp_sc = f"mapped.sc.vz.{ds_name}.{vz_sample['name']}.{version}.tangram.h5ad"
    file_sc = f"sc.{ds_name}.{file_sc_adata}"
    file_sp = f"VZ.{ds_name}.{vz_sample['name']}.{file_sc_adata}"

    # Creating ANNDATA for VZ data
    filename = ROOT_FOLDER_DIR + d_sample_to_dirname_vz[vz_sample['name']] + "/outs/filtered_feature_bc_matrix"
    
    # Loading VZ dataset
    adata_vz_lcl = preprocess_data_vz(filename, vz_sample['name'], vz_sample['file_morph_mapping'])
    adata_vz_list.append(
        add_coords(
            adata_vz_lcl, d_sample_to_dirname_vz[vz_sample['name']], vz_sample['flipped'], vz_sample['rotation_angle']
        )
    )

    adata_sp = sc.concat(adata_vz_list)

    use_annotations = file_annotations is not None
    cluster_type = "cell_type" if use_annotations else "leiden"
    
    if not os.path.exists(file_mapped_sp_sc):
        #logging.info("Saved mapped spatial vs single-cell data not found... It will be created from scratch.")
        
        if not os.path.exists(file_sc):
            #logging.info("Saved single-cell data not found... It will be created from scratch.")
            
            # Read sc data
            adata_sc_list = []
            for sample_code in sc_sample_ids:
                #logging.info(f"[SC] Adding {sample_code['library_name']}")
                filename = ROOT_FOLDER_DIR + d_sample_to_dirname_sc[sample_code["library_name"]] + "filtered_feature_bc_matrix.h5"
                adata_sc_list.append(preprocess_data_sc(filename, sample_code, use_annotations, file_annotations))

            adata_sc = sc.concat(adata_sc_list)
            adata_sc.obs["Condition"] = adata_sc.obs.Sample.map(d_sample_to_categ_sc)

            sc.pp.filter_genes(adata_sc, min_cells=50)
            #logging.info(f"Size after removing genes in low number of cells: {adata_sc.shape}")
            #logging.info(f"Observation attributes: {adata_sc.obs_names}")
            adata_sc.obs.groupby("Sample").count()
            adata_sc.raw_unnormalized = adata_sc.copy()
            sc.pp.normalize_total(adata_sc)
            sc.pp.log1p(adata_sc)
            adata_sc.raw = adata_sc.copy()
            sc.pp.highly_variable_genes(adata_sc, n_top_genes=top_genes_count)

            # Run MODEL SCVI & leiden
            if not use_annotations:
                cluster_type = "leiden"
                sc.settings.n_jobs = 100
                scvi.model.SCVI.setup_anndata(
                    adata_sc,
                    categorical_covariate_keys=["Sample"],
                    continuous_covariate_keys=["pct_counts_mt"]
                )
                model = scvi.model.SCVI(adata_sc)
                model.train(max_epochs=50)
                model.get_latent_representation().shape
                adata_sc.obsm["X_scVI"] = model.get_latent_representation()
                sc.pp.neighbors(adata_sc, use_rep="X_scVI")

                # Clustering
                #logging.info("Clustering sc data...")
                sc.tl.leiden(adata_sc, resolution=2.5)
            else:
                cluster_type = "cell_type"

            #logging.info("Saving sc-data h5ad")
            adata_sc.write_h5ad(file_sc)
        else:
            #logging.info("Loading the saved single-cell h5ad data ...")
            adata_sc = sc.read_h5ad(file_sc)

        hvg_genes = adata_sc.var.index[adata_sc.var['highly_variable']].tolist()
        tg.pp_adatas(adata_sc, adata_sp, genes=hvg_genes)

        #logging.info(f"#training genes in sc: {adata_sc.uns['training_genes'][:100]} ...")
        #logging.info(f"#training genes in sp: {adata_sp.uns['training_genes'][:100]} ...")

        adata_sc.write_h5ad(file_sc)
        adata_sp.write_h5ad(file_sp)

        #logging.info("Running mapping...")
        ad_map = tg.map_cells_to_space(
            adata_sc=adata_sc,
            adata_sp=adata_sp,
            device='cpu',
            num_epochs=epochs,
        )
        
        ad_map.write_h5ad(file_mapped_sp_sc)
        adata_sc.write_h5ad(file_sc)
        adata_sp.write_h5ad(file_sp)
        #logging.info("Saved ad_map, ad_sc, and ad_sp.")
    else:
        #logging.info("Found saved ad_map, ad_sc, and ad_sp. Loading now...")
        adata_sc = sc.read_h5ad(file_sc)
        ad_map = sc.read_h5ad(file_mapped_sp_sc)
        adata_sp = sc.read_h5ad(file_sp)

    try:
        #logging.info("Running spatial cell projections...")
        tg.project_cell_annotations(ad_map, adata_sp, annotation=cluster_type)
        annotation_list = sorted(list(pd.unique(adata_sc.obs[cluster_type])))
        tg.plot_cell_annotation_sc(
            adata_sp,
            annotation_list,
            x='x',
            y='y',
            spot_size=200,
            scale_factor=0.1,
            perc=0.001,
            filename=f".tg.{ds_name}.{version}.{cluster_type}.scvi.pdf"
        )

        ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
        logging.info(f"Gene expression on overlap data: {ad_ge}")

        # Printing figures for some regions/cell types
        genes_astrocyte = ['gfap', 'slc1a2', 'aqp4']
        genes_OPC = ['pdgfra', 'itpr2', 'tcf7l2']
        genes_oligodendrocyte = ['olig2', 'bmp4']

        fig = tg.plot_genes_sc(
            genes_astrocyte, 
            adata_measured=adata_sp, 
            adata_predicted=ad_ge, 
            spot_size=200, 
            scale_factor=0.1, 
            perc=0.001, 
            return_figure=True
        )
        fig.savefig(f".{ds_name}.{version}_genes.astr.pdf")
        fig = tg.plot_genes_sc(
            genes_OPC, 
            adata_measured=adata_sp, 
            adata_predicted=ad_ge, 
            spot_size=200, 
            scale_factor=0.1, 
            perc=0.001, 
            return_figure=True
        )
        fig.savefig(f".{ds_name}.{version}_genes.opc.pdf")
        fig = tg.plot_genes_sc(
            genes_oligodendrocyte, 
            adata_measured=adata_sp, 
            adata_predicted=ad_ge, 
            spot_size=200, 
            scale_factor=0.1, 
            perc=0.001, 
            return_figure=True)
        fig.savefig(f".{ds_name}.{version}_genes.olig.pdf")

        df_all_genes = tg.compare_spatial_geneexp(
            ad_ge, adata_sp, adata_sc
        )
        #logging.info(f"Similarity scores on all genes: {df_all_genes}")

        #logging.info("Plotting AUC graph...")
        tg.plot_auc(df_all_genes, f".{ds_name}.{version}.auc.pdf")

    except Exception as e:
        print(f"{e}")
        #logging.error(f"Error occurred: {e}")
    finally:
        print("continue.")
        #logging.info(f"Process completed for {ds_name}.")



