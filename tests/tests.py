import pytest
from data_integration.preprocessing import preprocess_data_sc, preprocess_data_vz
from data_integration.config import (
    FILE_ATAC_ANNOTATIONS_MC1, 
    FILE_ANNOTATIONS_MA
    ROOT_FOLDER_DIR, 
    d_sample_to_dirname_sc
)


def test_preprocess_data_sc():
    # Test the single-cell preprocessing pipeline
    sample = {
        "library_name": "3886",
        "atac_filename": FILE_ATAC_ANNOTATIONS_MC1,
        "sample_name": "MC1"
    }
    filename_sample_metadata = FILE_ANNOTATIONS_MA 
    filename = ROOT_FOLDER_DIR + d_sample_to_dirname_sc[sample_code["library_name"]] + "filtered_feature_bc_matrix.h5"
    adata = preprocess_data_sc(filename, sample, True, filename_sample_metadata)
    assert adata is not None
    assert 'Sample' in adata.obs

