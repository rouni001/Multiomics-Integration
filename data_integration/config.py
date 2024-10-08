from .filenames import (
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
    "P83P1":"p83_p1_run2",
    "P95P1":"p95_p1_run2",
    "P100P1":"p100_p1_run2",
    "P112P1":"p112_p1_run2",
    "P113P1":"p113_p1_run2",
    "P81P2":"p81_p2_run2",
    "P81P1":"p81_p1_run2",
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
    "P95P1":"CM",
    "P81P2":"EcigM",
    "P81P1":"EcigM",
    "P83P1":"EcigM",
    "P112P1":"EcigM",
    "P113P1":"EcigM",
}

# Define datasets
datasets =[
    {
        "name": "MAC_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P100A1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/MA100.cor.csv",  "rotation_angle": 180, "flipped": True}, #P94A1
        "params_sc_sample_ids": [
            {"library_name":"3886", "sample_name":"MC1", "atac_filename": FILE_ATAC_ANNOTATIONS_MC1},
            {"library_name":"3887", "sample_name":"MC2", "atac_filename": FILE_ATAC_ANNOTATIONS_MC2}
        ],
        "sc_annotation": FILE_ANNOTATIONS_MA,
    },
    {
        "name":"MAE_SC_VZ_ATAC",
        "params_vz_sample_ids": {"name":"P112A1", "file_morph_mapping": FILE_ANNOTATIONS_VZ_R + "/MA112.cor.csv",  "rotation_angle": 180, "flipped": True}, #P80A1
        "params_sc_sample_ids": [
            {"library_name":"3888", "sample_name":"ME1", "atac_filename": FILE_ATAC_ANNOTATIONS_ME1},
            {"library_name":"3889", "sample_name":"ME2", "atac_filename": FILE_ATAC_ANNOTATIONS_ME2}
        ],
        "sc_annotation": FILE_ANNOTATIONS_MA,
    }
]

