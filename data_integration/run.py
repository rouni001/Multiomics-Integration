from .integration import run_integration
from .config import datasets


def run_dataset_integration(dataset_name):
    """
    Run the integration for a specific dataset by its name.

    Parameters:
    - dataset_name (str): Name of the dataset to integrate.

    Returns:
    - None
    """
    # Find the dataset in the config file
    dataset = next((ds for ds in datasets if ds['name'] == dataset_name), None)

    if dataset:
        print(f"Running integration for dataset: {dataset['name']}")
        run_integration(
            dataset["name"],
            dataset["params_vz_sample_ids"],
            dataset["params_sc_sample_ids"],
            dataset["sc_annotation"],
            600,
            5000
        )
    else:
        print(f"Dataset {dataset_name} not found in configuration.")

