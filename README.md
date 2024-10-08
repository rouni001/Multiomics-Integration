# Data Integration Package

This Python package facilitates the integration and processing of scRNA-seq, ATAC-seq and spatial (10x/Visium) datasets using AnnData and Tangram. It allows for efficient mapping of single-cell data to spatial transcriptomics data, enabling downstream analysis and visualization of the spatial distribution of cell types and gene expression.


## Author & Date
Dr. Rachid Ounit.

September 2024.

## Version
v0.1.0

## **Features**

- **Preprocessing Pipelines**: Handles both single-cell and VZ (spatial transcriptomics) data. It preprocesses data by adding ATAC-seq information, filtering cells, detecting doublets, and performing quality control.
  
- **Tangram Mapping**: Efficiently maps single-cell data onto spatial transcriptomics data using the Tangram framework, allowing spatial projection of gene expression and cell-type annotations.

- **Modular Configuration**: Easily add or modify datasets by updating the configuration file (`config.py`). This allows flexible dataset handling and integration.

- **Command-line Interface**: Run dataset integration with simple command-line arguments to process different datasets without modifying the code.

- **Testing Suite**: Includes unit tests to validate core functions, ensuring reliable and consistent performance of the package.

---

## **Table of Contents**

- [Features](#features)
- [Installation](#installation)
- [Building Dependencies](#building-dependencies)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Configuration](#configuration)
- [Running Tests](#running-tests)
- [Contributing](#contributing)
- [License](#license)

---

## **Installation**

To install the package and set up the necessary environment, follow the steps below:

1. **Clone the repository**:

```bash
git clone https://github.com/your-repo/data-integration-package.git
cd data_integration_package

2. **Install the required dependencies**:

The package relies on several dependencies, including scanpy, tangram, scvi-tools, and others. To install all dependencies, run:

```pip install -r requirements.txt```

## **Building Dependencies**
scanpy
tangram
scvi-tools
pandas
seaborn
matplotlib
numpy
torch
pytest

To install these dependencies, simply run:
``pip install -r requirements.txt``


## **Configuration**
The package allows you to add and configure datasets easily. The config.py file contains all relevant information regarding file paths, sample mappings, and datasets.
You can extend this configuration by adding new datasets, updating file paths, or modifying the parameters used for mapping and preprocessing.

## **Running Tests**
The package includes unit tests to ensure that core functionality works as expected. The tests are written using pytest, and you can run them by executing:

pytest


## **Contributing**

We welcome contributions to this project! To contribute:

- Fork the repository.

- Create a new branch with your feature or bugfix: git checkout -b feature-name.

- Commit your changes: git commit -am 'Add new feature'.

- Push to the branch: git push origin feature-name.

- Open a Pull Request.

- Make sure to run the test suite (pytest) before submitting your changes to ensure everything works as expected.


