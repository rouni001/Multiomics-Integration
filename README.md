# Multiomics-Integration

This Python repo presents key functions facilitating the data integration and processing of scRNA-seq, ATAC-seq and spatial (10x/Visium) datasets using AnnData and Tangram. It allows for efficient mapping of single-cell data to spatial transcriptomics data, enabling downstream analysis and visualization of the spatial distribution of cell types and gene expression.

## Introduction

The Multiomics-Integration repository is a Python package designed to integrate multiple types of omics data, such as single-cell RNA sequencing (scRNA-seq), spatial transcriptomics (e.g., 10x Visium), and single-cell ATAC-seq (scATAC-seq). The integration is performed using the AnnData object format and the Tangram framework, enabling efficient mapping of single-cell data to spatial transcriptomics data. This integration facilitates the analysis and visualization of spatial distributions of gene expression and cell types, contributing significantly to multi-omics research.

Importance of Multi-Omics Integration
In bioinformatics, combining multiple types of omics data (e.g., RNA expression, chromatin accessibility, and spatial distribution) is essential for understanding complex biological systems. This integration allows for a more comprehensive view of cellular behaviors, regulatory mechanisms, and disease processes. Using spatial transcriptomics with single-cell RNA and ATAC-seq data enables researchers to investigate the spatial context of gene expression, how genes are regulated, and how they interact within specific tissue environments. This approach is crucial for advancing personalized medicine, improving drug discovery, and exploring tissue heterogeneity.

## Author & Date
Rachid Ounit, Ph.D

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


