from setuptools import setup, find_packages


setup(
    name='data_integration_package',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'scanpy', 
        'tangram', 
        'scvi-tools', 
        'pandas', 
        'numpy', 
        'matplotlib', 
        'seaborn',
        'torch'
    ],
    description='Package for scRNA, scATAC and spatial data integration.',
    author='Rachid Ounit, Ph.D.',
    author_email='rachid.ounit@gmail.com',
    url='https://github.com/rouni001/Multiomics-Integration',
)

