# IFNBoost

This repository contains the code accompanying the paper [**"IFNBoost: An interpretable computational model for identifying IFNγ-inducing peptides"**](https://www.biorxiv.org/content/10.1101/2025.01.21.634172v1), currently available on *bioRxiv (2025)*.

## Installation

1. Install [Anaconda](https://www.anaconda.com/products/distribution)
2. Create a virtual environment (example name: IFNBoost) with Python 3.11.8 using "conda create --name IFNBoost python=3.11.8"
3. Within the environment install the dependencies provided in file "requirements.txt" using "pip install -r requirements.txt"
4. The main code is in *MODEL.py* file and it uses the *allhost.xlsx* file as input and processes it. Make sure the input file is in the same folder, else specify the path to file in the script. The *functions.py* also needs to be in the same folder as *MODEL.py*.
5. The figures showing the performance of IFNBoost using 1000 bootsraps are generated in *Figures.py*.
6. Comparisons with other methods are in *Validation.py*, which utlises the 2024 data from IEDB contained in *Tcell2024.xlsx* 

## Web application

An interactive version of IFNBoost is available on streamlit: [IFNBoost Webserver](ifnboost.streamlit.app)

## Citation

If you use this code or model, please cite:

I. Azad et al., "IFNBoost: An interpretable computational model for identifying IFNγ inducing peptides," bioRxiv, 2025.\
[https://doi.org/10.1101/2025.01.21.634172]

## Troubleshooting

For any issues or feedback, please contact inam.ulhaqazad@student.unimelb.edu.au
