\# ER\_FA Masks Analysis



This repository contains the jupyter notebook and environment setup used for the ER\_FA project. 



Output from the jupyter notebook is then quantitatively analysed with the R scripts contained in this repo.



\## Repository Contents



\- `FA-ER masks analysis 1.ipynb` — mask analysis script 

\- `environment.yml` — conda/mamba environment specification

\- Other files as needed



\### Setting up the environment



The project uses Python 3.9 and `devbio-napari`. To recreate the environment, you can either use the exported environment file (recommended) or create it manually. Both options are shown below:



```bash

\# Recommended: create from the exported environment file

mamba env create -f environment.yml

\# or with conda

conda env create -f environment.yml



\# Alternative: manually create the environment

mamba create --name ERFA-devbio-napari-env python=3.9 devbio-napari -c conda-forge



\# Activate the environment

conda activate ERFA-devbio-napari-env



\# Open Jupyter Lab

jupyter lab



\### Storing ER and FA masks for reading into the script'



Binary masks should be generated (e.g. in FIJI) and stored as TIF files. Each file should retain the name of the original image followed by "\_threshold-ER" or "\_threshold-vinc" for ER and FA respectively. All masks are then stored together in the working directory, which is read into the jupyter notebook at the folder\_dir line (replace the path with the path to your working directory).

