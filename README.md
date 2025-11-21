\# ER\_FA Masks Analysis



This repository contains the jupyter notebook and environment setup used for the ER\_FA project. 



Output from the jupyter notebook is then quantitatively analysed with the R scripts contained in this repo.



\## Repository Contents



\- `FA-ER masks analysis 1.ipynb` — read in binary masks of ER and FA and produce csv files quantifying properties of FAs and the ER overlap/colocalisation

\- `environment.yml` — conda/mamba environment specification for `FA-ER masks analysis 1.ipynb`

\- `ERFA masks only analysis.R` — produce plots of the csv output from `FA-ER masks analysis 1.ipynb`

\- `README.md`

\- `.gitignore`

\- `LICENCE`


\### Setting up the environment



The project uses Python 3.9 and `devbio-napari`. To recreate the environment, you can either use the exported environment file (recommended) or create it manually. Both options are shown below:



```bash

# Recommended: create from the exported environment file

mamba env create -f environment.yml

# or with conda

conda env create -f environment.yml



# Alternative: manually create the environment

mamba create --name ERFA-devbio-napari-env python=3.9 devbio-napari -c conda-forge



# Activate the environment

conda activate ERFA-devbio-napari-env



# Open Jupyter Lab

jupyter lab

```

\### Storing ER and FA masks for reading into the script'



Binary masks should be generated (e.g. in FIJI) and stored as TIF files. Each file should retain the name of the original image followed by "\_threshold-ER" or "\_threshold-vinc" for ER and FA respectively. All masks are then stored together in the working directory, which is read into the jupyter notebook at the folder\_dir line (replace the path with the path to your working directory).

\### Analysing output from Jupyter notebook in R

The Jupyter notebook produces csv files which are loaded into R using the R script. Packages are listed at the start of the script. Load csv files in and use the plotting functions to produce different outputs.

