#!/bin/bash

# This script reproduces Python (3.6, 3.8) and R (4.2.1) environments
# used in our analyses. This script was developed to run on a Compute Canada cluster.
# It should trivially be portable to any multi-user compute cluster having a system-wide
# installation of Python 3.6, Python 3.8, and R 4.2.1.  

# python 3.6 environment
module load python/3.6
python -m venv p36
source p36/bin/activate
pip install -r p36_requirements.txt

# python 3.8 environment
module load python/3.8
python -m venv p38
source p38/bin/activate
pip install --upgrade pip
pip install -r p38_requirements.txt

# R 4.2.1 environment
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
module load r/4.2.1
Rscript -e 'install.packages("renv", repos="https://cran.r-project.org")'
Rscript ./init_R.R


echo "Script finished"
