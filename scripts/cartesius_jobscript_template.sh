#!/bin/bash
#SBATCH -t 80:00:00
#SBATCH -N 1
#SBATCH -n 23
#SBATCH -p normal

source activate i2

cd $HOME/retino_HCP/retino_HCP

python prf_fit.py ---SUBJECT--- ---HEMI---
