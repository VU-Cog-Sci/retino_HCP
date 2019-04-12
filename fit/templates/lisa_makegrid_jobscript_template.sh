#!/bin/bash
#SBATCH -t ---job_dur---
#SBATCH -N 1

cd $HOME/projects/retino_HCP/
python ---fun_script--- ---fit_model--- ---save_file--- ---base_dir---