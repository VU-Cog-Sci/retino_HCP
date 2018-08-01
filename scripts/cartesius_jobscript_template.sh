#!/bin/bash
#SBATCH -t ---job_dur---
#SBATCH -N 1
#SBATCH -n 23
#SBATCH -p normal

source activate i36
python ---fit_file--- ---subject--- ---start_idx--- ---end_idx--- ---data_file--- ---base_dir---