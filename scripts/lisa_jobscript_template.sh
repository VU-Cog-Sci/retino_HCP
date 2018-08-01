#!/bin/bash
#PBS -lwalltime=---job_dur---
#PBS -lnodes=1

source activate i36
python -m ---fit_file--- ---subject--- ---start_idx--- ---end_idx--- ---data_file--- ---base_dir---