import re
import os
import glob
import json
import sys

#############################################################################################
#############################################################################################
##### FMRIPREP
#############################################################################################
#############################################################################################

subjects = [1,2]

batch_string = """#!/bin/bash
#SBATCH -p normal -t 20:00:00 -N 1
# job requires at most 100 hours, 0 minutes
#     and 0 seconds wallclock time

# call the programs
echo "Job $SLURM_JOB_ID started at `date`" | mail $USER -s "Job $SLURM_JOB_ID"

PYTHONPATH="" singularity run -B /nfs/k_lab,/scratch \
/nfs/k_lab/software/fmriprep-tfpatch.simg \
/nfs/k_lab/data/cerebellum_prf/sourcedata/ /nfs/k_lab/data/cerebellum_prf/derivatives/ participant \
--participant-label sub-$SJ_NR --output-space T1w template fsaverage6 \
--nthreads 15 --omp-nthreads 15 --low-mem --fs-license-file /home/knapen/bin/freesurfer/license \
--ignore slicetiming -w /scratch --skip_bids_validation

wait          # wait until programs are finished

echo "Job $SLURM_JOB_ID finished at `date`" | mail $USER -s "Job $SLURM_JOB_ID"
"""

basedir = '/home/knapen/batch/'

os.chdir(basedir)

for subject in subjects:

    working_string = batch_string.replace('$SJ_NR', str(subject).zfill(2))

    js_name = os.path.join(basedir, str(subject).zfill(2) + '_cerebellar_retinotopy.sh')
    of = open(js_name, 'w')
    of.write(working_string)
    of.close()

    print('submitting ' + js_name + ' to queue')
    print(working_string)
    os.system('sbatch ' + js_name)