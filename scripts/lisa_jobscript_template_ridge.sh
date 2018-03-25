#PBS -S /bin/bash
#PBS -lwalltime=18:00:00 -lnodes=1

# call the programs
echo "Job $PBS_JOBID started at `date`" | mail $USER -s "Job $PBS_JOBID"

source activate i2

cd $HOME/retino_HCP/retino_HCP

python ridge.py ---SUBJECT--- ---HEMI---

wait          # wait until programs are finished

echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"