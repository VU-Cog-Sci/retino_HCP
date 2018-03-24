#PBS -S /bin/bash
#PBS -lwalltime=36:00:00 -lnodes=1:mem64gb
# job requires at most 36 hours, 0 minutes
#     and 0 seconds wallclock time

# call the programs
echo "Job $PBS_JOBID started at `date`" | mail $USER -s "Job $PBS_JOBID"

source activate i2

cd $HOME/retino_HCP/retino_HCP

python prf_fit.py ---SUBJECT---

wait          # wait until programs are finished

echo "Job $PBS_JOBID finished at `date`" | mail $USER -s "Job $PBS_JOBID"