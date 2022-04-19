#!/bin/bash
#SBATCH -J PLT_cluster_%a  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-60:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=100000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outs/PLT_cluster_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outs/PLT_cluster_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate milo_py37

python3 PLT_cluster.py ../../output/ "${SLURM_ARRAY_TASK_ID}" 

python3 PLT_remove_chimeras.py ../../output/ "${SLURM_ARRAY_TASK_ID}"
