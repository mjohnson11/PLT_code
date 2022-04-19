#!/bin/bash
#SBATCH -J PLT_cleanup_%a  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-00:05              # Runtime in D-HH:MM
#SBATCH -p shared       # Partition to submit to
#SBATCH --mem=5000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outs/PLT_cleanup_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outs/PLT_cleanup_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

source activate milo_py37


python3 PLT_remove_chimeras.py ../../output/ "${SLURM_ARRAY_TASK_ID}"
