#!/bin/bash
#SBATCH -J PLT_parse_%a  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-60:00              # Runtime in D-HH:MM
#SBATCH -p shared       # Partition to submit to
#SBATCH --mem=60000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../../shell_outs/PLT_parse_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../../shell_outs/PLT_parse_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent


python3 PLT_parse_predemultiplexed.py 

