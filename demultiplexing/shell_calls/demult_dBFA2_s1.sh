#!/bin/bash
#SBATCH -J demult_dBFA2_s1  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-60:00              # Runtime in D-HH:MM
#SBATCH -p shared       # Partition to submit to
#SBATCH --mem=60000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../../shell_outs/demult_dBFA2_s1.out      # File to which STDOUT will be written
#SBATCH -e ../../../shell_outs/demult_dBFA2_s1.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

python3 ../PLT_demult.py ../../file_info/demult_maps/BFAs/dBFA2_Stanford_1_indices.csv ../../../all_reads_demultiplexed/ /n/holyscratch01/desai_lab/mjohnson/PLT/raw_reads/BFAs/BFA_raw_reads/dBFA2_Stanford1/dBFA2_Stanford1_R1.fastq.gz dBFA2_s1
