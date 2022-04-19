import pandas as pd
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('out_base', help='basic output directory')
parser.add_argument('assay_number', help='number that identifies which assay to do')
args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
out_base = args.out_base
assay_number = int(args.assay_number)
with open('../file_info/All_assays.txt') as infile:
    a_list = infile.readlines()
    run_name = a_list[assay_number].strip()


full_out_base = out_base + run_name + '/' + run_name
bc_in = full_out_base + '_bc_counts_clustered_clean.csv'
bc_out = full_out_base + '_bc_counts_final_w_envs.csv'

bcd = pd.read_csv(bc_in)
subpool_cols = [i for i in bcd if 'Subpool' in i]
bcd['Subpool_total'] = np.sum(bcd[subpool_cols], axis=1)


# hBFA2-08M_NaCl_alpha-R2-Subpool
# Function to call the environment based on subpool sequencing, criteria: 
# 95% of the total subpool counts have to be from one environment
# This has to be at least 3 reads
def call_env(row, env_list):
    env_sums = {env: sum([row[c] for c in subpool_cols if env in c.split('-')]) for env in env_list}
    max_counts = max(list(env_sums.values()))
    if row['Subpool_total'] == 0:
        return 'none'
    if max_counts/row['Subpool_total'] > 0.95 and max_counts > 2:
        max_env = [env for env in env_list if env_sums[env] == max_counts]
        assert len(max_env) == 1
        return max_env[0]
    else:
        return 'none'
    
# Simple function to determine if the bc is present in the R1 subpool, the R2 subpool, or both
def call_env_reps(row):
    if row['Subpool_environment'] == 'none':
        return 'none'
    else:
        env = row['Subpool_environment'] 
        env_reps = [c for c in subpool_cols if env in c.split('-') if row[c]/row['Subpool_total'] >= 0.05]
        return ';'.join([e[e.index(env)+len(env):e.index('-Subpool')] for e in env_reps])

envs = set([i.split('-')[1] for i in subpool_cols])
bcd['Subpool_environment'] = bcd.apply(lambda r: call_env(r, envs), axis=1)
bcd['Subpool_replicate'] = bcd.apply(lambda r: call_env_reps(r), axis=1)

def column_sorter(c):
    if 'Time' in c:
        return '-'.join(c.split('-')[:-1])+c.split('-')[-1][4:].zfill(3)
    else:
        return 'z'+c

info_cols = ['Full_BC', 'Diverse_BC', 'Environment_BC', 'Total_Counts']
other_cols = sorted([i for i in bcd if i not in info_cols], key=column_sorter)

bcd[info_cols+other_cols].to_csv(bc_out, index=False)
