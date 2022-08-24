import pandas as pd
import numpy as np
import csv
from collections import defaultdict, Counter
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('env_file', help='csv with environment bc info')
parser.add_argument('input_base', help='directory with comma separated bc count files')
parser.add_argument('ebc_exclude_include_file', help='file with environment barcode exclusion / inclusion info')
parser.add_argument('tp_exclude_file', help='file with timepoint exclusion info')
parser.add_argument('output_base', help='directory for output png files')
parser.add_argument('-count_thresh', type=int, default=100000, help='count threshold to exclude timepoints by')
args = parser.parse_args()

envs = [
    '02M_NaCl_2N-R1',
    '02M_NaCl_alpha-R1',
    '08M_NaCl_2N-R1',
    '08M_NaCl_alpha-R1',
    'CLM_2N-R1',
    'CLM_alpha-R1',
    'FLC4_2N-R1',
    'FLC4_alpha-R1',
    'FLC32_2N-R1',
    'FLC32_alpha-R1',
    'M3_2N-R1',
    'M3_alpha-R1',
    '21C_2N-R1',
    '21C_alpha-R1',
    'GlyEtOH_2N-R1',
    'GlyEtOH_alpha-R1',
    '37C_2N-R1',
    '37C_2N-R2',
    'YPD_2N-R1',
    'YPD_2N-R2',
    'SC_2N-R1',
    'SC_2N-R2',
    '48Hr_2N-R1',
    '48Hr_2N-R2',
    'pH7_3_2N-R1',
    'pH7_3_2N-R2',
    'pH3_8_2N-R1',
    'pH3_8_2N-R2'
]

env_file = args.env_file
input_base = args.input_base
output_base = args.output_base
count_thresh = int(args.count_thresh)
ebc_exclude_include_file = args.ebc_exclude_include_file
tp_exclude_file = args.tp_exclude_file
ed = {i[0]: i[1] for i in pd.read_csv(env_file).as_matrix(['Environment.BC', 'Putative.Environment'])}
tp_exclude = defaultdict(list)
for i in pd.read_csv(tp_exclude_file).as_matrix(['Evolution', 'TP']):
    tp_exclude[i[0]].append(i[1])
ei_dict = defaultdict(lambda: defaultdict(list))
with open(ebc_exclude_include_file, 'r') as infile:
    reader = csv.reader(infile)
    rc = 0
    for row in reader:
        rc += 1
        if rc > 1:
            for evo in row[0].split(';'):
                ei_dict[evo][row[1]].append(row[2:])

print('Excluded timepoints', tp_exclude)

for env_name in envs:
    short_env_name = env_name.split('-')[0]
    if 'FLC32' in short_env_name or 'M3' in short_env_name:
        short_env_name = short_env_name.replace('FLC32', 'YPD').replace('M3', 'SC')
    d = pd.read_csv(input_base + env_name + '_lineage_tracking_reads.csv')
    all_tps = [i for i in d.columns if 'Time' in i and sum(d[i]) > 100000 and
               int(i[i.index('Time')+4:]) not in tp_exclude[env_name]]
    print(env_name, [int(i[i.index('Time')+4:]) for i in all_tps])
    # making separate barcode columns - backwards compatible with this old bc exclusion stuff
    d['Environment.BC'] = d['Barcode'].apply(lambda b: b.split('_')[1])
    d['Diverse.BC'] = d['Barcode'].apply(lambda b: b.split('_')[0])
    d['Full.BC'] = d['Environment.BC']+d['Diverse.BC']
    all_ebcs = set(d['Environment.BC'])
    use_ebcs = [i for i in all_ebcs if short_env_name in ed.setdefault(i, 'unknown')]
    excluded_full_bcs = []
    for ebc in ei_dict[env_name]:
        for tmp in ei_dict[env_name][ebc]:
            if tmp[0] == 'all' and tmp[1] == 'include':
                use_ebcs.append(ebc)
            elif tmp[1] == 'exclude':
                for dbc in tmp[0].split(';'):
                    excluded_full_bcs.append(ebc+dbc)

    d['Putative.Environment'] = d.apply(lambda row: ed.setdefault(row['Environment.BC'], 'unknown'), axis=1)
    d_use = d.loc[d['Environment.BC'].isin(use_ebcs)].loc[~d['Full.BC'].isin(excluded_full_bcs)]

    d_use[['Barcode']+all_tps].to_csv(output_base + env_name + '_clean.csv', index=False)
