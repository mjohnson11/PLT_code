import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

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
bc_in = full_out_base + '_bc_counts_clustered.csv'
bc_out = full_out_base + '_bc_counts_clustered_clean.csv'
chimera_output = full_out_base + '_chimeric.csv'

dat = pd.read_csv(bc_in)
dat['Full_BC'] = dat['Diverse_BC'] + '_' + dat['Environment_BC']

bc_cols = ['Full_BC', 'Diverse_BC', 'Environment_BC']

dat = dat.sort_values('Total_Counts', ascending=False)
# excluding chimeras
dd = dict(dat['Diverse_BC'].value_counts())
dbc_with_mult_hits = set([d for d in dat['Diverse_BC'] if dd[d] > 1])
fullbc_chimeras = []
dbc_top_counts = dict()
for entry in dat.as_matrix(['Full_BC', 'Diverse_BC', 'Total_Counts']):
    if entry[1] in dbc_with_mult_hits:
        if entry[1] in dbc_top_counts:
            if (entry[2] / dbc_top_counts[entry[1]]) < 0.01:
                fullbc_chimeras.append(entry[0])
        else:
            dbc_top_counts[entry[1]] = entry[2]
dat['is_chimeric'] = dat.apply(lambda drow: drow['Full_BC'] in fullbc_chimeras, axis=1)
d_no_chimeras = dat.loc[~dat['is_chimeric']]
d_yes_chimeras = dat.loc[dat['is_chimeric']]

d_yes_chimeras[[c for c in dat if c!='is_chimeric']].to_csv(chimera_output, index=False)

print(len(d_no_chimeras), 'bcs after combining and excluding chimeras.')

d_no_chimeras[[c for c in dat if c!='is_chimeric']].to_csv(bc_out, index=False)
