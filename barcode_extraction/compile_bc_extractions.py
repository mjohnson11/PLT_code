import csv
import pandas as pd
from collections import Counter
from glob import glob

def get_deletion_neighborhood(stringer):
    # returns a set of all single character deletions of a string (includes the string)
    return set([stringer] + [stringer[:x] + stringer[x+1:] for x in range(len(stringer))])

def find_bc_match(b, dd):
    # finds the correct barcode by looking for overlap with the deletion network
    # this accounts for all single errors - substitutions, insertions, and deletions
    if b in dd:
        return dd[b]
    else:
        hits = set()
        for b_edit in get_deletion_neighborhood(b):
            if b_edit in dd:
                hits.add(dd[b_edit])
        if len(hits) == 1:
            return hits.pop()
        else:
            return None

def call_bc_match(f):
    try: 
        td = pd.read_csv(f)
        dbc_counter = Counter()
        for entry in np.array(td.loc[td['Type'].isin(['DBC', 'unknown'])][['BC', 'Reads']]):
            check = find_bc_match(entry[0], del_dict)
            if check:
                dbc_counter[check] += entry[1]
        ebc_counter = Counter()
        for entry in np.array(td.loc[td['Type'].isin(['EBC', 'unknown'])][['BC', 'Reads']]):
            check = find_bc_match(entry[0], edel_dict)
            if check:
                ebc_counter[check] += entry[1]
        dbcs, ebcs = sorted([d for d in dbc_counter], key=lambda x: dbc_counter[x]*-1), sorted([e for e in ebc_counter], key=lambda x: ebc_counter[x]*-1)
        return [';'.join(dbcs), ';'.join([str(dbc_counter[d]) for d in dbcs]), ';'.join(ebcs), ';'.join([str(ebc_counter[e]) for e in ebcs])]
    except (FileNotFoundError, IOError):
        print("Wrong file or file path", f)
        return ['', '', '', '']

output_file = '../../output/BC_extraction/bcs_extracted_compiled.csv'
bfa_files = ['../../output/'+b+'/'+b+'_bc_counts_final_w_envs.csv' for b in ['hBFA1', 'hBFA2', 'dBFA2']]
# Reading in known barcodes
bcs_in_bfa = pd.concat([pd.read_csv(bfa_file)[['Diverse_BC', 'Environment_BC']] for bfa_file in bfa_files])

# Generating single-deletiion neighborhoods for each barcode for error correction
del_dict = dict()
edel_dict = dict()
for bc in set(bcs_in_bfa['Diverse_BC']):
    del_dict.update({d: bc for d in get_deletion_neighborhood(bc)})
for bc in set(bcs_in_bfa['Environment_BC']):
    edel_dict.update({d: bc for d in get_deletion_neighborhood(bc)})

bc_extract_files = glob('../../output/BC_extraction/bc_counts/*.csv')

with open(output_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['File', 'dbcs', 'dbc_counts', 'ebcs', 'ebc_counts'])
    c = 0
    for f in bc_extract_files:
        m = call_bc_match(f)
        if m[0] != '':
            c += 1
        writer.writerow([f.split('/')[-1].split('.csv')[0]] + m)
    print('Found a dbc for', c, 'out of', len(bc_extract_files))

    
# Final steps (tranferred over from jupyter notebook)
    
def get_barcode(row):
    dbc_counts = [int(i) for i in str(row['dbc_counts']).split(';')]
    ebc_counts = [int(i) for i in str(row['ebc_counts']).split(';')]
    if len(dbc_counts) > 1:
        dbc_check = (dbc_counts[0] >= 2 and dbc_counts[0]/dbc_counts[1] >= 3)
    else:
        dbc_check = dbc_counts[0] >= 2
    if len(ebc_counts) > 1:
        ebc_check = (ebc_counts[0] >= 2 and ebc_counts[0]/ebc_counts[1] >= 3)
    else:
        ebc_check = ebc_counts[0] >= 2
    if dbc_check and ebc_check:
        return str(row['dbcs']).split(';')[0]+'_'+str(row['ebcs']).split(';')[0]
    else:
        return np.nan

bcd = pd.read_csv(output_file)
bcd = bcd[(bcd['dbc_counts'].notnull()) & (bcd['ebc_counts'].notnull())]
bcd['Barcode'] = bcd.apply(get_barcode, axis=1)
print(len(bcd[bcd['Barcode'].notnull()]), 'successfully called BC')

fitds = []
bc_recs = dict()
for bfa_name in ['dBFA2', 'hBFA1', 'hBFA2']:
    tmp = pd.read_csv('../Final_data_sets/'+bfa_name+'_all_fitness.csv')[['Barcode', 'Home_Environment']].drop_duplicates()
    fitds.append(tmp)
    bc_recs[bfa_name] = list(tmp['Barcode'])
    
def get_bfas(bc):
    return ';'.join([bfa_name for bfa_name in ['dBFA2', 'hBFA1', 'hBFA2'] if bc in bc_recs[bfa_name]])
    
fitd = pd.concat(fitds)
fitd['BFAs'] = fitd['Barcode'].apply(get_bfas)

td = bcd.merge(fitd[['Barcode', 'BFAs', 'Home_Environment']].drop_duplicates(), on='Barcode', how='left')
td = td[td['Barcode'].notnull()]

# Adding file and folder info
folder_names = [
    'Dip_clone_Env',
    'CLM_R2_clones_1N_',
    'CLM_2N_Batch2',
    'GlyEtOH_2N_20170616',
    'GlyEtOH_1N',
    'CLM_2N',
    'FLC_R1_clones_1N',
    'FLC_R2_clones_1N',
    'FLC_2N_Plate1', 
    'FLC_2N_Plate3',
    'FLC_2N_Plate2',
    'GlyEtOH_1N',
    'CLM_clones_1N'
]
# have to sort longest to shortest to get the right one out of the filename
snv = sorted(folder_names, key=lambda x: -1*len(x))

def get_sn(f):
    for sn in snv:
        if f[:len(sn)]==sn:
            return sn
    return 'None'

def get_st(f):
    for sn in snv:
        if f[:len(sn)]==sn:
            return f[len(sn)+1:]
    return 'None'

td['Folder'] = td['File'].apply(get_sn)
td['Strain'] = td['File'].apply(get_st)
td.to_csv('../Final_data_sets/bc_assignments_plus.csv', index=False)