import pandas as pd
from glob import glob
import numpy as np
flist = glob('../../output/LT_*/*_bc_counts_clustered_clean.csv')
for f in flist:
    print(f)
    fs = f.split('/')[-1].split('_bc_counts_clustered_clean')[0]
    td = pd.read_csv(f).rename(columns={'Full_BC': 'Barcode'})
    time_cols = sorted([i for i in td if 'Time' in i], key=lambda x: '-'.join(x.split('-')[:-1])+x.split('Time')[-1].zfill(3))
    td[['Barcode']+time_cols].to_csv('../../LT_reads_final/'+fs+'_lineage_tracking_reads.csv', index=False)
    suf = '_clipped_log10_freq'
    for tp in time_cols:
        tmp_sum = np.sum(td[tp])
        td[tp+suf] = np.log10(np.clip(td[tp]/tmp_sum, 10**(-6), 1))
    td[['Barcode']+[i+suf for i in time_cols]].to_csv('../../LT_reads_final/'+fs+'_lineage_tracking_clipped_freqs.csv', index=False)

