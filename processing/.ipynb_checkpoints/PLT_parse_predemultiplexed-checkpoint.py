"""
A program to demultiplex reads and count barcodes for the PLT
Milo Johnson
Started 3_6_16
made simple 2_6_17
"""

import time
import gzip
import pandas as pd
import re
import csv
import numpy as np
import os
import subprocess
from PLT_htmlReport import make_html_report
from PLT_tools import FourLineFastq, reverse_transcribe


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('assay_number', help='number that identifies which assay to do')
parser.add_argument('out_base', help='basic output directory')
parser.add_argument("-quality_cutoff", type=int, default=30, help='quality threshold for bc region')
args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
assay_number = args.assay_number
out_base = args.out_base
QUALITY_CUTOFF = args.quality_cutoff
with open('../file_info/All_assays.txt') as infile:
    a_list = infile.readlines()
    run_name = a_list[assay_number].strip()
    
assay_file = pd.read_csv('../file_info/assay_files/'+run_name+'_assay.csv')

output_dir = out_base + run_name + '/'
output_file = output_dir + run_name + '_bc_counts.csv'
stats_out_base = output_dir + 'run_statistics/'
lib_stats_out = stats_out_base + run_name + '_library_statistics.csv'
umi_fam_size_out = stats_out_base + run_name + '_umi_family_sizes.csv'
html_report_out = stats_out_base + run_name + '_html_report'

all_primer_info = {row['Filename']: row for j, row in pd.read_csv('../file_info/All_file_primer_info.csv').iterrows()}
read_file_base = '../../all_reads_demultiplexed/'

if not os.path.isdir(output_dir):
    print('Making main output directory:', output_dir)
    subprocess.call(['mkdir', output_dir])
if not os.path.isdir(stats_out_base):
    print('Making stats output directory:', stats_out_base)
    subprocess.call(['mkdir', stats_out_base])

MY_REGEX = re.compile('\D*?(GTACC|GGACC|GGTCC|G.TACC|GG.ACC|GGT.CC|GGTA.C|GGTAC.)(\D{24,28})(.TAACT|A.AACT|AT.ACT|ATA.CT|ATAA.T|ATAAC|AAACT|ATACT|ATAAT)\D*')


class BcCounter:

    """
    BcCounter counts barcodes and corrects counts based on the unique molecular indices (UMIs) they are paired with.
    """

    def __init__(self, assay_file):
        # this bc_dict has entries like: [bc div, bc env, total counts, [counts_in_lib_1, counts_in_lib_2, etc.]]
        self.bc_dict = dict()

        self.dem_df = assay_file

        self.libraries = list(self.dem_df['Library'])
        self.num_libraries = len(self.libraries)

        # this dict has keys for each library, and entries like
        # [total reads, failed on quality, failed on regex, failed on UMI, primer dimer likely]
        self.lib_stats = {l: [0, 0, 0, 0, 0] for l in self.libraries}
        self.other_reads = 0  # reads not in any of the libraries
        # this dict has keys that are the libraries and values that are dictionaries
        # these dictionaries have concatenated UMIs and entries that are counts in that UMI family
        self.umi_dict = {l: dict() for l in self.libraries}

    def add_count(self, bc_div, bc_env, library_index):
        bc = bc_div + bc_env
        if bc in self.bc_dict:
            tmp_entry = self.bc_dict[bc]
            tmp_entry[2] += 1
            tmp_entry[3][library_index] += 1
        else:
            tmp_entry = self.bc_dict[bc] = [bc_div, bc_env, 1, [0]*self.num_libraries]
            tmp_entry[3][library_index] += 1

    def read_files(self, R1_files, tmp_lib, paired):
        for R1in in R1_files:
            file_info = all_primer_info[R1in]
            dbc_start = file_info['R1_bp_to_BC']
            if paired:
                ebc_start = file_info['R2_bp_to_BC']
            else:
                ebc_start = dbc_start + 60
            reads1 = gzip.open(read_file_base+R1in, 'rt')
            if paired:
                reads2 = gzip.open(read_file_base+R1in.replace('R1.fastq.gz', 'R2.fastq.gz'), 'rt')
                R2_iterator = FourLineFastq(reads2)
            rc = 0
            for R1_title, R1_seq, R1_qual in FourLineFastq(reads1):
                rc += 1
                if paired:
                    R2_title, R2_seq, R2_qual = next(R2_iterator)
                else: # for unpaired sequencing (subpools) - we are filtering for inline index here
                    # this doesn't happen for the others because they have already been demultiplexed
                    if R1_seq[:8] != file_info['R1_index']:
                        break
                tmp_lib_stats = self.lib_stats[tmp_lib]
                tmp_lib_stats[0] += 1   # count for total reads
                # Quality check
                if paired:
                    qual_failed = np.mean([ord(c)-33 for c in (R1_qual[dbc_start:dbc_start+26] +
                                                    R2_qual[ebc_start:ebc_start+26])]) < QUALITY_CUTOFF
                else:
                    # for the unpaired subpool sequencing, just looking on R1
                    qual_failed = np.mean([ord(c)-33 for c in (R1_qual[dbc_start:dbc_start+26] +
                                                    R1_qual[ebc_start:ebc_start+26])]) < QUALITY_CUTOFF
                # quality check
                if qual_failed:
                    tmp_lib_stats[1] += 1
                else:
                    # regex check
                    reghit1 = MY_REGEX.match(R1_seq[dbc_start-10:dbc_start+36])
                    if paired:
                        reghit2 = MY_REGEX.match(R2_seq[ebc_start-10:ebc_start+36])
                    else:
                        reghit2 = MY_REGEX.match(reverse_transcribe(R1_seq[ebc_start-10:ebc_start+36]))
                    if (not reghit1) or (not reghit2):
                        # regex failed
                        tmp_lib_stats[2] += 1
                        # checks if regex failed because this was a primer-dimer fragment
                        if paired:
                            primer_dimer_fail = 'TTGAATTCGA' in R1_seq
                        else:
                            primer_dimer_fail = 'CTGTCTCTT' in R1_seq
                        if primer_dimer_fail:
                            tmp_lib_stats[4] += 1
                    else:
                        bc_div = reghit1.group(2)
                        bc_env = reghit2.group(2)
                        # umi check
                        if paired:
                            umi_combined = R1_seq[:8] + R2_seq[:8]
                            if umi_combined in self.umi_dict[tmp_lib]:
                                self.umi_dict[tmp_lib][umi_combined] += 1
                                tmp_lib_stats[3] += 1
                            else:
                                self.umi_dict[tmp_lib][umi_combined] = 1
                                # all checks passed, group 2 in the regex match is the barcode region
                                self.add_count(bc_div, bc_env, tmp_lib_ind)
                        else: # for unpaired subpool sequencing we don't worry about UMIs
                            self.add_count(bc_div, bc_env, tmp_lib_ind)

        reads1.close()
        if paired:
            reads2.close()

    def write_output(self, fout):

        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Diverse.BC', 'Environment.BC', 'Total.Counts'] + self.libraries)
            sorted_bcs = sorted(self.bc_dict, key=lambda x: self.bc_dict[x][2], reverse=True)
            for bc in sorted_bcs:
                entry = self.bc_dict[bc]
                writer.writerow(entry[:3]+entry[3])

    def write_lib_stats(self, fout):

        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Library', 'Total.Reads', 'Quality.Failed', 'Regex.Failed', 'UMI.Repeats', 'Primer.Dimer',
                             'Usable.Reads'])
            for lib in self.libraries:
                entry = self.lib_stats[lib]
                writer.writerow([lib] + entry + [entry[0]-entry[1]-entry[2]-entry[3]])

    def write_umi_fam_sizes(self, fout):
        # outputs rows with umi family size distributions for all libraries in the run
        biggest_fam = max([max(self.umi_dict[l].values()) for l in self.libraries if len(self.umi_dict[l]) > 0])
        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Library'] + [str(i) for i in range(1, biggest_fam + 1)])
            for lib in self.libraries:
                tmp_row = [lib]
                lib_fam_sizes = list(self.umi_dict[lib].values())
                tmp_row += [lib_fam_sizes.count(i) for i in range(1, biggest_fam + 1)]
                writer.writerow(tmp_row)


otime = time.time()
bcc = BcCounter(assay_file)

for j, row in assay_file.iterrows():
    paired = '001.fastq' in row['R1_files'] # all the pre-demultiplexed, unpaired subpool sequencing has this format
    bcc.read_files(row['R1_files'].split(';'), row['Library'], paired)

bcc.write_lib_stats(lib_stats_out)

bcc.write_output(output_file)

bcc.write_umi_fam_sizes(umi_fam_size_out)

make_html_report(lib_stats_out, umi_fam_size_out, html_report_out, run_name)

print('Time:', time.time()-otime)