#!/usr/bin/env python3
# This Python file uses the following encoding: utf-8
# -*- coding: utf-8 -*-
"""
Demultiplex paired end samples by barcodes on one read.
"""

import argparse
import os.path
import gzip



def parse_arguments():
    """Argument parser"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-file1',
        help='File containing sequences',
        required=True
    )
    parser.add_argument(
        '-file2',
        help='File containing reverse end sequences',
        required=True
    )
    parser.add_argument(
        '-barcodes',
        help='File containing barcodes',
        required=True
    )
    parser.add_argument(
        '-bc_pos',
        help='Start and end index of barcodes (default: 7 14)',
        type=int,
        nargs=2,
        default=[7, 14]
    )
    parser.add_argument(
        '-umi_pos',
        help='Start and end index of barcodes (default: 7 14)',
        type=int,
        nargs=2,
    )
    parser.add_argument(
        '-prefix',
        help='Filename prefixes',
        required=True
    )
    parser.add_argument(
        '-r2bc',
        help='Flag that barcode is on read 2',
        action='store_true'
    )
    parser.add_argument(
        '-ignr1',
        help='do not output R1',
        action='store_true'
    )
    parser.add_argument(
        '-count_sum',
        help='Flag for printing count summary',
        action='store_true'
    )
    args = parser.parse_args()
    return args

def mean(num_vec):
    """Take mean of a numeric vector."""
    return float(sum(num_vec)) / len(num_vec)


class Read:
    """Representation of a single read of a FASTQ file"""
    def __init__(self, conn):
        """Takes a file connection and returns a FASTQ structured read"""
        self.anno = conn.readline().strip()
        if not self.anno:
            raise EOFError
        self.seq = conn.readline().strip()
        self.strand = conn.readline().strip()
        self.quality_str = conn.readline().strip()

    def data(self):
        """Returns the read in FASTQ format"""
        return self.anno + '\n' + self.seq + '\n' + self.strand + '\n' + self.quality_str + '\n'

    def gc_content(self):
        """Returns the GC percentage in read"""
        cnt = {'N': 0, 'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for c_ind in self.seq:
            cnt[c_ind] += 1
        return 100.0 * (cnt['G'] + cnt['C']) / sum(cnt.values())

    def base_content(self, base):
        """Returns the percentage of particular base"""
        cnt = {'N': 0, 'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for c_ind in self.seq:
            cnt[c_ind] += 1
        return 100.0 * cnt[base] / sum(cnt.values())

    def base_count(self, base):
        """Returns the count of particular base"""
        cnt = {'N': 0, 'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for c_ind in self.seq:
            cnt[c_ind] += 1
        return cnt[base]

    def valid(self):
        """Check if a read is valid"""
        # check strand info
        if self.strand not in ['+', '-']:
            print(('Invalid strand info: ' + self.strand))
            return False
        # check length of sequence is equal to length of quality
        elif len(self.seq) != len(self.quality_str):
            print('Invalid read: sequence and quality of different length')
            return False
        else:
            return True

    def quality(self, start=None, end=None):
        """Returns the quality string either from 'start' onwards or from 'start' to 'end'"""
        if start:
            if end:
                return [ord(x) - 33 for x in self.quality_str[(start - 1):end]]
            else:
                return [ord(x) - 33 for x in self.quality_str[:start]]
        else:
            return [ord(x) - 33 for x in self.quality_str[:start]]

    def avg_quality(self, start=None, end=None):
        """Returns the average quality either from 'start' onwards or from 'start' to 'end'"""
        if start:
            if end:
                return mean([ord(x) - 33 for x in self.quality_str[(start - 1):end]])
            else:
                return mean([ord(x) - 33 for x in self.quality_str[:start]])
        else:
            return mean([ord(x) - 33 for x in self.quality_str[:start]])


class Paired:
    """Representation of paired end read"""
    def __init__(self, read1, read2):
        self.read1 = Read(read1)
        self.read2 = Read(read2)



def main():
    """Main"""
    args = parse_arguments()

    files = (args.file1, args.file2)
    barcodes = args.barcodes
    bc_pos = args.bc_pos
    umi_pos = args.umi_pos
    prefix = args.prefix
    flags = {
        'r2bc': args.r2bc,
        'count_sum': args.count_sum,
        'ignore_r1': args.ignr1
    }

    fastq_split(files, barcodes, bc_pos, umi_pos, prefix, flags)


def fastq_split(files, barcodes, bc_pos, umi_pos, prefix, flags):
    """Core function for splitting barcodes"""
    file1 = files[0]
    file2 = files[1]
    
    bc_start = bc_pos[0] - 1
    bc_end = bc_pos[1] - 1
    umi_start = umi_pos[0] - 1
    umi_end = umi_pos[1] - 1    

    read2_barcode = flags['r2bc']
    count_summary = flags['count_sum']

    with open(barcodes) as barcode_file:
        # read in barcodes
        bc_list = []
        for line in barcode_file:
            bc_list.append(line.strip().split())

        for barcode in bc_list:
            barcode.append(barcode[0])
            barcode[0] = prefix + '_' + barcode[0]

    # generate dictionary of tuples for barcode output files
    if not flags["ignore_r1"]:
        targets = {
            x[1]: {
                'output_files': (
                    open_fastq(x[0] + '-1_R1.fastq', 'w'),
                    open_fastq(x[0] + '-1_R2.fastq', 'w')
                ),
                'counter': 0
            }
            for x in bc_list
        }
    else:
        targets = {
            x[1]: {
                'output_files': (
                    None,
                    open_fastq(x[0] + '-1.demultiplexed.fq', 'w')
                ),
                'counter': 0
            }
            for x in bc_list
        } 
    # generate tuple for files for undetermined barcodes
    und = {
        'output_files': (
            open_fastq(prefix + '_Undetermined_R1.fastq', 'w'),
            open_fastq(prefix + '_Undetermined_R2.fastq', 'w')
        ),
        'counter': 0
    }

    with open_fastq(file1, 'r') as readfile1, open_fastq(file2, 'r') as readfile2:
        while True:
            try:
                read_pair = Paired(readfile1, readfile2)
                # get barcode portion
                if read2_barcode:
                    read_bc = read_pair.read2.seq[bc_start:(bc_end + 1)]
                    read_umi = read_pair.read2.seq[umi_start:(umi_end + 1)]
                else:
                    read_bc = read_pair.read1.seq[bc_start:(bc_end + 1)]
                    read_umi = read_pair.read1.seq[umi_start:(umi_end + 1)]

                # identify barcode, default undetermined
                output = targets.get(read_bc, und)
                # increase counter and write to files
                output['counter'] += 1
                if not flags["ignore_r1"]:
                    output['output_files'][0].write(read_pair.read1.data())
                #read_pair.read2.anno = ':'.join(read_pair.read2.anno.split('+')[0].split(":")[:-1]) + ':' + read_bc + '+' + read_umi
                output['output_files'][1].write(read_pair.read2.data())

            except EOFError:
                break

        if count_summary:
            print_count_summary(targets, und)
    

def print_count_summary(targets, und):
    """Print count summaries for each barcode, undetermined and total"""
    total_count = 0
    for barcode in targets:
        bc_counter = targets[barcode]['counter']
        print(f'{barcode}:  {bc_counter}')

        total_count += bc_counter

    undef_count = und['counter']
    print(f'Undefined:  {undef_count}')

    total_count += undef_count
    print(f'Total:  {total_count}')

def open_fastq(file_path, mode):
    """Generic opener for fastq files"""
    _, file_extension = os.path.splitext(file_path)
    if file_extension == '.gz':
        return gzip.open(file_path, mode + "t")
    else:
        return open(file_path, mode)

if __name__ == '__main__':
    main()
