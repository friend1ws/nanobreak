#! /usr/bin/env python

import os, csv, subprocess
import pysam


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


def add_sequence(input_file, bam_file, output_file, bp_margin = 300):

    readid2info = {}
    with open(input_file, 'r') as hin:
        for row in csv.reader(hin, delimiter = '\t'):
            cluster_info = row[:5]
            read_info = row[6:]
            readid = row[5]
    
            if readid not in readid2info: readid2info[readid] = []
            readid2info[readid].append((cluster_info, read_info))

 
    with pysam.AlignmentFile(bam_file, 'rb') as hbam, open(output_file + ".tmp.unsorted", 'w') as hout:
        for read in hbam.fetch():

            if read.query_name in readid2info and not read.is_secondary and not read.is_supplementary:

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                for infos in readid2info[read.query_name]:
                    cluster_info, read_info = infos
                    cstrand = cluster_info[4]
                    qlen, qstrand, qstart, qend = int(read_info[0]), read_info[1], int(read_info[2]), int(read_info[3])

                    if cstrand != qstrand:
                        part_seq = read_seq[(qstart - 1):min(qlen, qend + bp_margin)]
                        part_seq = reverse_complement(part_seq)
                    else:
                        part_seq = read_seq[max(0, qstart - 1 - bp_margin):qend]

                    print('\t'.join(cluster_info) + '\t' + read.query_name + '\t' + '\t'.join(read_info) + '\t' + part_seq, file = hout)

    with open(output_file, 'w') as hout:
        subprocess.check_call(["sort", "-k1,1n", output_file + ".tmp.unsorted"], stdout = hout)

    os.remove(output_file + ".tmp.unsorted")


add_sequence(snakemake.input[0], snakemake.input[1], snakemake.output[0])

