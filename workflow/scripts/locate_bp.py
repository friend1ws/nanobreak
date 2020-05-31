#! /usr/bin/env python

import csv, shutil, tempfile
import pysam
import nanomonsv

from utils import reverse_complement

def locate_breakpoint(input_file, output_file, reference_file, margin = 200):

    class Breakpoint_locator(object):

        def __init__(self):

            self.consensus = None
            self.seq_around_bp = None
            self.seq_start = None
            self.seq_end = None
            self.seq_dir = None
            self.bp_pos_consensus = None
            self.bp_pos_reference = None
            self.tmp_dir = tempfile.mkdtemp()


        def __del__(self):
            shutil.rmtree(self.tmp_dir)            

        def initialize(self, cluster_id, consensus, seq_around_bp, seq_start, seq_end, seq_dir):

            self.cluster_id = cluster_id
            self.consensus = consensus
            self.seq_around_bp = seq_around_bp
            self.seq_start = seq_start
            self.seq_end = seq_end
            self.seq_dir = seq_dir
            self.bp_pos_consensus = None
            self.bp_pos_reference = None
            if len(self.consensus) >= 1000: self.consensus = self.consensus[:1000]

        def locate_by_alignment(self):
    
            with open(self.tmp_dir + '/' + self.cluster_id + ".query.fa", 'w') as hout:
                print(">query_%s\n%s" % (self.cluster_id, self.seq_around_bp), file = hout)

            with open(self.tmp_dir + '/' + self.cluster_id + ".target.fa", 'w') as hout:
                print(">target_%s\n%s" % (self.cluster_id, self.consensus), file = hout)

            alignment_info = nanomonsv.long_read_validate.ssw_check(
                self.tmp_dir + '/' + self.cluster_id + ".target.fa", 
                self.tmp_dir + '/' + self.cluster_id + ".query.fa")

            # print(self.seq_start, self.seq_end, self.seq_dir)
            # print(alignment_info["query_" + self.cluster_id])

            if "query_" + self.cluster_id not in alignment_info: return
            _, tstart_a, tend_a, qstart_a, qend_a, strand_a = alignment_info["query_" + self.cluster_id]
            if strand_a != '+': return

            self.bp_pos_consensus = tend_a
            if self.seq_dir == '+':
                self.bp_pos_reference = self.seq_end - (len(self.seq_around_bp) - qend_a) 
            else:
                self.bp_pos_reference = self.seq_start + (len(self.seq_around_bp) - qend_a)

            # print(self.bp_pos_reference, self.bp_pos_consensus)


    bp_loc = Breakpoint_locator()
    fasta_file = pysam.FastaFile(reference_file)

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        for row in csv.reader(hin,delimiter = '\t'):
            if row[4] == '+':
                seq_around_bp = fasta_file.fetch(row[1], int(row[2]) - margin, int(row[3]))
                bp_loc.initialize(row[0], row[5], seq_around_bp, int(row[2]) - margin + 1, int(row[3]), '+')
            else:
                seq_around_bp = fasta_file.fetch(row[1], int(row[2]), int(row[3]) + margin)
                seq_around_bp = reverse_complement(seq_around_bp)
                bp_loc.initialize(row[0], row[5], seq_around_bp, int(row[2]) + 1, int(row[3]) + margin, '-')

            bp_loc.locate_by_alignment()
    
            if bp_loc.bp_pos_reference is not None:
                print("%s\t%s\t%d\t%s\t%s" % (row[0], row[1], bp_loc.bp_pos_reference, row[4], row[5][bp_loc.bp_pos_consensus:]), file = hout)
  
    del bp_loc
    fasta_file.close()


locate_breakpoint(snakemake.input[0], snakemake.output[0], snakemake.input[1])

