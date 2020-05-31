#! /usr/bin/env python

import csv, gzip, statistics
import pysam


def cluster_breakpoints(tumor_bp_file, normal_bp_file, output_file, min_tumor_read_num = 3, max_normal_read_num = 2, min_median_mapq = 40.0, margin = 30):

    class Cluster_flusher(object):

        def __init__(self, output_h, normal_bp_tb, min_tumor_read_num, max_normal_read_num, min_median_mapq):
            self.chr = ''
            self.pos_plus = 0
            self.pos_minus = 0
            self.line_plus = []
            self.line_minus = []
            self.output_h = output_h
            self.normal_bp_tb = normal_bp_tb
            self.min_tumor_read_num = min_tumor_read_num
            self.min_normal_read_num = max_normal_read_num
            self.min_median_mapq = min_median_mapq
            self.cluster_id = 0

        def flush_lines(self, strand):

            tlines = self.line_plus if strand == '+' else self.line_minus

            if len(tlines) < self.min_tumor_read_num: return
            if statistics.median([int(x[9]) for x in tlines]) < self.min_median_mapq: return

            tstart = min([int(x[2]) for x in tlines]) - 1
            tend = max([int(x[2]) for x in tlines])
 
            normal_count = 0
            tabix_error_flag = False
            try:            
                records = self.normal_bp_tb.fetch(self.chr, tstart - 5, tend + 5)
            except:
                tabix_error_flag = True

            if not tabix_error_flag:
                for rec in records:
                    normal_count = normal_count + 1
        
            if normal_count >= self.min_normal_read_num: return

            for tl in tlines:
                print("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
                    (self.cluster_id, self.chr, tstart, tend, strand, tl[4], tl[5], tl[6], tl[7], tl[8], tl[9], tl[10]), 
                    file = self.output_h)  

            self.cluster_id = self.cluster_id + 1
    
    with gzip.open(tumor_bp_file, 'rt') as hin, open(output_file, 'w') as hout, pysam.TabixFile(normal_bp_file) as normal_bp_tb:

        cluster_flusher = Cluster_flusher(hout, normal_bp_tb, min_tumor_read_num, max_normal_read_num, min_median_mapq)
        for row in csv.reader(hin, delimiter = '\t'):

            # skip secondary alignment
            if row[11] == "True": continue
            
            # skip decoy
            if "decoy" in row[0]: continue
    
            if row[0] != cluster_flusher.chr:
                cluster_flusher.flush_lines('+')
                cluster_flusher.flush_lines('-')
                cluster_flusher.chr = row[0]
                cluster_flusher.line_plus, cluster_flusher.line_minus = [], []

            elif row[3] == '+' and int(row[2]) > cluster_flusher.pos_plus + margin:
                cluster_flusher.flush_lines('+')
                cluster_flusher.line_plus = []
            
            elif row[3] == '-' and int(row[2]) > cluster_flusher.pos_minus + margin:
                cluster_flusher.flush_lines('-')
                cluster_flusher.line_minus = []

            if row[3] == '+':
                cluster_flusher.line_plus.append(row)
                cluster_flusher.pos_plus = int(row[2])

            elif row[3] == '-':
                cluster_flusher.line_minus.append(row)
                cluster_flusher.pos_minus = int(row[2])

        cluster_flusher.flush_lines('+')
        cluster_flusher.flush_lines('-')


cluster_breakpoints(snakemake.input[0], snakemake.input[1], snakemake.output[0])

