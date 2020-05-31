#! /usr/bin/env python

import os, shutil, csv, tempfile, subprocess
import pysam, edlib

from utils import reverse_complement


def count_read_by_alignment(input_file, bam_file, reference, output_file, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8):

    class Alignment_counter(object):

        def __init__(self, score_ratio_thres, start_pos_thres, end_pos_thres):
            self.key = ''
            # self.hout = open(output_file, 'w')
            self.query_seq = None
            self.target_seq_list = []
            self.score_ratio_thres = score_ratio_thres
            self.start_pos_thres = start_pos_thres
            self.end_pos_thres = end_pos_thres
            self.tmp_dir = tempfile.mkdtemp()
            

        def __del__(self):
            # self.hout.close()
            shutil.rmtree(self.tmp_dir)

        def initialize(self, key, query_seq):
            self.key = key
            self.query_seq = query_seq
            self.target_seq_list = []

        def add_query_seq(self, read_name, read_seq):
            self.target_seq_list.append((read_name, read_seq))

        def count_alignment(self):

            if len(self.target_seq_list) == 0: return None
            cluster_id, _, _, _ = key.split('\t')

            """
            with open(self.tmp_dir + '/' + cluster_id + ".target.fa", 'w') as hout_ta:
                for read in self.target_seq_list:
                    print(">%s\n%s" % (read[0], read[1]), file = hout_ta)

            with open(self.tmp_dir + '/' + cluster_id + ".query.fa", 'w') as hout_qu:
                print(">query_%s\n%s" % (cluster_id, self.query_seq), file = hout_qu)

            alignment_info = nanomonsv.long_read_validate.ssw_check(
                self.tmp_dir + '/' + cluster_id + ".query.fa",
                self.tmp_dir + '/' + cluster_id + ".target.fa")

            all_rnames = list(set(alignment_info.keys()))
            supporting_reads = [rname for rname in all_rnames if \
                alignment_info[rname][0] > self.score_ratio_thres * len(self.query_seq) and \
                alignment_info[rname][1] < self.start_pos_thres * len(self.query_seq) and \
                alignment_info[rname][2] > self.end_pos_thres * len(self.query_seq)]

            if "15033196" in self.key:
                for a in alignment_info: print(a, alignment_info[a], a in supporting_reads)
                print(' ')
            """
            # print(cluster_id)

            all_rnames = []
            supporting_reads = []
            align_res = []
            for target_read in self.target_seq_list:
                all_rnames.append(target_read[0])                                                          
                tres = edlib.align(self.query_seq, target_read[1], mode = "HW", task = "locations")
                # print(target_read[0])
                # print(tres) 
                if tres["editDistance"] < len(self.query_seq) * 0.25: # and \
                    # tres["locations"][0][0] < 0.2 * len(self.query_seq) and \
                    # tres["locations"][0][1] > 0.8 * len(self.query_seq):
                        supporting_reads.append(target_read[0])

            return (len(all_rnames), len(supporting_reads))


    bam_hin = pysam.AlignmentFile(bam_file, 'rb')
    reference_fasta = pysam.FastaFile(reference)
    rname2key = {}
    key2contig = {}
    with open(input_file, 'r') as hin:
        for row in csv.reader(hin, delimiter = '\t'):
            key = '\t'.join(row[:4])
            for read in bam_hin.fetch(row[1], max(int(row[2]) - 100, 0), int(row[2]) + 100):
                if read.is_secondary: continue
                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

            if row[3] == '+':
                contig = reference_fasta.fetch(row[1], max(int(row[2]) - validate_sequence_length - 1, 0), int(row[2]))
            else:
                contig = reference_fasta.fetch(row[1], int(row[2]) - 1, int(row[2]) + validate_sequence_length - 1)
                contig = reverse_complement(contig)
            contig = contig + row[4][:validate_sequence_length]
            key2contig[key] = contig

    for rname in rname2key:
        keys = list(set(rname2key[rname]))
        rname2key[rname] = keys

    with open(output_file + ".tmp.long_read_seq.unsorted", 'w') as hout:
        for read in bam_hin.fetch():

            if read.is_secondary or read.is_supplementary: continue
            if read.qname in rname2key:
                read_seq = read.query_sequence
                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                for key in rname2key[read.qname]:
                    print("%s\t%s\t%s" % (key, read.qname, read_seq), file = hout)

    bam_hin.close()

    with open(output_file + ".tmp.long_read_seq.sorted", 'w') as hout:
        subprocess.check_call(["sort", "-k1,1", output_file + ".tmp.long_read_seq.unsorted"], stdout = hout)
    os.remove(output_file + ".tmp.long_read_seq.unsorted")

    # key2count = {}
    alignment_counter = Alignment_counter(score_ratio_thres, start_pos_thres, end_pos_thres)
    with open(output_file + ".tmp.long_read_seq.sorted", 'r') as hin, open(output_file, 'w') as hout:
        for row in csv.reader(hin, delimiter = '\t'):
            key = '\t'.join(row[:4])
            if key != alignment_counter.key:
                acount = alignment_counter.count_alignment()
                if acount is not None:
                    print("%s\t%d\t%d" % (alignment_counter.key, acount[0], acount[1]), file = hout)
                    # key2count[alignment_counter.key] = acount

                alignment_counter.initialize(key, key2contig[key])
            alignment_counter.add_query_seq(row[4], row[5])

        acount = alignment_counter.count_alignment()
        if acount is not None:
            print("%s\t%d\t%d" % (alignment_counter.key, acount[0], acount[1]), file = hout)
            # key2count[alignment_counter.key] = acount

    del alignment_counter
    os.remove(output_file + ".tmp.long_read_seq.sorted")


count_read_by_alignment(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.output[0])

