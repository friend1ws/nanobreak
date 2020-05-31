#! /usr/bin/env python

import csv, subprocess, tempfile

def add_consensus_seq(input_file, output_file, start_margin = 120, min_inclusion_ratio = 0.9, min_inclusion_count = 3):

    class Seq_gatherer(object):

        def __init__(self, output_file, start_margin, min_inclusion_ratio, min_inclusion_count):

            self.cluster_id = ''
            self.clustser_info = ''
            self.readid2seq = {}
            self.tmp_dir = tempfile.mkdtemp()
            self.start_margin = start_margin
            self.min_inclusion_ratio = min_inclusion_ratio
            self.min_inclusion_count = min_inclusion_count
            self.hout = open(output_file, 'w')
            
        def __del__(self):
            self.hout.close()
            # shutil.rmtree(self.tmp_dir)

        def initialize(self, cluster_id, cluster_info):
            self.cluster_id = cluster_id
            self.cluster_info = cluster_info
            self.readid2seq = {}

        def print_consensus(self):

            if len(self.readid2seq) == 0: return

            with open(self.tmp_dir + '/' + self.cluster_id + "_reads.fa", 'w') as hout_fa: 
                for readid in self.readid2seq:
                    print(">%s\n%s" % (readid, self.readid2seq[readid]), file = hout_fa)

            with open(self.tmp_dir + '/' + self.cluster_id + "_ava_minimap2.paf", 'w') as hout_paf:
                subprocess.check_call(["minimap2", "-x", "ava-ont", self.tmp_dir + '/' + self.cluster_id + "_reads.fa", self.tmp_dir + '/' + self.cluster_id + "_reads.fa"], 
                                      stderr = subprocess.DEVNULL, stdout = hout_paf)

            readid2inclusion_count = {}
            with open(self.tmp_dir + '/' + self.cluster_id + "_ava_minimap2.paf") as hin:
                for row in csv.reader(hin, delimiter = '\t'):
                    if int(row[2]) <= self.start_margin and (float(row[3]) - float(row[2])) / float(row[1]) >= self.min_inclusion_ratio:
                        if row[5] not in readid2inclusion_count: readid2inclusion_count[row[5]] = 0
                        readid2inclusion_count[row[5]] = readid2inclusion_count[row[5]] + 1
                    if int(row[7]) <= self.start_margin and (float(row[8]) - float(row[7])) / float(row[6]) >= self.min_inclusion_ratio:
                        if row[0] not in readid2inclusion_count: readid2inclusion_count[row[0]] = 0
                        readid2inclusion_count[row[0]] = readid2inclusion_count[row[0]] + 1

            if len(readid2inclusion_count) == 0: return
            readid_max = max(readid2inclusion_count, key = readid2inclusion_count.get)
            if readid2inclusion_count[readid_max] < self.min_inclusion_count: return

            with open(self.tmp_dir + '/' + self.cluster_id + "_ref.fa", 'w') as hout_ref:
                for readid in self.readid2seq:
                    if readid == readid_max:
                        print(">%s\n%s" % (readid, self.readid2seq[readid]), file = hout_ref)

            with open(self.tmp_dir + '/' + self.cluster_id + "_ova_minimap2.paf", 'w') as hout_paf2:
                subprocess.check_call(["minimap2", "-x", "map-ont", self.tmp_dir + '/' + self.cluster_id + "_ref.fa",
                                       self.tmp_dir + '/' + self.cluster_id + "_reads.fa"], stderr = subprocess.DEVNULL, stdout = hout_paf2)

            consensus = ''
            try:
                with open(self.tmp_dir + '/' + self.cluster_id + "_ref_polished.fa", 'w') as hout_ref2:
                    subprocess.check_call(["racon", self.tmp_dir + '/' + self.cluster_id + "_reads.fa",
                                           self.tmp_dir + '/' + self.cluster_id + "_ova_minimap2.paf",
                                           self.tmp_dir + '/' + self.cluster_id + "_ref.fa"], stderr = subprocess.DEVNULL, stdout = hout_ref2)

                with open(self.tmp_dir + '/' + self.cluster_id + "_ref_polished.fa", 'r') as hin_ref2:
                    for line in hin_ref2:
                        if line.startswith('>'): continue
                        consensus = consensus + line.rstrip('\n')
            except:
                consensus = ''

            if len(consensus) < 1000: return
                        
            print("%s\t%s\t%s" % (self.cluster_id, self.cluster_info, consensus), file = self.hout)

    seq_gatherer = Seq_gatherer(output_file, start_margin, min_inclusion_ratio, min_inclusion_count)

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        for row in csv.reader(hin, delimiter = '\t'):
            if row[0] != seq_gatherer.cluster_id:
                seq_gatherer.print_consensus()
                seq_gatherer.initialize(row[0], '\t'.join(row[1:5]))

            seq_gatherer.readid2seq[row[5]] = row[12]

        seq_gatherer.print_consensus()

    del seq_gatherer


add_consensus_seq(snakemake.input[0], snakemake.output[0])

