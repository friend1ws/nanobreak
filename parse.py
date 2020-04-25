#! /usr/bin/env python

import sys, os, gzip, subprocess, csv, shutil, statistics, tempfile
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


def parse_alignment_info(input_bam, output_file, min_clipping_size = 200):

    hout = open(output_file, 'w') 

    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch():

        query_name = read.query_name
        query_strand = '-' if read.is_reverse else '+'
        query_length = read.infer_read_length()

        reference_name = read.reference_name
        reference_start = read.reference_start + 1
        reference_end = read.reference_end
        mapping_quality = read.mapping_quality
        is_secondary = read.is_secondary
        is_supplementary = read.is_supplementary


        cigartuples = read.cigartuples
        # left_hard_clipping_size, right_hard_clipping_size = 0, 0
        left_clipping_size, right_clipping_size = 0, 0
        if cigartuples[0][0] in (4, 5): left_hard_clipping_size = cigartuples[0][1]
        if cigartuples[-1][0] in (4, 5): right_hard_clipping_size = cigartuples[-1][1]
        
        if left_hard_clipping_size < min_clipping_size and right_hard_clipping_size < min_clipping_size: continue
 
        if not is_supplementary:
            if query_strand == '+':
                query_start = read.query_alignment_start + 1 
                query_end = read.query_alignment_end
            else:
                query_start = query_length - read.query_alignment_end + 1
                query_end = query_length - read.query_alignment_start
        else:
            if query_strand == '+':
                query_start = left_hard_clipping_size + 1
                query_end = query_length - right_hard_clipping_size
            else:
                query_start = right_hard_clipping_size + 1
                query_end = query_length - left_hard_clipping_size


        if left_hard_clipping_size >= min_clipping_size:

            if query_strand == '+':
                print("%s\t%d\t%d\t-\t%s\t%d\t+\t%d\t%d\t%d\t%s\t%s" % 
                    (reference_name, reference_start - 1, reference_start, query_name, query_length,
                     1, query_start - 1, mapping_quality, is_supplementary, is_secondary), file = hout)
            else:
                print("%s\t%d\t%d\t-\t%s\t%d\t-\t%d\t%d\t%d\t%s\t%s" % 
                    (reference_name, reference_start - 1, reference_start, query_name, query_length,
                    query_end + 1, query_length,  mapping_quality, is_supplementary, is_secondary), file = hout)

        elif right_hard_clipping_size >= min_clipping_size:
            
            if query_strand == '+':
                print("%s\t%d\t%d\t+\t%s\t%d\t+\t%d\t%d\t%d\t%s\t%s" % 
                    (reference_name, reference_end - 1, reference_end, query_name, query_length,
                    query_end + 1, query_length,  mapping_quality, is_supplementary, is_secondary), file = hout)
            else:
                print("%s\t%d\t%d\t+\t%s\t%d\t-\t%d\t%d\t%d\t%s\t%s" %
                    (reference_name, reference_end - 1, reference_end, query_name, query_length,
                    1, query_start - 1, mapping_quality, is_supplementary, is_secondary), file = hout)


    hout.close()
    bamfile.close()



def cluster_breakpoints(input_file, output_file, min_read_num = 3, min_median_mapq = 40.0, margin = 30):

    class Cluster_flusher(object):

        def __init__(self, output_h, min_read_num, min_median_mapq):
            self.chr = ''
            self.pos_plus = 0
            self.pos_minus = 0
            self.line_plus = []
            self.line_minus = []
            self.output_h = output_h
            self.min_read_num = min_read_num
            self.min_median_mapq = min_median_mapq
            self.cluster_id = 0

        def flush_lines(self, strand):

            tlines = self.line_plus if strand == '+' else self.line_minus

            if len(tlines) < self.min_read_num: return
            if statistics.median([int(x[9]) for x in tlines]) < self.min_median_mapq: return
    
            tstart = min([int(x[2]) for x in tlines]) - 1
            tend = max([int(x[2]) for x in tlines])
            # print("%s\t%s\t%d\t%d\t%s" % (self.cluster_id, self.chr, tstart, tend, strand), file = self.cluster_info_file_h)

            for tl in tlines:
                print("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % 
                    (self.cluster_id, self.chr, tstart, tend, strand, tl[4], tl[5], tl[6], tl[7], tl[8], tl[9], tl[10]), 
                    file = self.output_h)  

            self.cluster_id = self.cluster_id + 1
    

    with gzip.open(input_file, 'rt') as hin, open(output_file, 'w') as hout:

        cluster_flusher = Cluster_flusher(hout, min_read_num, min_median_mapq)
        for row in csv.reader(hin, delimiter = '\t'):

            # skip secondary alignment
            if row[11] == "True": continue
        
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
            shutil.rmtree(self.tmp_dir)

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

            with open(self.tmp_dir + '/' + self.cluster_id + "_ref_polished.fa", 'w') as hout_ref2:
                subprocess.check_call(["racon", self.tmp_dir + '/' + self.cluster_id + "_reads.fa",
                                       self.tmp_dir + '/' + self.cluster_id + "_ova_minimap2.paf",
                                       self.tmp_dir + '/' + self.cluster_id + "_ref.fa"], stderr = subprocess.DEVNULL, stdout = hout_ref2)

            consensus = ''
            with open(self.tmp_dir + '/' + self.cluster_id + "_ref_polished.fa", 'r') as hin_ref2:
                for line in hin_ref2:
                    if line.startswith('>'): continue
                    consensus = consensus + line.rstrip('\n')
                        
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

 



if __name__ == "__main__":

    import sys

    input_bam = sys.argv[1]
    output_prefix = sys.argv[2]

    """
    parse_alignment_info(input_bam, output_prefix + ".tmp.nbp_parse.bed")

    hout = open(output_prefix + ".tmp.nbp_parse.sorted.bed", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", output_prefix + ".tmp.nbp_parse.bed"], stdout = hout)
    hout.close()

    hout = open(output_prefix + ".tmp.nbp_parse.sorted.bed.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", output_prefix + ".tmp.nbp_parse.sorted.bed"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", output_prefix + ".tmp.nbp_parse.sorted.bed.gz"])

    os.remove(output_prefix + ".tmp.nbp_parse.bed")
    os.remove(output_prefix + ".tmp.nbp_parse.sorted.bed")

    cluster_breakpoints(output_prefix + ".tmp.nbp_parse.sorted.bed.gz", 
                        output_prefix + ".tmp.nbp_parse.cluster_info.txt")
    
    add_sequence(output_prefix + ".tmp.nbp_parse.cluster_info.txt", input_bam, output_prefix + ".tmp.nbp_parse.cluster_info_seq.txt")
    """

    add_consensus_seq(output_prefix + ".tmp.nbp_parse.cluster_info_seq.txt", output_prefix + ".tmp.nbp_parse.cluster_info_consensus.txt")

