#! /usr/bin/env python

import sys, os, gzip, subprocess, csv, shutil, statistics, tempfile
import pysam
import nanomonsv

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

            print(self.seq_start, self.seq_end, self.seq_dir)
            print(alignment_info["query_" + self.cluster_id])

            if "query_" + self.cluster_id not in alignment_info: return
            _, tstart_a, tend_a, qstart_a, qend_a, strand_a = alignment_info["query_" + self.cluster_id]
            if strand_a != '+': return

            self.bp_pos_consensus = tend_a
            if self.seq_dir == '+':
                self.bp_pos_reference = self.seq_end - (len(self.seq_around_bp) - qend_a) 
            else:
                self.bp_pos_reference = self.seq_start + (len(self.seq_around_bp) - qend_a)

            print(self.bp_pos_reference, self.bp_pos_consensus)
            """
            60001 60203 -
            {'query_0': [374, 67, 265, 1, 203, '+']}
            
            """


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



def long_read_validate(input_file, tumor_bam, normal_bam, output_file, reference, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8,
    min_tumor_variant_read_num = 3, min_tumor_VAF = 0.05, max_control_variant_read_num = 1, max_control_VAF = 0.03):

    def count_read_by_alignment(input_file, bam_file, output_file, reference, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8):

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
    
        key2count = {}
        alignment_counter = Alignment_counter(score_ratio_thres, start_pos_thres, end_pos_thres)
        with open(output_file + ".tmp.long_read_seq.sorted", 'r') as hin:
            for row in csv.reader(hin, delimiter = '\t'):
                key = '\t'.join(row[:4])
                if key != alignment_counter.key:
                    acount = alignment_counter.count_alignment()
                    if acount is not None:
                        key2count[alignment_counter.key] = acount

                    alignment_counter.initialize(key, key2contig[key])
                alignment_counter.add_query_seq(row[4], row[5])

            acount = alignment_counter.count_alignment()
            if acount is not None:
                key2count[alignment_counter.key] = acount

        del alignment_counter
        os.remove(output_file + ".tmp.long_read_seq.sorted")


        return key2count


    key2count_tumor = count_read_by_alignment(input_file, tumor_bam, output_file, reference, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8)
    key2count_normal = count_read_by_alignment(input_file, normal_bam, output_file, reference, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8)

    key2contig = {}
    with open(input_file, 'r') as hin:
        for row in csv.reader(hin, delimiter = '\t'):
            key = '\t'.join(row[:4])
            key2contig[key] = row[4]

    with open(output_file, 'w') as hout:
        for key in key2count_tumor:
            count_tumor = key2count_tumor[key]
            count_normal = key2count_normal[key] if key in key2count_normal else (0, 0)

            if count_tumor[0] == 0: continue
            if count_tumor[1] < min_tumor_variant_read_num: continue
            if float(count_tumor[1]) / float(count_tumor[0]) < min_tumor_VAF: continue

            if count_normal[0] == 0: continue
            if count_normal[1] > max_control_variant_read_num: continue
            if float(count_normal[1]) / float(count_normal[0]) > max_control_VAF: continue

            print("%s\t%d\t%d\t%d\t%d\t%s" % (key, count_tumor[0], count_tumor[1], count_normal[0], count_normal[1], key2contig[key]), file = hout)


 
if __name__ == "__main__":

    import sys

    tumor_bam = sys.argv[1]
    normal_bam = sys.argv[2]
    output_prefix = sys.argv[3]
    reference_file = sys.argv[4]

    ##########
    # tumor parse
    parse_alignment_info(tumor_bam, output_prefix + ".tmp.tumor.bp_parse.bed")

    hout = open(output_prefix + ".tmp.tumor.bp_parse.sorted.bed", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", output_prefix + ".tmp.tumor.bp_parse.bed"], stdout = hout)
    hout.close()

    hout = open(output_prefix + ".tmp.tumor.bp_parse.sorted.bed.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", output_prefix + ".tmp.tumor.bp_parse.sorted.bed"], stdout = hout)
    hout.close()
    
    subprocess.check_call(["tabix", "-p", "bed", output_prefix + ".tmp.tumor.bp_parse.sorted.bed.gz"])
    
    os.remove(output_prefix + ".tmp.tumor.bp_parse.bed")
    os.remove(output_prefix + ".tmp.tumor.bp_parse.sorted.bed")
    ##########

    ##########
    # normal parse
    parse_alignment_info(normal_bam, output_prefix + ".tmp.normal.bp_parse.bed")

    hout = open(output_prefix + ".tmp.normal.bp_parse.sorted.bed", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", output_prefix + ".tmp.normal.bp_parse.bed"], stdout = hout)
    hout.close()
    
    hout = open(output_prefix + ".tmp.normal.bp_parse.sorted.bed.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", output_prefix + ".tmp.normal.bp_parse.sorted.bed"], stdout = hout)
    hout.close()
    
    subprocess.check_call(["tabix", "-p", "bed", output_prefix + ".tmp.normal.bp_parse.sorted.bed.gz"])
        
    os.remove(output_prefix + ".tmp.normal.bp_parse.bed")
    os.remove(output_prefix + ".tmp.normal.bp_parse.sorted.bed")
    ##########

    cluster_breakpoints(output_prefix + ".tmp.tumor.bp_parse.sorted.bed.gz", 
                        output_prefix + ".tmp.normal.bp_parse.sorted.bed.gz",
                        output_prefix + ".tmp.nbp_parse.cluster_info.txt")
  
    add_sequence(output_prefix + ".tmp.nbp_parse.cluster_info.txt", tumor_bam, output_prefix + ".tmp.nbp_parse.cluster_info_seq.txt")

    add_consensus_seq(output_prefix + ".tmp.nbp_parse.cluster_info_seq.txt", output_prefix + ".tmp.nbp_parse.cluster_info_consensus.txt")


    locate_breakpoint(output_prefix + ".tmp.nbp_parse.cluster_info_consensus.txt", 
        output_prefix + ".tmp.nbp_parse.bp_info_consensus.txt",
        reference_file)

    long_read_validate(output_prefix + ".tmp.nbp_parse.bp_info_consensus.txt", tumor_bam, normal_bam, output_prefix + ".tmp.bp_cluster.alignment_count.txt", reference_file)
