#! /usr/bin/env python

import pysam


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



parse_alignment_info(snakemake.input[0], snakemake.output[0])


