#! /usr/bin/env python

import csv

def join_count_result(consensus_file, tumor_count_file, normal_count_file, output_file, min_tumor_variant_read_num = 3, min_tumor_VAF = 0.05, max_control_variant_read_num = 1, max_control_VAF = 0.03):

    key2tumor_count = {}
    with open(tumor_count_file, 'r') as hin:
        for row in csv.reader(hin, delimiter = '\t'):
            key = '\t'.join(row[:4])
            key2tumor_count[key] = (int(row[4]), int(row[5]))

    key2normal_count = {}
    with open(normal_count_file, 'r') as hin:
        for row in csv.reader(hin, delimiter = '\t'):
            key = '\t'.join(row[:4])
            key2normal_count[key] = (int(row[4]), int(row[5]))

    with open(consensus_file, 'r') as hin, open(output_file, 'w') as hout:
        for row in csv.reader(hin, delimiter = '\t'):
            key = '\t'.join(row[:4])
            consensus = row[4]
            if key not in key2tumor_count or key not in key2normal_count: continue
            tumor_count = key2tumor_count[key]
            normal_count = key2normal_count[key]

            if tumor_count[0] == 0: continue
            if normal_count[0] == 0: continue
            if tumor_count[1] < min_tumor_variant_read_num: continue
            if normal_count[1] > max_control_variant_read_num: continue
            if float(tumor_count[1]) / tumor_count[0] < min_tumor_VAF: continue
            if float(normal_count[1]) / normal_count[0] > max_control_VAF: continue

            print("%s\t%d\t%d\t%d\t%d\t%s" % (key, tumor_count[0], tumor_count[1], normal_count[0], normal_count[1], consensus), file = hout)


join_count_result(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.output[0])


