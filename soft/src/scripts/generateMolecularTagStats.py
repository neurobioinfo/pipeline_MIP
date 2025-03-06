#!/usr/bin/env python

import os
import sys
import pysam
import subprocess
import getopt
import argparse


#parse script arguments
def usage():
    print sys.argv[0] + ' -m|--bam-list <comma_seperated_list_of_bamfiles> (required) -b|--bed <bedfile> (optional)'

parser = argparse.ArgumentParser(description='Generate molecular tag stats for somatic mip-captured bamfiles')
parser.add_argument('-m', '--bam-file', dest='bams', action='append', help='/path/to/input_bamfile (can be specified multiple times)', required=True)
parser.add_argument('-b', '--bed', dest='bed', help='/path/to/input_bedfile')
args   = parser.parse_args()
bams   = args.bams
bed    = args.bed
#if --bed argument given: make sure file exists and is non-empty
if bed is not None:
    try:
        s = os.stat(bed)
        if s.st_size == 0:
            print "The bed file file {} is empty".format(bed)
            sys.exit(42)
    except OSError as e:
        print e
        sys.exit(42)

#parse input data into global variables
mips, data = ([], {})
for bam in bams:
    reads, reads_processed = ([], {})
    samfile = pysam.AlignmentFile(bam, 'rb')
    if bam not in data: data[bam]={}
    if bed is None:
        #if no bedfile, process entire bam
        for read in samfile.fetch(until_eof=True): reads.append(read)
    else:
        #if bedfile, process only reads in bed region(s)
        regions	= [region.split("\t") for region in subprocess.check_output(['/RQexec/spiegelm/soft/src/bedops-2.4.2/bedops', '--merge', bed]).strip().split("\n")]
        #sys.exit(regions)
        for region in regions:
            for read in samfile.fetch(region[0], int(region[1]), int(region[2]), until_eof=True): reads.append(read)
    for read in reads:
        #little fix to only read each read once (may appear more than once if it overlaps >1 bed region)
        if read.query_name + str(read.is_read1) in reads_processed: continue
        reads_processed[read.query_name + str(read.is_read1)]=1
        #now go on to process mip data for the read
        mip, tag = ("_".join(read.query_name.split("_")[0:2]), read.query_name.split("_")[2])
        if mip not in mips: mips.append(mip)
        if mip not in data[bam]: data[bam][mip]={}
        if tag not in data[bam][mip]: data[bam][mip][tag]=0
        data[bam][mip][tag]+=1

#print header lines: 1 = list of bamfiles, 3 times; 2 = ['#tags', '#molecules'], len(list of bamfiles) times
outputCategories = ['#tags', '#reads']
print "\t".join(["mip", "\t".join(["\t".join([bam for bam in bams]) for cat in outputCategories])])
print "\t".join(["mip", "\t".join(["\t".join([cat for bam in bams]) for cat in outputCategories])])

#print data, mip by mip
for mip in sorted(mips, key=lambda mip: (int(mip.partition('_')[0]))):
    tagData = map(str, [len(data[bam][mip]) if mip in data[bam] else 0 for bam in bams])
    molData = map(str, [sum([data[bam][mip][tag] for tag in data[bam][mip].keys()]) if mip in data[bam] else 0 for bam in bams])
    print "\t".join([mip, "\t".join(tagData), "\t".join(molData)])
