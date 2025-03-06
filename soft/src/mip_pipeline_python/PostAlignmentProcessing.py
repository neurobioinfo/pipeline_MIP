#!/usr/bin/python

import os
import sys
import subprocess
import gzip
import Design
import ReadPair
import Mip
#import Alignment
import HomologySearch
import swalign
import pysam
import printStruct
import random
from MipPipelineHelperFunctions import *
from swalign import *


########
# MAIN #
########
# process inputs; assign + populate global variables
bam, designFile, designSource = sys.argv[1:]
design = Design.Design(designFile, designSource)
design.parse()

#add list of AlignedReads to each mip, from aligned bam file; for each mip then count how many reads are on vs offtarget
samfile = pysam.AlignmentFile(bam, 'rb')
#for read in samfile.fetch(until_eof=True): mip.Get("_".join(read.query_name.split("_")[0:len(read.query_name.split("_"))-1]))._AlignedSegments.append(read)	# this breaks when doing somatic mips, so:
#for read in samfile.fetch(until_eof=True): mip.Get("_".join(read.query_name.split("_")[0:2]))._AlignedSegments.append(read)					# and this breaks if mip name looks like "1_3_5.0", i.e. if there is more than one underscore-separated digit before the [5,3,-1,5.0,3.0] of the mip name, so:
#for read in samfile.fetch(until_eof=True): design._mips[0].Get("_".join([x for x in read.qname.split(':')[0].split('_') if x.replace('-','').replace('.','').isdigit()]))._AlignedSegments.append(read)	# this breaks if the mip name doesn't look like \d+_\d+(_[\d.-]+); i.e mipgen names like 'MNG_coding_0895'
#for read in samfile.fetch(until_eof=True): design._mips[0].Get('_'.join([x for i,x in enumerate(read.qname.split(':')[0].split('_'))][0:[i for i,x in enumerate(read.qname.split(':')[0].split('_')) if str(x).isdigit()][-1]+1]))._AlignedSegments.append(read)	# this could break if the machine_id part of the read name contains a string like '_\d+_'
mips_strings = map(str, design._mips)
for read in samfile.fetch(until_eof=True):		# this way actually parses the longest possible string from the beginning of the read name, which matches a mip
  mip_string = ['_'.join(read.qname.split(':')[0].split('_')[0:i+1]) for i,x in enumerate(read.qname.split(':')[0].split('_')) if '_'.join(read.qname.split(':')[0].split('_')[0:i+1]) in mips_strings][-1]
  design._mips[0].Get(mip_string)._AlignedSegments.append(read)


#for each mip, count # of reads aligning on vs offtarget vs unmapped
# output file name: [sample bamfile].mip-by-mip_counts
# output file format: col1=mip_name; col2=ontarget_count; col3=offtarget_count; col4=unmapped_count
# NOTE: these stats only apply to reads with a mip assigned during PreAlignmentProcessing.py
for mip in design._mips:
  onT, offT, unM, mMQ0, onTMQ10, onTMQ60 = 0,0,0,0,0,0
  chr, start, end = mip._onTargetBed
  for read in mip._AlignedSegments:
    if read.reference_id == -1:
        unM+=1
        continue
    if read.mapping_quality == 0: mMQ0+=1
    if samfile.getrname(read.reference_id) == chr and read.reference_start >= start and read.reference_end <= end:
      onT+=1
      if read.mapping_quality >= 10:
        onTMQ10+=1
        if read.mapping_quality == 60: onTMQ60+=1
    else: offT+=1
  mip.setOnTargetDepth(onT)
  mip.setOffTargetDepth(offT)
  mip.setUnmappedDepth(unM)
  mip.setMappedMQ0Count(mMQ0)
  mip.setOnTargetMQ10Count(onTMQ10)
  mip.setOnTargetMQ60Count(onTMQ60)

#generate global on/off/unmapped stats
# output file name: [sample bamfile].on_off_unmapped
# output file format: col1=total_ontarget; col2=total_offtaget; col3=total_unmapped
# NOTE: these stats only apply to reads with a mip assigned during PreAlignmentProcessing.py
ON, OFF, UN = 0,0,0
gen = (mip._on_target_depth for mip in design._mips)
for i in gen: ON+=i
gen = (mip._off_target_depth for mip in design._mips)
for i in gen: OFF+=i
gen = (mip._unmapped_depth for mip in design._mips)
for i in gen: UN+=i

#write stats to output files: 
with open(bam + ".mip-by-mip_counts", 'w') as byMipStats:
  for mip in design._mips:
    byMipStats.write("\t".join(map(str, [mip._name, mip._on_target_depth, mip._off_target_depth, mip._unmapped_depth, mip._mapped_MQ0_count, mip._on_target_MQ10_depth, mip._on_target_MQ60_depth])) + "\n")
with open(bam + ".on_off_unmapped", 'w') as bySampleStats:
  bySampleStats.write("\t".join(map(str, [ON,OFF,UN])) + "\n")
