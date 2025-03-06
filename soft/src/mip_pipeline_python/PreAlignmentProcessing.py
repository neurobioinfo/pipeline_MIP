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
import itertools
from MipPipelineHelperFunctions import *
from swalign import *
#from itertools import izip, islice
import getopt
import argparse

PS = printStruct.printStruct

####################
# HELPER FUNCTIONS #
####################
def processLine(line):
  out = line.rstrip().split(' ')[0]
  return out

def updateTaggedReads(read1, read2, tag_length, tagged_read_number):	# removes molecular tag sequence from start of read(s) with tag; adds molecular tag to read names; note: arg tagged_read_index is a string with possible values of 1 or 2
  dict         = {1:read1, 2:read2}							# dict to choose which read to modify seq, from tagged_read_number
  read         = dict[tagged_read_number]						# identifies read to parse (and modify) for tag sequence
  tag          = read['seq'][0:tag_length]						# parses tag string from read
  read['seq']  = read['seq'][tag_length:]						# removes tag bases from read sequence
  read['qual'] = read['qual'][tag_length:]						# removes tag quals from read sequence
  for read in [read1,read2]: read['name'] = read['name'].replace('@', '@' + tag + '_')	# adds molecular tag sequence to read name

def trimArmsFromReads(mip):             				# removes mip arm sequences from read
  [e, l] = map(lambda x:len(x.getSeq()), mip.getArms())
  (r1start, r1end) = (l, l+112)
  (r2start, r2end) = (e, e+112)
  (read1['seq'], read1['qual']) = (read1['seq'][r1start:r1end], read1['qual'][r1start:r1end])
  (read2['seq'], read2['qual']) = (read2['seq'][r2start:r2end], read2['qual'][r2start:r2end])

def pickRandomElementFromList(tab):
  return tab[random.randint(0,len(tab)-1)]

def opposite(arm,read):
  ARM  = 'ext' if arm=='lig' else 'lig'
  READ = 2 if read==1 else 1
  return (ARM,READ)

def get_matching_mips():
  result=[]
  for mip in [mip for mip in test if test[mip].keys()>1]:
    for arm in test[mip]:
      for read in test[mip][arm]:
        (ARM,READ) = opposite(arm,read)
        if ARM in test[mip] and read in test[mip][arm] and READ in test[mip][ARM]: result.append(mip)
# TODO: UPDATE following properly from somatic pipeline
#  if len(result)>1: #we have an ambigous match, trying to resolve it
#    result = []
#    for IOC_mip in result:
#      STR_ligArmSeq = IOC_mip.getArm('lig').getSeq()
#      STR_extArmSeq = IOC_mip.getArm('ext').getSeq()
#      STR_read1LigSeq = reverseComplement(IOC_read1.getSeq()[0:len(STR_ligArmSeq)])
#      STR_read2ExtSeq = IOC_read2.getSeq()[0:len(STR_extArmSeq)]
#      if ( STR_ligArmSeq == STR_read1LigSeq and STR_extArmSeq == STR_read2ExtSeq ):
#        result.append(IOC_mip)

  return list(set(result))

def printReadWithoutMipFastq(r1, r2, FH1, FH2):
  FH1.write("\n".join(["/".join([r1['name'], str(r1['member'])]), r1['seq'], "+", r1['qual']]) + '\n')
  FH2.write("\n".join(["/".join([r2['name'], str(r2['member'])]), r2['seq'], "+", r2['qual']]) + '\n')

def printReadWithMipFastq (mip, r1, r2, FH1, FH2):
  FH1.write("\n".join(["/".join([r1['name'].replace('@', '@' + mip.getId() + "_"), str(r1['member'])]),r1['seq'],"+",r1['qual']]) + '\n')
  FH2.write("\n".join(["/".join([r2['name'].replace('@', '@' + mip.getId() + "_"), str(r2['member'])]),r2['seq'],"+",r2['qual']]) + '\n')

###############
# MAIN SCRIPT #
###############
##process inputs - test data
#read1File, read2File, designFile, designSource = ['/alloc/grouleau/COMMUN/runs/dan/MIP/rawfiles/ALS_Somatic/ALSSOM_MIP36/S15641/Illumina_HiSeq4000_Paired-Macrogen-MIP-2016_08_25-Macrogen_Run_1607KHF_0048/RAW/L1_S15641_MIP36_1.fastq.gz', '/alloc/grouleau/COMMUN/runs/dan/MIP/rawfiles/ALS_Somatic/ALSSOM_MIP36/S15641/Illumina_HiSeq4000_Paired-Macrogen-MIP-2016_08_25-Macrogen_Run_1607KHF_0048/RAW/L1_S15641_MIP36_2.fastq.gz', '/RQexec/spiegelm/data/pipeline_MIPs.svn/data/targets/design/ALSSOMv3.design', 'breakpoint_resolution_wgs_mips']
#(molecular_tag_length, tagged_read_number) = (12,2)

##process inputs
#read1File, read2File, designFile, designSource = sys.argv[1:5]
#if len(sys.argv[1:]) != 6:
#  print "PreAlignmentProcessing.py (designFile, read1File, read2File, molecular_tag_length, tagged_read_number) missing args; using default values of molecular_tag_length=0 and tagged_read_number=2"
#  (molecular_tag_length, tagged_read_number) = (0,2)
#else:
#  molecular_tag_length, tagged_read_number = map(int,sys.argv[5:])
#print "script args: " + str(sys.argv)

#parse script arguments
parser = argparse.ArgumentParser(description='Preprocessing of raw MIP fastq files: adds MIP names to reads')
parser.add_argument('-r1', '--read1_fastq_file', dest='read1File', help='/path/to/read1_fastq.gz_file', required=True)
parser.add_argument('-r2', '--read2_fastq_file', dest='read2File', help='/path/to/read2_fastq.gz_file', required=True)
parser.add_argument('-d', '--design_file', dest='designFile', help='/path/to/mip_design_file', required=True)
parser.add_argument('-s', '--design_source', dest='designSource', help='mip design source (preferably \"generic\", or else \"mipgen\" or \"breakpoint_resolution_wgs_mips\")', required=True)
parser.add_argument('-l', '--molecular_tag_length', dest='molecular_tag_length', help='bp length of molecular tags (integer)', default=0)
parser.add_argument('-t', '--tagged_read_number', dest='tagged_read_number', help='read pair member containing molecular tags (integer)', default=2)
parser.add_argument('-k', '--mip_key_length', dest='mip_key_length', help='first N bases of arm sequence to use for lookup (integer)', default=6)
args   = parser.parse_args()
(read1File, read2File, designFile, designSource, molecular_tag_length, tagged_read_number, mip_key_length) = (args.read1File, args.read2File, args.designFile, args.designSource, int(args.molecular_tag_length), int(args.tagged_read_number), int(args.mip_key_length))

#check that all input variables are good: files must exist, integers must be integers and within allowed parameters
for file in (read1File, read2File, designFile):
    try:
        s = os.stat(file)
        if s.st_size == 0:
            print("The input file {} is empty".format(file))
            sys.exit(42)
    except OSError as e:
        print("The input file {} does not exist".format(file))
        sys.exit(42)
if isinstance(molecular_tag_length, int) == False: raise TypeError("--molecular_tag_length must be an integer")
if isinstance(tagged_read_number, int) == False:   raise TypeError("--tagged_read_number must be an integer")
if tagged_read_number not in [1,2]: raise ValueError("--tagged_read_number must be either 1 or 2")

print("script args: " + str({'--read1_fastq_file': read1File, '--read2_fastq_file': read2File, '--design_file': designFile, '--design_source': designSource, '--molecular_tag_length': str(molecular_tag_length), '--tagged_read_number': str(tagged_read_number), '--mip_key_length': str(mip_key_length)}))

#create/populate global variables
outputFile1 = read1File + ".withMips.fastq.gz"
outputFile2 = read2File + ".withMips.fastq.gz"
scrapFile1 = read1File + ".withoutMips.fastq.gz"
scrapFile2 = read2File + ".withoutMips.fastq.gz"
design = Design.Design(designFile, designSource)
design.parse()
passed_reads_counter = 0
total_reads_counter = 0
output_prefix = read1File.replace('_read1','').replace('.fastq.gz','')
inputFH1 = gzip.open(read1File)
inputFH2 = gzip.open(read2File)
for file in [outputFile1, outputFile2, scrapFile1, scrapFile2]:
  if os.path.isfile(file): os.remove(file)
outputFH1 = gzip.open(outputFile1, 'ab')
outputFH2 = gzip.open(outputFile2, 'ab')
scrapFH1 = gzip.open(scrapFile1, 'ab')
scrapFH2 = gzip.open(scrapFile2, 'ab')
lookup = {}						# a dict to quickly identify mips by the first 6 bases of their arm sequences
redundant_mips={}					# a dict to store info about redundant mips; to be used in later reporting

#speed this process up with a dictionary linking the first 16 bases of each mip arm (normal + revcomp) to the mip
for mip in design._mips:
  for arm in mip.getArms():
    key=arm.getSeq()[0:mip_key_length]
    if key not in lookup: lookup[key]=[]
    lookup[key].append(arm)
    try:
        key=revcomp(arm.getSeq())[0:mip_key_length]
    except: raise RuntimeError(mip, arm.getSeq(), mip_key_length)
    if key not in lookup: lookup[key]=[]
    lookup[key].append(arm)

#read both read files simultaneously; assign a mip if possible, print to outputs, increment counter for stats, without storing previous reads in memory
for ((num1,line1),(num2,line2)) in zip(enumerate(inputFH1), enumerate(inputFH2)):
    if num1%4==0: content = {0:('',''),1:('',''),2:('',''),3:('','')}
    content[num1%4] = (processLine(line1), processLine(line2))
    if num1%4==3:
        [names, seqs, strands, quals] = content.values()
        read1 = {'name':names[0], 'member':1, 'seq':seqs[0], 'qual':quals[0], 'mips':[]}
        read2 = {'name':names[1], 'member':2, 'seq':seqs[1], 'qual':quals[1], 'mips':[]}
        total_reads_counter+=1
        if molecular_tag_length > 0: updateTaggedReads(read1, read2, molecular_tag_length, tagged_read_number)
        for read in [read1, read2]:
            key = read['seq'][0:mip_key_length]
            if key in lookup:
                for mip in lookup[key]: read['mips'].append(mip)
        test={}
        for read in [read1,read2]:
            for arm in read['mips']:
              mip=arm.getMip()
              if mip not in test: test[mip]={}
              if arm._type not in test[mip]: test[mip][arm._type]=[]
              test[mip][arm._type].append(read['member'])
              test[mip][arm._type] = list(set(test[mip][arm._type]))
        hits = get_matching_mips()
        if len(hits)==0:
            printReadWithoutMipFastq(read1,read2,scrapFH1,scrapFH2)
            continue
        MIP = pickRandomElementFromList(hits)
        trimArmsFromReads(MIP)
        if len(hits)>1:
            mips_string = str(sorted(hits))
            if mips_string not in redundant_mips: redundant_mips[mips_string]=0
            redundant_mips[mips_string]+=1
        printReadWithMipFastq(MIP, read1, read2, outputFH1, outputFH2)
        passed_reads_counter +=1

# be a responsible programmer and close file handles that are no longer needed
inputFH1.close(); inputFH2.close()
outputFH1.close(); outputFH2.close()
scrapFH1.close(); scrapFH2.close()

# report redundant mips (i.e. that match indistinguishably for some read pairs)
redundantMip_report_filename = output_prefix + "_redundant_mips"	# have to pass a filename and not a file handle to printStruct function
PS(redundant_mips, file=redundantMip_report_filename)			# need to run this function twice for some bizzare reason; otherwise output file is empty
PS(redundant_mips, file=redundantMip_report_filename)

#print to file a summary of how many reads were included (and by implication, how many excluded) for further analysis
readPair_report = open(output_prefix + ".preprocessing_report.txt", 'wb')
readPair_report.write("assigned a mip to " + str(passed_reads_counter) + " out of " + str(total_reads_counter) + " total read pairs" + "\n")
readPair_report.close()
