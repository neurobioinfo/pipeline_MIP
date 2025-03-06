#!/usr/bin/python

########################################################################################
#
# Todo list
#
# DONE: Mask out bases from raw reads that are under a certain value (13?)
# Review the gapped multi-alignment mecanism
# DONE: Needs a minimum of reads to generate a consensus
# Use varscan somatic and MuTect2
# Test different parameters (kmer size, phread score assignment from multi-align, )
# DONE: Do variant calling using -mmq 13 and  -A TandemRepeatAnnotator
# DONE: Add option to define the output prefix
# DONE: chenge ambigous to redundant
# DONE: output the number of molecules in report
# DONE: Add input for mip arm length assignment 
# DONE: Change the ambigous mip assignments
#
########################################################################################

from __future__ import print_function

import Design
import ReadPair
import Mip
import Molecule
from MipPipelineHelperFunctions import reverseComplement

import os
import sys
import subprocess
import gzip
import pysam
import printStruct
#from itertools import izip
import argparse

PS = printStruct.printStruct

####################
# HELPER FUNCTIONS #
####################
def processLine(line):
  out = line.rstrip().split(' ')[0]
  return out

def get_matching_mips(IOC_readPair, DIC_lookup):
  #NOTE: read1 sequence must be matched against lig arm data, and read2 against ext arm data
  [IOC_read1, IOC_read2] = IOC_readPair.getReads()
  key1 = IOC_read1.getSeq()[0:mip_key_length]
  key2 = IOC_read2.getSeq()[0:mip_key_length]
  SET_matches1 = set(DIC_lookup['lig'][key1]) if key1 in DIC_lookup['lig'] else set()
  SET_matches2 = set(DIC_lookup['ext'][key2]) if key2 in DIC_lookup['ext'] else set()
  # true hits are any mips which appear in the intersection of read1-lig / read2-ext sequence hits
  LIS_matches = list(SET_matches1.intersection(SET_matches2))

  if len(LIS_matches)>1: #we have an ambigous match, trying to resolve it
    LIS_matches = []
    for IOC_mip in LIS_matches:
      STR_ligArmSeq = IOC_mip.getArm('lig').getSeq()
      STR_extArmSeq = IOC_mip.getArm('ext').getSeq()
      STR_read1LigSeq = reverseComplement(IOC_read1.getSeq()[0:len(STR_ligArmSeq)])
      STR_read2ExtSeq = IOC_read2.getSeq()[0:len(STR_extArmSeq)]
      if ( STR_ligArmSeq == STR_read1LigSeq and STR_extArmSeq == STR_read2ExtSeq ):
        LIS_matches.append(IOC_mip)
    
  return LIS_matches

def printReadWithoutMipFastq(IOC_read1, IOC_read2, FH1, FH2):
  FH1.write(IOC_read1.printFastq() + '\n')
  FH2.write(IOC_read2.printFastq() + '\n')

###############
# MAIN SCRIPT #
###############
#parse script arguments
parser = argparse.ArgumentParser(description='Preprocessing of raw MIP fastq files: adds MIP names to reads')
parser.add_argument('-r1', '--read1_fastq_file', dest='read1File', help='/path/to/read1_fastq.gz_file', required=True)
parser.add_argument('-r2', '--read2_fastq_file', dest='read2File', help='/path/to/read2_fastq.gz_file', required=True)
parser.add_argument('-d', '--design_file', dest='designFile', help='/path/to/mip_design_file', required=True)
parser.add_argument('-s', '--design_source', dest='designSource', help='mip design source (preferably \"generic\", or else \"mipgen\" or \"breakpoint_resolution_wgs_mips\")', required=True)
parser.add_argument('-t', '--tagged_read_number', dest='tagged_read_number', help='read pair member containing molecular tags (integer)', default=2)
parser.add_argument('-l', '--molecular_tag_length', dest='molecular_tag_length', help='bp length of molecular tags (integer)', required=True)
parser.add_argument('-k', '--mip_key_length', dest='mip_key_length', help='first N bases of arm sequence to use for lookup (integer)', default=6)
parser.add_argument('-m', '--min_reads_per_molecule', dest='min_reads_per_mol', help='define the minimum of reads to consider a molecule (integer)', default=6)
parser.add_argument('-q', '--min_qual_per_base', dest='min_qual_per_base', help='define the minimal quality value of the raw base to be considered in the molecule (integer)', default=13)
parser.add_argument('-n', '--threads', dest='threads', help='define the number of threads to use (integer)', default=16)
parser.add_argument('-o', '--output', dest='output_prefix', help='specify the prefix for all output files (string)', required=True)
args   = parser.parse_args()
(read1File, read2File, designFile, designSource, molecular_tag_length, tagged_read_number, mip_key_length, min_reads_per_mol, min_qual_per_base, threads, output_prefix) = (args.read1File, args.read2File, args.designFile, args.designSource, int(args.molecular_tag_length), int(args.tagged_read_number), int(args.mip_key_length), int(args.min_reads_per_mol), int(args.min_qual_per_base), int(args.threads), args.output_prefix)

#check that all input variables are good: files must exist, integers must be integers and within allowed parameters
for file in (read1File, read2File, designFile):
    try:
        s = os.stat(file)
        if s.st_size == 0:
            print ("The input file {} is empty" + format(file))
            sys.exit(42)
    except OSError as e:
        print ("The input file {} does not exist" + format(file))
        sys.exit(42)
if isinstance(molecular_tag_length, int) == False: raise TypeError("--molecular_tag_length must be an integer")
if isinstance(tagged_read_number, int) == False:   raise TypeError("--tagged_read_number must be an integer")
if molecular_tag_length <= 0: raise ValueError("--molecular_tag_length must be greater than zero")
if tagged_read_number not in [1,2]: raise ValueError("--tagged_read_number must be either 1 or 2")

#print "script args: " + str({'--read1_fastq_file': read1File, '--read2_fastq_file': read2File, '--design_file': designFile, '--design_source': designSource, '--molecular_tag_length': str(molecular_tag_length), '--tagged_read_number': str(tagged_read_number), '--mip_key_length': str(mip_key_length), '--threads': str(threads)})
print ("Arguments:")
print ("--read1_fastq_file " + read1File)
print ("--read2_fastq_file " + read2File)
print ("--design_file " + designFile)
print ("--design_source " + designSource)
print ("--molecular_tag_length " + str(molecular_tag_length))
print ("--tagged_read_number " + str(tagged_read_number))
print ("--mip_key_length " + str(mip_key_length))
print ("--min_reads_per_molecule " + str(min_reads_per_mol))
print ("--min_qual_per_base " + str(min_qual_per_base))
print ("--threads " + str(threads))
print ("--output " + output_prefix)
print


#create/populate global variables
tmp_outputFile = output_prefix + '.molecular_consensus.fastq'
outputFile = output_prefix + '.molecular_consensus.fastq.gz'
scrapFile1 = output_prefix + ".withoutMips.read1.fastq.gz"
scrapFile2 = output_prefix + ".withoutMips.read2.fastq.gz"
redundantFile1 = output_prefix + ".redundantMips.read1.fastq.gz"
redundantFile2 = output_prefix + ".redundantMips.read2.fastq.gz"
design = Design.Design(designFile, designSource)
design.parse()
passed_reads_counter = 0
total_reads_counter = 0
inputFH1 = gzip.open(read1File)
inputFH2 = gzip.open(read2File)
for file in [tmp_outputFile, outputFile, scrapFile1, scrapFile2, redundantFile1, redundantFile2]:
  if os.path.isfile(file): os.remove(file)
try:
    outputFH = open(tmp_outputFile, 'wb')
except IOError as ioex:
    print ("The output  file {} cannot be created" + format(tmp_outputFile))
    sys.exit(42)
scrapFH1 = gzip.open(scrapFile1, 'wb')
scrapFH2 = gzip.open(scrapFile2, 'wb')
redundantFH1 = gzip.open(redundantFile1, 'wb')
redundantFH2 = gzip.open(redundantFile2, 'wb')
DIC_lookup = {'ext':{}, 'lig':{}}				# a dict to quickly identify mips by the first 6 bases of their arm sequences
redundant_mips={}					# a dict to store info about redundant mips; to be used in later reporting

#speed this process up with a dictionary linking the first N bases (default 6) of each mip arm (ext + reverseComplement(lig)) to the mip
for mip in design._mips:
  extKey = mip.getArm('ext').getSeq()[0:mip_key_length]
  ligKey = reverseComplement(mip.getArm('lig').getSeq())[0:mip_key_length]
  if extKey not in DIC_lookup['ext']: DIC_lookup['ext'][extKey]=[]
  if ligKey not in DIC_lookup['lig']: DIC_lookup['lig'][ligKey]=[]
  DIC_lookup['ext'][extKey].append(mip)
  DIC_lookup['lig'][ligKey].append(mip)

print ("Parsing reads")
#read both read files simultaneously; link readPairs, Mip and Molecule objects if possible, print to designated outputs, track stats for later reporting
line_num1=0
for ((line_num1,line_content1),(line_num2,line_content2)) in zip(enumerate(inputFH1), enumerate(inputFH2)):
    if line_num1%4==0: content = {0:('',''),1:('',''),2:('',''),3:('','')}
    #verbosity
    if line_num1%4000==0 and line_num1!=0:
        sys.stdout.write('.')
        sys.stdout.flush()
        if line_num1%400000==0:
          sys.stdout.write( "  [" + str(line_num1/4) + ']\n')

    content[line_num1%4] = (line_content1.rstrip(), line_content2.rstrip())
    if line_num1%4==3:
        total_reads_counter+=1
        # process last 4 lines of data into Read and ReadPair objects
        [names, seqs, strands, quals] = content.values()
        IOC_read1 = ReadPair.Read(names[0].replace('@',''), seqs[0], quals[0])
        IOC_read2 = ReadPair.Read(names[1].replace('@',''), seqs[1], quals[1])
        IOC_readPair = ReadPair.ReadPair(IOC_read1, IOC_read2)
        # remove molecular tag sequence from appropriate read
        STR_tag = IOC_readPair.updateMolecularTaggedReads(molecular_tag_length, tagged_read_number)
        # identify the MIP sequenced by the ReadPair
        LIS_matchingMips = get_matching_mips(IOC_readPair, DIC_lookup)
        # for ReadPairs with no matching mips: print to scrapFiles
        if len(LIS_matchingMips) == 0: printReadWithoutMipFastq(IOC_read1,IOC_read2,scrapFH1,scrapFH2)
        # for ReadPairs with 1 matching mip: create Molecule object that associates Molecule, ReadPair and MIP, remove MIP sequence from reads; and store for later processing and printing
        elif len(LIS_matchingMips) == 1:
            passed_reads_counter +=1
            IOC_mip = LIS_matchingMips[0]
            IOC_readPair.trimMipArmsFromReads(IOC_mip)
            IOC_readPair.maskReadBasesFromQual(min_qual_per_base)
            STR_molID = IOC_mip.getId() + '_' + STR_tag
            IOC_molecule = Molecule.Molecule.Get(STR_molID) if Molecule.Molecule.Get(STR_molID) else Molecule.Molecule(STR_tag, IOC_mip)
            IOC_molecule.addReadPair(IOC_readPair)
        # for ReadPairs with >1 matching mips: print to redundantFiles, track offenders
        elif len(LIS_matchingMips)>1:
            # pass redundant ReadPairs to designated output file
            redundantFH1.write(IOC_read1.printFastq())
            redundantFH2.write(IOC_read2.printFastq())
            # track counter of redundant mip groups for report
            mips_string = str(sorted(LIS_matchingMips))
            if mips_string not in redundant_mips: redundant_mips[mips_string]=0
            redundant_mips[mips_string]+=1

#verbosity
spaces=((100000-((line_num1/4)%100000))//1000)+1
for x in range(0, spaces):
    sys.stdout.write( ' ')
print ( "  [" + str(line_num1/4) + "]\n")

#iterate over all molecules, create consensus reads, and print to outfile

#the interative way, without threads
#for IOC_molecule in Molecule.Molecule.GetInstances():
#    IOC_consensusRead = IOC_molecule.createConsensusRead()
#    outputFH.write(IOC_consensusRead.printFastq() + '\n')

total_molecules=sum(1 for _ in Molecule.Molecule.GetInstances())
print ("Finding consensus for " + str(total_molecules) + " molecules")

#threaded way:
import pymp
import inspect
print("DEBUG: using pymp module found here: " + inspect.getfile(pymp))

#molecules = [i for i in Molecule.Molecule.GetInstances()]
used_molecules_array = pymp.shared.array((1,), dtype='uint32')
min_read_ok=True
#print molecules
with pymp.Parallel(threads) as p:
    for IOC_molecule in p.iterate(Molecule.Molecule.GetInstances()):
        if IOC_molecule.getNumberReads() < min_reads_per_mol: 
          min_read_ok=False
          continue
        IOC_consensusRead = IOC_molecule.createConsensusRead()
        with p.lock:
          used_molecules_array[0] += 1
          #verbosity
          if used_molecules_array[0]%10==0:
            if min_read_ok: p.print('.',end="") 
            else: p.print('x',end="")
            min_read_ok=True
            if used_molecules_array[0] % 1000==0:
              val = "  [" + str(used_molecules_array[0]) + "]"
              p.print (val)
          #append consensus read to output file
          p.print (IOC_consensusRead.printFastq(), file=outputFH)
          outputFH.flush()

#verbosity
spaces=((1000-((used_molecules_array[0])%1000))//10)+1
for x in range(0, spaces):
    sys.stdout.write( ' ')
print ( "  [" + str(used_molecules_array[0]) + "]\n")
valid_molecules = str(used_molecules_array[0])

# be a responsible programmer and close file handles that are no longer needed
inputFH1.close(); inputFH2.close()
outputFH.close()
scrapFH1.close(); scrapFH2.close()
redundantFH1.close(); redundantFH2.close()

#compress the output file
import shutil
with open(tmp_outputFile, 'rb') as f_in, gzip.open(outputFile , 'wb') as f_out:
    shutil.copyfileobj(f_in, f_out)
os.remove(tmp_outputFile)

# report redundant mips (i.e. that match indistinguishably for some read pairs)
redundantMip_report_filename = output_prefix + "_redundant_mips"	# have to pass a filename and not a file handle to printStruct function
PS(redundant_mips, file=redundantMip_report_filename)			# need to run this function twice for some bizzare reason; otherwise output file is empty
PS(redundant_mips, file=redundantMip_report_filename)

#print to file a summary of how many reads were included (and by implication, how many excluded) for further analysis
readPair_report = open(output_prefix + ".preprocessing_report.txt", 'wb')
readPair_report.write("assigned a mip to " + str(passed_reads_counter) + " out of " + str(total_reads_counter) + " total read pairs" + "\n")
readPair_report.write("usable molecules is " + valid_molecules + " out of a total of " + str(total_molecules) + "\n")
readPair_report.close()
