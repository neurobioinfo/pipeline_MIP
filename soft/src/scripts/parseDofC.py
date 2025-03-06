#!/usr/bin/env python

# import modules
import sys
import subprocess
import os
import getopt
import argparse
PIPELINE_HOME=os.environ.get('PIPELINE_HOME')
try: sys.path.append('/'.join([PIPELINE_HOME, 'soft/src/mip_pipeline_python']))
except: raise RuntimeError('\n\n\nPIPELINE_HOME environment variable does not exist. Did you forget to set it?\n\n\n')
from printStruct import printStruct as ps

# process inputs
parser = argparse.ArgumentParser(description='generate somatic calls from GATK DepthOfCoverage data')
parser.add_argument('-i', '--input_file', dest='inputData', help='either /path/to/DepthOfCoverage_file or else \"-\" (if piped)', required=True)
parser.add_argument('-b', '--bed_file', dest='bed', help='/path/to/bed_file', required=True)
parser.add_argument('-r', '--min_supporting_reads', dest='min_supporting_reads', help='minimum # of reads per molecule, supporting variant', default=6)
parser.add_argument('-f', '--perfect_fraction', dest='perfectFraction', help='fraction of reads (per molecule) above which the call is considered \"perfect\"', default=0.9)
args   = parser.parse_args()
inputData, bed, min_supporting_reads, perfectFraction = (args.inputData, args.bed, int(args.min_supporting_reads), float(args.perfectFraction))
inputFile = sys.stdin if inputData == '-' else inputData
sys.stderr.write('input args: ' + str({'--input_file':inputData, '--bed_file':bed, '--min_supporting_reads':str(min_supporting_reads), '--perfect_fraction':str(perfectFraction)}) + '\n')

# error-check inputs
if inputFile != sys.stdin and (os.path.isfile(inputFile) == False or os.path.getsize(inputFile) == 0): raise ValueError('input file ' + inputFile + 'empty or does not exist')
if os.path.isfile(bed) == False or os.path.getsize(bed) == 0: raise ValueError('bed file ' + bed + 'empty or does not exist')
if isinstance(perfectFraction, int) and (perfectFraction < 0 or perfectFraction > 1): raise ValueError('--perfect_fraction must be between 0 and 1')

# define helper functions
def isPerfectCall(sampleDict):
    return True if max(sampleDict[i] for i in ACGT)/sampleDict['depth'] >= perfectFraction else False

def getAltAllelesCov(sampleDict, referenceAllele):
    return sum([sampleDict[i] for i in ACGT if i != referenceAllele])

def passesThresholds(sampleDict, threshold=None):
    if threshold is None: threshold = min_supporting_reads
    return True if sampleDict['depth'] >= threshold else False

# get bases for all bed positions
ref='/'.join([PIPELINE_HOME, 'data/reference/human_g1k_v37.fasta'])
lbed , dbed, pos2seq = ([],{},{})
p1 = subprocess.Popen(['bedtools', 'merge', '-i', bed], stdout=subprocess.PIPE)
p2 = subprocess.Popen(['bedtools', 'getfasta', '-fi', ref, '-bed', '-', '-fo', os.path.basename(bed) + '.fa'], stdin=p1.stdout, stdout=subprocess.PIPE)
p2.communicate()

with open(os.path.basename(bed) + '.fa') as f:
    for line in f.read().splitlines():
        lbed.append(line.replace('>', '').replace('-',':'))

dbed = dict(zip(lbed[0::2], lbed[1::2]))
for pos in dbed:
    chr, start, end = pos.split(":"); start, end = map(int, [start,end])
    start=start+1
    seq = dbed[pos]
    for i in xrange(0, len(seq)):
        pos2seq[chr + ":" + str(start + i)] = list(seq)[i]

# parse + process DofC data one variant at a time
ACGT = ['A','C','G','T']
outputNames = [	'var',
		'ref_allele',
		'tags_all',
		'cov_all',
		'tags_all_perfect',
		'cov_all_perfect',
		'tags_all_perfect_ref',
		'cov_all_perfect_ref',
		'tags_all_perfect_alt',
		'cov_all_perfect_alt',
		'tags_all_imperfect',
		'cov_all_imperfect',
		'tags_all_imperfect_ref',
		'cov_all_imperfect_ref',
		'tags_all_imperfect_alt',
		'cov_all_imperfect_alt',
		'tags_filtered',
		'cov_filtered',
		'tags_filtered_perfect',
		'cov_filtered_perfect',
		'tags_filtered_perfect_ref',
		'cov_filtered_perfect_ref',
		'tags_filtered_perfect_alt',
		'cov_filtered_perfect_alt',
		'tags_filtered_imperfect',
		'cov_filtered_imperfect',
		'tags_filtered_imperfect_ref',
		'cov_filtered_imperfect_ref',
		'tags_filtered_imperfect_alt',
		'cov_filtered_imperfect_alt',
		'alt_alleles_all',
		'alt_alleles_all_perfect',
		'alt_alleles_all_imperfect',
		'alt_alleles_filtered',
		'alt_alleles_filtered_perfect',
		'alt_alleles_filtered_imperfect']

inputFH = sys.stdin if inputFile == sys.stdin else open(inputFile, 'r')
for i, line in enumerate(inputFH):
    # header line: associate column numbers with sample names; print header
    lineData=line.split("\t")
    if not lineData: continue
    if i==0:
       col2sample = {x:lineData[x][10:] for x in xrange(3,len(lineData),2)}
       print "\t".join(outputNames)
       continue
    # variant lines
    if lineData[0][0] == "-" or lineData[0][0:4] == "Done": continue
    var, total_depth = lineData[0:2]
    ref_allele = pos2seq[var]
    sampleData={}
    for col in xrange(3,len(lineData),2):
        try: sample = col2sample[col]
	except:
		sys.stderr.write('parsing problem at position ' + str(var) + ' for sample ' + str(sample) + ' (probably a truncated/fused input line; it happens sometimes with huge DepthOfCoverage runs. skipping the rest of the line...\n')
		break
        try: depth = int(lineData[col])
	except:
		sys.stderr.write('parsing problem at position ' + str(var) + ' for sample ' + str(sample) + ' (probably a truncated/fused input line; it happens sometimes with huge DepthOfCoverage runs. skipping the rest of the line...\n')
		break
        if depth == 0: continue
        sampleData[sample] = {key:0 for key in ['depth', 'A', 'C', 'G', 'T', 'N']}
        sampleData[sample]['depth'] = depth
        try: 
            baseData = dict([k, int(v)] for k,v in dict([val.split(':') for val in lineData[col+1].strip().split()]).iteritems())
        except ValueError:
            print 'script shits a brick at line ' + str(i) + ' in file ' + str(inputFile)
            print 'dump of line data: ' + str(ps(lineData))
            break
        sampleData[sample].update(baseData)
    # apply various filters to variant data as appropriate and print output
    d					= sampleData		# a little shorthand to save some typing
    tags_all				= len([ sample	for sample in d																	])
    tags_all_perfect			= len([ sample	for sample in d	if											        isPerfectCall(d[sample])	])
    tags_all_perfect_ref		= len([ sample	for sample in d	if					d[sample][ref_allele]			!=0	and     isPerfectCall(d[sample])	])
    tags_all_perfect_alt		= len([ sample	for sample in d	if					getAltAllelesCov(d[sample], ref_allele)	!=0	and     isPerfectCall(d[sample])	])
    tags_all_imperfect			= len([ sample	for sample in d	if											    not isPerfectCall(d[sample])	])
    tags_all_imperfect_ref		= len([ sample	for sample in d	if					d[sample][ref_allele]			!=0	and not isPerfectCall(d[sample])	])
    tags_all_imperfect_alt		= len([ sample	for sample in d	if					getAltAllelesCov(d[sample], ref_allele)	!=0	and not isPerfectCall(d[sample])	])
    tags_filtered			= len([ sample	for sample in d	if passesThresholds(d[sample])													])
    tags_filtered_perfect		= len([ sample	for sample in d	if passesThresholds(d[sample])								and     isPerfectCall(d[sample])	])
    tags_filtered_perfect_ref		= len([ sample	for sample in d	if passesThresholds(d[sample]) and	d[sample][ref_allele]			!=0	and     isPerfectCall(d[sample])	])
    tags_filtered_perfect_alt		= len([ sample	for sample in d	if passesThresholds(d[sample]) and	getAltAllelesCov(d[sample], ref_allele)	!=0	and     isPerfectCall(d[sample])	])
    tags_filtered_imperfect		= len([ sample	for sample in d	if passesThresholds(d[sample]) 								and not isPerfectCall(d[sample])	])
    tags_filtered_imperfect_ref		= len([ sample	for sample in d	if passesThresholds(d[sample]) and	d[sample][ref_allele]			!=0	and not isPerfectCall(d[sample])	])
    tags_filtered_imperfect_alt		= len([ sample	for sample in d	if passesThresholds(d[sample]) and	getAltAllelesCov(d[sample], ref_allele)	!=0	and not isPerfectCall(d[sample])	])
    cov_all				= sum([ d[sample]['depth']			for sample in d										])
    cov_all_perfect			= sum([ d[sample]['depth']			for sample in d	if				            isPerfectCall(d[sample])	])
    cov_all_perfect_ref			= sum([ d[sample][ref_allele]			for sample in d	if				            isPerfectCall(d[sample])	])
    cov_all_perfect_alt			= sum([ getAltAllelesCov(d[sample], ref_allele)	for sample in d	if				            isPerfectCall(d[sample])	])
    cov_all_imperfect			= sum([ d[sample]['depth']			for sample in d	if				        not isPerfectCall(d[sample])	])
    cov_all_imperfect_ref		= sum([ d[sample][ref_allele]			for sample in d	if				        not isPerfectCall(d[sample])	])
    cov_all_imperfect_alt		= sum([ getAltAllelesCov(d[sample], ref_allele)	for sample in d	if				        not isPerfectCall(d[sample])	])
    cov_filtered			= sum([ d[sample]['depth']			for sample in d	if	passesThresholds(d[sample])					])
    cov_filtered_perfect		= sum([ d[sample]['depth']			for sample in d	if	passesThresholds(d[sample]) and     isPerfectCall(d[sample])	])
    cov_filtered_perfect_ref		= sum([ d[sample][ref_allele]			for sample in d	if	passesThresholds(d[sample]) and     isPerfectCall(d[sample])	])
    cov_filtered_perfect_alt		= sum([ getAltAllelesCov(d[sample], ref_allele)	for sample in d	if	passesThresholds(d[sample]) and     isPerfectCall(d[sample])	])
    cov_filtered_imperfect		= sum([ sum([d[sample][i] for i in ACGT])	for sample in d	if	passesThresholds(d[sample]) and not isPerfectCall(d[sample])	])
    cov_filtered_imperfect_ref		= sum([ d[sample][ref_allele]			for sample in d	if	passesThresholds(d[sample]) and not isPerfectCall(d[sample])	])
    cov_filtered_imperfect_alt		= sum([ getAltAllelesCov(d[sample], ref_allele)	for sample in d	if	passesThresholds(d[sample]) and not isPerfectCall(d[sample])	])
    alt_alleles_all			= ':'.join(list(set([i for i in ACGT for sample in d if d[sample][i]!=0 and i!=ref_allele									])))
    alt_alleles_all_perfect		= ':'.join(list(set([i for i in ACGT for sample in d if d[sample][i]!=0 and i!=ref_allele                                 and     isPerfectCall(d[sample])	])))
    alt_alleles_all_imperfect		= ':'.join(list(set([i for i in ACGT for sample in d if d[sample][i]!=0 and i!=ref_allele                                 and not isPerfectCall(d[sample])	])))
    alt_alleles_filtered		= ':'.join(list(set([i for i in ACGT for sample in d if d[sample][i]!=0 and i!=ref_allele and passesThresholds(d[sample])					])))
    alt_alleles_filtered_perfect	= ':'.join(list(set([i for i in ACGT for sample in d if d[sample][i]!=0 and i!=ref_allele and passesThresholds(d[sample]) and     isPerfectCall(d[sample])	])))
    alt_alleles_filtered_imperfect	= ':'.join(list(set([i for i in ACGT for sample in d if d[sample][i]!=0 and i!=ref_allele and passesThresholds(d[sample]) and not isPerfectCall(d[sample])	])))
    outputValues = [var,
    		ref_allele,
    		tags_all,
    		cov_all,
    		tags_all_perfect,
    		cov_all_perfect,
    		tags_all_perfect_ref,
    		cov_all_perfect_ref,
    		tags_all_perfect_alt,
    		cov_all_perfect_alt,
    		tags_all_imperfect,
    		cov_all_imperfect,
    		tags_all_imperfect_ref,
    		cov_all_imperfect_ref,
    		tags_all_imperfect_alt,
    		cov_all_imperfect_alt,
    		tags_filtered,
    		cov_filtered,
    		tags_filtered_perfect,
    		cov_filtered_perfect,
    		tags_filtered_perfect_ref,
    		cov_filtered_perfect_ref,
    		tags_filtered_perfect_alt,
    		cov_filtered_perfect_alt,
    		tags_filtered_imperfect,
    		cov_filtered_imperfect,
    		tags_filtered_imperfect_ref,
    		cov_filtered_imperfect_ref,
    		tags_filtered_imperfect_alt,
    		cov_filtered_imperfect_alt,
    		alt_alleles_all,
    		alt_alleles_all_perfect,
    		alt_alleles_all_imperfect,
    		alt_alleles_filtered,
    		alt_alleles_filtered_perfect,
    		alt_alleles_filtered_imperfect]
    print "\t".join(map(str, outputValues))

if inputFile != sys.stdin: inputFH.close()
