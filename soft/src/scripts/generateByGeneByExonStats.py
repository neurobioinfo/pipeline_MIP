#!/usr/bin/env python

import sys

INFILE = sys.argv[1]
OUTPUT_PREFIX = sys.argv[2]

GENES = {}
EXONS = {}
RESULTS = {'by_gene':{},'by_exon':{}}

with open(sys.argv[1], 'r') as f:
  for line in f.readlines():
    g_chr, g_start, g_stop, g_gene, m_chr, m_start, m_stop, mip, m_gene, type, av_cov, x20, x50, x100, overlap = line.rstrip().split('\t')
    exon = '\t'.join([g_chr, g_start, g_stop, m_gene])
    gene = m_gene
    if gene not in GENES: GENES[gene]={}
    if exon not in EXONS: EXONS[exon]={}
    for pos in xrange(max(map(int, [g_start, m_start])), min(map(int, [g_stop, m_stop]))+1):
      if pos not in GENES[gene]: GENES[gene][pos]=0
      if pos not in EXONS[exon]: EXONS[exon][pos]=0
      GENES[gene][pos]+=float(av_cov)
      EXONS[exon][pos]+=float(av_cov)

for gene in GENES:
  size = len(GENES[gene].keys())
  av_num, x20, x50, x100 = (float(), float(), float(), float())
  for pos in GENES[gene]:
    cov = GENES[gene][pos]
    av_num += cov
    if cov >=20: x20+=1
    if cov >=50: x50+=1
    if cov >=100: x100+=1
  RESULTS['by_gene'][gene] = {'gene' : gene, 'size': size, 'av_cov' : av_num/size, '%20x' : x20/size, '%50x' : x50/size, '%100x' : x100/size}

for exon in EXONS:
  size = len(EXONS[exon].keys())
  av_num, x20, x50, x100 = (float(), float(), float(), float())
  for pos in EXONS[exon]:
    cov = EXONS[exon][pos]
    av_num += cov
    if cov >=20: x20+=1
    if cov >=50: x50+=1
    if cov >=100: x100+=1
  RESULTS['by_exon'][exon] = {'exon' : exon, 'size': size, 'av_cov' : av_num/size, '%20x' : x20/size, '%50x' : x50/size, '%100x' : x100/size}

#print gene outputs to OUTPUT_PREFIX.by_gene.coverage_data
with open(OUTPUT_PREFIX + '.by_gene.coverage_data', 'w') as f:
  LIS_header = ['gene','av_cov','%20x','%50x','%100x']
  f.write('\t'.join(LIS_header) + '\n')
  for gene in RESULTS['by_gene']:
    LIS_data = [RESULTS['by_gene'][gene][x] for x in LIS_header]
    f.write('\t'.join(map(str, LIS_data)) + '\n')

#print exon outputs to OUTPUT_PREFIX.by_exon.coverage_data
with open(OUTPUT_PREFIX + '.by_exon.coverage_data', 'w') as f:
  LIS_header = ['chrom','start','end','gene','av_cov','%20x','%50x','%100x']
  f.write('\t'.join(LIS_header) + '\n')
  for exon in RESULTS['by_exon']:
    LIS_data = [RESULTS['by_exon'][exon][x] for x in ['exon','av_cov','%20x','%50x','%100x']]
    f.write('\t'.join(map(str, LIS_data)) + '\n')
