import sys
import subprocess
import gzip
import ReadPair
import Mip
#import Alignment
import HomologySearch
import swalign
from swalign import *
import pysam
import printStruct
import random
import string

def parseReadPairFilesIntoReadPairs(read1File, read2File):
  readPairs = list()
  with gzip.open(read1File) as file1:
    content1 = [line.rstrip().split(" ")[0] for line in file1.readlines()]
  with gzip.open(read2File) as file2:
    content2 = [line.rstrip().split(" ")[0] for line in file2.readlines()]
  for i in xrange(len(content1)/4):
    name = content1[4*i]
    read1 = ReadPair.Read(name, '1', content1[4*i+1], content1[4*i+3])
    read2 = ReadPair.Read(name, '2', content2[4*i+1], content2[4*i+3])
    rp = ReadPair.ReadPair(name, read1, read2)
    read1.setReadPair(rp)
    read2.setReadPair(rp)
    readPairs.append(rp)
  return readPairs

def parseMipGenDesignIntoMips_breakpoint_resolution_wgs_mips(designFile):
  mips = list()
  with open(designFile) as design:
    design_contents = [line.rstrip() for line in design.readlines()]
  for i in xrange(1, len(design_contents)):
    line = design_contents[i].split("\t")
    (name, chr) = ("_".join(line[0:2]), line[2])
    gene = line[-1] if len(line) == 21 else "dummy_gene"
    ext = Mip.Arm(name, 'ext', chr, line[3], line[4], line[5])
    lig = Mip.Arm(name, 'lig', chr, line[7], line[8], line[9])
    mip = Mip.Mip(name, ext, lig, chr, line[11], line[12], line[13], line[17], gene)
    ext.setMip(mip)
    lig.setMip(mip)
    mips.append(mip)
    mip.setOnTargetBed(mip.getOnTargetBed())
  return mips

def parseMipGenDesignIntoMips_mipgen(designFile):
  mips = list()
  with open(designFile) as design:
    design_contents = [line.rstrip() for line in design.readlines()]
  for i in xrange(1, len(design_contents)):
    line = design_contents[i].split("\t")
    (name, chr) = ("_".join(line[0:2]), line[2])
    gene = line[-1] if len(line) == 21 else "dummy_gene"
    ext = Mip.Arm(name, 'ext', chr, line[3], line[4], line[5])
    lig = Mip.Arm(name, 'lig', chr, line[7], line[8], line[9])
    mip = Mip.Mip(name, ext, lig, chr, line[11], line[12], line[13], line[17], gene)
    ext.setMip(mip)
    lig.setMip(mip)
    mips.append(mip)
    mip.setOnTargetBed(mip.getOnTargetBed())
  return mips

def basicAlign(read, arm):
  search = HomologySearch.BasicSearch(read, arm)
  search.alignForward()
  search.alignRevComp()

def pairwise(readPairAlignments):
  out = []
  a = readPairAlignments[0]
  b = readPairAlignments[1]
  if len(a)>0 and len(b)>0:
    for x in xrange(0,len(a)):
      for y in xrange(0,len(b)):
        out.append((a[x],b[y]))
    return out

def findMipInReadPair(matrix):
  mipHits = []
  if matrix is not None:
    for test in matrix:
      if test[0]._arm.getId().split("_")[0:2] == test[1]._arm.getId().split("_")[0:2] and test[0]._arm.getId().split("_")[-1] != test[1]._arm.getId().split("_")[-1]: mipHits.append(test[0]._arm.getMip())
    return mipHits

def linkReadAndMips(readPair, hits):    # 'hits' is an array of mips
  if hits is not None:
    for mip in hits:
      #mip.addReadPairToMip(readPair)
      readPair.addMipToReadPair(mip)

def reverseComplement(STR_sequence):
  DIC_trans = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
  LIS_seq = list(STR_sequence)
  LIS_comp = [ DIC_trans[i] for i in LIS_seq ]
  LIS_comp.reverse()
  STR_revComp = ''.join(LIS_comp)
  return STR_revComp
