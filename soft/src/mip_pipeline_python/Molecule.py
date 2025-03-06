import Mip
import ReadPair
from MipPipelineHelperFunctions import reverseComplement
import subprocess
from pysam import *
import cStringIO
import itertools
from collections import Counter
from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO
from Bio import AlignIO
from tempfile import NamedTemporaryFile
from printStruct import printStruct
import sys

PS=printStruct

class Molecule:

  # class variables and methods

  Instances = {}

  @classmethod
  def InstanceCount(cls):
    return len(cls.Instances)
  @classmethod
  def GetInstances(cls):
    return (value for value in cls.Instances.values())
  @classmethod
  def Get(cls, target):
    """Return the instance whose GID is target"""
    return cls.Instances.get(target, None)

  # fundamental methods

  def __init__(self, STR_tagSequence, IOC_mip):
    self._tagSequence = STR_tagSequence
    self._readPairs = []
    self._mip = IOC_mip
    Molecule.Instances[self._mip.getId() + '_' + self._tagSequence] = self

  def __repr__(self):
    return self._mip.getId() + '_' + self._tagSequence

  # Setters

  def addReadPair(self, IOC_readPair):
    self._readPairs.append(IOC_readPair)

  # Getters

  def getNumberReads(self):
    return len (self._readPairs)


  # Doers

  def createConsensusRead(self):
    # first create fasta from all reads associated with this molecule
    STR_fasta = ''
    LIS_reads = []
    LIS_fasta = []
    for IOC_readPair in self._readPairs: LIS_fasta.append(IOC_readPair.printMafftInput())
    STR_fasta = '\n'.join(LIS_fasta)
    # run custom mafft (externally via subprocess) to create multiple sequence alignment
    mafft_exe="/RQexec/dionnela/soft/bin/one_to_rule_them_all"
    mafft_args = " _ -u 0.0 -g -0.10 -f -2.00 -h 0.1 -A _ -+ 16 -W 0.00001 -V -1.53 -s 0.0 -C 0.0 -b 62  -Q 100.0 -F -B -l 2.7 -X 0.1 -i "
    mafft_in = NamedTemporaryFile(dir="/dev/shm")
    mafft_cmd = mafft_exe + mafft_args + mafft_in.name
    with open(mafft_in.name, 'wb') as f: f.write(STR_fasta)
    P = subprocess.Popen(mafft_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = P.communicate()
    align = AlignIO.read(StringIO(stdout), "fasta")
    # create consensus sequence, assign phred score and output a Read object
    # phred score: based on fraction of reads with most common base (scaled using an exponential function that passes through the x-y points (0.5,3), (0.9,13) and (1,41)
    # rationale: normal min/max phred values in fastq are 3/41; also GATK has hidden filter to ignore bases with phred < 13. So we designed an exponential curve that passed through these 3 points
    # great resource: https://www.symbolab.com/solver/solve-for-equation-calculator ... couldn't have solved the system of equations without it
    # its formula: general y=ac^x + k; solving the 3 variables (a, c, k) for the 3 desired points: y = 6.28x10^-5 x 605673.9624^x + 2.95111
    STR_seq=""; STR_qual=""
    for i in xrange(0,align.get_alignment_length()-1):
      c = Counter(list(align[:,i]))
      c_mc_nogaps = filter(lambda x:x[0]!="-", c.most_common())
      mc = c_mc_nogaps[0]
      mc_count = mc[1]
      tot = sum(c.values()) - c['-']
      base = mc[0]
      freq = float(mc_count/float(tot))
      phred = int(6.28*10**-5 * 605673.9624**freq + 2.95111) + 33
      STR_seq = STR_seq + base.upper()
      STR_qual = STR_qual + chr(phred)
    STR_moleculeName = self._mip.getId() + '_' + self._tagSequence
    IOC_read = ReadPair.Read(STR_moleculeName, STR_seq, STR_qual)
    return IOC_read
