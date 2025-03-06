import Mip
from MipPipelineHelperFunctions import reverseComplement

class Read:

# Fundamental Methods

  def __init__(self, id, seq, qual):
    self._id	= id
    self._seq	= seq
    self._qual	= qual
  def __repr__(self):
    return self._id

# Getters

  def getSeq(self):
    return self._seq
  def getQual(self):
    return self._qual

# Setters

  def updateSeq(self, STR_seq):
    self._seq = STR_seq
  def updateQual(self, STR_qual):
    self._qual = STR_qual

# Doers

  def printFastq(self):
    return '\n'.join(['@' + self._id, self.getSeq(), '+', self.getQual()])
  def printFasta(self):
    return '\n'.join(['>' + self._id, self.getSeq()])
  def printFastaRevComp(self):
    return '\n'.join(['>' + self._id, reverseComplement(self.getSeq())])

class ReadPair:

# Fundamental Methods

  def __init__(self, read1, read2):
    self._read1	= read1
    self._read2	= read2
  def __repr__(self):
    return str([self._read1, self._read2])

# Getters

  def getRead1(self):
    return self._read1
  def getRead2(self):
    return self._read2
  def getReads(self):
    return [self._read1, self._read2]

# Doers

  def printFasta(self):
    return '\n'.join([IOC_read.printFasta() for IOC_read in self.getReads()])

  def printMafftInput(self):
    STR_fastaR1 = self._read1.printFasta()
    STR_fastaR2 = self._read2.printFastaRevComp()
    return '\n'.join([STR_fastaR1, STR_fastaR2])

  def updateMolecularTaggedReads(self, INT_tagLength, INT_taggedReadNumber):
    """ parses (and returns) molecular tag sequence from a read pair, and removes the tag sequence from the appropriate read """
    IOC_read = self.getRead2() if INT_taggedReadNumber == 2 else self.getRead1()
    STR_tag = IOC_read.getSeq()[0:INT_tagLength]
    STR_newSeq = IOC_read.getSeq()[INT_tagLength:]
    STR_newQual = IOC_read.getQual()[INT_tagLength:]
    IOC_read.updateSeq(STR_newSeq)
    IOC_read.updateQual(STR_newQual)
    return STR_tag

  def trimMipArmsFromReads(self, IOC_mip):
    [e, l] = map(lambda x:len(x.getSeq()), IOC_mip.getArms())
    (r1start, r1end) = (l, l+112)
    (r2start, r2end) = (e, e+112)
    self._read1.updateSeq(self._read1.getSeq()[r1start:r1end])
    self._read1.updateQual(self._read1.getQual()[r1start:r1end])
    self._read2.updateSeq(self._read2.getSeq()[r2start:r2end])
    self._read2.updateQual(self._read2.getQual()[r2start:r2end])

  def maskReadBasesFromQual(self, INT_minQual):
    index=0
    seq1 = self._read1.getSeq()
    qual1 = self._read1.getQual()
    for qual_value in list(qual1):
      if ( ord(qual_value) < INT_minQual+33 ):
        seq1 = seq1[:index] + 'N' + seq1[index + 1:]
      index+=1
    self._read1.updateSeq(seq1)

    index=0
    seq2 = self._read2.getSeq()
    qual2 = self._read2.getQual()
    for qual_value in list(qual2):
      if ( ord(qual_value) < INT_minQual+33 ):
        seq2 = seq2[:index] + 'N' + seq2[index + 1:]
      index+=1
    self._read2.updateSeq(seq2)
