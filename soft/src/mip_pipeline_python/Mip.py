import itertools

class Arm:

# class methods

  Instances                            			 = {}
  @classmethod
  def Instancecount(self):
    return len(self.Instances)
  @classmethod
  def getInstances(self):
    return (value for value in self.Instances.keys())
  @classmethod
  def Get(self, target):
    """return the Instance whose gid is target"""
    return self.Instances.get(target, None)

# fundamental methods

  def __init__(self, name, type, chr, start, stop, seq):
    self._name						= name				# same as mip._name
    self._type						= type
    self._chr						= chr
    self._start						= start
    self._stop						= stop
    self._seq						= seq
    self._mip						= None				# will eventually point to mip object
    self._alignments                                    = list()			# will eventually have list of alignment objects (corresponding to successful alignments against read)
    self.Instances["_".join([self._name,self._type])]	= self				# unique identifier to access later; combines mip._name and arm._type
  def __str__(self):
    return self._name
  def __repr__(self):
    return self._name

# predicate methods
# access methods ("getters")

  def getId(self):
    return "_".join([self._name,self._type])
  def getMip(self):
    return self._mip
  def getStartPos(self):
    return self._start
  def getSeq(self):
    return self._seq
  def getType(self):
    return self._type

# modification methods ("setters")

  def setMip(self, mip):
    self._mip = mip
  def addAlignment(self, alignment):
    self._alignments.append(alignment)

# action methods ("doers")
# helper methods


class Mip:

# class methods

  Instances			= {}
  @classmethod
  def Instancecount(self):
    return len(self.Instances)
  @classmethod
  def getInstances(self):
    return (value for value in self.Instances.keys())
  @classmethod
  def Get(self, target):
    """return the Instance whose gid is target"""
    return self.Instances.get(target, None)

# fundamental methods

  def __init__(self, name, ext, lig, target_chr, probe_strand):
    self._name 			= name
    self._ext			= ext				# arm object
    self._lig			= lig				# arm object
    self._target_chr		= target_chr
    self._probe_strand		= probe_strand
    self._read_pairs		= list()			# list of readpair objects
    self._AlignedSegments	= list()			# will eventually hold a list of pysam.AlignedSegment objects - one for each aligned read corresponding to this mip
    self._onTargetBed		= list()			# will eventually hold a list of the corresponding [chr, start, end] for the mip
    self._on_target_depth	= None
    self._off_target_depth	= None
    self._unmapped_depth	= None
    self._mapped_MQ0_count	= None
    self._on_target_MQ10_depth	= None
    self._on_target_MQ60_depth	= None
    self.Instances[self._name]	= self				# adds key:value pair of name:object of current Instance of mip to class dictionnary of all mips
  def __str__(self):
    return self._name
  def __repr__(self):
    return self._name

# Predicates
# Access Methods ("getters")

  def getId(self):
    return self._name
  def getArm(self, type):
    return self._ext if type == 'ext' else self._lig
  def getArms(self):
    return [self._ext, self._lig]
  def getOnTargetBed(self):					# NOTE: transforms the start position to 0-based, keeps the stop position as 1-based, i.e. proper BED format
    return [self._target_chr, min(map(int, [self._ext._start, self._ext._stop, self._lig._start, self._lig._stop]))-1, max(map(int, [self._ext._start, self._ext._stop, self._lig._start, self._lig._stop]))]
  def getFasta(self):
    return "\n".join(list(itertools.chain(*[["_".join([">"+arm._name, arm._type]), arm.getSeq()] for arm in self.getArms()])))
  def getTargetLength(self):
    chr, start, stop = self.getOnTargetBed()
    return stop-start

# Modification Methods ("setters")

  def addReadPairToMip(self, readPair):				# appends one ReadPair object to the list variable self._read_pairs
    self._read_pairs.append(readPair)
  def setOnTargetBed(self, onBed):
    self._onTargetBed = onBed
  def setOnTargetDepth(self, depth):
    self._on_target_depth = depth
  def setOffTargetDepth(self, depth):
    self._off_target_depth = depth
  def setUnmappedDepth(self, depth):
    self._unmapped_depth = depth
  def setMappedMQ0Count(self, depth):
    self._mapped_MQ0_count = depth
  def setOnTargetMQ10Count(self, depth):
    self._on_target_MQ10_depth = depth
  def setOnTargetMQ60Count(self, depth):
    self._on_target_MQ60_depth = depth

# Action Methods ("doers")

# Helper Methods ("helpers")


