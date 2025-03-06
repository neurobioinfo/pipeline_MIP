from swalign import *

class ReadArmAlignment(Alignment):

  def __init__(self, read, arm, q_pos = None, r_pos = None, cigar = (), score = None, rc = False):
    Alignment.__init__(self, read.getSeq(), arm.getSeq(), q_pos, r_pos, cigar, score, read.getId(), arm.getId(), rc)
    self._read  = read
    self._arm	= arm
    self._rc	= rc

  def __str__(self):
    return ":".join([str(self._read), str(self._arm)])

  def getDetails(self):
    return "\n".join(["read =" + self._read.getId(), "arm = " + self._arm.getId(), "revcomp = " + str(self._rc)])
