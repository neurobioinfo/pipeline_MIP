from swalign import *
import Alignment

class BasicSearch:

  def __init__(self, read, arm):
    self._read	= read
    self._arm	= arm

  def alignForward(self):
    read_seq	= self._read.getSeq()
    arm_seq	= self._arm.getSeq()
    if read_seq[0:len(arm_seq)] == arm_seq:
      alignment = Alignment.ReadArmAlignment(self._read, self._arm, 0, 0, [(len(self._arm.getSeq()), "M")], 1, False)
      self._read.addAlignment(alignment)
      self._arm.addAlignment(alignment)


  def alignRevComp(self):
    read_seq	= self._read.getSeq()
    arm_seq	= revcomp(self._arm.getSeq())
    if read_seq[0:len(arm_seq)] == arm_seq:
      alignment = Alignment.ReadArmAlignment(self._read, self._arm, 0, 0, [(len(self._arm.getSeq()), "M")], 1, True)
      self._read.addAlignment(alignment)
      self._arm.addAlignment(alignment)

#class SmithWatermanSearch:
