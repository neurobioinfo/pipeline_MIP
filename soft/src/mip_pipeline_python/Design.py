import os
import sys
import re
import Mip
#from itertools import izip

class DesignException(Exception):

  def __init__(self, errmsg):
    self._errmsg	= errmsg

  def __str__(self):
    return repr(self._errmsg)

class Design:

# fundamental methods

  def __init__(self, file, source):
    self._file			= file
    self._source		= source
    self._mips			= list()	# empty list of mip objects, to be parsed from design file using function 'parse'
    # try to open input file, and throw exception if it fails
    try:
      self._filestream	= open(file, 'r')
    except IOError as e:
      print("I/O Error({0}): {1} (while trying to open file {2})").format(e.errno, e.strerror, file)
    # parse header line: ...
    line = self._filestream.readline().rstrip()
    self._header_dict		= dict(zip(line.split('\t'), [line.split('\t').index(i) for i in line.split('\t')]))
    self._mandatory_columns 	= ['chr', 'lig_probe_start', 'lig_probe_stop', 'lig_probe_sequence', 'ext_probe_start', 'ext_probe_stop', 'ext_probe_sequence', 'probe_strand']
    identifier_columns	= 	{	'mipgen':['mip_name'],
					'breakpoint_resolution_wgs_mips':['>mip_pick_count','rank_score'],
					'generic':['mip_id']
				}
    mip_target_columns	= 	{	'mipgen':['mip_scan_start_position','mip_scan_stop_position'],
					'breakpoint_resolution_wgs_mips':['mip_target_start_position','mip_target_stop_position'],
					'generic':['target_start','target_stop']
				}
    # ... define which column(s) will supply mip names, depending on input source (generic vs mipgen vs breakpoint_resolution_wgs_mips) ...
    if self._source not in identifier_columns:
      print('Input Argument Error: given source \"{0}\" not among accepted options [\"generic\", \"mipgen\",\"breakpoint_resolution_wgs_mips\"]'.format(self._source))
      raise ValueError
    try:
      self._name_cols		= [self._header_dict[e] for i,e in enumerate(identifier_columns[self._source])]
      self._mip_target_cols	= [self._header_dict[e] for i,e in enumerate(mip_target_columns[self._source])]
      self._mandatory_columns.extend(identifier_columns[self._source])
      self._mandatory_columns.extend(mip_target_columns[self._source])
    except KeyError:
      print("Input Argument Error: given source \"{0}\" does not match the header line").format(self._source)
      raise RuntimeError(self._header_dict)
    # ... check if any mandatory columns are missing; die + list missing columns if they are
    missing_mandatory	= [col for col in self._mandatory_columns if col not in self._header_dict]
    if len(missing_mandatory) != 0: raise DesignException('missing mandatory columns: ' + str(missing_mandatory))

  def parse(self):
    self._filestream.seek(0); self._filestream.readline()	#to set filestream to 2nd line of input file (i.e. 1st non-header line)
    lines=[line.rstrip() for line in self._filestream.readlines()]
    for line in lines:
      mip_name = '_'.join([line.split('\t')[i] for i in self._name_cols]).rstrip()
      try: line_dict = dict(zip([x for x in self._mandatory_columns], [line.split('\t')[self._header_dict[x]] for x in self._mandatory_columns]))
      except: raise RuntimeError(line)
      ext = Mip.Arm(mip_name, 'ext', line_dict['chr'], line_dict['ext_probe_start'], line_dict['ext_probe_stop'], line_dict['ext_probe_sequence'])
      lig = Mip.Arm(mip_name, 'lig', line_dict['chr'], line_dict['lig_probe_start'], line_dict['lig_probe_stop'], line_dict['lig_probe_sequence'])
      mip = Mip.Mip(mip_name, ext, lig, line_dict['chr'], line_dict['probe_strand'])
      ext.setMip(mip)
      lig.setMip(mip)
      mip.setOnTargetBed(mip.getOnTargetBed())
      self._mips.append(mip)

  def generateBedFiles(self, outFile):
    self._filestream.seek(0); self._filestream.readline()	#to set filestream to 2nd line of input file (i.e. 1st non-header line)
    lines=[line.rstrip() for line in self._filestream.readlines()]
    with open(outFile, 'w') as BED, open('.'.join([outFile, 'targets-only.bed']), 'w') as TARG:
      for line in lines:
        mip_name = '_'.join([line.split('\t')[i] for i in self._name_cols]).rstrip()
        mip = self._mips[0].Get(mip_name)
        line_dict = dict(zip([x for x in self._mandatory_columns], [line.split('\t')[self._header_dict[x]] for x in self._mandatory_columns]))
        BED.write('\t'.join(['\t'.join(map(str, mip.getOnTargetBed())), mip._name]) + '\n')
        TARG.write('\t'.join([line_dict['chr'], '\t'.join(sorted([line_dict[x] for x in line_dict if re.match('(mip_)*(target|scan)_(start|stop)(_position)*', x)])), mip_name]) + '\n')
