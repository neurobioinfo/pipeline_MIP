#!/usr/bin/python

import sys

def printStruct(struc, indent=0, file=None):
  if file is not None:
    sys.stdout = open(file, 'wb')
  if isinstance(struc, dict):
    print('  '*indent+'{')
    for key,val in struc.iteritems():
      if isinstance(val, (dict, list, tuple)):
        print('  '*(indent+1) + str(key) + '=> ')
        printStruct(val, indent+2)
      else:
        print('  '*(indent+1) + str(key) + '=> ' + str(val))
    print('  '*indent+'}')
  elif isinstance(struc, list):
    print('  '*indent + '[')
    for item in struc:
      printStruct(item, indent+1)
    print('  '*indent + ']')
  elif isinstance(struc, tuple):
    print('  '*indent + '(')
    for item in struc:
      printStruct(item, indent+1)
    print('  '*indent + ')')
  else: print('  '*indent + str(struc))
