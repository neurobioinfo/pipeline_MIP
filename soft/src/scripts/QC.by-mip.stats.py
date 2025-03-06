import sys
sys.path.append('/RQexec/spiegelm/data/pipeline_MIPs.svn/soft/src/mip_pipeline_python')
#from printStruct import printStruct as ps

files=sys.argv[1:]
mipData={}

#stuff all data from all files into one nice hash of hash of hashes
for file in files:
  with open(file) as f:
    for line in f.read().splitlines():
      try: mip,on,off,un,MQ0,MQ10,MQ60 = line.split("\t")
      except: raise RuntimeError(line.split("\t"))
      if mip not in mipData: mipData[mip]={file: {key: 0 for key in ['on','off','un','MQ0','MQ10','MQ60']}}
      mipData[mip][file]={'on': on,'off': off,'un': un, 'MQ0':MQ0, 'MQ10':MQ10, 'MQ60':MQ60}

#print out all data in nicely ordered format
print "\t".join(['', "\t".join(["\t".join([i for file in mipData[mipData.keys()[1]]]) for i in ['on','off','un','MQ0','MQ10','MQ60']])])	#header line 1 ('on', 'off', 'un' sections)
print "\t".join(['mip', "\t".join(["\t".join([file for file in mipData[mipData.keys()[1]]]) for i in xrange(0,3)])])	#header line 2 (file names)
try:
  for mip in sorted(mipData, key=int):
    print "\t".join([mip,"\t".join(["\t".join([str(mipData[mip][file][key]) for file in mipData[mip]]) for key in ['on','off','un','MQ0','MQ10','MQ60']])])	#body
except:
  for mip in sorted(mipData):
    print "\t".join([mip,"\t".join(["\t".join([str(mipData[mip][file][key]) for file in mipData[mip]]) for key in ['on','off','un','MQ0','MQ10','MQ60']])])     #body
