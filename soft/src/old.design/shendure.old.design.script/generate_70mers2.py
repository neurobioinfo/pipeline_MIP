import sys

if len(sys.argv)!=3:
	print "usage: python generate_70mers.py <input_file> <output_file>"

infile = sys.argv[1]
outfile = sys.argv[2]


#left_flank = 'AGGACCGGATCAACT'
#right_flank = 'CATTGCGTGAACCGA'
backbone_sequence = 'CTTCAGCTTCCCGATATCCGACGGTAGTGT'

fout = open(outfile, 'w')
fin = open(infile)
#fin.readline()
for line in fin:
	if line.startswith('>'):
		line=line.rstrip('\n')
		fout.write(line+ '\t'+"70mer"+ '\n')
	else:
		values = line.split('\t')
		if (len(values) > 9 and line[0]!=">mip_pick_count"):
			ligation_sequence = values[9]
			extension_sequence = values[5]
			if ligation_sequence != 'NA':
				#mip_sequence = left_flank.upper() + ligation_sequence.lower() + backbone_sequence.upper() + extension_sequence.lower() + right_flank.upper()
				mip_sequence = ligation_sequence.lower() + backbone_sequence.upper() + extension_sequence.lower()
				fout.write(line.strip() + '\t' + mip_sequence + '\n')
fin.close()
fout.close()
