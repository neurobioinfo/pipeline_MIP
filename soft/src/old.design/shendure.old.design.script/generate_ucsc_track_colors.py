#generate_ucsc_track_colors.py
#tweaked by Beth Martin, lot of stuff in here isn't used, but I just left it.
#usage: python generate_ucsc_track.py <input_file> <output_file> <mipname> <plus color> <minus color> color choices are: purple/green/red/gold/blue/orange/pink/brown

import sys

if len(sys.argv)!=6:
	print "usage: python generate_ucsc_track.py <input_file> <output_file> <mipname> <plus color> <minus color> color choices are: purple/green/red/gold/blue/orange/pink/brown"

infile = sys.argv[1]
outfile = sys.argv[2]
mipname = sys.argv[3]
pluscolor = sys.argv[4]
minuscolor = sys.argv[5]

color = {'purple':'106,90,205', 'green':'85,107,47', 'red':'205,92,92', 'gold':'218,165,32', 'blue':'65,105,225', 'orange':'255,69,0', 'pink':'255,105,180', 'brown':'139,69,19'}

plus_strand_lines = []
minus_strand_lines = []


fin = open(infile)
fin.readline()
for line in fin:
	values = line.strip().split('\t')
	if (len(values) > 9 and values[5] != 'NA' and values[4].isdigit()):
		name = values[0] + "_" + values[1]
		score = str((float(values[1]) * 3.0481996570775385) + 500)
		chromosome = values[2]
		ext_probe_start = int(values[3])
		ext_probe_stop = int(values[4])
		lig_probe_start = int(values[7])
		lig_probe_stop = int(values[8])
		strand = values[17]
		start_stop = [ext_probe_start, ext_probe_stop, lig_probe_start, lig_probe_stop]
		start_stop.sort()
		start_val = start_stop[0]-1
		thick_start = start_stop[1]
		thick_end = start_stop[2]-1
		end_val = start_stop[3]
		if strand==('+'):
			plus_strand_lines.append('chr' + chromosome + '\t' + str(start_val) + '\t' + str(end_val) + '\t' + name + '\t' + score + '\t' + strand + '\t' + str(thick_start) + '\t' + 
str(thick_end) + '\t' + color[pluscolor] + '\n')
		else:
			minus_strand_lines.append('chr' + chromosome + '\t' + str(start_val) + '\t' + str(end_val) + '\t' + name + '\t' + score + '\t' + strand + '\t' + str(thick_start) + '\t' + 
str(thick_end) + '\t' + color[minuscolor] + '\n')
fin.close()
fout = open(outfile, 'w')
fout.write(
'track name=' + mipname + '_plus' + ' itemRgb=on\n')
for entry in plus_strand_lines:
	fout.write(entry)

fout.write('\n')
fout.write(
'track name=' + mipname + '_minus' + ' itemRgb=on\n')
for entry in minus_strand_lines:
	fout.write(entry)

fout.close()

