import os
import sys
import ConfigParser
import subprocess

# parse input variables 
config_file, gene_list, experiment_name, flanking_interval, input_bed = sys.argv[1:]

# read config file and update global variables
# note: config file must be in python configParser format - i.e. section names surrounded by []; and variables defined as variable: value, one per line; this script will only load variables from the first section; 
#     Variables in this config file must include: 
#	scripts_path (path to MIPgen scripts)
#	genome_dir (directory of genome fa files)
#	fasta_list (lists full paths of all files in genome_dir)
#	snp_file (UCSC data dump of your favourite dbsnp version)
#	genome_reference (fasta.fai of entire genome)
#     These are global variables to be used in subsequent steps
cfg = ConfigParser.RawConfigParser()
cfg.read(config_file)
globals().update(cfg.items(cfg.sections()[0]))

# create intermediate output file names
experiment_bed			= experiment_name + '.bed' if input_bed == '' else os.path.basename(input_bed)
mip_design_script		= experiment_name + ".mip_design.sh"
mip_design_script_outfile	= experiment_bed + '.picked_mip_probe_arms.' + experiment_bed + '.scan112.all_mips.copy_counts.ranked.ranked_list.txt'
redesign_w_snps_outfile		= mip_design_script_outfile + ".fixed_snps"
generateUCSCtrack_outfile	= experiment_name + '_ucsc.bed'
generate70mers_outfile		= experiment_name + '_70mers.txt'

# create various system commands to be used by subprocess.Popen in subsequent steps
body_awk                   = ['awk', """BEGIN{FS=OFS="\t"}{if (FNR==NR) {genes[$1]++; next}; if (substr(FILENAME,length(FILENAME)-2,3)!="bed") {if ($13 in genes) isoforms[$2]++; next}; if (substr(FILENAME,length(FILENAME)-2,3)=="bed") {iso=$4; gsub(/_cds.*/,"",iso); if (iso in isoforms) print $1,$2,$3}}""", gene_list, refseq_flat, refseq_bed]
body_bedtoolsSlop          = ['bedtools', 'slop', '-i', '-', '-g', genome_reference, '-b', flanking_interval]
body_bedtoolsMerge         = ['bedtools', 'merge', '-i',  '-']
body_createMipDesignScript = ['python', scripts_path + '/makeRunMipDesignShellScript2.py', experiment_bed, '112', '112', '1', experiment_name, scripts_path, genome_dir, fasta_list, 'y', snp_file, '0']
body_runMipDesignScript    = ['sh', mip_design_script]
body_redesignMipsWithSnps  = ['perl', scripts_path + '/redesign_mips_with_snps.pl', '-mips_file', mip_design_script_outfile, '-genome_dir', genome_dir, '-snp_file', snp_file]
body_generateUCSCtrack     = ['python', scripts_path + '/generate_ucsc_track_colors.py', redesign_w_snps_outfile, generateUCSCtrack_outfile, experiment_name, 'blue', 'red']
body_generate70mers        = ['python', scripts_path + '/generate_70mers2.py', redesign_w_snps_outfile, generate70mers_outfile]

# generate bed file from gene list
if input_bed == '':
  cmd_awk           = subprocess.Popen(body_awk, stdout=subprocess.PIPE)
  cmd_bedtoolsSlop  = subprocess.Popen(body_bedtoolsSlop, stdin=cmd_awk.stdout, stdout=subprocess.PIPE)
  cmd_bedtoolsMerge = subprocess.Popen(body_bedtoolsMerge, stdin=cmd_bedtoolsSlop.stdout, stdout=open(experiment_bed, 'w'))
  cmd_bedtoolsMerge.communicate()

# initial mip design
cmd_createMipDesignScript = subprocess.Popen(body_createMipDesignScript, stdout=open(mip_design_script, 'wb'))
cmd_createMipDesignScript.wait()
cmd_runMipDesignScript    = subprocess.Popen(body_runMipDesignScript, stderr=subprocess.STDOUT, stdout=open(mip_design_script + '.log', 'w'))
cmd_runMipDesignScript.wait()

# redesign mips with snps in them
cmd_redesignMipsWithSnps  = subprocess.Popen(body_redesignMipsWithSnps, stdout=open(redesign_w_snps_outfile + ".log", 'w'))
cmd_redesignMipsWithSnps.wait()

## generate UCSC custom track for visualization
cmd_generateUCSCtrack     = subprocess.Popen(body_generateUCSCtrack, stdout=open(generateUCSCtrack_outfile + '.log', 'w'))

# make the 70mers file (== the designFile for PreAlignmentProcessing.py) - name will be [experiment_name]_70mers.txt
cmd_generate70mers        = subprocess.Popen(body_generate70mers)

for process in [cmd_generateUCSCtrack, cmd_generate70mers]:
  process.wait()
#cmd_generate70mers.wait()
