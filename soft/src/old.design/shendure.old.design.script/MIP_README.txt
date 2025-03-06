#Beth Martin's MIP Notes. martin91@uw.edu
#General mip design pipeline for tiling mips which alternate strands, 112 bp inserts.

#2013-03-25 update - I'm not sure that the scripts work if you put "n" under the "check for snps" option. So just put "y".

#2012-12-13 update
Clarifying that you do NOT want to use the fasta_list_hg19 or _hg18 that is on our server. Just use it as an example fasta list and make your own with the FULL
PATHS to the fasta files in your local genome directory. 

#2012-12-04 PLEASE READ
The script makeRunMipDesignShellScript2.py has been updated to fix a path.
Some are having issues with the binary file genome_compare. If you don't get files that look like: *.all_mips.copy_counts then this is likely the problem. 
We've put the uncompiled version in the folder: genome_compare_distrib
Also, we've added some test input and output files in the folder: testfiles
Make sure to check out the README file in that directory, and it'll walk you through everything. 


#MIP DESIGN
#for all these scripts, make sure you are running the latest version of python
#Need a start/stops bed file for design regions
#Bed files are base 0 start, base 1 stop

#An option for refseq genes
Go to table browser at UCSC, choose:
	Genes and Gene Predictions track
	track RefSeq genes
	table refGene
	paste in your gene names
	output format BED
	coding exons
	output to Galaxy
		in Galaxy:
			Operate on Genomic Intervals - Merge
In Excel, add (-5) to starts and (+5) to stops. Save as a tab-delimited BED file. The amount of flanking sequence is optional.

#Run the shell creation script:
python makeRunMipDesignShellScript2.py <bedfile> 112 112 1 <mip set name> <path to scripts> <genome directory> <fasta list> y <snp file> 0 > <outputfilename.sh>

the arguments in this command line are:
<bedfile> = the name of the bedfile with the starts and stops
112 = insert size
112 = feature flank #i.e. how much of the flanking sequence to scan for possible mip probes
1 = mip overlap #i.e. the minimum tiling mips need to overlap
<mip set name> = log prefix - the name that all the resulting files will start with
<path to scripts> = directory where all the scripts are kept, example = /net/shendure/vol6/mip_pipeline
<genome directory> = the directory with your reference sequence, we use hg19, example = /net/shendure/vol6/mip_pipeline/hg19
<fasta list> = a text file with all the fasta file names, needs the full path in the file names, example = /net/shendure/vol6/mip_pipeline/fasta_list_hg19
note: do NOT use the fasta_list_hg19/hg18 from our server, those have our paths in them and they won't work fo you. You need to make your own file.
y = check for snps - y or n (just put "y", I don't think it works if you put "n")
<snp file> = file with all the snps in it. specify a file, even if not using check snp option. You can use a dummy file, one line, if you wish to ignore SNPs. This is the ROD file format see:https://cgwb.nci.nih.gov/cgi-bin/hgTables?db=hg18&hgta_group=varRep&hgta_track=snp129&hgta_table=snp129&hgta_doSchema=describe+table+schema, example = /net/shendure/vol6/mip_pipeline/snp132.txt
0 = capture flank #i.e. add any additional capture bases to the flanks of the bed file coordinates 
> <outputfilename.sh> = output to a file

#Then you need to change the permissions of this shell file to make executable:
chmod 777 <outputfilename.sh>

#then execute it:
./<outputfilename.sh>

#it will take a while to run. 
The final output file at the end in "...all_mips.copy_counts.ranked.ranked_list.txt"

#If there are snps in the mips (and there usually are), run the following (replacing the genome directory and the snpfile with your path/filename):
perl redesign_mips_with_snps.pl -mips_file <the output file with long mip filename> -genome_dir /net/shendure/vol6/mip_pipeline/hg19 -snp_file /net/shendure/vol6/mip_pipeline/snp132.txt

#To make a UCSC track:
python generate_ucsc_track_colors.py <huge filename> <yourprefix_ucsc.bed> <mipname> <plus color> <minus color>
#color choices are: purple/green/red/gold/blue/orange/pink/brown

#Everything look okay? Make the 70mers:
python generate_70mers2.py <output file with even bigger filename> <prefix_70mers.txt>
If importing into excel, choice the tab delimited option, with treat consecutive delimiters as one option. 

#Best practices
You should remove from the final list any mips with very high arms counts, greater than 10. Especially if both arms are repetitive. Mips with scores of 5 or 3 generally work well. Those with -1 scores, often require re-balancing, at higher concentration, which can increase the capture coverage. 
