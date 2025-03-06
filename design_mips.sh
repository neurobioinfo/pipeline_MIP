#!/bin/bash

unset PIPELINE_HOME dbsnp_file genome_directory gene_list project_name exon_type flank

# function to print help
Usage() {
        echo
        echo -e "Usage:\t$0 [arguments]"
        echo -e "\tmandatory arguments:\n" \
          "\t\t-l (--gene_list)        = gene list file\n" \
          "\t\t-p (--project_name)     = project name\n" \
          "\t\t-x (--exon_type)        = exon type (options are \"coding\" or \"all\")\n"
        echo -e "\toptional arguments:\n " \
          "\t\t-h (--help)             = get the program options and exit\n" \
          "\t\t-P (--pipeline_home)    = pipeline home (default: /RQexec/spiegelm/data/pipeline_MIPs.svn); hint: use the default!)\n" \
          "\t\t-s (--dbsnp_file)       = dbsnp file (default: /RQexec/spiegelm/data/pipeline_MIPs.svn/soft/src/mipgen/common_all.SNP.vcf.gz)\n" \
          "\t\t-g (--genome_directory) = genome directory (default: /RQexec/spiegelm/data/pipeline_MIPs.svn/data/reference/mipgen.hg19)\n" \
          "\t\t-f (--flanking_bases)   = number of flanking bases to add around each target exon\n" \
          "\t\t-i (--min_capture_size) = min_capture_size to use with mipgen design. Default is 152, minimum is 120.\n" \
          "\t\t-a (--max_capture_size) = max_capture_size to use with mipgen design. Default is 152, maximum is 250.\n"
        echo
}

# ===============================================
# PARSING ARGUMENTS
# ===============================================

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name design --alternative --unquoted --options hP:s:g:l:p:x:f:i:a: --longoptions help,pipeline_home:,dbsnp_file:,genome_directory:,gene_list:,project_name:,exon_type:,flanking_bases:,min_capture_size:,max_capture_size: -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi
set -- $options
while [ $# -gt 0 ]
do
    case $1 in
    -h| --help) Usage; exit 0;;
    # for options with required arguments, an additional shift is required
    -P| --pipeline_home)    PIPELINE_HOME=$2; shift;;
    -s| --dbsnp_file)       dbsnp_file=$2; shift;;
    -g| --genome_directory) genome_directory=$2; shift;;
    -l| --gene_list)        gene_list=$2; shift;;
    -p| --project_name)     project_name=$2; shift;;
    -x| --exon_type)        exon_type=$2; shift;;
    -f| --flanking_bases)   flank=$2; shift;;
    -i| --min_capture_size) min_capture_size=$2; shift;;
    -a| --max_capture_size) max_capture_size=$2; shift;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done

# error checking
unset failure
if [[ ! -s $gene_list ]]; then echo "error: gene_list file $gene_list (-g or --gene_list) empty or not specified"; failure=1; fi
if [[ -z $project_name ]]; then echo "error: project_name (-p or --project_name) not specified"; failure=1; fi
if [[ -z $exon_type ]]; then echo "error: exon_type (-x or --exon_type) not specified"; failure=1; fi
if [[ $exon_type != "coding" && $exon_type != "all" ]]; then echo "error: exon_type (-x or --exon_type) must be one of \"coding\" or \"all\""; failure=1; fi
if [[ $failure -eq 1 ]]; then echo "$(basename $0) error. exiting"; exit 42; fi
# set default variables, if empty
if [[ -z $PIPELINE_HOME    ]]; then PIPELINE_HOME=/RQexec/spiegelm/data/pipeline_MIPs.svn; fi
if [[ -z $dbsnp_file       ]]; then dbsnp_file=/RQexec/spiegelm/data/pipeline_MIPs.svn/soft/src/mipgen/common_all.SNP.vcf.gz; fi
if [[ -z $genome_directory ]]; then genome_directory=/RQexec/spiegelm/data/pipeline_MIPs.svn/data/reference/mipgen.hg19; fi
if [[ -z $min_capture_size ]]; then min_capture_size=152; fi
if [[ -z $max_capture_size ]]; then max_capture_size=152; fi

# init pipeline (initializes other variables, e.g. $REF)
. $PIPELINE_HOME/data/templates/init_pipeline.sh $PIPELINE_HOME

#check if all genes in $gene_list are in refSeq data
awk -F"\t" '{if (NR==FNR) {genes[$1]++; next}; if ($13 in genes) good[$13]++}END{n=0; d=length(genes); missing=""; print length(good)"/"length(genes)" genes have refseq data"; if (length(good)!=length(genes)) {for (i in genes) if (i in good) {continue} else {missing=i","missing}; gsub(/,$/,"",missing); print "missing genes: "missing}}' $gene_list $GENE_LIST > $project_name.design_mips.out

#run mip design
#choose between two different ways to run design: coding exons only, or coding + UTR
if [[ $exon_type == "coding" ]]; then extract_script=extract_coding_gene_exons.sh; else extract_script=extract_gene_exons_plus_utrs.sh; fi
#create bedfile to input to mipgen
$extract_script $gene_list $GENE_LIST|sort -k2n|sortByRef.pl - $PIPELINE_HOME/data/reference/gatk.ref.dict |awk -v flank=$flank 'BEGIN{FS=OFS="\t"}{print $1,$2-flank,$3+flank}'|bedtools merge -i - > $project_name.bed
# run mipgen
mipgen -genome_dir $genome_directory -project_name $project_name -bwa_genome_index $REF -regions_to_scan $project_name.bed -snp_file $dbsnp_file -min_capture_size $min_capture_size -max_capture_size $max_capture_size -tag_sizes 0,0 2>&1|tee -a $project_name.design_mips.out
# create UCSC track
python $PIPELINE_HOME/soft/src/mipgen/tools/generate_ucsc_track.py $project_name.picked_mips.txt ${project_name}_ucsc_track blue red
# if any mips have SNPs in the arms, add mips for allele "b"; create final output file named $project.design
cat <(grep -v SNP_a $project_name.picked_mips.txt) <(grep SNP_a $project_name.picked_mips.txt) <(for i in $(grep SNP_a $project_name.picked_mips.txt|awk '{print $NF}'|sed 's/_a/_b/g'); do grep $i $project_name.snp_mips.txt; done)|sort -t$'\t' -k20 > $project.design
