#!/bin/bash

#script usage
echo "usage: $(dirname $0)/$(basename $0) [gene_list file] [config file] [project_name]"

#define scripts and variables used by this program
. $PIPELINE_HOME/data/templates/init_pipeline.sh $PIPELINE_HOME
GENE_LIST=$PIPELINE_HOME/data/targets/RefSeq.Genes.v37.refGene
old_design_script=$PIPELINE_HOME/soft/src/old.design/shendure.old.design.script/runMipDesignPipeline_v1.py

#read input variables
read gene_list cfg_file project_name <<< $@

#error checking
unset failure
if [[ ! -s $gene_list ]]; then echo "error: gene_list file $gene_list (1st argument) empty or not specified"; failure=1; fi
if [[ -z $project_name ]]; then echo "error: project_name (2nd argument) not specified"; failure=1; fi
if [[ $failure -eq 1 ]]; then echo "$(basename $0) error. exiting"; exit 42; fi

#check if all genes in $gene_list are in refSeq data
awk -F"\t" '{if (NR==FNR) {genes[$1]++; next}; if ($13 in genes) good[$13]++}END{n=0; d=length(genes); missing=""; print length(good)"/"length(genes)" genes have refseq data"; if (length(good)!=length(genes)) {for (i in genes) if (i in good) {continue} else {missing=i","missing}; gsub(/,$/,"",missing); print "missing genes: "missing}}' $gene_list $GENE_LIST > $project_name.design_mips.out

#helper function in case you need to run script multiple times to cover hard-to-design regions
checkIfMissingAfterRedo() {
  read previous_bed previous_70mers counter <<< $@
  new_bed=${project_name}.redo$[$counter+1].bed
  grep -P "($(diff <(awk '{print $2+1,$3}' $previous_bed|sort -n) <(awk '{if (NR!=1) print $15,$16}' $previous_70mers |sort -n|uniq)|awk '{print $3}'|sort|grep -v ^$|tr '\n' '|'|sed 's/|$//g'))$" $previous_bed > $new_bed
}

#run mip design:
#step 1: make bed file
#extract_coding_gene_exons.sh $gene_list $GENE_LIST|sort -k2n|sortByRef.pl - /RQusagers/spiegelm/references/gatk.ref.dict |bedtools merge -i - > $project_name.bed
extract_gene_exons_plus_utrs.sh $gene_list $GENE_LIST|sort -k2n|sortByRef.pl - /RQusagers/spiegelm/references/gatk.ref.dict |bedtools merge -i - > $project_name.bed
#step 2: run design script once and check if any pieces are missing from target bed
python $old_design_script $cfg_file $gene_list $project_name 20 $project_name.bed
c=1	#create redo counter; will increment if necessary later with ((c++))
checkIfMissingAfterRedo $project_name.bed ${project_name}_70mers.txt $c

#step 3: if any regions missing, keep running design script on missing pieces until it's all done
while [[ -s ${project_name}.redo$c.bed ]]; do
  python $old_design_script $cfg_file $gene_list ${project_name}.redo$c 20 ${project_name}.redo$c.bed
  checkIfMissingAfterRedo ${project_name}.redo$c.bed $c
  ((c++))
done

#NOTE: I NEVER FINISHED THIS SCRIPT! THERE STILL NEEDS TO BE A STEP WHICH COMBINES THE *_70mers FILES FROM ALL SUB-STEPS INTO A SINGLE FILE
