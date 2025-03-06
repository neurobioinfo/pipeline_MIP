#!/bin/bash

#howto
echo "usage: bash $0 [project_name] [gene_list]" >/dev/stderr

#parse input variables
read project_name gene_list <<< $@

#check input args
unset failure
if [[ -z ${project_name} ]]; then echo "error: project_name (1st argument) not specified"; failure=1; fi
if [[ ! -s $gene_list ]]; then echo "error: gene list (2nd argument) empty or not specified"; failure=1; fi
if [[ $failure -eq 1 ]]; then echo "encountered errors; exiting"; exit 42; fi

#set dependent variables
new_design_script=/alloc/grouleau/COMMUN/spiegelm/runs/MIP/scripts/design_mips.sh
old_design_script=/alloc/grouleau/COMMUN/spiegelm/runs/MIP/design/runMipDesignPipeline_v1.py
cfg_file=/alloc/grouleau/COMMUN/spiegelm/runs/MIP/design/runMipDesignPipeline_v1.cfg

#run mipgen (current Schendure lab pipeline)
bash $new_design_script $gene_list ${project_name}

#run old Shendure mip design pipeline iteratively, for any missing mips
	#helper function
checkIfMissingAfterRedo() {
  read previous_bed previous_70mers counter <<< $@
  new_bed=${project_name}.redo$[$counter+1].bed
  grep -P "($(diff <(awk '{print $2+1,$3}' $previous_bed|sort -n) <(awk '{if (NR!=1) print $15,$16}' $previous_70mers |sort -n|uniq)|awk '{print $3}'|sort|grep -v ^$|tr '\n' '|'|sed 's/|$//g'))$" $previous_bed > $new_bed
}
	#if any regions missing, redesign once; if any regions still missing, launch iterative rounds of redesign until all regions covered
if [[ -s ${project_name}.coverage_failed.bed ]]; then
  c=1	#create redo counter; will increment if necessary later with ((c++))
  python $old_design_script $cfg_file $gene_list ${project_name}.redo$c 20 ${project_name}.coverage_failed.bed
  checkIfMissingAfterRedo ${project_name}.coverage_failed.bed ${project_name}_70mers.txt $c
  while [[ -s ${project_name}.redo$[$c+1].bed ]]; do
    python $old_design_script $cfg_file $gene_list ${project_name}.redo$[$c+1] 20 ${project_name}.redo$[$c+1].bed
    checkIfMissingAfterRedo ${project_name}.redo$c.bed $[$c+1]
    ((c++))
  done
fi
