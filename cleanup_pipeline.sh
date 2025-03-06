#!/bin/bash

# NOTE: this script MUST be run in the sample's analysis directory - i.e. S?????/Illumina_etc/GATK_BWA_v37

rm -rf temp
rm -v *trimmomaticPE.*paired.gz .*.done *.TMP_MERGE* *sam *fastq.gz *aligned.merged.sorted.bam* *intervals
find -perm 600 -exec chmod 664 \{\} \;;
mkdir LOGS;
mv *.o[0-9]* *.out *.stderr *DataProcessingPipeline.*.txt launch_record* launch.txt *jobreport.pdf LOGS;
mkdir STATS;
mv DepthOfCoverage.* *.metrics *.idxstats *.flagstat *.vs.*.csv *.verifyBamID.* fastqc.* *on_off_unmapped *mip-by-mip_counts *preprocessing_report.txt *redundant_mips *Callability.bed STATS;
tar czvf SMALL_FILES.tar.gz $(echo  $(find . -maxdepth 1 -type f -size -1024k)    $(find -name "MT.*.ba*")    $(find -name "*non_pri_aln.ba*")    LOGS STATS | tr ' ' '\n' |  sed -e '/gz$/d' -e '/tbi$/d'  -e '/clean.ba/d' -e '/HC.ba/d' -e '/preprocessing_report/d' -e '/redundant_mips/d'| sort | uniq);
rm -rfv $(echo    $(find . -maxdepth 1 -type f -size -1024k)    $(find -name "MT.*.ba*")    $(find -name "*non_pri_aln.ba*")    LOGS STATS temp | tr ' ' '\n' |  sed -e 's/ /\n/g' -e '/gz$/d' -e '/tbi$/d' -e '/clean.ba/d' -e '/HC.ba/d' -e '/preprocessing_report/d' -e '/redundant_mips/d'| sort | uniq);
