#-----------------------------------------------------#
# NEEDS ONE ARGUMENT : PIPELINE_HOME                  #
#-----------------------------------------------------#
if [[ $# -eq 1 ]]; then
  export PIPELINE_HOME=$1
else
  echo "ERROR INITIALIZING PIPELINE: PROVIDE VALID PIPELINE_HOME"
  exit 1
fi
#-----------------------------------------------------#
# System Defaults
umask 002
export PATH=${PIPELINE_HOME}/soft/bin/:$PATH
export PYTHONPATH=${PIPELINE_HOME}/soft/src/python/custom_modules:$PYTHONPATH
export JAVA_MEM=4

#-----------------------------------------------------#
# Analysis files and paths
export GATK_TORQUE_DIR=${PIPELINE_HOME}/soft/packages/gatk07_24_12_torque
#export GATK_DIR=${PIPELINE_HOME}/soft/packages/GATK-git-2.6-4
export GATK_DIR=${PIPELINE_HOME}/soft/packages/GenomeAnalysisTK-3.5
export PICARD_DIR=${PIPELINE_HOME}/soft/packages/picard
export BWA_DIR=${PIPELINE_HOME}/soft/src/bwa
export BWA_MEM_DIR=${PIPELINE_HOME}/soft/src/bwa-0.7.5a
export SAMTOOLS_1_15_DIR=${PIPELINE_HOME}/soft/src/samtools-0.1.15
export SAMTOOLS_DIR=${PIPELINE_HOME}/soft/src/samtools-0.1.18
export TRIMMOMATIC_JAR=${PIPELINE_HOME}/soft/packages/trimmomatic/trimmomatic-0.27.jar
export DINDEL_DIR=${PIPELINE_HOME}/soft/packages/dindel-1.01
#-----------------------------------------------------#
export VERIFYBAMID_GENOTYPES=${PIPELINE_HOME}/soft/packages/verifyBamID/data/Omni25_genotypes.1525_samples_v2.b37.PASS.ALL_sites.in_concat_38mb_50mb_V4.vcf.gz
export DBSNP=${PIPELINE_HOME}/data/reference/dbsnp_132.b37.vcf
export INDELS=${PIPELINE_HOME}/data/reference/1000G_indels_for_realignment.b37.vcf
export HAPMAP=${PIPELINE_HOME}/data/reference/hapmap_3.3.b37.sites.vcf
export OMNI=${PIPELINE_HOME}/data/reference/1000G_omni2.5.b37.sites.vcf
export MILLS_KG=${PIPELINE_HOME}/data/reference/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
export REF=${PIPELINE_HOME}/data/reference/human_g1k_v37.fasta
export REF_MT=${PIPELINE_HOME}/data/reference/human_g1k_v37_MT.fasta
export REF_MEM=${PIPELINE_HOME}/data/reference/human_g1k_v37.bwa-0_7_5a.fasta
export REF_MEM_MT=${PIPELINE_HOME}/data/reference/human_g1k_v37_MT.bwa-0_7_5a.fasta
export GENE_BED=${PIPELINE_HOME}/data/targets/RefSeq.Genes.v37.coding_exons.no_flanks.sort.bed
export GENE_LIST=${PIPELINE_HOME}/data/targets/RefSeq.Genes.v37.refGene
#-----------------------------------------------------#
export PERL5LIB=${PIPELINE_HOME}/soft/packages/perl/share/perl5
export CMHTOOLS=${PIPELINE_HOME}/soft/packages/cmhtools/
#-----------------------------------------------------#
# GATK Scala files
export MIP_POST_PROCESSING_SCALA=${PIPELINE_HOME}/data/templates/DataProcessingPipeline.noTorque.MIPs.fastqs.scala
export GATK_KEY=${PIPELINE_HOME}/data/templates/dan.spiegelman_crchum.qc.ca.key
#-----------------------------------------------------#
# GATK results files
export GATK_MIP_PREFIX=MIP
export GATK_MT_PREFIX=MT
export GATK_MIP_BAMFILE=${GATK_MIP_PREFIX}.${SAMPLE}.clean.bam
export GATK_MIP_TOTAL_RECALL_VCFFILE=${GATK_MIP_PREFIX}.${SAMPLE}.GATK_UG.vcf.gz
export GATK_MIP_TOTAL_RECALL_VARIANT_SITES_VCFFILE=${GATK_MIP_PREFIX}.${SAMPLE}.GATK_UG.variant_sites.vcf.gz
export GATK_MIP_HC_GVCFFILE=${GATK_MIP_PREFIX}.${SAMPLE}.GATK_HC.g.vcf.gz
export GATK_MIP_HC_VCFFILE=${GATK_MIP_PREFIX}.${SAMPLE}.GATK_HC.variant_sites.vcf.gz
export GATK_HC_BAMFILE=${GATK_MIP_PREFIX}.${SAMPLE}.GATK_HC.bam
export GATK_MT_BAMFILE=${GATK_MT_PREFIX}.${SAMPLE}.clean.dedup.recal.bam
export GATK_BAMFILE_MT_FILTERED=${GATK_MT_PREFIX}.${SAMPLE}.ref_MT.MT_L30.bam
#-----------------------------------------------------#
