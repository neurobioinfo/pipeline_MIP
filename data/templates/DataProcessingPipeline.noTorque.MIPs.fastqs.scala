package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.tools.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.gatk.utils.baq.BAQ.CalculationMode

import collection.JavaConversions._
import htsjdk.samtools._
//import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMFileHeader.SortOrder
//import htsjdk.samtools.SamReaderFactory

import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.function.ListWriterFunction
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.commandline

// ALEX ADD
//import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import org.broadinstitute.gatk.engine.phonehome.GATKRunReport

class DataProcessingPipeline extends QScript {
  qscript =>

  /****************************************************************************
  * Required Parameters
  ****************************************************************************/

  @Input(doc="input BAM file - or list of BAM files", fullName="input", shortName="i", required=true)
  var input: File = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = Seq()

// ALEX ADD
  @Argument(doc="Max processors/node", fullName="procs_per_node", shortName="ppn", required=true)
  var ppn: Int = _

// ALEX ADD
  @Argument(doc="Max Memory Limit", fullName="max_memory_limit", shortName="ml", required=true)
  var maxMemoryLimit: Int = _

  /****************************************************************************
  * Optional Parameters
  ****************************************************************************/

// ALEX ADD
  @Input(doc="Key for GATK ET", fullName="key", shortName="key", required=false)
  var key: File = ""

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: Seq[File] = Seq()

  @Input(doc="The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)", fullName="path_to_bwa", shortName="bwa", required=false)
  var bwaPath: File = "bwa"

  @Argument(doc="Run the alignment using bwa", fullName = "realign", shortName = "realign", required=false)
  var realign: Boolean = false

// ALEX remove "Project"
  @Argument(doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", fullName="project", shortName="p", required=false)
  var projectName: String = ""

  @Argument(doc="Output path for the processed BAM files.", fullName="output_directory", shortName="outputDir", required=false)
  var outputDir: String = ""

  @Argument(doc="the -L interval string to be used by GATK - output bams at interval only", fullName="gatk_interval_string", shortName="L", required=false)
  var intervalString: String = ""

  @Input(doc="an intervals file to be used by GATK - output bams at intervals only", fullName="gatk_interval_file", shortName="intervals", required=false)
  var intervals: File = _

  @Argument(doc="Cleaning model: KNOWNS_ONLY, USE_READS or USE_SW", fullName="clean_model", shortName="cm", required=false)
  var cleaningModel: String = "USE_READS"

  @Argument(doc="Perform validation on the BAM files", fullName="validation", shortName="vs", required=false)
  var validation: Boolean = false

  /****************************************************************************
  * Hidden Parameters
  ****************************************************************************/
// ALEX MOD nContigs to scatterThreads
  @Hidden
  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  //  var nContigs: Int = -1
  var scatterThreads: Int = 1

  @Hidden
  @Argument(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""

  @Hidden
  @Argument(doc="Run the pipeline in test mode only", fullName = "test_mode", shortName = "test", required=false)
  var testMode: Boolean = false


  /****************************************************************************
  * Global Variables
  ****************************************************************************/

  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output
  var cleanModelEnum: ConsensusDeterminationModel = ConsensusDeterminationModel.USE_READS

  /****************************************************************************
  * Helper classes and methods
  ****************************************************************************/

  class ReadGroup (val id: String,
                   val lb: String,
                   val pl: String,
                   val pu: String,
                   val sm: String,
                   val cn: String,
                   val ds: String)
  {}

  // Utility function to merge all bam files of similar samples. Generates one BAM file per sample.
  // It uses the sample information on the header of the input BAM files.
  //
  // Because the realignment only happens after these scripts are executed, in case you are using
  // bwa realignment, this function will operate over the original bam files and output over the
  // (to be realigned) bam files.
  def createSampleFiles(bamFiles: Seq[File], realignedBamFiles: Seq[File]): Map[String, Seq[File]] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, Seq[File]]
    val realignedIterator = realignedBamFiles.iterator
    for (bam <- bamFiles) {
      val rBam = realignedIterator.next()  // advance to next element in the realignedBam list so they're in sync.

      // val samReader = new SAMFileReader(bam)
      val samReader = SamReaderFactory.makeDefault().open(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if (!sampleTable.contains(sample))
          sampleTable(sample) = Seq(rBam)
        else if ( !sampleTable(sample).contains(rBam))
          sampleTable(sample) :+= rBam
      }
    }
    sampleTable.toMap
  }

  // Rebuilds the Read Group string to give BWA
  def addReadGroups(inBam: File, outBam: File, samReader: SamReader) {
    val readGroups = samReader.getFileHeader.getReadGroups
    var index: Int = readGroups.length
    for (rg <- readGroups) {
      val intermediateInBam: File = if (index == readGroups.length) { inBam } else { swapExt(outBam, ".bam", index+1 + "-rg.bam") }
      val intermediateOutBam: File = if (index > 1) {swapExt(outBam, ".bam", index + "-rg.bam") } else { outBam}
      val readGroup = new ReadGroup(rg.getReadGroupId, rg.getLibrary, rg.getPlatform, rg.getPlatformUnit, rg.getSample, rg.getSequencingCenter, rg.getDescription)
      add(addReadGroup(intermediateInBam, intermediateOutBam, readGroup))
      index = index - 1
    }
  }

  // Takes a list of processed BAM files and realign them using the BWA option requested  (bwase or bwape).
  // Returns a list of realigned BAM files.

  def performAlignment(bams: Seq[File]): Seq[File] = {
    var realignedBams: Seq[File] = Seq()
    var index = 1
    for (bam <- bams) {
      // ALEX add 'qscript.projectName' to the file names
      val saiFile1 = swapExt(bam, ".bam", "." + qscript.projectName + "." + index + ".1.sai")
      val saiFile2 = swapExt(bam, ".bam", "." + qscript.projectName + "." + index + ".2.sai")
      val realignedSamFile = swapExt(bam, ".bam", "." + qscript.projectName + "." + index + ".realigned.sam")
      val realignedSortedBamFile = swapExt(bam, ".bam", "." + qscript.projectName + "." + index + ".realigned.sorted.bam")
      val rgRealignedSortedBamFile = swapExt(bam, ".bam", "." + qscript.projectName + "." + index + ".realigned.sorted.rg.bam")

      // first revert the BAM file to the original qualities
      // ALEX we do not revert BAMs anymore
      // ALEX we use useBWApe by default
      add(bwa_aln_pe(bam, saiFile1, 1),
          bwa_aln_pe(bam, saiFile2, 2),
          bwa_sam_pe(bam, saiFile1, saiFile2, realignedSamFile))

      add(sortSam(realignedSamFile, realignedSortedBamFile, SortOrder.coordinate))
      addReadGroups(realignedSortedBamFile, rgRealignedSortedBamFile, SamReaderFactory.makeDefault().open(bam))
      realignedBams :+= rgRealignedSortedBamFile
      index = index + 1
    }
    realignedBams
  }

  def getIndelCleaningModel: ConsensusDeterminationModel = {
    if (cleaningModel == "KNOWNS_ONLY")
      ConsensusDeterminationModel.KNOWNS_ONLY
    else if (cleaningModel == "USE_SW")
      ConsensusDeterminationModel.USE_SW
    else
      ConsensusDeterminationModel.USE_READS
  }

  /****************************************************************************
  * Main script
  ****************************************************************************/
  def script() {
/** not needed anymore
    // final output list of processed bam files
    var cohortList: Seq[File] = Seq()
*/

    // sets the model for the Indel Realigner
    cleanModelEnum = getIndelCleaningModel

    // keep a record of the number of contigs in the first bam file in the list
    val bams = QScriptUtils.createSeqFromFile(input)
// ALEX from nContigs to scatterThreads
    if (scatterThreads < 0)
     scatterThreads = QScriptUtils.getNumberOfContigs(bams(0))

// ALEX no useBWApe option anymore
    // val realignedBAMs = if (useBWApe || useBWAse  || useBWAsw) {performAlignment(bams)} else {revertBams(bams, false)}
    var realignedBAMs = bams; // already aligned by default
    if ( realign ) { realignedBAMs = performAlignment(bams) }
//    else { realignedBAMs = bams }

    // generate a BAM file per sample joining all per lane files if necessary
    val sampleBAMFiles: Map[String, Seq[File]] = createSampleFiles(bams, realignedBAMs)

    // if this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir + projectName + ".intervals")
    if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY)
      add(target(null, globalIntervals))

    // put each sample through the pipeline
    for ((sample, bamList) <- sampleBAMFiles) {

      // BAM files generated by the pipeline
      val bam        = new File(qscript.projectName + "." + sample + ".bam")
      val cleanedBam = swapExt(bam, ".bam", ".clean.bam")
//      val dedupedBam = swapExt(bam, ".bam", ".clean.dedup.bam")
//      val recalBam   = swapExt(bam, ".bam", ".clean.recal.bam")

      // Accessory files
      val targetIntervals = if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
// ALEX 'MarkDuplicates' to the metrics file
      val metricsFile     = swapExt(bam, ".bam", ".MarkDuplicates.metrics")
/* REMOVED ALL STEPS RELATED TO RECAL, SINCE WE DON'T WANT TO DO THIS FOR MIPS     val preRecalFile    = swapExt(bam, ".bam", ".pre_recal.table")
      val postRecalFile   = swapExt(bam, ".bam", ".post_recal.table")
      val preOutPath      = swapExt(bam, ".bam", ".pre")
      val postOutPath     = swapExt(bam, ".bam", ".post")
      val preValidateLog  = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")

      // Validation is an optional step for the BAM file generated after
      // alignment and the final bam file of the pipeline.
      if (validation) {
        for (sampleFile <- bamList)
        add(validate(sampleFile, preValidateLog),
            validate(recalBam, postValidateLog))
      }
*/

      if (cleaningModel != ConsensusDeterminationModel.KNOWNS_ONLY)
        add(target(bamList, targetIntervals))

      add(clean(bamList, targetIntervals, cleanedBam))
//          dedup(cleanedBam, dedupedBam, metricsFile),
//          cov(cleanedBam, preRecalFile),
//          recal(cleanedBam, preRecalFile, recalBam),
//          cov(recalBam, postRecalFile))

//      cohortList :+= recalBam
    }

/** Not needed anymore 
    // output a BAM list with all the processed per sample files
    val cohortFile = new File(qscript.outputDir + qscript.projectName + ".cohort.list")
    add(writeList(cohortList, cohortFile))
*/
  }



  /****************************************************************************
  * Classes (GATK Walkers)
  ****************************************************************************/
  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
// ALEX memory from 4 to maxMemoryLimit
    this.memoryLimit = qscript.maxMemoryLimit
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
// ALEX ADD the key file to disable phone home
    if ( qscript.key != null ) {
      this.et = GATKRunReport.PhoneHomeOption.NO_ET // Disable phone home
      this.K = qscript.key
    }
  }

  trait PicardBamFunctionArgs extends PicardBamFunction with ExternalCommonArgs {
// ALEX REMOVE and replace
//      this.maxRecordsInRam = 100000
    this.maxRecordsInRam = (qscript.maxMemoryLimit * 250000 * 0.8).toInt                // rule of thumb for 100bp reads: 250,000 records / GB assigned to java
                                                                                        // (with -Xmx, equivalent to ExternalCommonArgs.memoryLimit);
                                                                                        // we add a 20% overhead tax to be sure not to exceed the memory limit
  }

  case class target (inBams: Seq[File], outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (cleanModelEnum != ConsensusDeterminationModel.KNOWNS_ONLY)
      this.input_file = inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= qscript.dbSNP
    if (indels != null)
      this.known ++= qscript.indels
    this.scatterCount = scatterThreads
    this.nt=qscript.ppn
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class clean (inBams: Seq[File], tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known ++= qscript.dbSNP
    if (qscript.indels != null)
      this.known ++= qscript.indels
    this.consensusDeterminationModel = cleanModelEnum
    this.compress = 9
    this.isIntermediate = false
    this.noPGTag = qscript.testMode;
    this.scatterCount = scatterThreads
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class cov (inBam: File, outRecalFile: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.knownSites ++= qscript.dbSNP
    this.covariate ++= Seq("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "ContextCovariate")
    this.input_file :+= inBam
    this.disable_indel_quals = true
    this.out = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    if (!qscript.intervalString.isEmpty) this.intervalsString ++= Seq(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.scatterCount = scatterThreads
    this.nct=qscript.ppn
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.BQSR = inRecalFile
// ALEX REMOVE baq
    //  this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.baq = CalculationMode.OFF
//    this.solid_nocall_strategy = SOLID_NOCALL_STRATEGY.PURGE_READ
    this.out = outBam
    if (!qscript.intervalString.isEmpty) this.intervalsString ++= Seq(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.scatterCount = scatterThreads
    this.isIntermediate = false
// ALEX produce the OQ flag in bams
    this.emit_original_quals=true
// ALEX compress to max
    this.compress = 9
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }


  /****************************************************************************
   * Classes (non-GATK programs)
   ****************************************************************************/

// ALEX for all picard functions, modify the ExternalCommonArgs to PicardBamFunctionArgs

/*
  case class dedup (inBam: File, outBam: File, metricsFile: File) extends MarkDuplicates with PicardBamFunctionArgs {
    this.input :+= inBam
    this.output = outBam
    this.metrics = metricsFile
// ALEX assumeSorted = true
    this.assumeSorted = true
// ALEX memoryLimit set by ExternalCommonArgs
    this.memoryLimit = 16
    this.analysisName = queueLogDir + outBam + ".dedup"
    this.jobName = queueLogDir + outBam + ".dedup"
    override def commandLine = super.commandLine + " 'PROGRAM_RECORD_ID=null'"
  }
*/
/*
  case class joinBams (inBams: Seq[File], outBam: File) extends MergeSamFiles with PicardBamFunctionArgs {
    this.input = inBams
    this.output = outBam
// ALEX this.assumeSorted = true
    this.assumeSorted = true
    this.analysisName = queueLogDir + outBam + ".joinBams"
    this.jobName = queueLogDir + outBam + ".joinBams"
  }
*/

  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam with PicardBamFunctionArgs {
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class validate (inBam: File, outLog: File) extends ValidateSamFile with PicardBamFunctionArgs {
    this.input :+= inBam
    this.output = outLog
// ALEX this.assumeSorted = true
    this.assumeSorted = true
    this.REFERENCE_SEQUENCE = qscript.reference
    this.isIntermediate = false
    this.analysisName = queueLogDir + outLog + ".validate"
    this.jobName = queueLogDir + outLog + ".validate"
  }

  case class addReadGroup (inBam: File, outBam: File, readGroup: ReadGroup) extends AddOrReplaceReadGroups with PicardBamFunctionArgs {
    this.input :+= inBam
    this.output = outBam
    this.RGID = readGroup.id
    this.RGCN = readGroup.cn
    this.RGDS = readGroup.ds
    this.RGLB = readGroup.lb
    this.RGPL = readGroup.pl
// ALEX modify the RGPU
//    this.RGPU = readGroup.pu
    this.RGPU = readGroup.id + "_" + readGroup.sm
    this.RGSM = readGroup.sm
    this.analysisName = queueLogDir + outBam + ".rg"
    this.jobName = queueLogDir + outBam + ".rg"
  }

  case class revert (inBam: File, outBam: File, removeAlignmentInfo: Boolean) extends RevertSam with PicardBamFunctionArgs {
    this.output = outBam
    this.input :+= inBam
// ALEX this.assumeSorted = true
    this.assumeSorted = true
    this.removeAlignmentInformation = removeAlignmentInfo;
// ALEX for bioscope 
//    this.restoreOriginalQualities=false
    this.sortOrder = if (removeAlignmentInfo) {SortOrder.queryname} else {SortOrder.coordinate}
    this.analysisName = queueLogDir + outBam + "revert"
    this.jobName = queueLogDir + outBam + ".revert"
  }

/*
  case class convertToFastQ (inBam: File, outFQ: File) extends SamToFastq with PicardBamFunctionArgs {
    this.input :+= inBam
    this.fastq = outFQ
    this.analysisName = queueLogDir + outFQ + "convert_to_fastq"
    this.jobName = queueLogDir + outFQ + ".convert_to_fastq"
  }
*/

  case class bwa_aln_pe (inBam: File, outSai1: File, index: Int) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for 1st mating pair") var sai = outSai1
// ALEX bwaThreads -> ppn, added '-R 1000'
    def commandLine = bwaPath + " aln -t " + ppn + " -q 5 -R 1000 " + reference + " -b" + index + " " + bam + " > " + sai
    this.analysisName = queueLogDir + outSai1 + ".bwa_aln_pe1"
    this.jobName = queueLogDir + outSai1 + ".bwa_aln_pe1"
  }

  case class bwa_sam_pe (inBam: File, inSai1: File, inSai2:File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc="bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc="output aligned bam file") var alignedBam = outBam
// ALEX add ppn, add -a 700
    def commandLine = bwaPath + " sampe -t " + ppn + " -a 700 " + reference + " " + sai1 + " " + sai2 + " " + bam + " " + bam + " > " + alignedBam
// ALEX this.memoryLimit 6 -> qscript.maxMemoryLimit
    this.memoryLimit = qscript.maxMemoryLimit
    this.analysisName = queueLogDir + outBam + ".bwa_sam_pe"
    this.jobName = queueLogDir + outBam + ".bwa_sam_pe"
  }

/** Not needed anymore
  case class writeList(inBams: Seq[File], outBamList: File) extends ListWriterFunction {
    this.inputFiles = inBams
    this.listFile = outBamList
    this.analysisName = queueLogDir + outBamList + ".bamList"
    this.jobName = queueLogDir + outBamList + ".bamList"
  }
*/
}
