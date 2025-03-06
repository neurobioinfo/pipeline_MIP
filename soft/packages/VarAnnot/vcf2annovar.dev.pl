#!/usr/bin/perl

use File::Basename;
use File::Path;
use warnings;
use strict;
use POSIX ":sys_wait_h";
use Pod::Usage;
use Getopt::Long;
BEGIN { unshift @INC, '/RQexec/dionnela/soft/packages/VarAnnot'; }
use Vcf;

# Command line parameters
my ($VERBOSE, $HELP, $MAN);
my ($INFILE, $ID, $BUILDVER, $VARIANT_CALLER, $HEADER, $DBSNP, $EXCLUDE, $INCLUDE, $SERIAL, $CPU);

# Global variables
my $ANNOVARPATH = "/RQexec/dionnela/soft/packages/VarAnnot";
my $BODYMAPPATH = "$ANNOVARPATH/humandb/IlluminaBodyMap2";
my $GENOMETRAXPATH = "$ANNOVARPATH/humandb/GenomeTrax";
my $SMPINPUTFILE;               # path of the smp-runner input file (set if parallel mode is used)
my $ANNOINFILE;                 # path of the annovar input file
my $OUTFILE;                    # path of the output file
my $OUTFILE_ANNO;               # path of the annovar --geneanno output file
my $OUTFILE_DBSNP;              # path of the annovar --filter dbSNP output file
my $OUTFILE_1000G;              # path of the annovar --filter 1000g output file
my $OUTFILE_TSI;                # path of the annovar --filter 1000g toscani allele frequency output file
my $OUTFILE_ASN;                # path of the annovar --filter 1000g east asian allele frequency output file
my $OUTFILE_AFR;                # path of the annovar --filter 1000g african allele frequency output file
my $OUTFILE_AMR;                # path of the annovar --filter 1000g ad mixed american allele frequency output file
my $OUTFILE_EUR;                # path of the annovar --filter 1000g european allele frequency output file
my $OUTFILE_CG;                 # path of the annovar --filter Complete Genomics output file
my $OUTFILE_SIFT;               # path of the annovar --filter avsift output file
my $OUTFILE_PP2;                # path of the annovar --filter ljb_pp2 output file
my $OUTFILE_MT;                 # path of the annovar --filter ljb_mt output file
my $OUTFILE_LRT;                # path of the annovar --filter ljb_lrt output file
my $OUTFILE_PHYLOP;             # path of the annovar --filter ljb_phylop output file
my $OUTFILE_GERP;               # path of the annovar --filter ljb_gerp++ output file
my $OUTFILE_ESP6500;            # path of the annovar --filter esp6500 output file
my $OUTFILE_MIRNA;              # path of the annovar --regionanno mirna output file
my $OUTFILE_MIRNATARGET;        # path of the annovar --regionanno mirnatarget output file
my $OUTFILE_RNASEQ_CUFF;        # path of the annovar --regionanno rnaseq_cuff output file
my $OUTFILE_RNASEQ_SCRI;        # path of the annovar --regionanno rnaseq_scri output file
my $OUTFILE_GTX_CHIPSEQ;        # path of the annovar --regionanno gtx_chipseq output file
my $OUTFILE_GTX_CLINVAR;        # path of the annovar --regionanno gtx_clinvar output file
my $OUTFILE_GTX_COMMON_SNP;     # path of the annovar --regionanno gtx_common_snp output file
my $OUTFILE_GTX_COSMIC;         # path of the annovar --regionanno gtx_cosmic output file
my $OUTFILE_GTX_CPG;            # path of the annovar --regionanno gtx_cpg output file
my $OUTFILE_GTX_DBNSFP;         # path of the annovar --regionanno gtx_dbnsfp output file
my $OUTFILE_GTX_DISEASE;        # path of the annovar --regionanno gtx_disease output file
my $OUTFILE_GTX_DRUG;           # path of the annovar --regionanno gtx_drug output file
my $OUTFILE_GTX_EVS;            # path of the annovar --regionanno gtx_evs output file
my $OUTFILE_GTX_GWAS;           # path of the annovar --regionanno gtx_gwas output file
my $OUTFILE_GTX_HGMD;           # path of the annovar --regionanno gtx_hgmd output file
my $OUTFILE_GTX_HGMD_DISGENE;   # path of the annovar --regionanno gtx_hgmd_disgene output file
my $OUTFILE_GTX_MICROSAT;       # path of the annovar --regionanno gtx_microsat output file
my $OUTFILE_GTX_MIRNA;          # path of the annovar --regionanno gtx_mirna output file
my $OUTFILE_GTX_OMIM;           # path of the annovar --regionanno gtx_omim output file
my $OUTFILE_GTX_PATH;           # path of the annovar --regionanno gtx_path output file
my $OUTFILE_GTX_PGX;            # path of the annovar --regionanno gtx_pgx output file
my $OUTFILE_GTX_PTMS;           # path of the annovar --regionanno gtx_ptms output file
my $OUTFILE_GTX_TRANSFAC;       # path of the annovar --regionanno gtx_transfac output file
my $OUTFILE_GTX_TSS;            # path of the annovar --regionanno gtx_tss output file
my $OUTPATH;                    # path of the results output file (if not provided, will be interpreted from $OUTFILE, or even $INFILE if $OUTFILE is not provided)
my $SUBOUTPATH;                 # path of the all the annovar output files (if not provided, will be interpreted from $OUTFILE, or even $INFILE if $OUTFILE is not provided)
my $USEDBUILD;                  # build to be use during the annovar analysis (can be different than $BUILDVER, which will bw the one printed in name of output files)
my @EXCLUDETAB = ();
my @INCLUDETAB = ();

&main;

sub main {
    # First process the provided arguments
    processArguments();

    # Parse input data into hash
    $VERBOSE and printerr("\n**********\n* Starting parsing $VARIANT_CALLER input\n**********\n\n");
    my ($rH_data, $rH_vcfVarIndex) = convertVcf();
    $VERBOSE and printerr("\n**********\n* Parsing $VARIANT_CALLER input done...\n**********\n\n");

    my $rO_vcf = Vcf->new(file => $INFILE);
    $rO_vcf->parse_header();

    if ($SERIAL) {
        # Proceedto the sequential (-geneanno, dbsnp, 1000g...) executions of Annovar
        launchAnnovarSerial();
    } else {
        # Proceed to the parallel (-geneanno, dbsnp, 1000g...) executions of Annovar
        launchAnnovarParallel();
    }
    
    # Now use all Annovar output files to build an object with all the fetched annotations
    $VERBOSE and printerr("\n**********\n* Starting building VCF object\n**********\n\n");
    buildVcf(
        -data => $rH_data,
        -vcf  => $rO_vcf
    );
    $VERBOSE and printerr("\n**********\n* Building VCF object done...\n**********\n\n");

    # Print the previously built object in a vcf file
    $VERBOSE and printerr("\n**********\n* Starting printing VCF output\n**********\n\n");
    printVcf(
        -data  => $rH_data,
        -vcf   => $rO_vcf,
        -index => $rH_vcfVarIndex
    );
    $VERBOSE and printerr("\n**********\n* Printing VCF output done...\n**********\n\n");

    print "\n\nAnnotation of $VARIANT_CALLER $BUILDVER variants for $ID is a great success !!\nThank you for using this program, see you next time !! :)\n\n\n";
}

sub processArguments {
    my @command_line = @ARGV;       # command line argument
    GetOptions('verbose|v'=>\$VERBOSE, 'help'=>\$HELP, 'man|m'=>\$MAN,
               'infile|i=s'=>\$INFILE, 'outfile|o=s'=>\$OUTFILE, 'identifier|id=s'=>\$ID, 'buildver|b=s'=>\$BUILDVER, 'vcaller|vc=s'=>\$VARIANT_CALLER,
               'dbsnp|d=s'=>\$DBSNP, 'exclude|e=s'=>\$EXCLUDE, 'include|inc=s'=>\$INCLUDE,
               'cpu=s'=>\$CPU, 'serial'=>\$SERIAL, 'header|h=s'=>\$HEADER) || pod2usage();

    $HELP and pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDOUT);
    $MAN  and pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDOUT);
    @ARGV == 0 || pod2usage("\n**************************************************************\n Syntax error\n\n\n");

    if ($INCLUDE && $EXCLUDE && ($EXCLUDE !~ /^all$/i && $INCLUDE !~ /^all$/)) {
        pod2usage("\n**************************************************************\nError in argument: you can only use both --include and --exclude arguments if one is set to 'all'!!\n\n\n");
    }

    if (not $ID) {
        pod2usage("\n**************************************************************\nError in argument: the --id has to be provided !!\n\n\n");
    }
    $ID = ucfirst $ID;

    if (($VARIANT_CALLER) && ($VARIANT_CALLER !~ /^(varscan|dibayes|smallindel|smallindels|mpileup|gatk|dindel|custom)$/i)) {
        pod2usage("\n**************************************************************\nError in argument: the --vcaller must be one of these : varscan, dibayes, smllindels, mpileup, gatk, dindel, custom\n\n\n");
    }
    $VARIANT_CALLER = lc $VARIANT_CALLER if ($VARIANT_CALLER);

    if ((not $BUILDVER) || ($BUILDVER !~ /^(hg18|18|hg19|19|b37|v37)$/i)) {
        pod2usage("\n**************************************************************\nError in argument: the --buildver must be one of these : hg18 (or 18), hg19 (or 19), b37, v37\n\n\n");
    } elsif (($BUILDVER !~ /^hg/i) && ($BUILDVER !~ /b37|v37/i)) {
        $BUILDVER = "hg" . $BUILDVER;
    }
    $BUILDVER = lc $BUILDVER;

    if ($DBSNP) {
        ($DBSNP !~ /[snp130|130|snp131|131|snp132|132|snp135|135]/i) and pod2usage("\n**************************************************************\nError in argument: the --dbsnp must be one of these : snp130(or 130), snp131 (or 131), snp132 (or 132), snp135 (or 135)\n\n\n");
        ($DBSNP !~ /^snp/) and $DBSNP = "snp" . $DBSNP;
    }

    if ($EXCLUDE) {
        ($EXCLUDE !~ /^((all|gene|snp|1kg|cg|sift|pp2|mt|lrt|phylop|gerp|esp6500|mirna|mirnatarget|tsi|asn|afr|amr|eur|rnaseq_cuff|rnaseq_scri|gtx_chipseq|gtx_common_snp|gtx_cosmic|gtx_cpg|gtx_dbnsfp|gtx_disease|gtx_drug|gtx_evs|gtx_gwas|gtx_hgmd|gtx_hgmd_disgene|gtx_microsat|gtx_mirna|gtx_omim|gtx_path|gtx_pgx|gtx_ptms|gtx_transfac|gtx_tss),?)+$/i) and pod2usage("\n**************************************************************\nError in argument: the --exclude must be one of these : gene,snp,1kg,cg,sift,pp2,mt,lrt,phylop,gerp,esp6500,mirna,mirnatarget,tsi,asn,afr,amr,eur,rnaseq_cuff,rnaseq_scri,gtx_chipseq,gtx_common_snp,gtx_cosmic,gtx_cpg,gtx_dbnsfp,gtx_disease,gtx_drug,gtx_evs,gtx_gwas,gtx_hgmd,gtx_hgmd_disgene,gtx_microsat,gtx_mirna,gtx_omim,gtx_path,gtx_pgx,gtx_ptms,gtx_transfac,gtx_tss\n\n\n");
        @EXCLUDETAB = split /,/, $EXCLUDE; # / just to avoid a bug in eclipse display...
    }

    if ($INCLUDE) {
        ($INCLUDE !~ /^((all|gene|snp|1kg|cg|sift|pp2|mt|lrt|phylop|gerp|esp6500|mirna|mirnatarget|tsi|asn|afr|amr|eur|rnaseq_cuff|rnaseq_scri|gtx_chipseq|gtx_common_snp|gtx_cosmic|gtx_cpg|gtx_dbnsfp|gtx_disease|gtx_drug|gtx_evs|gtx_gwas|gtx_hgmd|gtx_hgmd_disgene|gtx_microsat|gtx_mirna|gtx_omim|gtx_path|gtx_pgx|gtx_ptms|gtx_transfac|gtx_tss),?)+$/i) and pod2usage("\n**************************************************************\nError in argument: the --include must be one of these : gene,snp,1kg,cg,sift,pp2,mt,lrt,phylop,gerp,esp6500,mirna,mirnatarget,tsi,asn,afr,amr,eur,rnaseq_cuff,rnaseq_scri,gtx_chipseq,gtx_common_snp,gtx_cosmic,gtx_cpg,gtx_dbnsfp,gtx_disease,gtx_drug,gtx_evs,gtx_gwas,gtx_hgmd,gtx_hgmd_disgene,gtx_microsat,gtx_mirna,gtx_omim,gtx_path,gtx_pgx,gtx_ptms,gtx_transfac,gtx_tss\n\n\n");
        @INCLUDETAB = split /,/, $INCLUDE; # / just to avoid a bug in eclipse display...
    }

    $VERBOSE ||= "";
    $VERBOSE &&= "--verbose";

    my $headerNotice = undef;
    if ($HEADER) {
        $headerNotice = "NOTICE : --header argument is obsolete and no longer used ; it will hence be ignored...\n";
#        (not -e $HEADER) and $HEADER = undef;
#        (not -e $HEADER) and $headerNotice = "NOTICE : the provided --header $HEADER does not exists... no header file will be used during analysis...\n";
#        (-d $HEADER) and $HEADER = undef;
#        (-d $HEADER) and $headerNotice = "NOTICE : the provided --header is a directory... since it has to be a txt file, no header file will be used during analysis...\n";
#    } else {
#        $headerNotice = "NOTICE : no --header argument provided : no header file will be used during analysis...\n";
    }

    if (not $INFILE) {
        pod2usage("\n**************************************************************\nError in argument: the --infile has to be provided !!\n\n\n");
    } elsif (!(-e $INFILE)) {
        pod2usage("\n**************************************************************\nError in argument: the provided --infile $INFILE does not exists !!\n\n\n");
    }

    my $outFileNotice = undef;
    if (not $OUTFILE) {
        my ($directory, $filename) = $INFILE =~ m/(.*\/)(.*)$/;
        (not $directory) and $directory = "./";
        $OUTPATH = $directory;
    } else {
        my ($directory, $filename) = $OUTFILE =~ m/(.*\/)(.*)$/;
        (not $directory) and $directory = "./";
        if ($filename && -d $filename) {
            $directory .= $filename . "/";
            $filename = undef;
        }
        $OUTPATH = $directory;
        (not $filename) and $OUTFILE = "";
        if (not -e $OUTPATH) {
            my $print = "path $OUTPATH does not exists...";
            my ($directory, $filename) = $INFILE =~ m/(.*\/)(.*)$/;
            (not $directory) and $directory = "./";
            $OUTPATH = $directory;
            $outFileNotice = "NOTICE : $print $OUTPATH will be used instead...\n";
        }
    }
    if (not $OUTPATH) {
        my ($directory, $filename) = $INFILE =~ m/(.*\/)(.*)$/;
        (not $directory) and $directory = "./";
        $OUTPATH = $directory;
    }
    $OUTPATH .= "/" if ($OUTPATH !~ /\/$/);

    my $date = GetCurrentTime();
    $SUBOUTPATH = $OUTPATH . $ID . "_ANNOVAR_analysis_files_" . $date . "/";
    mkpath $SUBOUTPATH, 002;

    $ANNOINFILE                 = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_input." . $VARIANT_CALLER . "." . $BUILDVER : $SUBOUTPATH . $ID . "_ANNOVAR_input." . $BUILDVER;
    $OUTFILE_ANNO               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_geneanno_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_geneanno_output";
    $OUTFILE_DBSNP              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp_output";
    $OUTFILE_1000G              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_output";
    $OUTFILE_TSI                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_tsi_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_tsi_output";
    $OUTFILE_ASN                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_asn_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_asn_output";
    $OUTFILE_AFR                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_afr_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_afr_output";
    $OUTFILE_AMR                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_amr_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_amr_output";
    $OUTFILE_EUR                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_1000g_eur_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_1000g_eur_output";
    $OUTFILE_CG                 = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_cg_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_cg_output";
    $OUTFILE_SIFT               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_sift_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_sift_output";
    $OUTFILE_PP2                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_pp2_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_pp2_output";
    $OUTFILE_MT                 = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_mt_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_mt_output";
    $OUTFILE_LRT                = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_lrt_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_lrt_output";
    $OUTFILE_PHYLOP             = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_phylop_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_phylop_output";
    $OUTFILE_GERP               = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_gerp++_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_gerp++_output";
    $OUTFILE_ESP6500            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_esp6500_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_esp6500_output";
    $OUTFILE_MIRNA              = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_mirna_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_mirna_output";
    $OUTFILE_MIRNATARGET        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_mirnatarget_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_mirnatarget_output";
    $OUTFILE_RNASEQ_CUFF        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Cufflinks_allTissues_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Cufflinks_allTissues_output";
    $OUTFILE_RNASEQ_SCRI        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Scripture_allTissues_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_rnaseq_Scripture_allTissues_output";
    $OUTFILE_GTX_CHIPSEQ        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_chipseq_sites_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_chipseq_sites_output";
    $OUTFILE_GTX_CLINVAR        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_clinvar_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_clinvar_output";
    $OUTFILE_GTX_COMMON_SNP     = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_common_snp_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_common_snp_output";
    $OUTFILE_GTX_COSMIC         = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_cosmic_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_cosmic_output";
    $OUTFILE_GTX_CPG            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_cpg_islands_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_cpg_islands_output";
    $OUTFILE_GTX_DBNSFP         = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_dbnsfp_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_dbnsfp_output";
    $OUTFILE_GTX_DISEASE        = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_disease_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_disease_output";
    $OUTFILE_GTX_DRUG           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_drug_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_drug_output";
    $OUTFILE_GTX_EVS            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_evs_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_evm_output";
    $OUTFILE_GTX_GWAS           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_gwas_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_gwas_output";
    $OUTFILE_GTX_HGMD           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_output";
    $OUTFILE_GTX_HGMD_DISGENE   = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_disease_genes_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_hgmd_disease_genes_output";
    $OUTFILE_GTX_MICROSAT       = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_microsatellites_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_microsatellites_output";
    $OUTFILE_GTX_MIRNA          = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_mirna_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_mirna_output";
    $OUTFILE_GTX_OMIM           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_omim_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_omim_output";
    $OUTFILE_GTX_PATH           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pathway_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pathway_output";
    $OUTFILE_GTX_PGX            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pgx_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_pgx_output";
    $OUTFILE_GTX_PTMS           = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_ptms_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_ptms_output";
    $OUTFILE_GTX_TRANSFAC       = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_transfac_sites_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_transfac_sites_output";
    $OUTFILE_GTX_TSS            = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_tss_output." . $VARIANT_CALLER : $SUBOUTPATH . $ID . "_ANNOVAR_GenomeTrax_tss_output";

    my $logFile = $SUBOUTPATH . $ID . "_ANNOVAR_unresolved.log";
    open(LOG, ">$logFile") || die "Could not open log file $logFile...\n($!)\n";

    if ($SERIAL) {
        $VERBOSE and printerr("NOTICE : using SERIAL mode; ANNOVAR calls will be made in serial, one database at a time...\n");
        $CPU and printerr("NOTICE : --cpu argument will be ignored using SERIAL mode; ANNOVAR calls will be made in serial, one database at a time...\n");
    } else {
        $VERBOSE and printerr("NOTICE : using PARALLEL mode; ANNOVAR calls will be made in parallel on 12 CPUs with smp-runner...\n");
        $SMPINPUTFILE = ($VARIANT_CALLER) ? $SUBOUTPATH . $ID . "_SMP-RUNNER_input." . $VARIANT_CALLER . "." . $BUILDVER : $SUBOUTPATH . $ID . "_SMP-RUNNER_input" . $BUILDVER;
    }

    ($headerNotice  && $VERBOSE) && printerr($headerNotice);
    ($outFileNotice && $VERBOSE) && printerr($outFileNotice);
}

sub convertVcf {

    # Execute convert2annovar to prepare the annovar input file
    $VERBOSE and printerr("\n**********\n* Starting converting vcf input to annovar input file sing the following command line :\n* perl $ANNOVARPATH/convert2annovar.pl --format vcf4 --allallele --includeinfo --outfile $ANNOINFILE $VERBOSE $INFILE\n**********\n\n");
    system "perl $ANNOVARPATH/convert2annovar.pl --format vcf4 --allallele --includeinfo --outfile $ANNOINFILE $VERBOSE $INFILE";
    $VERBOSE and printerr("\n**********\n* Input conversion done...\n**********\n\n");

    open(ANNOINFILE, "<$ANNOINFILE") || die "cannot open file $ANNOINFILE:\n$!\n";

    my %vcfVariants;
    my %input;
    my @annoinfile;

    # Parse converted file to prepare VariantClass (VC), build a hash (%vcfVariants) to track back variants from original vcf input file
    # and also remove some triling columns to make it lighter.
    while (my $inputLine = <ANNOINFILE>) {
        chomp $inputLine;
        my ($chrom, $start, $stop, $ref, $var, $vcfChrom, $vcrPos, $vcfJunk, $vcrRef, $vcfVar, @junk) = split /\t+/, $inputLine;    # / just to avoid a bug in eclipse display...

        my $variantKey = $chrom . ":" . $start . "-" . $stop . "_" . $ref . "/" . $var; # " just to avoid a bug in eclipse display...

        # MULTIPLE CASES if start=stop then
        if ($start == $stop) {

            # INSERTION
            if ($ref eq '-') {
                $input{$variantKey}{VC} = 'insertion';

            # DELETION
            } elsif ($var eq '-') {
                $input{$variantKey}{VC} = 'deletion';

            # SUBSTITUTION
            } elsif (length ($var) > 1) {
                $input{$variantKey}{VC} = 'substitution';

            # SNP
            } else {
                $input{$variantKey}{VC} = 'SNP';
            }

        # DELETION if var='-'
        } elsif ($var eq '-') {
            $input{$variantKey}{VC} = 'deletion';

        # if SUBSTITUTION
        } else {
            $input{$variantKey}{VC} = 'substitution';
        }

        push @{$vcfVariants{$vcfChrom . ":" . $vcrPos . "_" . $vcrRef . "/" . $vcfVar}}, $variantKey;
        push @annoinfile, $chrom . "\t" . $start . "\t" . $stop . "\t" . $ref . "\t" . $var . "\t" . $vcfChrom . "\t" . $vcrPos . "\t" . $vcrRef . "\t" . $vcfVar . "\n";
    }
    close(ANNOINFILE);

    open(ANNOINFILE, ">$ANNOINFILE") || die "cannot open file $ANNOINFILE:\n$!\n";
    print ANNOINFILE @annoinfile;
    close(ANNOINFILE);

    return (\%input, \%vcfVariants);
}

sub launchAnnovarSerial {
    my %arg = @_;
    my $rO_vcf  = $arg{-vcf};
    my $rH_data = $arg{-data};

    $OUTFILE_ANNO .= "." . $BUILDVER;

    $USEDBUILD = $BUILDVER;
    ($USEDBUILD =~ /b37|v37/) and $USEDBUILD = "hg19";

    if (doIt(-elem => "gene")) {
        # run ANNOVAR -geneanno, first pass
        $VERBOSE and printerr("\n**********\n* Starting genome annotation using the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -geneanno -buildver $USEDBUILD -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $ANNOVARPATH/humandb/ -separate -exonicsplicing -transcript_function $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -geneanno -buildver $USEDBUILD -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $ANNOVARPATH/humandb/ -separate -exonicsplicing -transcript_function $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Genome annotation done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output dbSNP variants
    if (doIt(-elem => "snp")) {
        my $dbsnp;
        if ($DBSNP) {
            $dbsnp = $DBSNP;
        } elsif ($USEDBUILD eq "hg18") {
            $dbsnp = "snp130";
        } elsif ($USEDBUILD eq "hg19") {
            $dbsnp = "snp131";
        }

        # First run annovar with clinically associated dbSNP
        $dbsnp .= "_clinic";
        $VERBOSE and printerr("\n**********\n* Starting filtering against dbSNP with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system("perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE");
        $VERBOSE and printerr("\n**********\n* Filtering against dbSNP done...\n**********\n\n");

        # Then run annovar with non-clinically associated dbSNP
        $dbsnp =~ s/_clinic$/_non_clinic/;
        $VERBOSE and printerr("\n**********\n* Starting filtering against dbSNP with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system("perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE");
        $VERBOSE and printerr("\n**********\n* Filtering against dbSNP done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg")) {
        my $dbtype;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $dbtype = "1000g2010jul_ceu";
            $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_CEU.sites.2010_07_dropped";
        } elsif ($USEDBUILD eq "hg19") {
            $dbtype = "1000g2012feb_all";
            $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_ALL.sites.2012_02_dropped";
        }
        $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -outfile $OUTFILE_1000G $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_1000G $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes db done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "tsi")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_TSI_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_TSI.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'tsi' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_TSI_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_TSI.$USEDBUILD\_generic_dropped";

            $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_TSI $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_TSI $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes Toscani frequencies db done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "asn")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_ASN_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_ASN.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'asn' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_ASN_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_ASN.$USEDBUILD\_generic_dropped";

            $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ASN $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ASN $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes Toscani frequencies db done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "afr")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_AFR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AFR.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'afr' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_AFR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AFR.$USEDBUILD\_generic_dropped";

            $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AFR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AFR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes Toscani frequencies db done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "amr")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_AMR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AMR.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'amr' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_AMR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AMR.$USEDBUILD\_generic_dropped";

            $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AMR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AMR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes Toscani frequencies db done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "eur")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_EUR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_EUR.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'eur' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_EUR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_EUR.$USEDBUILD\_generic_dropped";

            $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EUR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EUR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes Toscani frequencies db done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -filter to output Complete Genomics variants
    if (doIt(-elem => "cg")) {
        my $dbfile = $USEDBUILD . "_cg69.txt";
        $VERBOSE and printerr("\n**********\n* Starting filtering against Complete Genomics db with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -outfile $OUTFILE_CG $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_CG $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Filtering against Complete Genomics db done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output SIFT scores
    if (doIt(-elem => "sift")) {
        $VERBOSE and printerr("\n**********\n* Starting fetching of SIFT prediction scores with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype avsift -sift_threshold 0 -buildver $USEDBUILD -outfile $OUTFILE_SIFT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype avsift -sift_threshold 0 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_SIFT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of SIFT prediction scores done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output PolyPhen.v2 scores
    if (doIt(-elem => "pp2")) {
        $VERBOSE and printerr("\n**********\n* Starting fetching of PolyPhen (v2) prediction scores with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_pp2 -buildver $USEDBUILD -outfile $OUTFILE_PP2 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_pp2 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PP2 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of PolyPhen (v2) prediction scores done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output LRT scores
    if (doIt(-elem => "lrt")) {
        $VERBOSE and printerr("\n**********\n* Starting fetching of Likelyhood-Ratio-Test scores with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_lrt -buildver $USEDBUILD -outfile $OUTFILE_LRT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_lrt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_LRT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of Likelyhood-Ratio-Test scores done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output Mutation Taster scores
    if (doIt(-elem => "mt")) {
        $VERBOSE and printerr("\n**********\n* Starting fetching of Mutation Taster predictions with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_mt -buildver $USEDBUILD -outfile $OUTFILE_MT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_mt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_MT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of Mutation Taster predictions done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output PhyloP scores
    if (doIt(-elem => "phylop")) {
        $VERBOSE and printerr("\n**********\n* Starting fetching of PhyloP predictions with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_phylop -buildver $USEDBUILD -outfile $OUTFILE_PHYLOP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_phylop -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PHYLOP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of PhyloP predictions done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output GERP++ scores
    if (doIt(-elem => "gerp")) {
        my $dbfile = $USEDBUILD . "_ljb_gerp++.txt";
        $VERBOSE and printerr("\n**********\n* Starting fetching of GERP++ predictions with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -outfile $OUTFILE_GERP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GERP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of GERP++ predictions done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -filter to output ESP6500 scores
    if (doIt(-elem => "esp6500")) {
        my $dbfile = $USEDBUILD . "_esp6500si_all.txt";
        $VERBOSE and printerr("\n**********\n* Starting fetching of ESP6500 frequencies with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -outfile $OUTFILE_ESP6500 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP6500 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of ESP6500 frequencies done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -regionanno to output MIRNA regions
    if (doIt(-elem => "mirna")) {
        my $dbfile = $USEDBUILD . "_wgRna.txt";
        $VERBOSE and printerr("\n**********\n* Starting fetching of  MIRNA annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirna -buildver $USEDBUILD -outfile $OUTFILE_MIRNA $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirna -buildver $USEDBUILD -outfile $OUTFILE_MIRNA $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of MIRNA annotations done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -regionanno to output MIRNATARGET regions
    if (doIt(-elem => "mirnatarget")) {
        my $dbfile = $USEDBUILD . "_targetScanS.txt";
        $VERBOSE and printerr("\n**********\n* Starting fetching of  MIRNATARGET annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirnatarget -buildver $USEDBUILD -outfile $OUTFILE_MIRNATARGET $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirnatarget -buildver $USEDBUILD -outfile $OUTFILE_MIRNATARGET $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of MIRNATARGET annotations done...\n**********\n\n");
    }

    # if not excluded, run ANNOVAR -regionanno to output RNAseq regions
    if (doIt(-elem => "rnaseq_cuff")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'rnaseq_cuff' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hg19_Cufflinks_allTissues.sorted.custom.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  RNASEQ_CUFF annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_CUFF $ANNOINFILE $BODYMAPPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_CUFF $ANNOINFILE $BODYMAPPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of RNASEQ_CUFF annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output RNAseq regions
    if (doIt(-elem => "rnaseq_scri")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'rnaseq_scri_brain_r' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hg19_Scripture_allTissues.sorted.custom.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  RNASEQ_SCRI annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_SCRI $ANNOINFILE $BODYMAPPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_SCRI $ANNOINFILE $BODYMAPPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of RNASEQ_SCRI annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_chipseq")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_chipseq' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "chip_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_CHIPSEQ annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CHIPSEQ $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CHIPSEQ $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_CHIPSEQ annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_clinvar")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_clinvar' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "clinvar_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_CLINVAR annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CLINVAR $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CLINVAR $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_CLINVAR annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_common_snp")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_common_snp' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "common_snp_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_COMMON_SNP annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COMMON_SNP $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COMMON_SNP $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_COMMON_SNP annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_cosmic")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_cosmic' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "cosmic_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_COSMIC annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COSMIC $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COSMIC $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_COSMIC annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_cpg")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_cpg' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "cpg_islands_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_CPG annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CPG $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CPG $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_CPG annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_dbnsfp")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_dbnsfp' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "dbnsfp_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_DBNSFP annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DBNSFP $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DBNSFP $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_DBNSFP annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_disease")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_disease' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "disease_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_DISEASE annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DISEASE $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DISEASE $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_DISEASE annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_drug")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_drug' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "drug_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_DRUG annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DRUG $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DRUG $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_DRUG annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_evs")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_evs' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "evs_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_EVS annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_EVS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_EVS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_EVS annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_gwas")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_gwas' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "gwas_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_GWAS annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_GWAS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_GWAS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_GWAS annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_hgmd")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_hgmd' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hgmd_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_HGMD annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_HGMD annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_hgmd_disgene")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_hgmd_disgene' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hgmd_disease_genes_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_HGMD_DISGENE annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD_DISGENE $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD_DISGENE $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_HGMD_DISGENE annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_microsat")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_microsat' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "microsatellites_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_MICROSAT annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MICROSAT $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MICROSAT $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_MICROSAT annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_mirna")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_mirna' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "miRNA_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_MIRNA annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MIRNA $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MIRNA $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_MIRNA annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_omim")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_omim' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "omim_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_OMIM annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_OMIM $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_OMIM $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_OMIM annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_path")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_path' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "pathway_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_PATH annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PATH $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PATH $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_PATH annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_pgx")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_pgx' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "pgx_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of  GTX_PGX annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PGX $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PGX $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_PGX annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_ptms")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_ptms' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "ptms_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_PTMS annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PTMS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PTMS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_PTMS annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_transfac")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_transfac' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "transfac_sites_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_TRANSFAC annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TRANSFAC $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TRANSFAC $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_TRANSFAC annotations done...\n**********\n\n");
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_tss")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_tss' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "tss_hg19.bed";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GTX_TSS annotations with the following command line :\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TSS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\n**********\n\n");
            system "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TSS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GTX_TSS annotations done...\n**********\n\n");
        }
    }
}

sub launchAnnovarParallel {
    my %arg = @_;
    my $rO_vcf  = $arg{-vcf};
    my $rH_data = $arg{-data};

    $OUTFILE_ANNO .= "." . $BUILDVER;

    $USEDBUILD = $BUILDVER;
    ($USEDBUILD =~ /b37|v37/) and $USEDBUILD = "hg19";

    # First prepare smp-runner input file
    prepareSmpInputFile();

    # Launch Annovar with smp-runner
    my $smpLogFolderPath = $SUBOUTPATH . "smp.log";
    system "smp-runner --command=\"\" --smp=$CPU --quiet --log=$smpLogFolderPath < $SMPINPUTFILE;";
    wait();
}

sub prepareSmpInputFile {
    my %arg = @_;

    $VERBOSE and printerr("\n***********\n* Preparing smp-sunner input file...\n**********\n\n");
    open(SMPINPUT, ">$SMPINPUTFILE") || die "cannot open file $SMPINPUTFILE:\n$!\n";

    if (doIt(-elem => "gene")) {
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting genome annotation using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -geneanno -buildver $USEDBUILD -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $ANNOVARPATH/humandb/ -separate -exonicsplicing -transcript_function $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -geneanno -buildver $USEDBUILD -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $ANNOVARPATH/humandb/ -separate -exonicsplicing -transcript_function $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Genome annotation done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "snp")) {
        my $dbsnp;
        if ($DBSNP) {
            $dbsnp = $DBSNP;
        } elsif ($USEDBUILD eq "hg18") {
            $dbsnp = "snp130";
        } elsif ($USEDBUILD eq "hg19") {
            $dbsnp = "snp131";
        }

        # First run annovar with clinically associated dbSNP
        $dbsnp .= "_clinic";
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against dbSNP with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against dbSNP done\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";

        # Then run annovar with non-clinically associated dbSNP
        $dbsnp =~ s/_clinic$/_non_clinic/;
        $smpArgs = "";
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against dbSNP with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against dbSNP done\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "1kg")) {
        my $dbtype;
        if ($USEDBUILD eq "hg18") {
            $dbtype = "1000g2010jul_ceu";
        } elsif ($USEDBUILD eq "hg19") {
            $dbtype = "1000g2012feb_all";
        }

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against 1000 Genomes db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_1000G $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_1000G $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against 1000 Genomes db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "tsi")) {
        my $genericdbfile;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_TSI_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'tsi' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_TSI_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
        }

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against 1000 Genomes Toscani frequencies db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_TSI $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_TSI $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against 1000 Genomes Toscani frequencies db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "asn")) {
        my $genericdbfile;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_ASN_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'asn' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_ASN_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
        }

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against 1000 Genomes Asian frequencies db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ASN $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ASN $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against 1000 Genomes Asian frequencies db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "afr")) {
        my $genericdbfile;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_AFR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'afr' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_AFR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
        }

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against 1000 Genomes African frequencies db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AFR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AFR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against 1000 Genomes African frequencies db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "amr")) {
        my $genericdbfile;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_AMR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'amr' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_AMR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
        }

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against 1000 Genomes American frequencies db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AMR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_AMR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against 1000 Genomes American frequencies db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "eur")) {
        my $genericdbfile;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_EUR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'eur' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_EUR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
        }

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against 1000 Genomes European frequencies db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EUR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $genericdbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_EUR $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against 1000 Genomes European frequencies db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "cg")) {
        my $dbfile = $USEDBUILD . "_cg69.txt";

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting filtering against Complete Genomics db with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_CG $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_CG $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Filtering against Complete Genomics db done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "sift")) {
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of SIFT prediction scores with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype avsift -sift_threshold 0 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_SIFT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype avsift -sift_threshold 0 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_SIFT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of SIFT prediction scores done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "pp2")) {
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of PolyPhen (v2) prediction scores with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_pp2 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PP2 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_pp2 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PP2 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of PolyPhen (v2) prediction scores done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "lrt")) {
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of Likelyhood-Ratio-Test scores with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_lrt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_LRT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_lrt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_LRT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of Likelyhood-Ratio-Test scores done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "mt")) {
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of Mutation Taster prediction scores with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_mt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_MT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_mt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_MT $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of Mutation Taster prediction scores done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "phylop")) {
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of PhyloP prediction scores with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_phylop -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PHYLOP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype ljb_phylop -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PHYLOP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of PhyloP  prediction scores done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "gerp")) {
        my $dbfile = $USEDBUILD . "_ljb_gerp++.txt";

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of GERP++ prediction scores with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GERP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GERP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GERP++ prediction scores done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "esp6500")) {
        my $dbfile = $USEDBUILD . "_esp6500si_all.txt";
        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of ESP6500 frequencies with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP6500 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP6500 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of ESP6500 frequencies done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "mirna")) {
        my $dbfile = $USEDBUILD . "_wgRna.txt";

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of MIRNA annotations with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirna -buildver $USEDBUILD -outfile $OUTFILE_MIRNA $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirna -buildver $USEDBUILD -outfile $OUTFILE_MIRNA $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of MIRNA annotations done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "mirnatarget")) {
        my $dbfile = $USEDBUILD . "_targetScanS.txt";

        my $smpArgs;
        $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of MIRNATARGET annotations with the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirnatarget -buildver $USEDBUILD -outfile $OUTFILE_MIRNATARGET $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\\n**********\\n\\n';";
        $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype mirnatarget -buildver $USEDBUILD -outfile $OUTFILE_MIRNATARGET $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE";
        $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of MIRNATARGET annotations done...\\n**********\\n\\n'";
        print SMPINPUT $smpArgs . "\n";
    }

    if (doIt(-elem => "rnaseq_cuff")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'rnaseq_cuff' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hg19_Cufflinks_allTissues.sorted.custom.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of RNASEQ_CUFF annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_CUFF $ANNOINFILE $BODYMAPPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_CUFF $ANNOINFILE $BODYMAPPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of RNASEQ_CUFF annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "rnaseq_scri")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'rnaseq_scri_brain_r' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hg19_Scripture_allTissues.sorted.custom.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting fetching of RNASEQ_SCRI annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_SCRI $ANNOINFILE $BODYMAPPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_RNASEQ_SCRI $ANNOINFILE $BODYMAPPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of RNASEQ_SCRI annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_chipseq")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_chipseq' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "chip_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_CHIPSEQ annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CHIPSEQ $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CHIPSEQ $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n*Fetching of GTX_CHIPSEQ annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_clinvar")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_clinvar' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "clinvar_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_CLINVAR annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CLINVAR $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CLINVAR $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n*Fetching of GTX_CLINVAR annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_common_snp")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_common_snp' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "common_snp_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_COMMON_SNP annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COMMON_SNP $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COMMON_SNP $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_COMMON_SNP annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_cosmic")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_cosmic' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "cosmic_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_COSMIC annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COSMIC $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_COSMIC $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_COSMIC annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_cpg")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_cpg' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "cpg_islands_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_CPG annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CPG $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_CPG $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_CPG annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_dbnsfp")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_dbnsfp' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "dbnsfp_hg19.bed";

           my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_DBNSFP annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DBNSFP $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DBNSFP $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_DBNSFP annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_disease")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_disease' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "disease_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_DISEASE annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DISEASE $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DISEASE $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_DISEASE annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_drug")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_drug' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "drug_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_DRUG annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DRUG $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_DRUG $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_DRUG annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_evs")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_evs' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "evs_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_EVS annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_EVS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_EVS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_EVS annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_gwas")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_gwas' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "gwas_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_GWAS annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_GWAS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_GWAS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_GWAS annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_hgmd")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_hgmd' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hgmd_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_HGMD annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_HGMD annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_hgmd_disgene")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_hgmd_disgene' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hgmd_disease_genes_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_HGMD_DISGENE annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD_DISGENE $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_HGMD_DISGENE $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_HGMD_DISGENE annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_microsat")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_microsat' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "microsatellites_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_MICROSAT annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MICROSAT $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MICROSAT $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_MICROSAT annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_mirna")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_mirna' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "miRNA_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_MIRNA annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MIRNA $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_MIRNA $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_MIRNA annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_omim")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_omim' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "omim_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_OMIM annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_OMIM $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_OMIM $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_OMIM annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_path")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_path' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "pathway_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_PATH annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PATH $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PATH $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_PATH annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_pgx")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_pgx' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "pgx_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_PGX annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PGX $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PGX $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_PGX annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_ptms")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_ptms' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "ptms_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_PTMS annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PTMS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_PTMS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_PTMS annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_transfac")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_transfac' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "transfac_sites_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_TRANSFAC annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TRANSFAC $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TRANSFAC $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_TRANSFAC annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }

    if (doIt(-elem => "gtx_tss")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_tss' regionanno is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "tss_hg19.bed";

            my $smpArgs;
            $VERBOSE and $smpArgs = "echo -e '\\n**********\\n* Starting Fetching of GTX_TSS annotations using the following command line :\\n* perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TSS $ANNOINFILE $GENOMETRAXPATH $VERBOSE\\n**********\\n\\n';";
            $smpArgs .= "perl $ANNOVARPATH/annotate_variation.pl -regionanno -dbtype bed -bedfile $dbfile -colsWanted 4 -buildver $USEDBUILD -outfile $OUTFILE_GTX_TSS $ANNOINFILE $GENOMETRAXPATH $VERBOSE";
            $VERBOSE and $smpArgs .= ";echo -e '\\n**********\\n* Fetching of GTX_TSS annotations done...\\n**********\\n\\n'";
            print SMPINPUT $smpArgs . "\n";
        }
    }
    close(SMPINPUT);
}

sub buildVcf {
    my %arg = @_;
    my $rO_vcf  = $arg{-vcf};
    my $rH_data = $arg{-data};

    $USEDBUILD = $BUILDVER;
    ($USEDBUILD =~ /b37|v37/) and $USEDBUILD = "hg19";

    if (doIt(-elem => "gene")) {
        # open ANNOVAR output files for processing
        open(ANNOVAR_MAIN_GENERAL, "<$OUTFILE_ANNO.variant_function") || die "cannot open file $OUTFILE_ANNO.variant_function:\n$!\n";
        open(ANNOVAR_MAIN_DETAILED, "<$OUTFILE_ANNO.exonic_variant_function") || die "cannot open file $OUTFILE_ANNO.exonic_variant_function:\n$!\n";

        $rO_vcf->add_header_line({key=>'INFO', ID=>'VC',   Number=>'A', Type=>'String', Description=>'Variant class'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'VFT',  Number=>'A', Type=>'String', Description=>'Variant function type'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'VT',   Number=>'A', Type=>'String', Description=>'Coding variant type'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'Gene', Number=>1, Type=>'String', Description=>'Gene symbol'});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'DA',   Number=>'A', Type=>'String', Description=>'Detailed annotation of the variant'});

        # annotate noncoding variants - get gene symbol (col 2) and general functional effect (col 1)
        while (my $line = <ANNOVAR_MAIN_GENERAL>) {
            chomp $line;  # $1    $2            $3                  $4        $5     $6       $7         $8
            if ($line =~ /^(.+?)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $varFuncType = $1;
                my $gene = $2;
                $gene =~ s/;/:/g;   # ABC;NM_blabla becomes ABC:NM_blabla
                push @{$rH_data->{$key}->{'VFT'}}, $varFuncType || "";
                if ($varFuncType eq "intergenic") {
                    $gene =~ s/dist=/dist:/g;
                    $gene =~ s/,/-/g;
                    $gene =~ s/(\w+)\:\w+(\(dist:\d+\))/$1$2/g;
                    push @{$rH_data->{$key}->{'DA'}}, "$varFuncType:$gene";
                } else {
                    $gene =~ s/,/|/g;
                    push @{$rH_data->{$key}->{'Gene'}}, $gene;
                    my @genes = split /\|/, $gene;  # / just to avoid a bug in eclipse display...
                    $gene = join "|$varFuncType:", sort @genes;
                    if ($varFuncType !~ /^exonic/) {
                        push @{$rH_data->{$key}->{'DA'}}, "$varFuncType:$gene"; # / just to avoid a bug in eclipse display...
                    }
                }
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MAIN_GENERAL);

        # further annotate coding variants - using line number (col 1) to compare to key of global input hash, annotate synonymous/nonsynonymous/frameshift/nonframeshift/insertion/deletion (col 2) and isoform/protein/cdna/AA info (col 3)
        while (my $line = <ANNOVAR_MAIN_DETAILED>) {
            chomp $line;    #        $1    $2            $3               $4        $5     $6       $7         $8
            if ($line =~ /^line\d+\t(.+?)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                my $variantType = $1;
                my $annotation = $2;
                $annotation =~ s/,\s*$//g;
                $annotation =~ s/,/||/g;
                $variantType =~ s/ /_/;
                my $annotString = "";
                foreach my $varFuncType (@{$rH_data->{$key}->{'VFT'}}) {
                    if ($varFuncType eq "exonic" || $varFuncType eq "exonic_splicing") {
                        my @annotations = split /\|\|/, $annotation;    # / just to avoid a bug in eclipse display...
                        $annotation = join "|$varFuncType:$variantType:", sort @annotations;
                        $annotString = "$varFuncType:$variantType:$annotation";
                    }
                }
                if ($annotString eq "") {
                    my @annotations = split /\|\|/, $annotation; # / just to avoid a bug in eclipse display...
                    $annotation = join "|exonic:$variantType:", sort @annotations;
                }
                push @{$rH_data->{$key}->{'DA'}}, $annotString;
                $rH_data->{$key}->{'VT'} = $variantType;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MAIN_DETAILED);
    }


    # if not excluded, run ANNOVAR -filter to output dbSNP variants
    if (doIt(-elem => "snp")) {
        my $dbsnp;
        if ($DBSNP) {
            $dbsnp = $DBSNP;
        } elsif ($USEDBUILD eq "hg18") {
            $dbsnp = "snp130";
        } elsif ($USEDBUILD eq "hg19") {
            $dbsnp = "snp131";
        }
        $rO_vcf->add_header_line({key=>'INFO', ID=>'CdbSNP_rs', Number=>1, Type=>'String', Description=>"dbSNP id of a clinically associated variation (version $dbsnp)"});
        $rO_vcf->add_header_line({key=>'INFO', ID=>'NCdbSNP_rs', Number=>1, Type=>'String', Description=>"dbSNP id of a non-clinically associated variation (version $dbsnp)"});

        # open ANNOVAR output file for processing
        $dbsnp .= "_clinic";
        open(ANNOVAR_DBSNP, "<$OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped") || die "cannot open file $OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped:\n$!\n";

        while (my $line = <ANNOVAR_DBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #     $1           $2            $3               $4        $5     $6       $7         $8
            if ($line =~ /(snp\d+_clinic)\t(rs\d+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{CdbSNP_rs} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP);

        # open ANNOVAR output file for processing
        $dbsnp =~ s/_clinic$/_non_clinic/;
        open(ANNOVAR_DBSNP, "<$OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped") || die "cannot open file $OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped:\n$!\n";

        while (my $line = <ANNOVAR_DBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #        $1           $2             $3               $4        $5     $6       $7         $8
            if ($line =~ /(snp\d+_non_clinic)\t(rs\d+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{NCdbSNP_rs} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP);
        $dbsnp =~ s/_non_clinic$//;
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "1kg")) {
        my $dbtype;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $dbtype = "1000g2010jul_ceu";
            $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_CEU.sites.2010_07_dropped";
        } elsif ($USEDBUILD eq "hg19") {
            $dbtype = "1000g2012feb_all";
            $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_ALL.sites.2012_02_dropped";
        }
        $rO_vcf->add_header_line({key=>'INFO', ID=>'1KG', Number=>'A', Type=>'Float', Description=>"Frequency in the 1000 Genome project (version $dbtype)"});

        # open ANNOVAR output file for processing
        open(ANNOVAR_1000g, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_1000g>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1          $2               $3               $4        $5     $6       $7         $8
            if ($line =~ /($dbtype)\t([01]\.*\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{'1KG'} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_1000g);
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "tsi")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_TSI_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_TSI.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'tsi' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_TSI_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_TSI.$USEDBUILD\_generic_dropped";

            $rO_vcf->add_header_line({key=>'INFO', ID=>'AF_TSI', Number=>'A', Type=>'Float', Description=>"Allele frequency of Toscani population of the 1000 Genome project (phase1_release_v3.20101123)"});

            # open ANNOVAR output file for processing
            open(ANNOVAR_TSI, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_TSI>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    # $1          $2               $3               $4        $5     $6       $7         $8
                if ($line =~ /(generic)\t([01]\.*\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{AF_TSI} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_TSI);
        }
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "asn")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_ASN_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_ASN.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'asn' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_ASN_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_ASN.$USEDBUILD\_generic_dropped";

            $rO_vcf->add_header_line({key=>'INFO', ID=>'AF_ASN', Number=>'A', Type=>'Float', Description=>"Allele frequency of East Asian population of the 1000 Genome project (phase1_release_v3.20101123)"});

            # open ANNOVAR output file for processing
            open(ANNOVAR_ASN, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_ASN>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    # $1          $2               $3               $4        $5     $6       $7         $8
                if ($line =~ /(generic)\t([01]\.*\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{AF_ASN} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_ASN);
        }
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "afr")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_AFR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AFR.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'afr' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_AFR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AFR.$USEDBUILD\_generic_dropped";

            $rO_vcf->add_header_line({key=>'INFO', ID=>'AF_AFR', Number=>'A', Type=>'Float', Description=>"Allele frequency of African population of the 1000 Genome project (phase1_release_v3.20101123)"});

            # open ANNOVAR output file for processing
            open(ANNOVAR_AFR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_AFR>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    # $1          $2               $3               $4        $5     $6       $7         $8
                if ($line =~ /(generic)\t([01]\.*\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{AF_AFR} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_AFR);
        }
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "amr")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_AMR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AMR.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'amr' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_AMR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_AMR.$USEDBUILD\_generic_dropped";

            $rO_vcf->add_header_line({key=>'INFO', ID=>'AF_AMR', Number=>'A', Type=>'Float', Description=>"Allele frequency of Ad Mixed American population of the 1000 Genome project (phase1_release_v3.20101123)"});

            # open ANNOVAR output file for processing
            open(ANNOVAR_AMR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_AMR>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    # $1          $2               $3               $4        $5     $6       $7         $8
                if ($line =~ /(generic)\t([01]\.*\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{AF_AMR} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_AMR);
        }
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if (doIt(-elem => "eur")) {
        my $genericdbfile;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $genericdbfile = "hg18_EUR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_EUR.$USEDBUILD\_generic_dropped";
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & filter arguments... 'eur' filter is only supported by hg19 (or v37, or b37) build...\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            $genericdbfile = "hg19_EUR_FR.ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites";
            $annovarOutput = "$OUTFILE_EUR.$USEDBUILD\_generic_dropped";

            $rO_vcf->add_header_line({key=>'INFO', ID=>'AF_EUR', Number=>'A', Type=>'Float', Description=>"Allele frequency of European population of the 1000 Genome project (phase1_release_v3.20101123)"});

            # open ANNOVAR output file for processing
            open(ANNOVAR_EUR, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_EUR>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    # $1          $2               $3               $4        $5     $6       $7         $8
                if ($line =~ /(generic)\t([01]\.*\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{AF_EUR} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_EUR);
        }
    }


    # if not excluded, run ANNOVAR -filter to output Complete Genomics variants
    if (doIt(-elem => "cg")) {
        my $dbfile = $USEDBUILD . "_cg69.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'CG', Number=>'A', Type=>'Float', Description=>'allele frequency in 69 human subjects sequenced by Complete Genomics (db updated on 2012Feb22)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_CG, "<$OUTFILE_CG.$USEDBUILD\_generic\_dropped") || die "cannot open file $OUTFILE_CG.$USEDBUILD\_generic\_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_CG>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1         $2               $3                $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t([01]\.\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{CG} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_CG);
    }


    # if not excluded, run ANNOVAR -filter to output SIFT scores
    if (doIt(-elem => "sift")) {
        $rO_vcf->add_header_line({key=>'INFO', ID=>'SIFT', Number=>'A', Type=>'Float', Description=>'whole-exome SIFT scores (db updated on 2011Mar01)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_SIFT, "<$OUTFILE_SIFT.$USEDBUILD\_avsift_dropped") || die "cannot open file $OUTFILE_SIFT.$USEDBUILD\_avsift_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_SIFT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1     $2            $3               $4        $5     $6       $7         $8
            if ($line =~ /(avsift)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{SIFT} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_SIFT);
    }


    # if not excluded, run ANNOVAR -filter to output PolyPhen.v2 scores
    if (doIt(-elem => "pp2")) {
        $rO_vcf->add_header_line({key=>'INFO', ID=>'PPv2', Number=>'A', Type=>'Float', Description=>'whole-exome PolyPhen version 2 scores (db updated 2011May11)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_PP2, "<$OUTFILE_PP2.$USEDBUILD\_ljb_pp2_dropped") || die "cannot open file $OUTFILE_PP2.$USEDBUILD\_ljb_pp2_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_PP2>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2            $3               $4        $5     $6       $7         $8
            if ($line =~ /(ljb_pp2)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{PPv2} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_PP2);
    }


    # if not excluded, run ANNOVAR -filter to output LRT scores
    if (doIt(-elem => "lrt")) {
        $rO_vcf->add_header_line({key=>'INFO', ID=>'LRT', Number=>'A', Type=>'Float', Description=>'whole-exome LRT scores (db updated 2011May11)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_LRT, "<$OUTFILE_LRT.$USEDBUILD\_ljb_lrt_dropped") || die "cannot open file $OUTFILE_LRT.$USEDBUILD\_ljb_lrt_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_LRT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1     $2            $3               $4        $5     $6       $7         $8
            if ($line =~ /(ljb_lrt)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{LRT} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_LRT);
    }


    # if not excluded, run ANNOVAR -filter to output Mutation Taster scores
    if (doIt(-elem => "mt")) {
        $rO_vcf->add_header_line({key=>'INFO', ID=>'MT', Number=>'A', Type=>'Float', Description=>'whole-exome MutationTaster scores (db updated 2011May11)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_MT, "<$OUTFILE_MT.$USEDBUILD\_ljb_mt_dropped") || die "cannot open file $OUTFILE_MT.$USEDBUILD\_ljb_mt_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_MT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    # $1      $2          $3                $4        $5     $6       $7         $8
            if ($line =~ /(ljb_mt)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{MT} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MT);
    }


    # if not excluded, run ANNOVAR -filter to output PhyloP scores
    if (doIt(-elem => "phylop")) {
        $rO_vcf->add_header_line({key=>'INFO', ID=>'PhyloP', Number=>'A', Type=>'Float', Description=>'whole-exome PhyloP scores (db updated 2011May11)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_PHYLOP, "<$OUTFILE_PHYLOP.$USEDBUILD\_ljb_phylop_dropped") || die "cannot open file $OUTFILE_PHYLOP.$USEDBUILD\_ljb_phylop_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_PHYLOP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
            if ($line =~ /(ljb_phylop)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{PhyloP} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_PHYLOP);
    }


    # if not excluded, run ANNOVAR -filter to output GERP++ scores
    if (doIt(-elem => "gerp")) {
        my $dbfile = $USEDBUILD . "_ljb_gerp++.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'GERP++', Number=>'A', Type=>'Float', Description=>'whole-exome GERP++ scores (db updated 2012Feb22)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_GERP, "<$OUTFILE_GERP.$USEDBUILD\_generic_dropped") || die "cannot open file $OUTFILE_GERP.$USEDBUILD\_ljb_gerp++_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_GERP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{GERP} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_GERP);
    }


    # if not excluded, run ANNOVAR -filter to output ESP6500 scores
    if (doIt(-elem => "esp6500")) {
        my $dbfile = $USEDBUILD . "_esp6500si_all.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'ESP6500', Number=>'A', Type=>'Float', Description=>'alternative allele frequency in all subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls (db updated 2013Jan22)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_ESP6500, "<$OUTFILE_ESP6500.$USEDBUILD\_generic_dropped") || die "cannot open file $OUTFILE_ESP6500.$USEDBUILD\_esp6500_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_ESP6500>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ESP6500} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ESP6500);
    }


    # if not excluded, run ANNOVAR -regionanno to output MIRNA regions
    if (doIt(-elem => "mirna")) {
        my $dbfile = $USEDBUILD . "_wgRna.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'MIRNA', Number=>1, Type=>'Float', Description=>'snoRNA and miRNA annotations (from UCSC wgRna table, updated 2012Apr03)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_MIRNA, "<$OUTFILE_MIRNA.$USEDBUILD\_wgRna") || die "cannot open file $OUTFILE_MIRNA.$USEDBUILD\_wgRna\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_MIRNA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
            if ($line =~ /(mirna)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{MIRNA} = $2;
                $rH_data->{$key}->{MIRNA} =~ s/^Name\=//;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MIRNA);
    }


    # if not excluded, run ANNOVAR -regionanno to output MIRNATARGET regions
    if (doIt(-elem => "mirnatarget")) {
        my $dbfile = $USEDBUILD . "_targetScanS.txt";
        $rO_vcf->add_header_line({key=>'INFO', ID=>'MIRNATARGET', Number=>1, Type=>'Float', Description=>'TargetScan generated miRNA target site predictions (from UCSC TargetScanS table, updated 2012Apr03)'});

        # open ANNOVAR output file for processing
        open(ANNOVAR_MIRNATARGET, "<$OUTFILE_MIRNATARGET.$USEDBUILD\_targetScanS") || die "cannot open file $OUTFILE_MIRNATARGET.$USEDBUILD\_targetScanS\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_MIRNATARGET>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
            if ($line =~ /(mirnatarget)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{MIRNATARGET} = $2;
                $rH_data->{$key}->{MIRNATARGET} =~ s/^(Score\=\d+;)?Name\=//;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_MIRNATARGET);
    }


    # if not excluded, run ANNOVAR -regionanno to output RNAseq regions
    if (doIt(-elem => "rnaseq_cuff")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'rnaseq_cuff' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hg19_Cufflinks_allTissues.sorted.custom.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'RNASEQ_CUFF', Number=>1, Type=>'String', Description=>'Broad Cufflinks RNASeq alignment of BodyMap v2 (all tissues)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_RNASEQ_CUFF, "<$OUTFILE_RNASEQ_CUFF.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_RNASEQ_CUFF.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_RNASEQ_CUFF>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{RNASEQ_CUFF} = $2;
                    $rH_data->{$key}->{RNASEQ_CUFF} =~ s/^Name\=//;
                    ($rH_data->{$key}->{RNASEQ_CUFF} eq "NA") and $rH_data->{$key}->{RNASEQ_CUFF} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_RNASEQ_CUFF);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output RNAseq regions
    if (doIt(-elem => "rnaseq_scri")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'rnaseq_scri_brain_r' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hg19_Scripture_allTissues.sorted.custom.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'RNASEQ_SCRI', Number=>1, Type=>'String', Description=>'Broad Scripture RNASeq alignment of BodyMap v2 (all tissues)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_RNASEQ_SCRI, "<$OUTFILE_RNASEQ_SCRI.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_RNASEQ_SCRI.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_RNASEQ_SCRI>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{RNASEQ_SCRI} = $2;
                    $rH_data->{$key}->{RNASEQ_SCRI} =~ s/^Name\=//;
                    ($rH_data->{$key}->{RNASEQ_SCRI} eq "NA") and $rH_data->{$key}->{RNASEQ_SCRI} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_RNASEQ_SCRI);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_chipseq")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_chipseq' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "chip_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_CHIPSEQ', Number=>1, Type=>'String', Description=>'Genome Trax Predicted ChIP-Seq TFBS description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_CHIPSEQ, "<$OUTFILE_GTX_CHIPSEQ.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_CHIPSEQ.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_CHIPSEQ>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_CHIPSEQ} = $2;
                    $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_CHIPSEQ} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_CHIPSEQ} eq "NA") and $rH_data->{$key}->{GTX_CHIPSEQ} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_CHIPSEQ);
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_clinvar")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_clinvar' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "chip_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_CLINVAR', Number=>1, Type=>'String', Description=>'Genome Trax clinically significant variants from NCBI\'s ClinVar resource (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_CLINVAR, "<$OUTFILE_GTX_CLINVAR.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_CLINVAR.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_CLINVAR>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_CLINVAR} = $2;
                    $rH_data->{$key}->{GTX_CLINVAR} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_CLINVAR} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_CLINVAR} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_CLINVAR} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_CLINVAR} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_CLINVAR} eq "NA") and $rH_data->{$key}->{GTX_CLINVAR} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stopped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_CLINVAR);
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_common_snp")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_common_snp' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "common_snp_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_COMMON_SNP', Number=>1, Type=>'String', Description=>'Genome Trax Predicted Common SNP description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_COMMON_SNP, "<$OUTFILE_GTX_COMMON_SNP.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_COMMON_SNP.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_COMMON_SNP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_COMMON_SNP} = $2;
                    $rH_data->{$key}->{GTX_COMMON_SNP} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_COMMON_SNP} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_COMMON_SNP} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_COMMON_SNP} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_COMMON_SNP} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_COMMON_SNP} eq "NA") and $rH_data->{$key}->{GTX_COMMON_SNP} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_COMMON_SNP);
        }
    }

    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_cosmic")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_cosmic' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "cosmic_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_COSMIC', Number=>1, Type=>'String', Description=>'Genome Trax COSMIC (Catalogue of Somatic Mutations in Cancer) somatic disease mutations description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_COSMIC, "<$OUTFILE_GTX_COSMIC.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_COSMIC.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_COSMIC>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_COSMIC} = $2;
                    $rH_data->{$key}->{GTX_COSMIC} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_COSMIC} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_COSMIC} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_COSMIC} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_COSMIC} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_COSMIC} eq "NA") and $rH_data->{$key}->{GTX_COSMIC} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_COSMIC);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_cpg")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_cpg' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "cpg_islands_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_CPG', Number=>1, Type=>'String', Description=>'Genome Trax CpG Islands description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_CPG, "<$OUTFILE_GTX_CPG.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_CPG.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_CPG>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_CPG} = $2;
                    $rH_data->{$key}->{GTX_CPG} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_CPG} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_CPG} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_CPG} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_CPG} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_CPG} eq "NA") and $rH_data->{$key}->{GTX_CPG} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_CPG);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_dbnsfp")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_dbnsfp' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "dbnsfp_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_DBNSFP', Number=>1, Type=>'String', Description=>'Genome Trax DBNSFP description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_DBNSFP, "<$OUTFILE_GTX_DBNSFP.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_DBNSFP.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_DBNSFP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_DBNSFP} = $2;
                    $rH_data->{$key}->{GTX_DBNSFP} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_DBNSFP} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_DBNSFP} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_DBNSFP} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_DBNSFP} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_DBNSFP} eq "NA") and $rH_data->{$key}->{GTX_DBNSFP} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_DBNSFP);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_disease")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_disease' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "disease_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_DISEASE', Number=>1, Type=>'String', Description=>'Genome Trax Disease associations description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_DISEASE, "<$OUTFILE_GTX_DISEASE.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_DISEASE.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_DISEASE>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_DISEASE} = $2;
                    $rH_data->{$key}->{GTX_DISEASE} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_DISEASE} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_DISEASE} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_DISEASE} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_DISEASE} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_DISEASE} eq "NA") and $rH_data->{$key}->{GTX_DISEASE} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_DISEASE);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_drug")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_drug' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "drug_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_DRUG', Number=>1, Type=>'String', Description=>'Genome Trax Drug Targets description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_DRUG, "<$OUTFILE_GTX_DRUG.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_DRUG.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_DRUG>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_DRUG} = $2;
                    $rH_data->{$key}->{GTX_DRUG} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_DRUG} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_DRUG} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_DRUG} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_DRUG} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_DRUG} eq "NA") and $rH_data->{$key}->{GTX_DRUG} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_DRUG);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_evs")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_evs' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "evs_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_EVS', Number=>1, Type=>'String', Description=>'Genome Trax EVS description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_EVS, "<$OUTFILE_GTX_EVS.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_EVS.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_EVS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_EVS} = $2;
                    $rH_data->{$key}->{GTX_EVS} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_EVS} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_EVS} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_EVS} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_EVS} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_EVS} eq "NA") and $rH_data->{$key}->{GTX_EVS} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_EVS);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_gwas")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_gwas' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "gwas_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_GWAS', Number=>1, Type=>'String', Description=>'Genome Trax GWAS Catalogue description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_GWAS, "<$OUTFILE_GTX_GWAS.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_GWAS.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_GWAS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_GWAS} = $2;
                    $rH_data->{$key}->{GTX_GWAS} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_GWAS} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_GWAS} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_GWAS} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_GWAS} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_GWAS} eq "NA") and $rH_data->{$key}->{GTX_GWAS} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_GWAS);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_hgmd")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_hgmd' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hgmd_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_HGMD', Number=>1, Type=>'String', Description=>'Genome Trax HGMD inherited (germ-line) disease mutations description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_HGMD, "<$OUTFILE_GTX_HGMD.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_HGMD.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_HGMD>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_HGMD} = $2;
                    $rH_data->{$key}->{GTX_HGMD} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_HGMD} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_HGMD} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_HGMD} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_HGMD} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_HGMD} eq "NA") and $rH_data->{$key}->{GTX_HGMD} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_HGMD);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_hgmd_disgene")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_hgmd_disgene' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "hgmd_disease_genes_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_HGMD_DISGENE', Number=>1, Type=>'String', Description=>'Genome Trax HGMD disease genes description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_HGMD_DISGENE, "<$OUTFILE_GTX_HGMD_DISGENE.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_HGMD_DISGENE.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_HGMD_DISGENE>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_HGMD_DISGENE} = $2;
                    $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_HGMD_DISGENE} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_HGMD_DISGENE} eq "NA") and $rH_data->{$key}->{GTX_HGMD_DISGENE} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_HGMD_DISGENE);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_microsat")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_microsat' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "microsatellites_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_MICROSAT', Number=>1, Type=>'String', Description=>'Genome Trax Microsatellite repeats description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_MICROSAT, "<$OUTFILE_GTX_MICROSAT.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_MICROSAT.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_MICROSAT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_MICROSAT} = ($2) ? 1 : undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_MICROSAT);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_mirna")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_mirna' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "miRNA_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_MIRNA', Number=>1, Type=>'String', Description=>'Genome Trax microRNA sequences description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_MIRNA, "<$OUTFILE_GTX_MIRNA.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_MIRNA.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_MIRNA>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_MIRNA} = $2;
                    $rH_data->{$key}->{GTX_MIRNA} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_MIRNA} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_MIRNA} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_MIRNA} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_MIRNA} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_MIRNA} eq "NA") and $rH_data->{$key}->{GTX_MIRNA} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_MIRNA);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_omim")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_omim' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "omim_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_OMIM', Number=>1, Type=>'String', Description=>'Genome Trax OMIM disorders (NCBI - 22 feb. 2011)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_OMIM, "<$OUTFILE_GTX_OMIM.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_OMIM.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_OMIM>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_OMIM} = $2;
                    $rH_data->{$key}->{GTX_OMIM} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_OMIM} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_OMIM} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_OMIM} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_OMIM} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_OMIM} eq "NA") and $rH_data->{$key}->{GTX_OMIM} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_OMIM);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_path")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_path' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "pathway_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_PATH', Number=>1, Type=>'String', Description=>'Genome Trax Pathways membership description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_PATH, "<$OUTFILE_GTX_PATH.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_PATH.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_PATH>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_PATH} = $2;
                    $rH_data->{$key}->{GTX_PATH} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_PATH} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_PATH} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_PATH} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_PATH} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_PATH} eq "NA") and $rH_data->{$key}->{GTX_PATH} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_PATH);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_pgx")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_pgx' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "pgx_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_PGX', Number=>1, Type=>'String', Description=>'Genome Trax Pharmacogenomics Variants description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_PGX, "<$OUTFILE_GTX_PGX.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_PGX.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_PGX>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_PGX} = $2;
                    $rH_data->{$key}->{GTX_PGX} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_PGX} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_PGX} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_PGX} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_PGX} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_PGX} eq "NA") and $rH_data->{$key}->{GTX_PGX} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_PGX);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_ptms")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_ptms' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "ptms_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_PTMS', Number=>1, Type=>'String', Description=>'Genome Trax Post translational modifications description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_PTMS, "<$OUTFILE_GTX_PTMS.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_PTMS.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_PTMS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_PTMS} = $2;
                    $rH_data->{$key}->{GTX_PTMS} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_PTMS} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_PTMS} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_PTMS} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_PTMS} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_PTMS} eq "NA") and $rH_data->{$key}->{GTX_PTMS} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_PTMS);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_transfac")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_transfac' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "transfac_sites_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_TRANSFAC', Number=>1, Type=>'String', Description=>'Genome Trax TRANSFAC experimentally verified TFBS description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_TRANSFAC, "<$OUTFILE_GTX_TRANSFAC.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_TRANSFAC.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_TRANSFAC>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_TRANSFAC} = $2;
                    $rH_data->{$key}->{GTX_TRANSFAC} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_TRANSFAC} =~ s/\s+;/\|/g;
                    $rH_data->{$key}->{GTX_TRANSFAC} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_TRANSFAC} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_TRANSFAC} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_TRANSFAC} eq "NA") and $rH_data->{$key}->{GTX_TRANSFAC} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_TRANSFAC);
        }
    }


    # if not excluded, run ANNOVAR -regionanno to output Genome Trax HGMD Disease Genes Description
    if (doIt(-elem => "gtx_tss")) {
        if ($USEDBUILD eq "hg18") {
            printerr("\n!!!!!!!!!!\n! Incompatible buildver & regioanno arguments... 'gtx_tss' regionanno is only supported by hg19 (or v37, or b37) build.\n!!!!!!!!!!\n\n");
        } elsif ($USEDBUILD eq "hg19") {
            my $dbfile = "tss_hg19.bed";
            $rO_vcf->add_header_line({key=>'INFO', ID=>'GTX_TSS', Number=>1, Type=>'String', Description=>'Genome Trax TSSs (transcription start sites) description (BIOBASE - GenomeTrax.2013_1)'});

            # open ANNOVAR output file for processing
            open(ANNOVAR_GTX_TSS, "<$OUTFILE_GTX_TSS.$USEDBUILD\_bed") || die "cannot open file $OUTFILE_GTX_TSS.$USEDBUILD\_bed\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GTX_TSS>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(bed)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GTX_TSS} = $2;
                    $rH_data->{$key}->{GTX_TSS} =~ s/^Name\=//;
                    $rH_data->{$key}->{GTX_TSS} =~ s/\s*;/\|/g;
                    $rH_data->{$key}->{GTX_TSS} =~ s/, /,/g;
                    $rH_data->{$key}->{GTX_TSS} =~ s/ /_/g;
                    $rH_data->{$key}->{GTX_TSS} =~ s/__/_/g;
                    ($rH_data->{$key}->{GTX_TSS} eq "NA") and $rH_data->{$key}->{GTX_TSS} = undef;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line " . __LINE__ . " : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GTX_TSS);
        }
    }

}

sub printVcf {
    my %arg = @_;
    my $rO_vcf      = $arg{-vcf};
    my $rH_data     = $arg{-data};
    my $rH_vcfIndex = $arg{-index};

    # output results to a results file
    my $resultFile = $OUTFILE;
    (not $resultFile) and $resultFile = ($VARIANT_CALLER) ? $OUTPATH . $ID . "_ANNOVAR_output." . $VARIANT_CALLER . "." . $BUILDVER . ".RESULTS.vcf" : $OUTPATH . $ID . "_ANNOVAR_output." . $BUILDVER . ".RESULTS.vcf";
    open(RESULTS, ">$resultFile") || die ("Could not open RESULTS file $resultFile :\n$!\n");

    print RESULTS $rO_vcf->format_header();

    my $i = 0;
    while (my $rH_x = $rO_vcf->next_data_hash()) {
        $i++;
        my ($VC, $VT, $VFT, $Gene, $DA, $CdbSNP_rs, $NCdbSNP_rs, $oneKG, $CG, $SIFT, $PPv2) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($LRT, $MT, $PhyloP, $GERP, $ESP6500, $MIRNA, $MIRNATARGET, $TSI, $ASN, $AFR, $AMR) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
        my ($EUR, $RNASEQ_CUFF, $RNASEQ_SCRI, $GTX_CHIPSEQ, $GTX_CLINVAR, $GTX_COMMON_SNP, $GTX_COSMIC, $GTX_CPG) = (undef, undef, undef, undef, undef, undef, undef, undef);
        my ($GTX_DBNSFP, $GTX_DISEASE, $GTX_DRUG, $GTX_EVS, $GTX_GWAS, $GTX_HGMD, $GTX_HGMD_DISGENE, $GTX_MICROSAT) = (undef, undef, undef, undef, undef, undef, undef, undef);
        my ($GTX_MIRNA, $GTX_OMIM, $GTX_PATH, $GTX_PGX, $GTX_PTMS, $GTX_TRNASFAC, $GTX_TSS) = (undef, undef, undef, undef, undef, undef, undef);

        foreach my $variantKey (@{$rH_vcfIndex->{$rH_x->{CHROM} . ":" . $rH_x->{POS} . "_" . $rH_x->{REF} . "/" . join ",", @{$rH_x->{ALT}}}}) {
            $VC      .= (defined $rH_data->{$variantKey}->{VC})      ? $rH_data->{$variantKey}->{VC} . "," : "NIL,";
            $VFT     .= (defined $rH_data->{$variantKey}->{VFT})     ? join("|", @{$rH_data->{$variantKey}->{VFT}}) . "," : "NIL,";
            $VT      .= (defined $rH_data->{$variantKey}->{VT})      ? $rH_data->{$variantKey}->{VT} . "," : "NIL,";
            $DA      .= (defined $rH_data->{$variantKey}->{DA})      ? join("|", @{$rH_data->{$variantKey}->{DA}}) . "," : "NIL,";
            $oneKG   .= (defined $rH_data->{$variantKey}->{'1KG'})   ? $rH_data->{$variantKey}->{'1KG'} . "," : "-1,";
            $TSI     .= (defined $rH_data->{$variantKey}->{AF_TSI})  ? $rH_data->{$variantKey}->{AF_TSI} . "," : "-1,";
            $ASN     .= (defined $rH_data->{$variantKey}->{AF_ASN})  ? $rH_data->{$variantKey}->{AF_ASN} . "," : "-1,";
            $AFR     .= (defined $rH_data->{$variantKey}->{AF_AFR})  ? $rH_data->{$variantKey}->{AF_AFR} . "," : "-1,";
            $AMR     .= (defined $rH_data->{$variantKey}->{AF_AMR})  ? $rH_data->{$variantKey}->{AF_AMR} . "," : "-1,";
            $EUR     .= (defined $rH_data->{$variantKey}->{AF_EUR})  ? $rH_data->{$variantKey}->{AF_EUR} . "," : "-1,";
            $CG      .= (defined $rH_data->{$variantKey}->{CG})      ? $rH_data->{$variantKey}->{CG} . "," : "-1,";
            $LRT     .= (defined $rH_data->{$variantKey}->{LRT})     ? $rH_data->{$variantKey}->{LRT} . "," : "-1,";
            $SIFT    .= (defined $rH_data->{$variantKey}->{SIFT})    ? $rH_data->{$variantKey}->{SIFT} . "," : "-1,";
            $PPv2    .= (defined $rH_data->{$variantKey}->{PPv2})    ? $rH_data->{$variantKey}->{PPv2} . "," : "-1,";
            $MT      .= (defined $rH_data->{$variantKey}->{MT})      ? $rH_data->{$variantKey}->{MT} . "," : "-1,";
            $PhyloP  .= (defined $rH_data->{$variantKey}->{PhyloP})  ? $rH_data->{$variantKey}->{PhyloP} . "," : "-1,";
            $GERP    .= (defined $rH_data->{$variantKey}->{GERP})    ? $rH_data->{$variantKey}->{GERP} . "," : "-1,";
            $ESP6500 .= (defined $rH_data->{$variantKey}->{ESP6500}) ? $rH_data->{$variantKey}->{ESP6500} . "," : "-1,";
            ((defined $rH_data->{$variantKey}->{Gene}) && (!defined($Gene)))                         and $Gene             .= join("|", @{$rH_data->{$variantKey}->{Gene}}) . ",";
            ((defined $rH_data->{$variantKey}->{CdbSNP_rs}) && (!defined($CdbSNP_rs)))               and $CdbSNP_rs        .= $rH_data->{$variantKey}->{CdbSNP_rs} . ",";
            ((defined $rH_data->{$variantKey}->{NCdbSNP_rs}) && (!defined($NCdbSNP_rs)))             and $NCdbSNP_rs       .= $rH_data->{$variantKey}->{NCdbSNP_rs} . ",";
            ((defined $rH_data->{$variantKey}->{MIRNA}) && (!defined($MIRNA)))                       and $MIRNA            .= $rH_data->{$variantKey}->{MIRNA} . ",";
            ((defined $rH_data->{$variantKey}->{MIRNATARGET}) && (!defined($MIRNATARGET)))           and $MIRNATARGET      .= $rH_data->{$variantKey}->{MIRNATARGET} . ",";
            ((defined $rH_data->{$variantKey}->{RNASEQ_CUFF}) && (!defined($RNASEQ_CUFF)))           and $RNASEQ_CUFF      .= $rH_data->{$variantKey}->{RNASEQ_CUFF} . ",";
            ((defined $rH_data->{$variantKey}->{RNASEQ_SCRI}) && (!defined($RNASEQ_SCRI)))           and $RNASEQ_SCRI      .= $rH_data->{$variantKey}->{RNASEQ_SCRI} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_CHIPSEQ}) && (!defined($GTX_CHIPSEQ)))           and $GTX_CHIPSEQ      .= $rH_data->{$variantKey}->{GTX_CHIPSEQ} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_CLINVAR}) && (!defined($GTX_CLINVAR)))           and $GTX_CLINVAR      .= $rH_data->{$variantKey}->{GTX_CLINVAR} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_COMMON_SNP}) && (!defined($GTX_COMMON_SNP)))     and $GTX_COMMON_SNP   .= $rH_data->{$variantKey}->{GTX_COMMON_SNP} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_COSMIC}) && (!defined($GTX_COSMIC)))             and $GTX_COSMIC       .= $rH_data->{$variantKey}->{GTX_COSMIC} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_CPG}) && (!defined($GTX_CPG)))                   and $GTX_CPG          .= $rH_data->{$variantKey}->{GTX_CPG} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_DBNSFP}) && (!defined($GTX_DBNSFP)))             and $GTX_DBNSFP       .= $rH_data->{$variantKey}->{GTX_DBNSFP} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_DISEASE}) && (!defined($GTX_DISEASE)))           and $GTX_DISEASE      .= $rH_data->{$variantKey}->{GTX_DISEASE} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_DRUG}) && (!defined($GTX_DRUG)))                 and $GTX_DRUG         .= $rH_data->{$variantKey}->{GTX_DRUG} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_EVS}) && (!defined($GTX_EVS)))                   and $GTX_EVS          .= $rH_data->{$variantKey}->{GTX_EVS} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_GWAS}) && (!defined($GTX_GWAS)))                 and $GTX_GWAS         .= $rH_data->{$variantKey}->{GTX_GWAS} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_HGMD}) && (!defined($GTX_HGMD)))                 and $GTX_HGMD         .= $rH_data->{$variantKey}->{GTX_HGMD} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_HGMD_DISGENE}) && (!defined($GTX_HGMD_DISGENE))) and $GTX_HGMD_DISGENE .= $rH_data->{$variantKey}->{GTX_HGMD_DISGENE} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_MICROSAT}) && (!defined($GTX_MICROSAT)))         and $GTX_MICROSAT     .= $rH_data->{$variantKey}->{GTX_MICROSAT} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_MIRNA}) && (!defined($GTX_MIRNA)))               and $GTX_MIRNA        .= $rH_data->{$variantKey}->{GTX_MIRNA} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_OMIM}) && (!defined($GTX_OMIM)))                 and $GTX_OMIM         .= $rH_data->{$variantKey}->{GTX_OMIM} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_PATH}) && (!defined($GTX_PATH)))                 and $GTX_PATH         .= $rH_data->{$variantKey}->{GTX_PATH} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_PGX}) && (!defined($GTX_PGX)))                   and $GTX_PGX          .= $rH_data->{$variantKey}->{GTX_PGX} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_PTMS}) && (!defined($GTX_PTMS)))                 and $GTX_PTMS         .= $rH_data->{$variantKey}->{GTX_PTMS} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_TRNASFAC}) && (!defined($GTX_TRNASFAC)))         and $GTX_TRNASFAC     .= $rH_data->{$variantKey}->{GTX_TRNASFAC} . ",";
            ((defined $rH_data->{$variantKey}->{GTX_TSS}) && (!defined($GTX_TSS)))                   and $GTX_TSS          .= $rH_data->{$variantKey}->{GTX_TSS} . ",";
        }

        # Remove all trailing comma of the annotation strings
        foreach my $annot ($VC, $VFT, $VT, $Gene, $DA, $CdbSNP_rs, $NCdbSNP_rs, $oneKG, $CG, $SIFT, $PPv2, $LRT, $MT, $PhyloP, $GERP, $ESP6500, $MIRNA, $MIRNATARGET, $TSI, $ASN,
                           $AFR, $AMR, $EUR, $RNASEQ_CUFF, $RNASEQ_SCRI, $GTX_CHIPSEQ, $GTX_CLINVAR, $GTX_COMMON_SNP, $GTX_COSMIC, $GTX_CPG, $GTX_DBNSFP, $GTX_DISEASE, $GTX_DRUG,
                           $GTX_EVS,$GTX_GWAS, $GTX_HGMD, $GTX_HGMD_DISGENE, $GTX_MICROSAT, $GTX_MIRNA, $GTX_OMIM, $GTX_PATH, $GTX_PGX, $GTX_PTMS, $GTX_TRNASFAC, $GTX_TSS) {
            ($annot) and $annot =~ s/,$//;
        }

        #  => when number of annotations has to reflects number of alleles, their empty values must be set to 'NIL' for string fields
        foreach my $annot ($VC, $VT) {
            my @values = split /,/, $annot; # / just to avoid a bug in eclipse display...
            if ($annot . "," eq "NIL," x scalar(@values)) {
                $annot = undef;
            } else {
                $annot =~ s/NIL//g;
            }
        }

        # Custom ranking of VFT and DA values
        foreach my $annot ($VFT, $DA) {
           next if (!$annot);
            my @values = split /,/, $annot; # / just to avoid a bug in eclipse display...
            foreach my $value (@values) {
               next if ($value !~ /\|/);
                my @unsortedStrings = split /\|/, $value; # / just to avoid a bug in eclipse display...
                my $rA_sortedStrings = customRanking( -input => \@unsortedStrings );
                $value = join "|", @{$rA_sortedStrings};
            }
            my $commaCount = ($annot =~ tr/,//);
            $annot = join ",", @values;
            $annot =~ s/NIL/unknown/g;

        }

        if ($Gene) {
            my @GenesByAllele = split /,/, $Gene; # / just to avoid a bug in eclipse display...
            foreach my $alleleGene (@GenesByAllele) {
                if ($alleleGene =~ /\|/) {
                    my @distinctGenes = split /\|/, $alleleGene; # / just to avoid a bug in eclipse display...
                    foreach my $distinctGene (@distinctGenes) {
                        if ($distinctGene =~ /.+\:.+/) {
                            my @geneTab = split /\:/, $distinctGene; # / just to avoid a bug in eclipse display...
                            $distinctGene = $geneTab[0];
                        }
                    }
                    my %distinctGenes = map { $_ => 1 } sort @distinctGenes;
                    $alleleGene = join "|", sort keys %distinctGenes;
                } else {
                    if ($alleleGene =~ /.+\:.+/) {
                        my @geneTab = split /\:/, $alleleGene; # / just to avoid a bug in eclipse display...
                        $alleleGene = $geneTab[0];
                    }
                }
            }
            $Gene = join ",", sort @GenesByAllele;
        }

        #  => when number of annotations has to reflects number of alleles, their empty values must be set to '-1' for float (or interger) fields
        foreach my $annot ($oneKG, $CG, $TSI, $ASN, $AFR, $AMR, $EUR, $LRT, $SIFT, $PPv2, $MT, $PhyloP, $GERP, $ESP6500) {
            my @values = split /,/, $annot; # / just to avoid a bug in eclipse display...
            if (-(scalar(@values)) == eval join '+', @values) {
                $annot = undef;
            } else {
                $annot =~ s/-1//g;
            }
        }
        my $infoString = "";
        foreach my $infoKeys (keys %{$rH_x->{INFO}}) {
            ($rH_x->{INFO}->{$infoKeys}) and $infoString .= $infoKeys . "=" . $rH_x->{INFO}->{$infoKeys} . ";"
        }
        $infoString =~ s/;$//;

        $infoString = $rO_vcf->add_info_field($infoString, 'VC'               => $VC || undef,
                                                           'VFT'              => $VFT || undef,
                                                           'VT'               => $VT || undef,
                                                           'Gene'             => $Gene || undef,
                                                           'DA'               => $DA || undef,
                                                           'CdbSNP_rs'        => $CdbSNP_rs || undef,
                                                           'NCdbSNP_rs'       => $NCdbSNP_rs || undef,
                                                           '1KG'              => $oneKG || undef,
                                                           'AF_TSI'           => $TSI || undef,
                                                           'AF_ASN'           => $ASN || undef,
                                                           'AF_AFR'           => $AFR || undef,
                                                           'AF_AMR'           => $AMR || undef,
                                                           'AF_EUR'           => $EUR || undef,
                                                           'CG'               => $CG || undef,
                                                           'SIFT'             => $SIFT || undef,
                                                           'PPv2'             => $PPv2 || undef,
                                                           'LRT'              => $LRT || undef,
                                                           'MT'               => $MT || undef,
                                                           'PhyloP'           => $PhyloP || undef,
                                                           'GERP++'           => $GERP || undef,
                                                           'ESP6500'          => $ESP6500 || undef,
                                                           'MIRNA'            => $MIRNA || undef,
                                                           'MIRNATARGET'      => $MIRNATARGET || undef,
                                                           'RNASEQ_CUFF'      => $RNASEQ_CUFF || undef,
                                                           'RNASEQ_SCRI'      => $RNASEQ_SCRI || undef,
                                                           'GTX_CHIPSEQ'      => $GTX_CHIPSEQ || undef,
                                                           'GTX_CLINVAR'      => $GTX_CLINVAR || undef,
                                                           'GTX_COMMON_SNP'   => $GTX_COMMON_SNP || undef,
                                                           'GTX_COSMIC'       => $GTX_COSMIC || undef,
                                                           'GTX_CPG'          => $GTX_CPG || undef,
                                                           'GTX_DBNSFP'       => $GTX_DBNSFP || undef,
                                                           'GTX_DISEASE'      => $GTX_DISEASE || undef,
                                                           'GTX_DRUG'         => $GTX_DRUG || undef,
                                                           'GTX_EVS'          => $GTX_EVS || undef,
                                                           'GTX_GWAS'         => $GTX_GWAS || undef,
                                                           'GTX_HGMD'         => $GTX_HGMD || undef,
                                                           'GTX_HGMD_DISGENE' => $GTX_HGMD_DISGENE || undef,
                                                           'GTX_MICROSAT'     => $GTX_MICROSAT || undef,
                                                           'GTX_MIRNA'        => $GTX_MIRNA || undef,
                                                           'GTX_OMIM'         => $GTX_OMIM || undef,
                                                           'GTX_PATH'         => $GTX_PATH || undef,
                                                           'GTX_PGX'          => $GTX_PGX || undef,
                                                           'GTX_PTMS'         => $GTX_PTMS || undef,
                                                           'GTX_TRNASFAC'     => $GTX_TRNASFAC || undef,
                                                           'GTX_TSS'          => $GTX_TSS || undef );
        %{$rH_x->{INFO}} = split /[=;]/, $infoString;

        print RESULTS $rO_vcf->format_line($rH_x);
    }
    $rO_vcf->close();
    close RESULTS;
}

sub doIt {
    my %arg = @_;
    my $element = $arg{-elem};
    (!$EXCLUDE && !$INCLUDE) and return 1;
    if ($INCLUDE) {
        if ($INCLUDE eq "all") {
            (inExcTab(-elem => $element)) ? return 0 : return 1;
        } elsif (inIncTab(-elem => $element)) {
             return 1;
        } else {
            return 0;
        }
    }
    if ($EXCLUDE) {
        if ($EXCLUDE eq "all") {
            ($element eq "gene") and return 1;
            (inIncTab(-elem => $element)) ? return 1 : return 0;
        } elsif (inExcTab(-elem => $element)) {
             return 0;
        } else {
            return 1;
        }
    }
}

sub customRanking {
    my %arg = @_;
    my $rA_input = $arg{-input};

    my @input = @{$rA_input};

    my @FUNCTION_RANKING = qw(exonic exonic_splicing intronic_splicing intronic ncRNA_exonic ncRNA_splicing ncRNA_intronic UTR3_splicing UTR5_splicing UTR3 UTR5 downstream upstream intergenic);
    my %FUNCTION_RANKING_MAP = map { $FUNCTION_RANKING[$_] => $_ } 0 .. $#FUNCTION_RANKING;
    my $RANKING_PATTERN = join '|', @FUNCTION_RANKING;

    my @sorted = reverse sort {
        my ($x, $y) = map /($RANKING_PATTERN)/, $a, $b;
        $FUNCTION_RANKING_MAP{$y} <=> $FUNCTION_RANKING_MAP{$x};
    } @input;
    return \@sorted;
}

sub inIncTab {
    my %arg = @_;
    my $element = $arg{-elem};

    my %hash = map { lc($_) => 1 } @INCLUDETAB;

    return $hash{$element};
}

sub inExcTab {
    my %arg = @_;
    my $element = $arg{-elem};

    my %hash = map { lc($_) => 1 } @EXCLUDETAB;

    return $hash{$element};
}

sub GetCurrentTime {
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year = 1900 + $yearOffset;
    $month++;
    return sprintf "%04d-%02d-%02d_%02dh%02dm%02ds", ($year, $month, $dayOfMonth, $hour, $minute, $second);
}

sub printHeader {
    my %arg = @_;
    my $relultFH = $arg{-fh};

    open(HEAD, "<$HEADER") || die "Could not open file $HEADER...\n$!\n";

    while (my $line = <HEAD>) {
        chomp $line;
        print $relultFH $line, "\n";
    }
}

sub printerr {
    print STDERR @_;
    print LOG @_;
}

=head1 SYNOPSIS

 vcf2Annovar.pl [arguments]

 Mandatory arguments to control input and output files
        -i,  --infile     <file_path>    path of the input file(s)
        -b,  --buildver   <string>       genome build version (hg18 (or 18), hg19 (or 19), b37 or v37)
        -vc, --vcaller    <string>       name of the software (algorithm) that did generate the input files (varscan, diBayes, smallIndel(s), mpileup, gatk, dindel, custom...)
        -id, --identifier <string>       specify the identificator (e.g. sample) on which the script will run (used to create proper working files...)

 Optional arguments:
        -o,  --outfile <file_path>       path of the output file(s) (default: current user directory, filename generated based on variant_caller, identifier...)
        -h,  --header  <file_path>       path of the header file that will be print at the beginning of the output result file (default: no header)
        -d,  --dbsnp   <string>          dbSNP version (snp130, snp131 or snp132)
        -inc,--include <string>          database filter(s) to include in analysis (e.g snp or snp,1kg,cg), cannot be combined with --exclude unless --exclude is set to 'all'
        -e,  --exclude <string>          database filter(s) to exclude from analysis (e.g snp or snp,1kg,cg), cannot be combined with --include unless --include is set to 'all'
             --cpu     <string>          number of CPUs (maximum: 12) that smp-runner will be allow to use when parallel (default) mode is launched (default is 4)
             --serial                    force to execute annovar in serial mode, i.e. one database at a time (by default, annovar calls are made in parallel with smp-runner)
             --help                      print help message
        -m,  --man                       print complete documentation
        -v,  --verbose                   use verbose output

 Function: annotate a vcf input with Annovar, then outputting result in a vcf containing new annotation tags

 Version: $LastChangedDate: 2013-04-04 14:37:31 -0500 (Thu, 04 Apr 2013) $

=head1 OPTIONS

=over 12

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--infile>

path of the input file (generated by varscan or dibayes or smallindels or etc...),
containing all the variants information.
no default value : it has to be provided !

=item B<--outfile>

path of the output results file. By default, the same path as the provided -infile will be used.
each other output (temporary, dropped, filtered, log, etc...) files will be written in the same path as the output file.

=item B<--header>

path of the file containing all the analysis parameters. The content of the file will be
printed as a header at the begining of the output result file.

=item B<--identifier>

identifier of the current analysis.
e.g. : name of the sample we want to analyse (if so, must correspond to S2D database sample).

=item B<--buildver>

genome build version to use. The build version will be used by ANNOVAR to identify
corresponding database files automatically, for example, when gene-based annotation
is used for hg18 build, ANNOVAR will search for the hg18_refGene.txt file, but if
the hg19 is used as --buildver, ANNOVAR will examine hg19_refGene.txt instead.

=item B<--dbsnp>

dbSNP build version to use. The dbSNP version will be used by ANNOVAR to identify
corresponding dbSNP files automatically.
supported choices are : snp130 or 130 (default for hg18 genome build version), snp131 or 131 (default for hg19 genome build version) & snp132 or 132.

=item B<--vcaller>

specifies wich variant caller has been used to genearate the input file.
until now, seven variant callers are supported : varscan, diBayes, smallIndel(s), mpileup, gatk, dindel, custom

=item B<--include>

specifies the database filters which the user wants to include in his analysis, cannot be used with --exclude unless --exclude is set to 'all'.
supported choices are : gene (exclude gene annotation), snp (exclude dbSNP filter), 1kg (exclude 1000 Genome filter), cg (exclude Complete Genomics filter), sift (exclude SIFT filter), pp2 (exclude PolyPhen.v2 filter), mt (exclude Mutation Taster filter), lrt (exclude LRT filter), phylop (exclude PhyloP filter), gerp (exclude gerp++ filter), esp6500 (exclude esp6500 filter), mirna (exclude mirna region annotation), mirnatarget (exclude mirnatarget region annotation), tsi (exclude filter of allele frequency of Toscani population of 1000genomes project), asn (exclude filter of allele frequency of East Asian population of 1000genomes project), afr (exclude filter of allele frequency of African population of 1000genomes project), amr (exclude filter of allele frequency of Ad Mixed American population of 1000genomes project), eur (exclude filter of allele frequency of European population of 1000genomes project), rnaseq_cuff (exclude Broad Cufflinks RNASeq alignment region annotation), rnaseq_scri (exclude Broad Scripture RNASeq alignment region annotation), gtx_chipseq (exclude GenomeTrax Predicted ChIP-Seq TFBS description region annotation), gtx_common_snp (exclude GenomeTrax Common SNPs annotation), gtx_cosmic (exclude GenomeTrax COSMIC (Catalogue of Somatic Mutations in Cancer) somatic disease mutations description region annotation), gtx_cpg (exclude GenomeTrax CpG Islands description region annotation), gtx_dbnsfp (exclude GenomeTrax integrated database of functional predictions from multiple algorithms for the comprehensive collection of human non-synonymous SNPs annotation), gtx_disease (exclude GenomeTrax Disease associations description region annotation), gtx_drug (exclude GenomeTrax Drug Targets description region annotation), gtx_evs (exclude GenomeTrax EVS annotation), gtx_gwas (exclude GenomeTrax GWAS Catalogue description region annotation), gtx_hgmd (exclude GenomeTrax HGMD inherited (germ-line) disease mutations description region annotation), gtx_hgmd_disgene (exclude GenomeTrax HGMD disease genes description region annotation), gtx_microsat (exclude GenomeTrax Microsatellite repeats description region annotation), gtx_mirna (exclude GenomeTrax microRNA sequences description region annotation), gtx_omim (exclude GenomeTrax OMIM description region annotation), gtx_path (exclude GenomeTrax Pathways membership description region annotation), gtx_pgx (exclude GenomeTrax Pharmacogenomics Variants description region annotation), gtx_ptms (exclude GenomeTrax Post translational modifications description region annotation), gtx_transfac (exclude GenomeTrax TRANSFAC experimentally verified TFBS description region annotation), gtx_tss (exclude GenomeTrax TSSs (transcription start sites) description region annotation)

=item B<--exclude>

specifies the database filters which the user wants to exclude from his analysis, cannot be used with --include unless --include is set to 'all'.
supported choices are : gene (exclude gene annotation), snp (exclude dbSNP filter), 1kg (exclude 1000 Genome filter), cg (exclude Complete Genomics filter), sift (exclude SIFT filter), pp2 (exclude PolyPhen.v2 filter), mt (exclude Mutation Taster filter), lrt (exclude LRT filter), phylop (exclude PhyloP filter), gerp (exclude gerp++ filter), esp6500 (exclude esp6500 filter), mirna (exclude mirna region annotation), mirnatarget (exclude mirnatarget region annotation), tsi (exclude filter of allele frequency of Toscani population of 1000genomes project), asn (exclude filter of allele frequency of East Asian population of 1000genomes project), afr (exclude filter of allele frequency of African population of 1000genomes project), amr (exclude filter of allele frequency of Ad Mixed American population of 1000genomes project), eur (exclude filter of allele frequency of European population of 1000genomes project), rnaseq_cuff (exclude Broad Cufflinks RNASeq alignment region annotation), rnaseq_scri (exclude Broad Scripture RNASeq alignment region annotation), gtx_chipseq (exclude GenomeTrax Predicted ChIP-Seq TFBS description region annotation), gtx_common_snp (exclude GenomeTrax Common SNPs annotation), gtx_cosmic (exclude GenomeTrax COSMIC (Catalogue of Somatic Mutations in Cancer) somatic disease mutations description region annotation), gtx_cpg (exclude GenomeTrax CpG Islands description region annotation), gtx_dbnsfp (exclude GenomeTrax integrated database of functional predictions from multiple algorithms for the comprehensive collection of human non-synonymous SNPs annotation), gtx_disease (exclude GenomeTrax Disease associations description region annotation), gtx_drug (exclude GenomeTrax Drug Targets description region annotation), gtx_evs (exclude GenomeTrax EVS annotation), gtx_gwas (exclude GenomeTrax GWAS Catalogue description region annotation), gtx_hgmd (exclude GenomeTrax HGMD inherited (germ-line) disease mutations description region annotation), gtx_hgmd_disgene (exclude GenomeTrax HGMD disease genes description region annotation), gtx_microsat (exclude GenomeTrax Microsatellite repeats description region annotation), gtx_mirna (exclude GenomeTrax microRNA sequences description region annotation), gtx_omim (exclude GenomeTrax OMIM description region annotation), gtx_path (exclude GenomeTrax Pathways membership description region annotation), gtx_pgx (exclude GenomeTrax Pharmacogenomics Variants description region annotation), gtx_ptms (exclude GenomeTrax Post translational modifications description region annotation), gtx_transfac (exclude GenomeTrax TRANSFAC experimentally verified TFBS description region annotation), gtx_tss (exclude GenomeTrax TSSs (transcription start sites) description region annotation)

=cut
