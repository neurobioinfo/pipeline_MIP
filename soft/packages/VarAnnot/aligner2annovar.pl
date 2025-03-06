#!/usr/bin/perl

use File::Basename;
use File::Path;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

#BEGIN { unshift @INC, dirname($0) }
#system "module load bioperl/1.6.1";
use Bio::Tools::GFF;
use Bio::Annotation::Collection;

# Command line parameters
my ($VERBOSE, $HELP, $MAN);
my ($INFILE, $OUTFILE, $ID, $BUILDVER, $VARIANT_CALLER, $OUTFORMAT, $HEADER, $DBSNP, $EXCLUDE);

# Global variables
my $ANNOINFILE;     # path of the annovar input file
my $OUTFILE_ANNO;   # path of the annovar --geneanno output file
my $OUTFILE_DBSNP;  # path of the annovar --filter dbSNP output file
my $OUTFILE_1000G;  # path of the annovar --filter 1000g output file
my $OUTFILE_CG69;   # path of the annovar --filter Complete Genomics output file
my $OUTFILE_SIFT;   # path of the annovar --filter avsift output file
my $OUTFILE_PP2;    # path of the annovar --filter ljb_pp2 output file
my $OUTFILE_MT;     # path of the annovar --filter ljb_mt output file
my $OUTFILE_LRT;    # path of the annovar --filter ljb_lrt output file
my $OUTFILE_PHYLOP; # path of the annovar --filter ljb_phylop output file
my $OUTFILE_GERP;   # path of the annovar --filter ljb_gerp++ output file
my $OUTFILE_ESP5400;# path of the annovar --filter ljb_esp5400 output file
my $OUTPATH;        # path of the results output file (if not provided, will be interpreted from $OUTFILE, or even $INFILE if $OUTFILE is not provided)
my $SUBOUTPATH;     # path of the all the annovar output files (if not provided, will be interpreted from $OUTFILE, or even $INFILE if $OUTFILE is not provided)
my $VARSCAN_TYPE;   # 'SNP', 'INDEL' or 'SNP_INDEL', used to know what kind of variants will be annotated and also appears in the output file name.
my $USEDBUILD;      # build to be use during the annovar analysis (can be different than $BUILDVER, which will bw the one printed in name of output files)
my %PARSE_INPUT = (
    "custom"      => \&parseCustom,
    "dibayes"     => \&parseDiBayes,
    "dindel"      => \&parseDindel,
    "gatk"        => \&parseGATK,
    "mpileup"     => \&parseMpileUp,
    "smallindels" => \&parseSmallIndels,
    "smallindel"  => \&parseSmallIndels,
    "varscan"     => \&parseVarscan,
);
my %OUTPUT_DATA = (
    "csv" => \&printCsv,
    "gff" => \&printGff,
    "vcf" => \&printVcf,
);
my %IUPAC = (
    R => 'AG', Y => 'CT',  S => 'GC',  W => 'AT',  K => 'GT',  M => 'AC', A => 'AA', C => 'CC', G => 'GG',
    T => 'TT', B => 'CGT', D => 'AGT', H => 'ACT', V => 'ACG', N => 'ACGT', '.' => '-', '-' => '-'
);
my $SOURCE;
my $COMMON_HEADER      = "";
my $CUSTOM_HEADER      = "";
my $DEFINITION_HEADER  = "";
my $DIBAYES_HEADER     = "";
my $DINDEL_HEADER      = "";
my $GATK_HEADER        = "";
my $MPILEUP_HEADER     = "";
my $SMALLINDELS_HEADER = "";
my $VARSCAN_HEADER     = "";

&main;

sub main {
    # First process the provided arguments
    processArguments();

    # Parse input data into hash
    $VERBOSE and printerr("\n**********\n* Starting parsing $VARIANT_CALLER input\n**********\n\n");
    my $rH_data = $PARSE_INPUT{$VARIANT_CALLER}->();
    my @header = keys %{$rH_data->{1}};
    $VERBOSE and printerr("\n**********\n* Parsing $VARIANT_CALLER input done...\n**********\n\n");

    # Now proceed to the sequential (-geneanno, dbsnp, 1000g...) executions of Annovar
    launchAnnovar(
        -data => $rH_data
    );

    # Build object to be printed out, specifically to the wanted output format
    $VERBOSE and printerr("\n**********\n* Starting printing $OUTFORMAT output\n**********\n\n");
    $OUTPUT_DATA{$OUTFORMAT}->(
        -data => $rH_data,
    );
    $VERBOSE and printerr("\n**********\n* Printing $OUTFORMAT output done...\n**********\n\n");

    print "\n\nAnnotation of $VARIANT_CALLER $BUILDVER variants for $ID is a great success !!\nThank you for using this program, see you next time !! :)\n\n\n";
}

sub processArguments {
    my @command_line = @ARGV;       # command line argument
    GetOptions('verbose|v'=>\$VERBOSE, 'help'=>\$HELP, 'man|m'=>\$MAN, 'infile|i=s'=>\$INFILE, 'outfile|o=s'=>\$OUTFILE, 'identifier|id=s'=>\$ID, 'buildver|b=s'=>\$BUILDVER, 'vcaller|vc=s'=>\$VARIANT_CALLER, 'outformat|of=s'=>\$OUTFORMAT, 'header|h=s'=>\$HEADER, 'dbsnp|d=s'=>\$DBSNP, 'exclude|e=s'=>\$EXCLUDE) || pod2usage();

    $HELP and pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDOUT);
    $MAN  and pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDOUT);
    @ARGV == 0 || pod2usage("\n**************************************************************\n Syntax error\n\n\n");

    if (not $ID) {
        pod2usage("\n**************************************************************\nError in argument: the --id has to be provided !!\n\n\n");
    }
    $ID = ucfirst $ID;

    if ((not $VARIANT_CALLER) || ($VARIANT_CALLER !~ /^(varscan|dibayes|smallindel|smallindels|mpileup|gatk|dindel|custom)$/i)) {
        pod2usage("\n**************************************************************\nError in argument: the --vcaller must be one of these : varscan, dibayes, smllindels, mpileup, gatk, dindel, custom\n\n\n");
    }
    $VARIANT_CALLER = lc $VARIANT_CALLER;

    if ((not $BUILDVER) || ($BUILDVER !~ /^(hg18|18|hg19|19|b37|v37)$/i)) {
        pod2usage("\n**************************************************************\nError in argument: the --buildver must be one of these : hg18 (or 18), hg19 (or 19), b37, v37\n\n\n");
    } elsif (($BUILDVER !~ /^hg/i) && ($BUILDVER !~ /b37|v37/i)) {
        $BUILDVER = "hg" . $BUILDVER;
    }
    $BUILDVER = lc $BUILDVER;
    $COMMON_HEADER .= "##Human Genome version   " . $BUILDVER . "\n";

    if ($DBSNP) {
        ($DBSNP !~ /[snp130|130|snp131|131|snp132|132]/i) and pod2usage("\n**************************************************************\nError in argument: the --dbsnp must be one of these : snp130(or 130), snp131 (or 131), snp132 (or 132)\n\n\n");
        ($DBSNP !~ /^snp/) and $DBSNP = "snp" . $DBSNP;
    }

    if ($EXCLUDE) {
        ($EXCLUDE !~ /[all|snp|1kg|cg|sift|pp2|mt|lrt|phylop|,]/i) and pod2usage("\n**************************************************************\nError in argument: the --dbsnp must be one of these : snp,1kg,cg,sift,pp2,mt,lrt,phylop\n\n\n");
    }

    my $outFormatNotice = undef;
    if (not $OUTFORMAT) {
        $OUTFORMAT = 'gff';
        $outFormatNotice = "NOTICE : The --outformat is set as 'gff' by default...\n";
    } elsif ($OUTFORMAT !~ /^(csv|gff)$/i) {
        pod2usage("\n**************************************************************\nError in argument: the --outformat must be one of these : csv, gff\n\n\n");
    }
    $OUTFORMAT = lc $OUTFORMAT;

    $VERBOSE ||= "";
    $VERBOSE &&= "--verbose";

    my $headerNotice = undef;
    if ($HEADER) {
        (not -e $HEADER) and $HEADER = undef;
        (not -e $HEADER) and $headerNotice = "NOTICE : the provided --header $HEADER does not exists... no header file will be used during analysis...\n";
        (-d $HEADER) and $HEADER = undef;
        (-d $HEADER) and $headerNotice = "NOTICE : the provided --header is a directory... since it has to be a txt file, no header file will be used during analysis...\n";
    } else {
        $headerNotice = "NOTICE : no --header argument provided : no header file will be used during analysis...\n";
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
        $OUTPATH = $directory;
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

    $ANNOINFILE     = $SUBOUTPATH . $ID . "_ANNOVAR_input." . $VARIANT_CALLER . "." . $BUILDVER;
    $OUTFILE_ANNO   = $SUBOUTPATH . $ID . "_ANNOVAR_geneanno_output." . $VARIANT_CALLER;
    $OUTFILE_DBSNP  = $SUBOUTPATH . $ID . "_ANNOVAR_dbsnp_output." . $VARIANT_CALLER;
    $OUTFILE_1000G  = $SUBOUTPATH . $ID . "_ANNOVAR_1000g_output." . $VARIANT_CALLER;
    $OUTFILE_CG69   = $SUBOUTPATH . $ID . "_ANNOVAR_cg69_output." . $VARIANT_CALLER;
    $OUTFILE_SIFT   = $SUBOUTPATH . $ID . "_ANNOVAR_sift_output." . $VARIANT_CALLER;
    $OUTFILE_PP2    = $SUBOUTPATH . $ID . "_ANNOVAR_pp2_output." . $VARIANT_CALLER;
    $OUTFILE_MT     = $SUBOUTPATH . $ID . "_ANNOVAR_mt_output." . $VARIANT_CALLER;
    $OUTFILE_LRT    = $SUBOUTPATH . $ID . "_ANNOVAR_lrt_output." . $VARIANT_CALLER;
    $OUTFILE_PHYLOP = $SUBOUTPATH . $ID . "_ANNOVAR_phylop_output." . $VARIANT_CALLER;

    my $logFile = $SUBOUTPATH . $ID . "_ANNOVAR_unresolved.log";
    open(LOG, ">$logFile") || die "Could not open log file $logFile...\n($!)\n";

    ($outFormatNotice && $VERBOSE) && printerr($outFormatNotice);
    ($headerNotice    && $VERBOSE) && printerr($headerNotice);
    ($outFileNotice   && $VERBOSE) && printerr($outFileNotice);
}

sub parseCustom {
    my %input;

    open(IN, "<$INFILE") || die "cannot open file $INFILE because: $!\n";

    $CUSTOM_HEADER .= "##INFO=<ID=Alleles,Number=1,Type=String,Description=\"Reference allele/Mutant allele\">\n";
    $CUSTOM_HEADER .= "##INFO=<ID=IL,Number=1,Type=Integer,Description=\"Length of insertion/deletion\">\n";

    $SOURCE = "custom input";

    # populate input file data into local and global hash tables
    my $i = 1;    # line counter, to be used as key in %input
    while (my $line = <IN>) {
        chomp $line;
        if ($line =~ /^(chr)?([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([-ACGT]+)\t([-ACGT]+)(\t)?(.*)?/i) {
            my %variant = (
                'Chrom'  => $2,
                'Start'  => $4,
                'Stop'   => $5,
                'Ref'    => $6,
                'Var'    => $7,
                'Info'   => $9       
            );
            my $key = $variant{Chrom} . ":" .  $variant{Start} . "-" .  $variant{Stop} . "_" .  $variant{Ref} . "/" .  $variant{Var};

            $variant{'Alleles'} = $variant{'Ref'} . "/" . $variant{'Var'};
            if ((length $variant{'Ref'} > 1) || (length $variant{'Var'} > 1) || ($variant{'Ref'} eq "-") || ($variant{'Var'} eq "-")) {
                $variant{'VC'} = 'INDEL';
                my $ref = $variant{'Ref'};
                my $var = $variant{'Var'};

                # if INSERTION
                if ((length $ref < length $var) || ($ref eq "-")) {
                    $variant{'VC'} = 'insertion';
                    $variant{'Ref'} = "-";

                    if ($variant{'Ref'} eq "-") {
                        $variant{'IL'} = length $variant{'Var'};
                    } else {
                        my $firstRef = substr($ref, 0, 1);
                        $ref =~ s/^.//;
    
                        $var =~ /^$firstRef(.*?)$ref$/;
                        my $pattern = $1;
                        $variant{'Var'} = $pattern;
                        $variant{'IL'} = length($pattern);
    
                        $var =~ /^(.*?)$pattern(.*)$/;
                        $variant{'Start'} = eval($variant{'Start'} + length($1));
                    }
                    $variant{'Stop'} = $variant{'Start'}

                # if DELETION
                } elsif ((length $ref > length $var) || ($var eq "-")) {
                    $variant{'VC'} = 'deletion';
                    $variant{'Var'} = "-";

                    if ($variant{'Var'} eq "-") {
                        $variant{'IL'} = length $variant{'Ref'};
                    } else {
                        my $firstVar = substr($var, 0, 1);
                        $var =~ s/^.//;
    
                        $ref =~ /^$firstVar(.*)$var$/;
                        my $pattern = $1;
                        $variant{'Ref'} = $pattern;
                        $variant{'IL'} = length($pattern);
    
                        $ref =~ /^(.*?)$pattern(.*)$/;
                        $variant{'Start'} = eval($variant{'Start'} + length($1));
                    }
                    $variant{'Stop'}  = eval($variant{'Start'} + $variant{'IL'} - 1);
                }
            } else {
                $variant{'VC'} = 'SNP';
            }

            $input{$key} = \%variant;
            $i++;
        } elsif ($line !~ /^#/) {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 280 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(IN);

    return \%input;
}

sub parseDiBayes {
    my %input;

    # open the gff input file
    open(GFFI, "<$INFILE") || die "cannot open file $INFILE because : $!\n";
    my $gffi = Bio::Tools::GFF->new(-fh => \*GFFI, -gff_version => 3);

    # loop over the input stream
    my $i = 1;    # line counter, to be used as key in %input
    while (my $feature = $gffi->next_feature()) {
        (!$SOURCE) and $SOURCE = $feature->source_tag;
        my $chrom = $feature->seq_id;
        $chrom =~ s/^chr//;
        my $reference         = ($feature->get_tag_values('reference')         && $feature->get_tag_values('reference') !~ /^\s+$/)         ? join(',', $feature->get_tag_values('reference')) : undef;
        my $call              = ($feature->get_tag_values('genotype')          && $feature->get_tag_values('genotype') !~ /^\s+$/)          ? join(',', $feature->get_tag_values('genotype')) : undef;
        my $coverage          = ($feature->get_tag_values('coverage')          && $feature->get_tag_values('coverage') !~ /^\s+$/)          ? join(',', $feature->get_tag_values('coverage')) : undef;
        my $refAlleleCounts   = ($feature->get_tag_values('refAlleleCounts')   && $feature->get_tag_values('refAlleleCounts') !~ /^\s+$/)   ? join(',', $feature->get_tag_values('refAlleleCounts')) : undef;
        my $refAlleleStarts   = ($feature->get_tag_values('refAlleleStarts')   && $feature->get_tag_values('refAlleleStarts') !~ /^\s+$/)   ? join(',', $feature->get_tag_values('refAlleleStarts')) : undef;
        my $refAlleleMeanQV   = ($feature->get_tag_values('refAlleleMeanQV')   && $feature->get_tag_values('refAlleleMeanQV') !~ /^\s+$/)   ? join(',', $feature->get_tag_values('refAlleleMeanQV')) : undef;
        my $novelAlleleCounts = ($feature->get_tag_values('novelAlleleCounts') && $feature->get_tag_values('novelAlleleCounts') !~ /^\s+$/) ? join(',', $feature->get_tag_values('novelAlleleCounts')) : undef;
        my $novelAlleleStarts = ($feature->get_tag_values('novelAlleleStarts') && $feature->get_tag_values('novelAlleleStarts') !~ /^\s+$/) ? join(',', $feature->get_tag_values('novelAlleleStarts')) : undef;
        my $novelAlleleMeanQV = ($feature->get_tag_values('novelAlleleMeanQV') && $feature->get_tag_values('novelAlleleMeanQV') !~ /^\s+$/) ? join(',', $feature->get_tag_values('novelAlleleMeanQV')) : undef;
        my $diColor1          = ($feature->get_tag_values('diColor1')          && $feature->get_tag_values('diColor1') !~ /^\s+$/)          ? join(',', $feature->get_tag_values('diColor1')) : undef;
        my $diColor2          = ($feature->get_tag_values('diColor2')          && $feature->get_tag_values('diColor2') !~ /^\s+$/)          ? join(',', $feature->get_tag_values('diColor2')) : undef;
        my $het               = ($feature->get_tag_values('het')               && $feature->get_tag_values('het') !~ /^\s+$/)               ? join(',', $feature->get_tag_values('het')) : undef;
        my $flag              = ($feature->get_tag_values('flag')              && $feature->get_tag_values('flag') !~ /^\s+$/)              ? join(',', $feature->get_tag_values('flag')) : undef;

        $flag =~ s/^\s+//;
        $flag =~ s/\s+$//;
        $flag =~ s/,$//;

        my $obs = $IUPAC{$call} || die "Error script line 319 : invalid best call in <$feature->gff_string>\n";
        my @obs = split (//, $obs);
        @obs == 2 || die "Error script line 321 : observed IUPAC allele $call should correspond to two nucleotide alleles : <$feature->gff_string>\n";

        my $genotype;
        if ($obs[0] eq $reference and $obs[1] eq $reference) {
           die "Error script line 325 : reference alleles are identical to observed alleles : <$feature->gff_string>\n";
        } elsif ($obs[0] eq $reference) {
            $genotype = $obs[1];
        } elsif ($obs[1] eq $reference) {
            $genotype = $obs[0];
        } elsif ($obs[1] ne $obs[0]) {
            $genotype = $obs[0];
            my %snp = (
                'Chrom'             => $chrom,
                'Start'             => $feature->start,
                'Stop'              => $feature->end,
                'Ref'               => $reference,
                'Var'               => $genotype,
                'score'             => $feature->score,
                'strand'            => $feature->strand,
                'frame'             => $feature->frame,
                'VC'                => $feature->primary_tag,
                'observedAlleles'   => $call,
                'coverage'          => $coverage,
                'refAlleleCounts'   => $refAlleleCounts,
                'refAlleleStarts'   => $refAlleleStarts,
                'refAlleleMeanQV'   => $refAlleleMeanQV,
                'novelAlleleCounts' => $novelAlleleCounts,
                'novelAlleleStarts' => $novelAlleleStarts,
                'novelAlleleMeanQV' => $novelAlleleMeanQV,
                'diColor1'          => $diColor1,
                'diColor2'          => $diColor2,
                'het'               => $het,
                'flag'              => $flag,
            );
            my $key = $snp{Chrom} . ":" .  $snp{Start} . "-" .  $snp{Stop} . "_" .  $snp{Ref} . "/" .  $snp{Var};
            $input{$key} = \%snp;
            $i++;
            $genotype = $obs[1];
        } else {
            $genotype = $obs[0];
        }
       die "Error script line 361 : how come there is no observed genotype ? no observed alleles ? : <$feature->gff_string>" if (!$genotype);

        my %snp = (
            'Chrom'             => $chrom,
            'Start'             => $feature->start,
            'Stop'              => $feature->end,
            'Ref'               => $reference,
            'Var'               => $genotype,
            'score'             => $feature->score,
            'strand'            => $feature->strand,
            'frame'             => $feature->frame,
            'VC'                => $feature->primary_tag,
            'observedAlleles'   => $call,
            'coverage'          => $coverage,
            'refAlleleCounts'   => $refAlleleCounts,
            'refAlleleStarts'   => $refAlleleStarts,
            'refAlleleMeanQV'   => $refAlleleMeanQV,
            'novelAlleleCounts' => $novelAlleleCounts,
            'novelAlleleStarts' => $novelAlleleStarts,
            'novelAlleleMeanQV' => $novelAlleleMeanQV,
            'diColor1'          => $diColor1,
            'diColor2'          => $diColor2,
            'het'               => $het,
            'flag'              => $flag,
        );
        my $key = $snp{Chrom} . ":" .  $snp{Start} . "-" .  $snp{Stop} . "_" .  $snp{Ref} . "/" .  $snp{Var};
        $input{$key} = \%snp;
        $i++;
    }
    close(GFFI);

    return \%input;
}

sub parseDindel {
    my %input;

    open(IN, "<$INFILE") || die "cannot open file $INFILE because: $!\n";

    $DINDEL_HEADER .= "##INFO=<ID=Alleles,Number=1,Type=String,Description=\"Reference allele/Mutant allele\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=IL,Number=1,Type=Integer,Description=\"Length of insertion/deletion\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total number of reads in haplotype window\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Reference homopolymer tract length\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=NF,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant on forward strand\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant on reverse strand\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=NFS,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant site on forward strand\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=NRS,Number=1,Type=Integer,Description=\"Number of reads covering non-ref variant site on reverse strand\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=ZY,Number=1,Type=String,Description=\"Zygocity of the variation\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=q20,Number=1,Type=String,Description=\"Quality below 20\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=hp10,Number=1,Type=String,Description=\"Reference homopolymer length was longer than 10\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=fr0,Number=1,Type=String,Description=\"Non-ref allele is not covered by at least one read on both strands\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=wv,Number=1,Type=String,Description=\"Other indel in window had higher likelihood\">\n";
    $DINDEL_HEADER .= "##INFO=<ID=PASS,Number=1,Type=String,Description=\"Pass all the filters\">\n";

    $SOURCE = "Dindel";

    # populate input file data into local and global hash tables
    my $i = 1;      # line counter
    my $count = 1;  # record counter, to be used as key in %input
    while (my $line = <IN>) {
        chomp $line;
        my %variant = ( );
        if ($line =~ /^(chr)?([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(.+)\t([ACGTN]+)\t([ACGTN]+)\t(\d+\.?\d*)\t(.+)\t(.+)\t(.+)\t(.+)/i) {
            my %variant = (
                'Chrom'     => $2,
                'Start'     => $4,
                'Stop'      => $4,
                'Ref'       => $6,
                'Var'       => $7,
                'strand'    => "",
                'score'     => $8,
                'frame'     => "",
                'filter'    => $9,
                'Info'      => $10,
                'GtyFormat' => $11,
                'Genotype'  => $12
            );
            my $key = $variant{Chrom} . ":" .  $variant{Start} . "-" .  $variant{Stop} . "_" .  $variant{Ref} . "/" .  $variant{Var};

           next if ($variant{'Var'} =~ /<DEL>/);
            ($variant{'Var'} =~ /,/) and $variant{'Var'} =~ s/(.+),.+/$1/;

            $variant{'Alleles'} = $variant{'Ref'} . "/" . $variant{'Var'};
            if ((length $variant{'Start'} > 1) || (length $variant{'Stop'} > 1)) {
                $variant{'VC'} = 'INDEL';

                my $ref = $variant{'Ref'};
                my $var = $variant{'Var'};

                # if INSERTION
                if (length $ref < length $var) {
                    $variant{'VC'} = 'insertion';
                    $variant{'Ref'} = "-";

                    my $firstRef = substr($ref, 0, 1);
                    $ref =~ s/^.//;

                    $var =~ /^$firstRef(.*?)$ref$/;
                    my $pattern = $1;
                    $variant{'Var'} = $pattern;
                    $variant{'IL'} = length($pattern);

                    $var =~ /^(.*?)$pattern(.*)$/;
                    $variant{'Start'} = eval($variant{'Start'} + length($1));
                    $variant{'Stop'}  = $variant{'Start'}

                # if DELETION
                } elsif (length $ref > length $var) {
                    $variant{'VC'} = 'deletion';
                    $variant{'Var'} = "-";

                    my $firstVar = substr($var, 0, 1);
                    $var =~ s/^.//;

                    $ref =~ /^$firstVar(.*)$var$/;
                    my $pattern = $1;
                    $variant{'Ref'} = $pattern;
                    $variant{'IL'} = length($pattern);

                    $ref =~ /^(.*?)$pattern(.*)$/;
                    $variant{'Start'} = eval($variant{'Start'} + length($1));
                    $variant{'Stop'}  = eval($variant{'Start'} + $variant{'IL'} - 1);
                }
            } else {
                $variant{'VC'} = 'SNP';
            }

            my %variantInfo = split(/[=;]/, $variant{'Info'}); delete $variant{'Info'};
            foreach my $infoKey (keys %variantInfo) {
                ($infoKey eq "DP")  and $variant{"DP"}  = $variantInfo{$infoKey};
                ($infoKey eq "HP")  and $variant{"HP"}  = $variantInfo{$infoKey};
                ($infoKey eq "NF")  and $variant{"NF"}  = $variantInfo{$infoKey};
                ($infoKey eq "NR")  and $variant{"NR"}  = $variantInfo{$infoKey};
                ($infoKey eq "NFS") and $variant{"NFS"} = $variantInfo{$infoKey};
                ($infoKey eq "NRS") and $variant{"NRS"} = $variantInfo{$infoKey};
            }

            my @variantFilter = split(/;/, $variant{'filter'}); delete $variant{'filter'};
            foreach my $filter (@variantFilter) {
                ($filter eq "q20")  and $variant{"q20"}  = 1;
                ($filter eq "hp10") and $variant{"hp10"} = 1;
                ($filter eq "fr0")  and $variant{"fr0"}  = 1;
                ($filter eq "wv")   and $variant{"wv"}   = 1;
                ($filter eq "PASS") and $variant{"PASS"} = 1;
            }

            my @gtyFormat = split(/:/, $variant{'GtyFormat'}); delete $variant{'GtyFormat'};
            my @gtyValues = split(/:/, $variant{'Genotype'});  delete $variant{'Genotype'};
            my %variantGty = map { $gtyFormat[$_] => $gtyValues[$_] } (0..$#gtyFormat);
            foreach my $gtyKey (keys %variantGty) {
                if ($gtyKey eq "GT") {
                    ($variantGty{$gtyKey} eq "1/2") and $variant{"ZY"} = "MULTIPLE_ALLELES";
                    ($variantGty{$gtyKey} eq "1/1") and $variant{"ZY"} = "hom";
                    ($variantGty{$gtyKey} eq "0/1") and $variant{"ZY"} = "het";
                }
                ($gtyKey eq "GQ") and $variant{"GQ"} = $variantGty{$gtyKey};
            }

           next if (($variant{"ZY"} eq "MULTIPLE_ALLELES") && (($variant{'Var'} !~ /,/)));

            $input{$key} = \%variant;
            $count++;
        } elsif ($line !~ /^#/) {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 523 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(IN);

    return \%input;
}

sub parseGATK {
    my %input;

    open(IN, "<$INFILE") || die "cannot open file $INFILE because: $!\n";

    $GATK_HEADER .= "##INFO=<ID=Alleles,Number=1,Type=String,Description=\"Reference allele/Mutant allele\">\n";
    $GATK_HEADER .= "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n";
    $GATK_HEADER .= "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n";
    $GATK_HEADER .= "##INFO=<ID=IL,Number=1,Type=Integer,Description=\"Length of insertion/deletion\">\n";
    $GATK_HEADER .= "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n";
    $GATK_HEADER .= "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n";
    $GATK_HEADER .= "##INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">\n";
    $GATK_HEADER .= "##INFO=<ID=Dels,Number=1,Type=Float,Description=\"Fraction of Reads Containing Spanning Deletions\">\n";
    $GATK_HEADER .= "##INFO=<ID=HRun,Number=1,Type=Integer,Description=\"Largest Contiguous Homopolymer Run of Variant Allele In Either Direction\">\n";
    $GATK_HEADER .= "##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with at most two segregating haplotypes\">\n";
    $GATK_HEADER .= "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n";
    $GATK_HEADER .= "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">\n";
    $GATK_HEADER .= "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n";
    $GATK_HEADER .= "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">\n";
    $GATK_HEADER .= "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n";
    $GATK_HEADER .= "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n";
    $GATK_HEADER .= "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
    $GATK_HEADER .= "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">\n";
    $GATK_HEADER .= "##INFO=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n";
    $GATK_HEADER .= "##INFO=<ID=ZY,Number=1,Type=String,Description=\"Zygocity of the variation\">\n";
    $GATK_HEADER .= "##INFO=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">\n";

    $SOURCE = "GATK UnifiedGenotyper";

    # populate input file data into local and global hash tables
    my $i = 1;      # line counter
    my $count = 1;  # record counter, to be used as key in %input
    while (my $line = <IN>) {
        chomp $line;
        my %variant = ( );
        ($line =~ /##UnifiedGenotyper=.+/) and $GATK_HEADER .= "##UnifiedGenotyper=\"analysis_type=UnifiedGenotyper input_file=[GATK_Pipeline_FixMateInfoFixMateInfo.S12386.clean.dedup.recal.bam] sample_metadata=[] read_buffer_size=null phone_home=STANDARD read_filter=[] intervals=[1:1-100000000] excludeIntervals=null reference_sequence=/RQexec/dionnela/NGS/GENOME_REF/GATK/human_g1k_v37.fasta rodBind=[] rodToIntervalTrackName=null BTI_merge_rule=UNION nonDeterministicRandomSeed=false DBSNP=null downsampling_type=null downsample_to_fraction=null downsample_to_coverage=null baq=OFF baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false defaultBaseQualities=-1 validation_strictness=SILENT unsafe=null num_threads=6 interval_merging=ALL read_group_black_list=null processingTracker=null restartProcessingTracker=false processingTrackerStatusFile=null processingTrackerID=-1 allow_intervals_with_unindexed_bam=false disable_experimental_low_memory_sharding=false logging_level=INFO log_to_file=null help=false genotype_likelihoods_model=BOTH p_nonref_model=EXACT heterozygosity=0.0010 pcr_error_rate=1.0E-4 genotyping_mode=DISCOVERY output_mode=EMIT_VARIANTS_ONLY standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=30.0 noSLOD=false assume_single_sample_reads=null abort_at_too_much_coverage=-1 min_base_quality_score=17 min_mapping_quality_score=20 max_deletion_fraction=0.05 min_indel_count_for_genotyping=5 indel_heterozygosity=1.25E-4 indelGapContinuationPenalty=10.0 indelGapOpenPenalty=45.0 indelHaplotypeSize=80 doContextDependentGapPenalties=true getGapPenaltiesFromData=false indel_recal_file=indel.recal_data.csv indelDebug=false dovit=false GSA_PRODUCTION_ONLY=false exactCalculation=LINEAR_EXPERIMENTAL output_all_callable_bases=false genotype=false out=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub NO_HEADER=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub debug_file=null metrics_file=null annotation=[]\"\n";
        if ($line =~ /^(chr)?([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(.+)\t([ACGTN]+)\t([ACGTN]+)\t(\d+\.?\d*)\t(.+)\t(.+)\t(.+)\t(.+)/i) {
            my %variant = (
                'Chrom'     => $2,
                'Start'     => $4,
                'Stop'      => $4,
                'Ref'       => $6,
                'Var'       => $7,
                'strand'    => "",
                'score'     => $8,
                'frame'     => "",
                'Info'      => $10,
                'GtyFormat' => $11,
                'Genotype'  => $12,
            );
            my $key = $variant{Chrom} . ":" .  $variant{Start} . "-" .  $variant{Stop} . "_" .  $variant{Ref} . "/" .  $variant{Var};

            $variant{'Alleles'} = $variant{'Ref'} . "/" . $variant{'Var'};
            if ((length $variant{'Ref'} > 1) || (length $variant{'Var'} > 1)) {
                $variant{'VC'} = 'INDEL';
                $variant{'Info'} =~ s/;DS;/;/;
                $variant{'Info'} =~ s/;DB;/;/;

                my $ref = $variant{'Ref'};
                my $var = $variant{'Var'};

                # if INSERTION
                if (length $ref < length $var) {
                    $variant{'VC'} = 'insertion';
                    $variant{'Ref'} = "-";

                    my $firstRef = substr($ref, 0, 1);
                    $ref =~ s/^.//;

                    $var =~ /^$firstRef(.*?)$ref$/;
                    my $pattern = $1;
                    $variant{'Var'} = $pattern;
                    $variant{'IL'} = length($pattern);

                    $var =~ /^(.*?)$pattern(.*)$/;
                    $variant{'Start'} = eval($variant{'Start'} + length($1));
                    $variant{'Stop'}  = $variant{'Start'}

                # if DELETION
                } elsif (length $ref > length $var) {
                    $variant{'VC'} = 'deletion';
                    $variant{'Var'} = "-";

                    my $firstVar = substr($var, 0, 1);
                    $var =~ s/^.//;

                    $ref =~ /^$firstVar(.*)$var$/;
                    my $pattern = $1;
                    $variant{'Ref'} = $pattern;
                    $variant{'IL'} = length($pattern);

                    $ref =~ /^(.*?)$pattern(.*)$/;
                    $variant{'Start'} = eval($variant{'Start'} + length($1));
                    $variant{'Stop'}  = eval($variant{'Start'} + $variant{'IL'} - 1);
                }
            } else {
                $variant{'VC'} = 'SNP';
            }

            if ($variant{'Info'} =~ /;DS;/) {
                $variant{'Info'} =~ s/;DS;/;/;
                $variant{'DS'} = 1;
            }
            if ($variant{'Info'} =~ /;DB;/) {
                $variant{'Info'} =~ s/;DB;/;/;
                $variant{'DB'} = 1;
            }

            my %variantInfo = split(/[=;]/, $variant{'Info'}); delete $variant{'Info'};
            foreach my $infoKey (keys %variantInfo) {
                ($infoKey eq "AC")             and $variant{"AC"}             = $variantInfo{$infoKey};
                ($infoKey eq "AF")             and $variant{"AF"}             = $variantInfo{$infoKey};
                ($infoKey eq "AN")             and $variant{"AN"}             = $variantInfo{$infoKey};
                ($infoKey eq "BaseQRankSum")   and $variant{"BaseQRankSum"}   = $variantInfo{$infoKey};
                ($infoKey eq "DS")             and $variant{"DS"}             = $variantInfo{$infoKey};
                ($infoKey eq "DB")             and $variant{"DB"}             = $variantInfo{$infoKey};
                ($infoKey eq "Dels")           and $variant{"Dels"}           = $variantInfo{$infoKey};
                ($infoKey eq "HRun")           and $variant{"HRun"}           = $variantInfo{$infoKey};
                ($infoKey eq "HaplotypeScore") and $variant{"HaplotypeScore"} = $variantInfo{$infoKey};
                ($infoKey eq "MQ")             and $variant{"MQ"}             = $variantInfo{$infoKey};
                ($infoKey eq "MQ0")            and $variant{"MQ0"}            = $variantInfo{$infoKey};
                ($infoKey eq "MQRankSum")      and $variant{"MQRankSum"}      = $variantInfo{$infoKey};
                ($infoKey eq "QD")             and $variant{"QD"}             = $variantInfo{$infoKey};
                ($infoKey eq "ReadPosRankSum") and $variant{"ReadPosRankSum"} = $variantInfo{$infoKey};
                ($infoKey eq "SB")             and $variant{"SB"}             = $variantInfo{$infoKey};
            }

            my @gtyFormat = split(/:/, $variant{'GtyFormat'}); delete $variant{'GtyFormat'};
            my @gtyValues = split(/:/, $variant{'Genotype'});  delete $variant{'Genotype'};
            my %variantGty = map { $gtyFormat[$_] => $gtyValues[$_] } (0..$#gtyFormat);
            foreach my $gtyKey (keys %variantGty) {
                ($gtyKey eq "GT") and $variant{"ZY"} = ($variantGty{$gtyKey} eq "1/1") ? "hom" : "het";
                ($gtyKey eq "AD") and $variant{"AD"} = $variantGty{$gtyKey};
                ($gtyKey eq "DP") and $variant{"DP"} = $variantGty{$gtyKey};
                ($gtyKey eq "GQ") and $variant{"GQ"} = $variantGty{$gtyKey};
                ($gtyKey eq "PL") and $variant{"PL"} = $variantGty{$gtyKey};    # comma separated gty likehoods : Ref/Ref, Ref/Alt, Alt/Alt
            }

            $input{$key} = \%variant;
            $count++;
        } elsif ($line !~ /^#/) {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 664 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(IN);

    return \%input;
}

sub parseMpileUp {
    my %input;

    open(IN, "<$INFILE") || die "cannot open file $INFILE because: $!\n";

    $MPILEUP_HEADER .= "##INFO=<ID=Alleles,Number=1,Type=String,Description=\"Reference allele/Mutant allele\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=IL,Number=1,Type=Integer,Description=\"Length of insertion/deletion\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=HQfor_ref,Number=1,Type=Integer,Description=\"# high-quality reference-forward bases\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=HQrev_ref,Number=1,Type=Integer,Description=\"# high-quality reference-reverse bases\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=HQfor_mut,Number=1,Type=Integer,Description=\"# high-quality mutant-forward bases\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=HQrev_mut,Number=1,Type=Integer,Description=\"# high-quality mutant-reverse bases\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=FQ,Number=1,Type=Float,Description=\"Phred probability of all samples being the same\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the site allele frequency of the first ALT allele\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=G3,Number=3,Type=Float,Description=\"ML estimate of genotype frequencies\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Chi^2 based HWE test P-value based on G3\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=CI95,Number=2,Type=Float,Description=\"Equal-tail Bayesian credible interval of the site allele frequency at the 95% level\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=strand_bias,Number=1,Type=Float,Description=\"P-values for strand bias\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=baseQ_bias,Number=1,Type=Float,Description=\"P-values for baseQ bias\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=mapQ_bias,Number=1,Type=Float,Description=\"P-values for mapQ bias\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=tail_dist_bias,Number=1,Type=Float,Description=\"P-values for tail distance bias\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=PC2,Number=2,Type=Integer,Description=\"Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=PCHI2,Number=1,Type=Float,Description=\"Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=QCHI2,Number=1,Type=Integer,Description=\"Phred scaled PCHI2.\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=PR,Number=1,Type=Integer,Description=\"# permutations yielding a smaller PCHI2.\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=ZY,Number=1,Type=String,Description=\"Zygocity of the variation\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    $MPILEUP_HEADER .= "##INFO=<ID=PL,Number=-1,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods, number of values is (#ALT+1)*(#ALT+2)/2\">\n";
#    $MPILEUP_HEADER .= "##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)\">\n";
#    $MPILEUP_HEADER .= "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"# high-quality bases\">\n";
#    $MPILEUP_HEADER .= "##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">\n";

    $SOURCE = "mpileup";

    # populate input file data into local and global hash tables
    my $i = 1;      # line counter
    my $count = 1;  # record counter, to be used as key in %input
    while (my $line = <IN>) {
        chomp $line;
        my %variant = ( );
        if ($line =~ /^(chr)?([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(.+)\t([ACGTN]+)\t([ACGTN]+)\t(\d+\.?\d*)\t(.+)\t(.+)\t(.+)\t(.+)/i) {
            my %variant = (
                'Chrom'     => $2,
                'Start'     => $4,
                'Stop'      => $4,
                'Ref'       => $6,
                'Var'       => $7,
                'strand'    => "",
                'score'     => $8,
                'frame'     => "",
                'Info'      => $10,
                'GtyFormat' => $11,
                'Genotype'  => $12,
            );
            my $key = $variant{Chrom} . ":" .  $variant{Start} . "-" .  $variant{Stop} . "_" .  $variant{Ref} . "/" .  $variant{Var};

            $variant{'Alleles'} = $variant{'Ref'} . "/" . $variant{'Var'};
            if ($variant{'Info'} =~ /^INDEL;/) {
                $variant{'VC'} = 'INDEL';
                $variant{'Info'} =~ s/^INDEL;//;

                my $ref = $variant{'Ref'};
                my $var = $variant{'Var'};

                # if INSERTION
                if (length $ref < length $var) {
                    $variant{'VC'} = 'insertion';
                    $variant{'Ref'} = "-";

                    my $firstRef = substr($ref, 0, 1);
                    $ref =~ s/^.//;

                    $var =~ /^$firstRef(.*?)$ref$/;
                    my $pattern = $1;
                    $variant{'Var'} = $pattern;
                    $variant{'IL'} = length($pattern);

                    $var =~ /^(.*?)$pattern(.*)$/;
                    $variant{'Start'} = eval($variant{'Start'} + length($1));
                    $variant{'Stop'}  = $variant{'Start'}

                # if DELETION
                } elsif (length $ref > length $var) {
                    $variant{'VC'} = 'deletion';
                    $variant{'Var'} = "-";

                    my $firstVar = substr($var, 0, 1);
                    $var =~ s/^.//;

                    $ref =~ /^$firstVar(.*)$var$/;
                    my $pattern = $1;
                    $variant{'Ref'} = $pattern;
                    $variant{'IL'} = length($pattern);

                    $ref =~ /^(.*?)$pattern(.*)$/;
                    $variant{'Start'} = eval($variant{'Start'} + length($1));
                    $variant{'Stop'}  = eval($variant{'Start'} + $variant{'IL'} - 1);
                }
            } else {
                $variant{'VC'} = 'SNP';
            }

            my %variantInfo = split(/[=;]/, $variant{'Info'}); delete $variant{'Info'};
            foreach my $infoKey (keys %variantInfo) {
                ($infoKey eq "DP")    and $variant{"DP"}    = $variantInfo{$infoKey};
                ($infoKey eq "MQ")    and $variant{"MQ"}    = $variantInfo{$infoKey};
                ($infoKey eq "FQ")    and $variant{"FQ"}    = $variantInfo{$infoKey};
                ($infoKey eq "AF1")   and $variant{"AF1"}   = $variantInfo{$infoKey};
                ($infoKey eq "G3")    and $variant{"G3"}    = $variantInfo{$infoKey};
                ($infoKey eq "HWE")   and $variant{"HWE"}   = $variantInfo{$infoKey};
                ($infoKey eq "CI95")  and $variant{"CI95"}  = $variantInfo{$infoKey};
                ($infoKey eq "PC2")   and $variant{"PC2"}   = $variantInfo{$infoKey};
                ($infoKey eq "CI95")  and $variant{"CI95"}  = $variantInfo{$infoKey};
                ($infoKey eq "PCHI2") and $variant{"PCHI2"} = $variantInfo{$infoKey};
                ($infoKey eq "PR")    and $variant{"PR"}    = $variantInfo{$infoKey};

#                  ##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">\n";
                ($infoKey eq "DP4") and ($variant{"HQfor_ref"}, $variant{"HQrev_ref"}, $variant{"HQfor_mut"}, $variant{"HQrev_mut"}) = split(/,/,  $variantInfo{$infoKey});

#                  ##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">\n";
                ($infoKey eq "PV4") and ($variant{"strand_bias"}, $variant{"baseQ_bias"}, $variant{"mapQ_bias"}, $variant{"tail_dist_bias"}) = split(/,/,  $variantInfo{$infoKey});
            }

            my @gtyFormat = split(/:/, $variant{'GtyFormat'}); delete $variant{'GtyFormat'};
            my @gtyValues = split(/:/, $variant{'Genotype'});  delete $variant{'Genotype'};
            my %variantGty = map { $gtyFormat[$_] => $gtyValues[$_] } (0..$#gtyFormat);
            foreach my $gtyKey (keys %variantGty) {
                ($gtyKey eq "GT") and $variant{"ZY"} = ($variantGty{$gtyKey} eq "1/1") ? "hom" : "het";
                ($gtyKey eq "GQ") and $variant{"GQ"} = $variantGty{$gtyKey};
                ($gtyKey eq "GL") and $variant{"GL"} = $variantGty{$gtyKey};
                ($gtyKey eq "DP") and $variant{"DP"} = $variantGty{$gtyKey};
                ($gtyKey eq "SP") and $variant{"SP"} = $variantGty{$gtyKey};
                ($gtyKey eq "PL") and $variant{"PL"} = $variantGty{$gtyKey};    # comma separated gty likehoods : Ref/Ref, Ref/Alt, Alt/Alt
            }

            $input{$key} = \%variant;
            $count++;
        } elsif ($line !~ /^#/) {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 810 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(IN);

    return \%input;
}

sub parseSmallIndels {
    my %input;

    # open the gff input file
    open(GFFI, "<$INFILE") || die "cannot open file $INFILE because : $!\n";
    my $gffio = Bio::Tools::GFF->new(-fh => \*GFFI, -gff_version => 3);

    # loop over the input stream
    my $i = 1;      # line counter
    my $count = 1;  # record counter, to be used as key in %input
    while (my $feature = $gffio->next_feature()) {
        (!$SOURCE) and $SOURCE = $feature->source_tag;
        my $chrom = $feature->seq_id;
        $chrom =~ s/^chr//;

        my $clustered_indel_sizes   = ($feature->has_tag('clustered-indel-sizes')   && $feature->get_tag_values('clustered-indel-sizes') !~ /^\s+$/)   ? join(',', $feature->get_tag_values('clustered-indel-sizes')) : undef;
        my $allele_call_pos         = ($feature->has_tag('allele-call-pos')         && $feature->get_tag_values('clustered-indel-sizes') !~ /^\s+$/)   ? join(',', $feature->get_tag_values('allele-call-pos')) : undef;
        my $allele_call             = ($feature->has_tag('allele-call')             && $feature->get_tag_values('allele-call') !~ /^\s+$/)             ? join(',', $feature->get_tag_values('allele-call')) : undef;
        my $allele_pos              = ($feature->has_tag('allele-pos')              && $feature->get_tag_values('allele-pos') !~ /^\s+$/)              ? join(',', $feature->get_tag_values('allele-pos')) : undef;
        my $alleles                 = ($feature->has_tag('alleles')                 && $feature->get_tag_values('alleles') !~ /^\s+$/)                 ? join(',', $feature->get_tag_values('alleles')) : undef;
        my $allele_counts           = ($feature->has_tag('allele-counts')           && $feature->get_tag_values('allele-counts') !~ /^\s+$/)           ? join(',', $feature->get_tag_values('allele-counts')) : undef;
        my $tight_chrom_pos         = ($feature->has_tag('tight_chrom_pos')         && $feature->get_tag_values('tight_chrom_pos') !~ /^\s+$/)         ? join(',', $feature->get_tag_values('tight_chrom_pos')) : undef;
        my $loose_chrom_pos         = ($feature->has_tag('loose_chrom_pos')         && $feature->get_tag_values('loose_chrom_pos') !~ /^\s+$/)         ? join(',', $feature->get_tag_values('loose_chrom_pos')) : undef;
        my $no_nonred_reads         = ($feature->has_tag('no_nonred_reads')         && $feature->get_tag_values('no_nonred_reads') !~ /^\s+$/)         ? join(',', $feature->get_tag_values('no_nonred_reads')) : undef;
        my $run_names               = ($feature->has_tag('run_names')               && $feature->get_tag_values('no_nonred_reads') !~ /^\s+$/)         ? join(',', $feature->get_tag_values('run_names')) : undef;
        my $bead_ids                = ($feature->has_tag('bead_ids')                && $feature->get_tag_values('bead_ids') !~ /^\s+$/)                ? join(',', $feature->get_tag_values('bead_ids')) : undef;
        my $overall_qvs             = ($feature->has_tag('overall_qvs')             && $feature->get_tag_values('bead_ids') !~ /^\s+$/)                ? join(',', $feature->get_tag_values('overall_qvs')) : undef;
        my $no_mismatches           = ($feature->has_tag('no_mismatches')           && $feature->get_tag_values('no_mismatches') !~ /^\s+$/)           ? join(',', $feature->get_tag_values('no_mismatches')) : undef;
        my $read_pos                = ($feature->has_tag('read_pos')                && $feature->get_tag_values('no_mismatches') !~ /^\s+$/)           ? join(',', $feature->get_tag_values('read_pos')) : undef;
        my $from_end_pos            = ($feature->has_tag('from_end_pos')            && $feature->get_tag_values('no_mismatches') !~ /^\s+$/)           ? join(',', $feature->get_tag_values('from_end_pos')) : undef;
        my $strands                 = ($feature->has_tag('strands')                 && $feature->get_tag_values('strands') !~ /^\s+$/)                 ? join(',', $feature->get_tag_values('strands')) : undef;
        my $tags                    = ($feature->has_tag('tags')                    && $feature->get_tag_values('tags') !~ /^\s+$/)                    ? join(',', $feature->get_tag_values('tags')) : undef;
        my $indel_sizes             = ($feature->has_tag('indel_sizes')             && $feature->get_tag_values('indel_sizes') !~ /^\s+$/)             ? join(',', $feature->get_tag_values('indel_sizes')) : undef;
        my $non_indel_no_mismatches = ($feature->has_tag('non_indel_no_mismatches') && $feature->get_tag_values('non_indel_no_mismatches') !~ /^\s+$/) ? join(',', $feature->get_tag_values('non_indel_no_mismatches')) : undef;
        my $unmatched_lengths       = ($feature->has_tag('unmatched-lengths')       && $feature->get_tag_values('unmatched-lengths') !~ /^\s+$/)       ? join(',', $feature->get_tag_values('unmatched-lengths')) : undef;
        my $ave_unmatched           = ($feature->has_tag('ave-unmatched')           && $feature->get_tag_values('unmatched-lengths') !~ /^\s+$/)       ? join(',', $feature->get_tag_values('ave-unmatched')) : undef;
        my $anchor_match_lengths    = ($feature->has_tag('anchor-match-lengths')    && $feature->get_tag_values('unmatched-lengths') !~ /^\s+$/)       ? join(',', $feature->get_tag_values('anchor-match-lengths')) : undef;
        my $ave_anchor_length       = ($feature->has_tag('ave-anchor-length')       && $feature->get_tag_values('ave-anchor-length') !~ /^\s+$/)       ? join(',', $feature->get_tag_values('ave-anchor-length')) : undef;
        my $read_seqs               = ($feature->has_tag('read_seqs')               && $feature->get_tag_values('read_seqs') !~ /^\s+$/)               ? join(',', $feature->get_tag_values('read_seqs')) : undef;
        my $base_qvs                = ($feature->has_tag('base_qvs')                && $feature->get_tag_values('base_qvs') !~ /^\s+$/)                ? join(',', $feature->get_tag_values('base_qvs')) : undef;
        my $non_indel_seqs          = ($feature->has_tag('non_indel_seqs')          && $feature->get_tag_values('non_indel_seqs') !~ /^\s+$/)          ? join(',', $feature->get_tag_values('non_indel_seqs')) : undef;
        my $non_indel_qvs           = ($feature->has_tag('non_indel_qvs')           && $feature->get_tag_values('non_indel_qvs') !~ /^\s+$/)           ? join(',', $feature->get_tag_values('non_indel_qvs')) : undef;

        $strands =~ s/ /+/g;

        my ($reference, $mutant);
        if ($allele_call =~ /^(.+)?\/(.+)?$/) {
            $reference = $1;
            $mutant = $2;
        } else {
            printerr("\nError script line 868 : file line " . $i . " : WRONG ALLELE_CALL SYNTAXE\n" . $feature->gff_string . "\n\n");
            $i++;
           next;
        }
        $reference = "-" if (!$reference );
        $mutant = "-" if (!$mutant);
        if ($reference ne "-" && $mutant ne "-") {
            # todo : to be determined

            printerr("\nError script line 877 : file line " . $i . " : CONFLICT BETWEEN ALLELES\n" . $feature->gff_string . "\n\n");
            $i++;
           next;
        }

        my $variantClass;
        my $indel_len = "";
        if ($feature->primary_tag =~ /^insertion/) {
            $variantClass = "insertion";
            $indel_len = ($feature->get_tag_values('ins_len')) ? join(',', $feature->get_tag_values('ins_len')) : "-";
        } elsif ($feature->primary_tag =~ /^deletion/) {
            $variantClass = "deletion";
            $indel_len = ($feature->get_tag_values('del_len')) ? join(',', $feature->get_tag_values('del_len')) : "-";
        } else {
            ($mutant eq "-" || length($reference) > length($mutant)) and $variantClass = "deletion";
            ($reference eq "-" || length($mutant) > length($reference)) and $variantClass = "insertion";
            if ($variantClass !~ /(insertion|deletion)/) {
                printerr("\nError script line 894 : file line " . $i . " : COULD NOT DETERMINE IF THIS IS INSERTION OR DELETION\n" . $feature->gff_string . "\n\n");
                $i++;
               next;
            }
        }

        my %indel = (
            'Chrom'                   => $chrom,
            'Start'                   => $feature->start,
            'Stop'                    => $feature->end,
            'Ref'                     => $reference,
            'Var'                     => $mutant,
            'score'                   => $feature->score,
            'strand'                  => $feature->strand,
            'frame'                   => $feature->frame,
            'indel_len'               => $indel_len,
            'clustered-indel-sizes'   => $clustered_indel_sizes,
            'allele-call-pos'         => $allele_call_pos,
            'allele-call'             => $allele_call,
            'allele-pos'              => $allele_pos,
            'alleles'                 => $alleles,
            'allele-counts'           => $allele_counts,
            'tight_chrom_pos'         => $tight_chrom_pos,
            'loose_chrom_pos'         => $loose_chrom_pos,
            'no_nonred_reads'         => $no_nonred_reads,
            'run_names'               => $run_names,
            'bead_ids'                => $bead_ids,
            'overall_qvs'             => $overall_qvs,
            'no_mismatches'           => $no_mismatches,
            'read_pos'                => $read_pos,
            'from_end_pos'            => $from_end_pos,
            'strands'                 => $strands,
            'tags'                    => $tags,
            'indel_sizes'             => $indel_sizes,
            'non_indel_no_mismatches' => $non_indel_no_mismatches,
            'unmatched-lengths'       => $unmatched_lengths,
            'ave-unmatched'           => $ave_unmatched,
            'anchor-match-lengths'    => $anchor_match_lengths,
            'ave-anchor-length'       => $ave_anchor_length,
            'read_seqs'               => $read_seqs,
            'base_qvs'                => $base_qvs,
            'non_indel_seqs'          => $non_indel_seqs,
            'non_indel_qvs'           => $non_indel_qvs,
            'VC'                      => $variantClass
        );
        my $key = $indel{Chrom} . ":" .  $indel{Start} . "-" .  $indel{Stop} . "_" .  $indel{Ref} . "/" .  $indel{Var};
        $input{$key} = \%indel;
        $count++;
        $i++;
    }
    close(GFFI);

    return \%input;
}

sub parseVarscan {
    my %input;

    open(IN, "<$INFILE") || die "cannot open file $INFILE because: $!\n";

    $SOURCE = "varscan";

    # populate input file data into local and global hash tables
    my $i = 1;    # line counter, to be used as key in %input
    while (my $line = <IN>) {
        chomp $line;
        if ($line =~ /^(chr)?([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t([ACGT])\t([\+\-]?[ACGT]+)\t(\d+)\t(\d+)\t([\.\d]+\%)\t([012])\t([12])\t(\d+)\t(\d+)\t([\.\d]+)\t([01])\t([01])\t(\d+)\t(\d+)\t(\d+)\t(\d+)/i) {
            my %variant = (
                'Chrom'        => $2,
                'Position'     => $4,
                'Ref'          => $5,
                'Var'          => $6,
                'strand'       => "",
                'score'        => "",
                'frame'        => "",
                'Reads1'       => $7,
                'Reads2'       => $8,
                'VarFreq'      => $9,
                'Strands1'     => $10,
                'Strands2'     => $11,
                'Qual1'        => $12,
                'Qual2'        => $13,
                'Pvalue'       => $14,
                'MapQual1'     => $15,
                'MapQual2'     => $16,
                'Reads1Plus'   => $27,
                'Reads1Minus'  => $28,
                'Reads2Plus'   => $29,
                'Reads2Minus'  => $20,
            );

            $variant{'Start'} = $variant{'Position'};
            $variant{'Stop'}  = $variant{'Position'};
            my $key = $variant{Chrom} . ":" .  $variant{Start} . "-" .  $variant{Stop} . "_" .  $variant{Ref} . "/" .  $variant{Var};

            # if INSERTION
            if ($variant{'Var'} =~ /^\+/) {
                $variant{'Ref'} = "-";
                $variant{'Var'} = substr($variant{'Var'}, 1, length($variant{'Var'}) - 1);
                $variant{'VC'} = 'INDEL';
                (!$VARSCAN_TYPE) and $VARSCAN_TYPE = 'INDEL';
                ($VARSCAN_TYPE eq 'SNP') and $VARSCAN_TYPE .= '_INDEL'

            # if DELETION
            } elsif ($variant{'Var'} =~ /^\-/) {
                $variant{'Start'} = eval($variant{'Position'} + 1);
                $variant{'Stop'}  = eval($variant{'Start'} + length($variant{'Var'}) - 1);
                $variant{'Ref'}   = substr($variant{'Var'}, 1, length($variant{'Var'}) - 1);
                $variant{'Var'}   = "-";
                $variant{'VC'} = 'INDEL';
                (!$VARSCAN_TYPE) and $VARSCAN_TYPE = 'INDEL';
                ($VARSCAN_TYPE eq 'SNP') and $VARSCAN_TYPE .= '_INDEL'
            }

            # if SNP
            (not $variant{'VC'}) and $variant{'VC'} = 'SNP';
            (not $VARSCAN_TYPE) and $VARSCAN_TYPE = 'SNP';

            # if both SNP & INDEL
            if ($VARSCAN_TYPE ne 'SNP_INDEL') {
                (($VARSCAN_TYPE eq 'INDEL') && ($variant{'Var'} !~ /^[\+\-]/)) and $variant{'VC'} = 'SNP_INDEL';
            }

            delete $variant{'Position'};

            $input{$key} = \%variant;
            $i++;
        } elsif ($line !~ /^#/) {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 1020 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(IN);

    return \%input;
}

sub launchAnnovar {
    my %arg = @_;
    my $rH_data = $arg{-data};

    my $annovarPath = "/RQexec/dionnela/soft/packages/VarAnnot";

    # ANNOVAR input consists of chr, start, stop, ref, var
    open(ANNOVAR_IN_MAIN, ">$ANNOINFILE") || die "Could not open input file $ANNOINFILE :\n$!\n";
    foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
        print ANNOVAR_IN_MAIN $rH_data->{$key}->{Chrom} . "\t" . $rH_data->{$key}->{Start} . "\t" . $rH_data->{$key}->{Stop} . "\t" . $rH_data->{$key}->{Ref} . "\t" . $rH_data->{$key}->{Var} . "\n";
    }
    close(ANNOVAR_IN_MAIN);

    if ($VARSCAN_TYPE) {
        $OUTFILE_ANNO   .= "." . $VARSCAN_TYPE;
        $OUTFILE_DBSNP  .= "." . $VARSCAN_TYPE;
        $OUTFILE_1000G  .= "." . $VARSCAN_TYPE;
        $OUTFILE_CG69   .= "." . $VARSCAN_TYPE;
        $OUTFILE_SIFT   .= "." . $VARSCAN_TYPE;
        $OUTFILE_PP2    .= "." . $VARSCAN_TYPE;
        $OUTFILE_MT     .= "." . $VARSCAN_TYPE;
        $OUTFILE_LRT    .= "." . $VARSCAN_TYPE;
        $OUTFILE_PHYLOP .= "." . $VARSCAN_TYPE;
    }
    $OUTFILE_ANNO .= "." . $BUILDVER;

    $USEDBUILD = $BUILDVER;
    ($USEDBUILD =~ /b37|v37/) and $USEDBUILD = "hg19";

    # run ANNOVAR -geneanno, first pass
    $VERBOSE and printerr("\n**********\n* Starting genome annotation using the following command line :\n* perl $annovarPath/annotate_variation.pl -geneanno -buildver $USEDBUILD -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
    system "perl $annovarPath/annotate_variation.pl -geneanno -buildver $USEDBUILD -splicing_threshold 6 -outfile $OUTFILE_ANNO $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
    $VERBOSE and printerr("\n**********\n* Genome annotation done...\n**********\n\n");

    # open ANNOVAR output files for processing
    open(ANNOVAR_MAIN_GENERAL, "<$OUTFILE_ANNO.variant_function") || die "cannot open file $OUTFILE_ANNO.variant_function:\n$!\n";
    open(ANNOVAR_MAIN_DETAILED, "<$OUTFILE_ANNO.exonic_variant_function") || die "cannot open file $OUTFILE_ANNO.exonic_variant_function:\n$!\n";

    $DEFINITION_HEADER .= "##INFO=<ID=VC,Number=1,Type=String,Description=\"Variant class e.g. SNP, INDEL...\">\n";
    $DEFINITION_HEADER .= "##INFO=<ID=VFT,Number=1,Type=String,Description=\"Variant function type e.g. intronic, exonic, intergenic...\">\n";
    $DEFINITION_HEADER .= "##INFO=<ID=Gene,Number=1,Type=String,Description=\"Gene symbol\">\n";

    # annotate noncoding variants - get gene symbol (col 2) and general functional effect (col 1)
    #line14273       synonymous SNV  BGN:NM_001711:exon2:c.G141A:p.S47S,     X       152770230       152770230       G       A

    while (my $line = <ANNOVAR_MAIN_GENERAL>) {
        chomp $line;
        if ($line =~ /^(.+?)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
            my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
            my $varFuncType = $2;
            my $gene = $3;
            $varFuncType =~ s/;/|/g;
            $gene =~ s/;/|/g;
            $gene =~ s/dist=/dist:/g;
            $gene =~ s/,/-/g;
            $rH_data->{$key}->{'VFT'} = $varFuncType || "";
            $rH_data->{$key}->{'Gene'} = $gene;
        } 
        else {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 1081 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(ANNOVAR_MAIN_GENERAL);

    $DEFINITION_HEADER .= "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type e.g. synonymous, nonsynonymous...\">\n";
    $DEFINITION_HEADER .= "##INFO=<ID=DA,Number=1,Type=String,Description=\"Detailed annotation of the variant\">\n";

    # further annotate coding variants - using line number (col 1) to compare to key of global input hash, annotate synonymous/nonsynonymous/frameshift/nonframeshift/insertion/deletion (col 2) and isoform/protein/cdna/AA info (col 3)
    while (my $line = <ANNOVAR_MAIN_DETAILED>) {
        chomp $line;
        if ($line =~ /^line\d+\t(.+?)\t(.+?)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
            my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
            my $variantType = $1;
            my $annotation = $2;
            $annotation =~ s/,\s*$/,/g;
            $rH_data->{$key}->{'VT'} = $variantType || "";
            $rH_data->{$key}->{'VT'} =~ s/ /_/;
            $rH_data->{$key}->{'DA'} = $annotation || "";
        } else {
           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script line 1100 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        }
    }
    close(ANNOVAR_MAIN_DETAILED);


    # if not excluded, run ANNOVAR -filter to output dbSNP variants
    if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /snp/i)) {

        my $dbsnp;
        if ($DBSNP) {
            $dbsnp = $DBSNP;
        } elsif ($USEDBUILD eq "hg18") {
            $dbsnp = "snp130";
        } elsif ($USEDBUILD eq "hg19") {
            $dbsnp = "snp131";
        }
        $COMMON_HEADER .= "##dbSNP version   " . $dbsnp . "\n";
        $DEFINITION_HEADER .= "##INFO=<ID=CdbSNP_rs,Number=1,Type=String,Description=\"dbSNP id of a clinically associated variation\">\n";
        $DEFINITION_HEADER .= "##INFO=<ID=NCdbSNP_rs,Number=1,Type=String,Description=\"dbSNP id of a non-clinically associated variation\">\n";

        # First run annovar with clinically associated dbSNP
        $dbsnp .= "_clinic";
        $VERBOSE and printerr("\n**********\n* Starting filtering against dbSNP with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -outfile $OUTFILE_DBSNP $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
        system("perl $annovarPath/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $annovarPath/humandb/ $VERBOSE");
        $VERBOSE and printerr("\n**********\n* Filtering against dbSNP done...\n**********\n\n");

        # open ANNOVAR output file for processing
        open(ANNOVAR_DBSNP, "<$OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped") || die "cannot open file $OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped:\n$!\n";

        while (my $line = <ANNOVAR_DBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;
            if ($line =~ /(snp\d+_clinic)\t(rs\d+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{'CdbSNP_rs'} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1143 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP);

        # Then run annovar with non-clinically associated dbSNP
        $dbsnp =~ s/_clinic$/_non_clinic/;
        system("perl $annovarPath/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_DBSNP $ANNOINFILE $annovarPath/humandb/ $VERBOSE");
        $VERBOSE and printerr("\n**********\n* Starting filtering against dbSNP with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype $dbsnp -buildver $USEDBUILD -outfile $OUTFILE_DBSNP $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
        $VERBOSE and printerr("\n**********\n* Filtering against dbSNP done...\n**********\n\n");

        # open ANNOVAR output file for processing
        open(ANNOVAR_DBSNP, "<$OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped") || die "cannot open file $OUTFILE_DBSNP.$USEDBUILD\_$dbsnp\_dropped:\n$!\n";

        while (my $line = <ANNOVAR_DBSNP>) {    # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;
            if ($line =~ /(snp\d+_non_clinic)\t(rs\d+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{'NCdbSNP_rs'} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1163 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_DBSNP);
        $dbsnp =~ s/_non_clinic$//;
    }


    # if not excluded, run ANNOVAR -filter to output 1000 genomes variants
    if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /1kg/i)) {
        my $dbtype;
        my $annovarOutput;
        if ($USEDBUILD eq "hg18") {
            $dbtype = "1000g2010jul_ceu";
            $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_CEU.sites.2010_07_dropped";
        } elsif ($USEDBUILD eq "hg19") {
            $dbtype = "1000g2010nov_all";
            $annovarOutput = "$OUTFILE_1000G.$USEDBUILD\_ALL.sites.2010_11_dropped";
        }
        $COMMON_HEADER .= "##1000 Genome release   " . $dbtype . "\n";
        $DEFINITION_HEADER .= "##INFO=<ID=1KG,Number=1,Type=Float,Description=\"Frequency in the 1000 Genome project\">\n";
        $VERBOSE and printerr("\n**********\n* Starting filtering against 1000 Genomes db with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -outfile $OUTFILE_1000G $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
        system "perl $annovarPath/annotate_variation.pl -filter -dbtype $dbtype -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_1000G $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Filtering against 1000 Genomes db done...\n**********\n\n");

        # open ANNOVAR output file for processing
        open(ANNOVAR_1000g, "<$annovarOutput") || die "cannot open file $annovarOutput\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_1000g>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;
            if ($line =~ /($dbtype)\t([01]\.\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{'1KG'} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1198 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_1000g);
    }


    # if not excluded, run ANNOVAR -filter to output Complete Genomics variants
    if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /cg/i)) {
        my $dbfile;
        if ($USEDBUILD eq "hg18") {
            $dbfile = "hg18_cg69.txt";
        } elsif ($USEDBUILD eq "hg19") {
            $dbfile = "hg19_cg69.txt";
        }
        $COMMON_HEADER .= "##Complete Genomics db_file   " . $dbfile . "\n";
        $DEFINITION_HEADER .= "##INFO=<ID=CG69,Number=1,Type=Float,Description=\"Frequency in the Complete Genomics (69) project\">\n";
        $VERBOSE and printerr("\n**********\n* Starting filtering against Complete Genomics db with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -outfile $OUTFILE_CG69 $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
        system "perl $annovarPath/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_CG69 $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Filtering against Complete Genomics db done...\n**********\n\n");

        # open ANNOVAR output file for processing
        open(ANNOVAR_CG69, "<$OUTFILE_CG69.$USEDBUILD\_generic\_dropped") || die "cannot open file $OUTFILE_CG69.$USEDBUILD\_generic\_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_CG69>) {      # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;
            if ($line =~ /(generic)\t([01]\.\d*)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{'CG69'} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1229 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_CG69);
    }


    if (($VARIANT_CALLER !~ /^(smallindel|smallindels)$/i) || ($VARSCAN_TYPE ne "INDEL")) {

        # if not excluded, run ANNOVAR -filter to output SIFT scores
        if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /sift/i)) {
            $COMMON_HEADER .= "##SIFT db_file   " . $USEDBUILD . "_avsift.txt\n";
            $DEFINITION_HEADER .= "##INFO=<ID=SIFT,Number=1,Type=Float,Description=\"SIFT scores\">\n";
            $VERBOSE and printerr("\n**********\n* Starting fetching of SIFT prediction scores with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype avsift -sift_threshold 0 -buildver $USEDBUILD -outfile $OUTFILE_SIFT $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
            system "perl $annovarPath/annotate_variation.pl -filter -dbtype avsift -sift_threshold 0 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_SIFT $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of SIFT prediction scores done...\n**********\n\n");

            # open ANNOVAR output file for processing
            open(ANNOVAR_SIFT, "<$OUTFILE_SIFT.$USEDBUILD\_avsift_dropped") || die "cannot open file $OUTFILE_SIFT.$USEDBUILD\_avsift_dropped\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_SIFT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;
                if ($line =~ /(avsift)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{'SIFT'} = $2;
                } else {
    	           die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1256 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
    	        }
            }
            close(ANNOVAR_SIFT);
        }


        # if not excluded, run ANNOVAR -filter to output PolyPhen.v2 scores
        if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /pp2/i)) {
            $COMMON_HEADER .= "##PolyPhen.v2 db_file   " . $USEDBUILD . "_ljb_pp2.txt\n";
            $DEFINITION_HEADER .= "##INFO=<ID=PPv2,Number=1,Type=Float,Description=\"PolyPhen v2 scores\">\n";
            $VERBOSE and printerr("\n**********\n* Starting fetching of PolyPhen (v2) prediction scores with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_pp2 -buildver $USEDBUILD -outfile $OUTFILE_PP2 $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
            system "perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_pp2 -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PP2 $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of PolyPhen (v2) prediction scores done...\n**********\n\n");

            # open ANNOVAR output file for processing
            open(ANNOVAR_PP2, "<$OUTFILE_PP2.$USEDBUILD\_ljb_pp2_dropped") || die "cannot open file $OUTFILE_PP2.$USEDBUILD\_ljb_pp2_dropped\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_PP2>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;
                if ($line =~ /(ljb_pp2)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{'PPv2'} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1281 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_PP2);
        }


        # if not excluded, run ANNOVAR -filter to output LRT scores
        if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /lrt/i)) {
            $COMMON_HEADER .= "##LRT db_file   " . $USEDBUILD . "_ljb_lrt.txt\n";
            $DEFINITION_HEADER .= "##INFO=<ID=LRT,Number=1,Type=Float,Description=\"LRT scores\">\n";
            $VERBOSE and printerr("\n**********\n* Starting fetching of Mutation Taster predictions with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_lrt -buildver $USEDBUILD -outfile $OUTFILE_LRT $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
            system "perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_lrt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_LRT $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of Mutation Taster predictions done...\n**********\n\n");

            # open ANNOVAR output file for processing
            open(ANNOVAR_LRT, "<$OUTFILE_LRT.$USEDBUILD\_ljb_lrt_dropped") || die "cannot open file $OUTFILE_LRT.$USEDBUILD\_ljb_lrt_dropped\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_LRT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;
                if ($line =~ /(ljb_lrt)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{'LRT'} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1306 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_LRT);
        }


        # if not excluded, run ANNOVAR -filter to output Mutation Taster scores
        if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /mt/i)) {
            $COMMON_HEADER .= "##Mutation Taster db_file   " . $USEDBUILD . "_ljb_mt.txt\n";
            $DEFINITION_HEADER .= "##INFO=<ID=MT,Number=1,Type=Float,Description=\"Mutation Taster scores\">\n";
            $VERBOSE and printerr("\n**********\n* Starting fetching of Mutation Taster predictions with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_mt -buildver $USEDBUILD -outfile $OUTFILE_MT $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
            system "perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_mt -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_MT $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of Mutation Taster predictions done...\n**********\n\n");

            # open ANNOVAR output file for processing
            open(ANNOVAR_MT, "<$OUTFILE_MT.$USEDBUILD\_ljb_mt_dropped") || die "cannot open file $OUTFILE_MT.$USEDBUILD\_ljb_mt_dropped\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_MT>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;
                if ($line =~ /(ljb_mt)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{'MT'} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1331 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_MT);
        }


        # if not excluded, run ANNOVAR -filter to output PhyloP scores
        if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /phylop/i)) {
            $COMMON_HEADER .= "##PhyloP db_file   " . $USEDBUILD . "_ljb_phylop.txt\n";
            $DEFINITION_HEADER .= "##INFO=<ID=PhyloP,Number=1,Type=Float,Description=\"PhyloP scores\">\n";
            $VERBOSE and printerr("\n**********\n* Starting fetching of Mutation Taster predictions with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_phylop -buildver $USEDBUILD -outfile $OUTFILE_PHYLOP $ANNOINFILE $annovarPath/humandb/ $VERBOSE\n**********\n\n");
            system "perl $annovarPath/annotate_variation.pl -filter -dbtype ljb_phylop -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_PHYLOP $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of Mutation Taster predictions done...\n**********\n\n");

            # open ANNOVAR output file for processing
            open(ANNOVAR_PHYLOP, "<$OUTFILE_PHYLOP.$USEDBUILD\_ljb_phylop_dropped") || die "cannot open file $OUTFILE_PHYLOP.$USEDBUILD\_ljb_phylop_dropped\n$!\n";

            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_PHYLOP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;
                if ($line =~ /(ljb_phylop)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{'PhyloP'} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1356 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_PHYLOP);
        }


        # if not excluded, run ANNOVAR -filter to output GERP++ scores
        if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /gerp/i)) {
            my $dbfile;
            if ($USEDBUILD eq "hg18") {
                $dbfile = "hg18_ljb_gerp++.txt";
            } elsif ($USEDBUILD eq "hg19") {
                $dbfile = "hg19_ljb_gerp++.txt";
            }
            $COMMON_HEADER .= "##GERP++ db_file   " . $USEDBUILD . "_ljb_gerp++.txt\n";
            $DEFINITION_HEADER .= "##INFO=<ID=GERP++,Number=1,Type=Float,Description=\"GERP++ scores\">\n";
            $VERBOSE and printerr("\n**********\n* Starting fetching of GERP++ predictions with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -outfile $OUTFILE_GERP $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
            system "perl $annovarPath/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_GERP $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
            $VERBOSE and printerr("\n**********\n* Fetching of GERP++ predictions done...\n**********\n\n");
    
            # open ANNOVAR output file for processing
            open(ANNOVAR_GERP, "<$OUTFILE_GERP.$USEDBUILD\_generic_dropped") || die "cannot open file $OUTFILE_GERP.$USEDBUILD\_ljb_gerp++_dropped\n$!\n";
    
            # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
            while (my $line = <ANNOVAR_GERP>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
                chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
                if ($line =~ /(generic)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                    my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                    $rH_data->{$key}->{GERP} = $2;
                } else {
                   die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1396 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
                }
            }
            close(ANNOVAR_GERP);
        }
    }
    
    
    # if not excluded, run ANNOVAR -filter to output ESP5400 scores
    if ((!$EXCLUDE) || ($EXCLUDE !~ /all/i && $EXCLUDE !~ /esp5400/i)) {
        my $dbfile;
        if ($USEDBUILD eq "hg18") {
            $dbfile = "hg18_esp5400_all.txt";
        } elsif ($USEDBUILD eq "hg19") {
            $dbfile = "hg19_esp5400_all.txt";
        }
        $COMMON_HEADER .= "##ESP5400 db_file   " . $USEDBUILD . "_esp5400.txt\n";
        $DEFINITION_HEADER .= "##INFO=<ID=ESP5400,Number=1,Type=Float,Description=\"ESP5400 frequency\">\n";
        $VERBOSE and printerr("\n**********\n* Starting fetching of ESP5400 frequencies with the following command line :\n* perl $annovarPath/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -outfile $OUTFILE_ESP5400 $ANNOINFILE $ANNOVARPATH/humandb/ $VERBOSE\n**********\n\n");
        system "perl $annovarPath/annotate_variation.pl -filter -dbtype generic -genericdbfile $dbfile -buildver $USEDBUILD -indexfilter_threshold 1 -outfile $OUTFILE_ESP5400 $ANNOINFILE $annovarPath/humandb/ $VERBOSE";
        $VERBOSE and printerr("\n**********\n* Fetching of ESP5400 predictions done...\n**********\n\n");

        # open ANNOVAR output file for processing
        open(ANNOVAR_ESP5400, "<$OUTFILE_ESP5400.$USEDBUILD\_generic_dropped") || die "cannot open file $OUTFILE_ESP5400.$USEDBUILD\_esp5400_dropped\n$!\n";

        # annotate noncoding variants, using the combination of chrom_start_ref in the results hash and in the *.dropped file
        while (my $line = <ANNOVAR_ESP5400>) {         # find value of Lookup in ANNOVAR dbSNP output, and annotate variant in %results
            chomp $line;    #  $1        $2           $3                $4        $5     $6       $7         $8
            if ($line =~ /(generic)\t(.+)\t([0-9XYM(MT)(Un)(GL\d+)]+([_\.].+)?)\t(\d+)\t(\d+)\t([ACGT\-]+)\t([ACGT\-]+)/i) {
                my $key = $3 . ":" . $5 . "-" . $6 . "_" . $7 . "/" . $8;
                $rH_data->{$key}->{ESP5400} = $2;
            } else {
               die "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!! script execution stoped at line 1427 : " . $line . "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            }
        }
        close(ANNOVAR_ESP5400);
    }
}

sub printCsv {
    my %arg = @_;
    my $rH_data = $arg{-data};

    # output results to a results file
    my $resultFile = $OUTFILE;
    (not $resultFile) and $resultFile = $OUTPATH . $ID . "_ANNOVAR_output." . $VARIANT_CALLER . "." . $BUILDVER . ".RESULTS." . $OUTFORMAT;
    open(RESULTS, ">$resultFile") || die "Could not open RESULTS file $resultFile :\n$!\n";

    if ($HEADER) {
        printHeader( -fh => \*RESULTS );
    }
    print RESULTS $COMMON_HEADER;
    print RESULTS $DEFINITION_HEADER;

    if ($VARIANT_CALLER eq "custom") {
        print RESULTS $CUSTOM_HEADER;
        my $resultsHeader = "Chrom\tStart\tStop\tRef\tVar\tAlleles\tIL\tsource\tInfo\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tSIFT\tPPv2\tLRT\tMT\tPhyloP\tGERP++\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;
            $printLine .= ($rH_data->{$key}->{'Chrom'})       ? "$rH_data->{$key}->{'Chrom'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})       ? "$rH_data->{$key}->{'Start'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Stop'})        ? "$rH_data->{$key}->{'Stop'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})         ? "$rH_data->{$key}->{'Ref'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})         ? "$rH_data->{$key}->{'Var'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Alleles'})     ? "$rH_data->{$key}->{'Alleles'}\t"                  : "\t";
            $printLine .= ($rH_data->{$key}->{'IL'})          ? "$rH_data->{$key}->{'IL'}\t"                       : "\t";
            $printLine .= ($SOURCE)                           ? "$SOURCE\t"                                        : "\t";
            $printLine .= ($rH_data->{$key}->{'Info'})        ? "$rH_data->{$key}->{'Info'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})          ? "$rH_data->{$key}->{'VC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})         ? "$rH_data->{$key}->{'VFT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})        ? "$rH_data->{$key}->{'Gene'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})          ? "$rH_data->{$key}->{'VT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})          ? "$rH_data->{$key}->{'DA'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})    ? "$rH_data->{$key}->{'dbSNP_rs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})   ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'})  ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})         ? "$rH_data->{$key}->{'1KG'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})        ? "$rH_data->{$key}->{'CG69'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'SIFT'})        ? "$rH_data->{$key}->{'SIFT'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'PPv2'})        ? "$rH_data->{$key}->{'PPv2'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'LRT'})         ? "$rH_data->{$key}->{'LRT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'MT'})          ? "$rH_data->{$key}->{'MT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PhyloP'})      ? "$rH_data->{$key}->{'PhyloP'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'GERP'})        ? "$rH_data->{$key}->{'GERP'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})     ? "$rH_data->{$key}->{'ESP5400'}\n"                  : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Ref | 4 - Var | 5 - Reads1 | 6 - Reads2 | 7 - VarFreq | 8 - Strands1 | 9 - Strands2 | 10 - Qual1 | 11 - Qual2 | 12 - Pvalue | 13 - MapQual1 | 14 - MapQual2 | 15 - Reads1Plus | 16 - Reads1Minus | 17 - Reads2Plus | 18 - Reads2Minus | 19 - VariantClass | 20 - Variant_function_type | 21 - Gene | 22 - Variant_type | 23 - Detailed_annotation_main | 24 - dbSNP_rs | 25 - 1KG | 26 - CG69 | 27 - SIFT | 28 - PPv2 | 29 - LRT | 30 - MT | 31 - PhyloP |
        }

    } elsif ($VARIANT_CALLER eq "dibayes") {
        print RESULTS $DIBAYES_HEADER;
        my $resultsHeader = "Chrom\tStart\tStop\tRef\tVar\tsource\tscore\tstrand\tframe\tobservedAlleles\tcoverage\trefAlleleCounts\trefAlleleStarts\trefAlleleMeanQV\tnovelAlleleCounts\tnovelAlleleStarts\tnovelAlleleMeanQV\tdiColor1\tdiColor2\thet\tflag\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tSIFT\tPPv2\tLRT\tMT\tPhyloP\tGERP++\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;
            $printLine .= ($rH_data->{$key}->{'Chrom'})             ? "$rH_data->{$key}->{'Chrom'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})             ? "$rH_data->{$key}->{'Start'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Stop'})              ? "$rH_data->{$key}->{'Stop'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})               ? "$rH_data->{$key}->{'Ref'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})               ? "$rH_data->{$key}->{'Var'}\t"                      : "\t";
            $printLine .= ($SOURCE)                                 ? "$SOURCE\t"                                        : "\t";
            $printLine .= ($rH_data->{$key}->{'score'})             ? "$rH_data->{$key}->{'score'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'strand'})            ? "$rH_data->{$key}->{'strand'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'frame'})             ? "$rH_data->{$key}->{'frame'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'observedAlleles'})   ? "$rH_data->{$key}->{'observedAlleles'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'coverage'})          ? "$rH_data->{$key}->{'coverage'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'refAlleleCounts'})   ? "$rH_data->{$key}->{'refAlleleCounts'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'refAlleleStarts'})   ? "$rH_data->{$key}->{'refAlleleStarts'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'refAlleleMeanQV'})   ? "$rH_data->{$key}->{'refAlleleMeanQV'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'novelAlleleCounts'}) ? "$rH_data->{$key}->{'novelAlleleCounts'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'novelAlleleStarts'}) ? "$rH_data->{$key}->{'novelAlleleStarts'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'novelAlleleMeanQV'}) ? "$rH_data->{$key}->{'novelAlleleMeanQV'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'diColor1'})          ? "$rH_data->{$key}->{'diColor1'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'diColor2'})          ? "$rH_data->{$key}->{'diColor2'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'het'})               ? "$rH_data->{$key}->{'het'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'flag'})              ? "$rH_data->{$key}->{'flag'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})                ? "$rH_data->{$key}->{'VC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})               ? "$rH_data->{$key}->{'VFT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})              ? "$rH_data->{$key}->{'Gene'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})                ? "$rH_data->{$key}->{'VT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})                ? "$rH_data->{$key}->{'DA'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})          ? "$rH_data->{$key}->{'dbSNP_rs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})         ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'})        ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})               ? "$rH_data->{$key}->{'1KG'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})              ? "$rH_data->{$key}->{'CG69'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'SIFT'})              ? "$rH_data->{$key}->{'SIFT'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'PPv2'})              ? "$rH_data->{$key}->{'PPv2'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'LRT'})               ? "$rH_data->{$key}->{'LRT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'MT'})                ? "$rH_data->{$key}->{'MT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PhyloP'})            ? "$rH_data->{$key}->{'PhyloP'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'GERP'})              ? "$rH_data->{$key}->{'GERP'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})           ? "$rH_data->{$key}->{'ESP5400'}\n"                  : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Ref | 4 - Var | 5 - source | 6 - score | 7 - strand | 8 - frame | 9 - observedAlleles | 10 - coverage | 11 - refAlleleCounts | 12 - refAlleleStarts | 13 - refAlleleMeanQV | 14 - novelAlleleCounts | 15 - novelAlleleStarts | 16 - novelAlleleMeanQV | 17 - diColor1 | 18 - diColor2 | 19 - het | 20 - flag | 21 - VariantClass | 22 - Variant_function_type | 23 - Gene | 24 - Variant_type | 25 - Detailed_annotation_main | 26 - dbSNP_rs | 27 - 1KG | 28 - CG69| 29 - SIFT| 30 - PPv2| 31 - LRT| 32 - MT| 33 - PhyloP |
        }

    } elsif ($VARIANT_CALLER eq "dindel") {
        print RESULTS $DINDEL_HEADER;
        my $resultsHeader = "Chrom\tStart\tStop\tRef\tVar\tSource\tScore\tStrand\tFrame\tAlleles\tIL\tZY\tDP\tGQ\tHP\tNF\tNR\tNFS\tNRS\tq20\thp10\tfr0\twv\tPASS\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tSIFT\tPPv2\tLRT\tMT\tPhyloP\tGERP++\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;

            $printLine .= ($rH_data->{$key}->{'Chrom'})      ? "$rH_data->{$key}->{'Chrom'}\t"      : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})      ? "$rH_data->{$key}->{'Start'}\t"      : "\t";
            $printLine .= ($rH_data->{$key}->{'Stop'})       ? "$rH_data->{$key}->{'Stop'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})        ? "$rH_data->{$key}->{'Ref'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})        ? "$rH_data->{$key}->{'Var'}\t"        : "\t";
            $printLine .= ($SOURCE)                          ? "$SOURCE\t"                          : "\t";
            $printLine .= ($rH_data->{$key}->{'score'})      ? "$rH_data->{$key}->{'score'}\t"      : "\t";
            $printLine .= ($rH_data->{$key}->{'strand'})     ? "$rH_data->{$key}->{'strand'}\t"     : "\t";
            $printLine .= ($rH_data->{$key}->{'frame'})      ? "$rH_data->{$key}->{'frame'}\t"      : "\t";
            $printLine .= ($rH_data->{$key}->{'Alleles'})    ? "$rH_data->{$key}->{'Alleles'}\t"    : "\t";
            $printLine .= ($rH_data->{$key}->{'IL'})         ? "$rH_data->{$key}->{'IL'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'ZY'})         ? "$rH_data->{$key}->{'ZY'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'DP'})         ? "$rH_data->{$key}->{'DP'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'GQ'})         ? "$rH_data->{$key}->{'GQ'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'HP'})         ? "$rH_data->{$key}->{'HP'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'NF'})         ? "$rH_data->{$key}->{'NF'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'NR'})         ? "$rH_data->{$key}->{'NR'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'NFS'})        ? "$rH_data->{$key}->{'NFS'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'NRS'})        ? "$rH_data->{$key}->{'NRS'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'q20'})        ? "$rH_data->{$key}->{'q20'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'hp10'})       ? "$rH_data->{$key}->{'hp10'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'fr0'})        ? "$rH_data->{$key}->{'fr0'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'wv'})         ? "$rH_data->{$key}->{'wv'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'PASS'})       ? "$rH_data->{$key}->{'PASS'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})         ? "$rH_data->{$key}->{'VC'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})        ? "$rH_data->{$key}->{'VFT'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})       ? "$rH_data->{$key}->{'Gene'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})         ? "$rH_data->{$key}->{'VT'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})         ? "$rH_data->{$key}->{'DA'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})   ? "$rH_data->{$key}->{'dbSNP_rs'}\t"   : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})  ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"  : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'}) ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t" : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})        ? "$rH_data->{$key}->{'1KG'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})       ? "$rH_data->{$key}->{'CG69'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'SIFT'})       ? "$rH_data->{$key}->{'SIFT'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'PPv2'})       ? "$rH_data->{$key}->{'PPv2'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'LRT'})        ? "$rH_data->{$key}->{'LRT'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'MT'})         ? "$rH_data->{$key}->{'MT'}\t"         : "\t";
            $printLine .= ($rH_data->{$key}->{'PhyloP'})     ? "$rH_data->{$key}->{'PhyloP'}\t"     : "\t";
            $printLine .= ($rH_data->{$key}->{'GERP'})       ? "$rH_data->{$key}->{'GERP'}\t"       : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})    ? "$rH_data->{$key}->{'ESP5400'}\n"    : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Stop | 4 - Ref | 5 - Var | 6 - source | 7 - score | 8 - strand | 9 - frame | 10 - Alleles | 11 - INDEL_length | 12 - Zygocity | 13 - Genotype_likelyhoods | 14 - Genotype_Quality | 15 - coverage | 16 - #HQ_forward_reference | 17 - #HQ_reverse_reference | 18 - #HQ_forward_mutant | 19 - #HQ_reverse_mutant | 20 - Mappin_Quality | 21 - FQ | 22 - AF1 | 23 - G3 | 24 - HWE | 25 - CI95 | 26 - Strand_bias | 27 - baseQ_bias | 28 - mapQ_dist_bias | 29 - tail_dist_bias | 30 - VariantClass | 31 - Variant_function_type | 32 - Gene | 33 - Variant_type | 34 - Detailed_annotation_main | 35 - dbSNP_rs | 36 - 1KG | 37 - CG69 | 38 - SIFT | 39 - PPv2 | 40 - LRT | 41 - MT | 42 - PhyloP |
        }

    } elsif ($VARIANT_CALLER eq "gatk") {
        print RESULTS $GATK_HEADER;
        my $resultsHeader = "Chrom\tStart\tStop\tRef\tVar\tSource\tScore\tStrand\tFrame\tAlleles\tIL\tZY\tAC\tAF\tDP\tAN\tBaseQRankSum\tDS\tDels\tHRun\tHaplotypeScore\tMQ\tMQ0\tMQRankSum\tQD\tReadPosRankSum\tSB\tAD\tGQ\tPL\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tSIFT\tPPv2\tLRT\tMT\tPhyloP\tGERP++\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;
            $printLine .= ($rH_data->{$key}->{'Chrom'})          ? "$rH_data->{$key}->{'Chrom'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})          ? "$rH_data->{$key}->{'Start'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Stop'})           ? "$rH_data->{$key}->{'Stop'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})            ? "$rH_data->{$key}->{'Ref'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})            ? "$rH_data->{$key}->{'Var'}\t"                      : "\t";
            $printLine .= ($SOURCE)                              ? "$SOURCE\t"                                        : "\t";
            $printLine .= ($rH_data->{$key}->{'score'})          ? "$rH_data->{$key}->{'score'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'strand'})         ? "$rH_data->{$key}->{'strand'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'frame'})          ? "$rH_data->{$key}->{'frame'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Alleles'})        ? "$rH_data->{$key}->{'Alleles'}\t"                  : "\t";
            $printLine .= ($rH_data->{$key}->{'IL'})             ? "$rH_data->{$key}->{'IL'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'ZY'})             ? "$rH_data->{$key}->{'ZY'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'AC'})             ? "$rH_data->{$key}->{'AC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'AF'})             ? "$rH_data->{$key}->{'AF'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DP'})             ? "$rH_data->{$key}->{'DP'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'AN'})             ? "$rH_data->{$key}->{'AN'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'BaseQRankSum'})   ? "$rH_data->{$key}->{'BaseQRankSum'}\t"             : "\t";
            $printLine .= ($rH_data->{$key}->{'DS'})             ? "$rH_data->{$key}->{'DS'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'Dels'})           ? "$rH_data->{$key}->{'Dels'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'HRun'})           ? "$rH_data->{$key}->{'HRun'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'HaplotypeScore'}) ? "$rH_data->{$key}->{'HaplotypeScore'}\t"           : "\t";
            $printLine .= ($rH_data->{$key}->{'MQ'})             ? "$rH_data->{$key}->{'MQ'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'MQ0'})            ? "$rH_data->{$key}->{'MQ0'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'MQRankSum'})      ? "$rH_data->{$key}->{'MQRankSum'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'QD'})             ? "$rH_data->{$key}->{'QD'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'ReadPosRankSum'}) ? "$rH_data->{$key}->{'ReadPosRankSum'}\t"           : "\t";
            $printLine .= ($rH_data->{$key}->{'SB'})             ? "$rH_data->{$key}->{'SB'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'AD'})             ? "$rH_data->{$key}->{'AD'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'GQ'})             ? "$rH_data->{$key}->{'GQ'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PL'})             ? "$rH_data->{$key}->{'PL'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})             ? "$rH_data->{$key}->{'VC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})            ? "$rH_data->{$key}->{'VFT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})           ? "$rH_data->{$key}->{'Gene'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})             ? "$rH_data->{$key}->{'VT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})             ? "$rH_data->{$key}->{'DA'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})       ? "$rH_data->{$key}->{'dbSNP_rs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})      ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'})     ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})            ? "$rH_data->{$key}->{'1KG'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})           ? "$rH_data->{$key}->{'CG69'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'SIFT'})           ? "$rH_data->{$key}->{'SIFT'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'PPv2'})           ? "$rH_data->{$key}->{'PPv2'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'LRT'})            ? "$rH_data->{$key}->{'LRT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'MT'})             ? "$rH_data->{$key}->{'MT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PhyloP'})         ? "$rH_data->{$key}->{'PhyloP'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'GERP'})           ? "$rH_data->{$key}->{'GERP'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})        ? "$rH_data->{$key}->{'ESP5400'}\n"                  : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Stop | 4 - Ref | 5 - Var | 6 - source | 7 - score | 8 - strand | 9 - frame | 10 - Alleles | 11 - INDEL_length | 12 - Zygocity | 13 - Genotype_likelyhoods | 14 - Genotype_Quality | 15 - coverage | 16 - #HQ_forward_reference | 17 - #HQ_reverse_reference | 18 - #HQ_forward_mutant | 19 - #HQ_reverse_mutant | 20 - Mappin_Quality | 21 - FQ | 22 - AF1 | 23 - G3 | 24 - HWE | 25 - CI95 | 26 - Strand_bias | 27 - baseQ_bias | 28 - mapQ_dist_bias | 29 - tail_dist_bias | 30 - VariantClass | 31 - Variant_function_type | 32 - Gene | 33 - Variant_type | 34 - Detailed_annotation_main | 35 - dbSNP_rs | 36 - 1KG | 37 - CG69 | 38 - SIFT | 39 - PPv2 | 40 - LRT | 41 - MT | 42 - PhyloP |
        }

    } elsif ($VARIANT_CALLER eq "mpileup") {
        print RESULTS $MPILEUP_HEADER;
        my $resultsHeader = "Chrom\tStart\tStop\tRef\tVar\tSource\tScore\tStrand\tFrame\tAlleles\tIL\tZY\tPL\tGQ\tDP\tHQfor_ref\tHQrev_ref\tHQfor_mut\tHQrev_mut\tMQ\tFQ\tAF1\tG3\tHWE\tCI95\tstrand_bias\tbaseQ_bias\tmapQ_bias\ttail_dist_bias\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tSIFT\tPPv2\tLRT\tMT\tPhyloP\tGERP++\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;
            $printLine .= ($rH_data->{$key}->{'Chrom'})          ? "$rH_data->{$key}->{'Chrom'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})          ? "$rH_data->{$key}->{'Start'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Stop'})           ? "$rH_data->{$key}->{'Stop'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})            ? "$rH_data->{$key}->{'Ref'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})            ? "$rH_data->{$key}->{'Var'}\t"                      : "\t";
            $printLine .= ($SOURCE)                              ? "$SOURCE\t"                                        : "\t";
            $printLine .= ($rH_data->{$key}->{'score'})          ? "$rH_data->{$key}->{'score'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'strand'})         ? "$rH_data->{$key}->{'strand'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'frame'})          ? "$rH_data->{$key}->{'frame'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Alleles'})        ? "$rH_data->{$key}->{'Alleles'}\t"                  : "\t";
            $printLine .= ($rH_data->{$key}->{'IL'})             ? "$rH_data->{$key}->{'IL'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'ZY'})             ? "$rH_data->{$key}->{'ZY'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PL'})             ? "$rH_data->{$key}->{'PL'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'GQ'})             ? "$rH_data->{$key}->{'GQ'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DP'})             ? "$rH_data->{$key}->{'DP'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'HQfor_ref'})      ? "$rH_data->{$key}->{'HQfor_ref'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'HQrev_ref'})      ? "$rH_data->{$key}->{'HQrev_ref'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'HQfor_mut'})      ? "$rH_data->{$key}->{'HQfor_mut'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'HQrev_mut'})      ? "$rH_data->{$key}->{'HQrev_mut'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'MQ'})             ? "$rH_data->{$key}->{'MQ'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'FQ'})             ? "$rH_data->{$key}->{'FQ'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'AF1'})            ? "$rH_data->{$key}->{'AF1'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'G3'})             ? "$rH_data->{$key}->{'G3'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'HWE'})            ? "$rH_data->{$key}->{'HWE'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CI95'})           ? "$rH_data->{$key}->{'CI95'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'strand_bias'})    ? "$rH_data->{$key}->{'strand_bias'}\t"              : "\t";
            $printLine .= ($rH_data->{$key}->{'baseQ_bias'})     ? "$rH_data->{$key}->{'baseQ_bias'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'mapQ_bias'})      ? "$rH_data->{$key}->{'mapQ_bias'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'tail_dist_bias'}) ? "$rH_data->{$key}->{'tail_dist_bias'}\t"           : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})             ? "$rH_data->{$key}->{'VC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})            ? "$rH_data->{$key}->{'VFT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})           ? "$rH_data->{$key}->{'Gene'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})             ? "$rH_data->{$key}->{'VT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})             ? "$rH_data->{$key}->{'DA'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})       ? "$rH_data->{$key}->{'dbSNP_rs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})      ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'})     ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})            ? "$rH_data->{$key}->{'1KG'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})           ? "$rH_data->{$key}->{'CG69'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'SIFT'})           ? "$rH_data->{$key}->{'SIFT'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'PPv2'})           ? "$rH_data->{$key}->{'PPv2'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'LRT'})            ? "$rH_data->{$key}->{'LRT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'MT'})             ? "$rH_data->{$key}->{'MT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PhyloP'})         ? "$rH_data->{$key}->{'PhyloP'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'GERP'})           ? "$rH_data->{$key}->{'GERP'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})        ? "$rH_data->{$key}->{'ESP5400'}\t"                  : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Stop | 4 - Ref | 5 - Var | 6 - source | 7 - score | 8 - strand | 9 - frame | 10 - Alleles | 11 - INDEL_length | 12 - Zygocity | 13 - Genotype_likelyhoods | 14 - Genotype_Quality | 15 - coverage | 16 - #HQ_forward_reference | 17 - #HQ_reverse_reference | 18 - #HQ_forward_mutant | 19 - #HQ_reverse_mutant | 20 - Mappin_Quality | 21 - FQ | 22 - AF1 | 23 - G3 | 24 - HWE | 25 - CI95 | 26 - Strand_bias | 27 - baseQ_bias | 28 - mapQ_dist_bias | 29 - tail_dist_bias | 30 - VariantClass | 31 - Variant_function_type | 32 - Gene | 33 - Variant_type | 34 - Detailed_annotation_main | 35 - dbSNP_rs | 36 - 1KG | 37 - CG69 | 38 - SIFT | 39 - PPv2 | 40 - LRT | 41 - MT | 42 - PhyloP |
        }

    } elsif ($VARIANT_CALLER =~ /^(smallindel|smallindels)$/i) {
        my $resultsHeader = "Chrom\tStart\tStop\tRef\tVar\tSource\tScore\tStrand\tFrame\tIndel_len\tClustered_indel_sizes\tAllele_call_pos\tAllele_call\tAllele_pos\tAlleles\tAllele_counts\tTight_chrom_pos\tLoose_chrom_pos\tNo_nonred_reads\tRun_names\tBead_ids\tOverall_qvs\tNo_mismatches\tRead_pos\tFrom_end_pos\tStrands\tTags\tIndel_sizes\tNon_indel_no_mismatches\tUnmatched_lengths\tAve_unmatched\tAnchor_match_lengths\tAve_anchor_length\tRead_seqs\tBase_qvs\tNon_indel_seqs\tNon_indel_qvs\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;
            $printLine .= ($rH_data->{$key}->{'Chrom'})                   ? "$rH_data->{$key}->{'Chrom'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})                   ? "$rH_data->{$key}->{'Start'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Stop'})                    ? "$rH_data->{$key}->{'Stop'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})                     ? "$rH_data->{$key}->{'Ref'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})                     ? "$rH_data->{$key}->{'Var'}\t"                      : "\t";
            $printLine .= ($SOURCE)                                       ? "$SOURCE\t"                                        : "\t";
            $printLine .= ($rH_data->{$key}->{'score'})                   ? "$rH_data->{$key}->{'score'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'strand'})                  ? "$rH_data->{$key}->{'strand'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'frame'})                   ? "$rH_data->{$key}->{'frame'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'indel_len'})               ? "$rH_data->{$key}->{'indel_len'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'clustered-indel-sizes'})   ? "$rH_data->{$key}->{'clustered-indel-sizes'}\t"    : "\t";
            $printLine .= ($rH_data->{$key}->{'allele-call-pos'})         ? "$rH_data->{$key}->{'allele-call-pos'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'allele-call'})             ? "$rH_data->{$key}->{'allele-call'}\t"              : "\t";
            $printLine .= ($rH_data->{$key}->{'allele-pos'})              ? "$rH_data->{$key}->{'allele-pos'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'alleles'})                 ? "$rH_data->{$key}->{'alleles'}\t"                  : "\t";
            $printLine .= ($rH_data->{$key}->{'allele-counts'})           ? "$rH_data->{$key}->{'allele-counts'}\t"            : "\t";
            $printLine .= ($rH_data->{$key}->{'tight_chrom_pos'})         ? "$rH_data->{$key}->{'tight_chrom_pos'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'loose_chrom_pos'})         ? "$rH_data->{$key}->{'loose_chrom_pos'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'no_nonred_reads'})         ? "$rH_data->{$key}->{'no_nonred_reads'}\t"          : "\t";
            $printLine .= ($rH_data->{$key}->{'run_names'})               ? "$rH_data->{$key}->{'run_names'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'bead_ids'})                ? "$rH_data->{$key}->{'bead_ids'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'overall_qvs'})             ? "$rH_data->{$key}->{'overall_qvs'}\t"              : "\t";
            $printLine .= ($rH_data->{$key}->{'no_mismatches'})           ? "$rH_data->{$key}->{'no_mismatches'}\t"            : "\t";
            $printLine .= ($rH_data->{$key}->{'read_pos'})                ? "$rH_data->{$key}->{'read_pos'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'from_end_pos'})            ? "$rH_data->{$key}->{'from_end_pos'}\t"             : "\t";
            $printLine .= ($rH_data->{$key}->{'strands'})                 ? "$rH_data->{$key}->{'strands'}\t"                  : "\t";
            $printLine .= ($rH_data->{$key}->{'tags'})                    ? "$rH_data->{$key}->{'tags'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'indel_sizes'})             ? "$rH_data->{$key}->{'indel_sizes'}\t"              : "\t";
            $printLine .= ($rH_data->{$key}->{'non_indel_no_mismatches'}) ? "$rH_data->{$key}->{'non_indel_no_mismatches'}\t"  : "\t";
            $printLine .= ($rH_data->{$key}->{'unmatched-lengths'})       ? "$rH_data->{$key}->{'unmatched-lengths'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'ave-unmatched'})           ? "$rH_data->{$key}->{'ave-unmatched'}\t"            : "\t";
            $printLine .= ($rH_data->{$key}->{'anchor-match-lengths'})    ? "$rH_data->{$key}->{'anchor-match-lengths'}\t"     : "\t";
            $printLine .= ($rH_data->{$key}->{'ave-anchor-length'})       ? "$rH_data->{$key}->{'ave-anchor-length'}\t"        : "\t";
            $printLine .= ($rH_data->{$key}->{'read_seqs'})               ? "$rH_data->{$key}->{'read_seqs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'base_qvs'})                ? "$rH_data->{$key}->{'base_qvs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'non_indel_seqs'})          ? "$rH_data->{$key}->{'non_indel_seqs'}\t"           : "\t";
            $printLine .= ($rH_data->{$key}->{'non_indel_qvs'})           ? "$rH_data->{$key}->{'non_indel_qvs'}\t"            : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})                      ? "$rH_data->{$key}->{'VC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})                     ? "$rH_data->{$key}->{'VFT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})                    ? "$rH_data->{$key}->{'Gene'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})                      ? "$rH_data->{$key}->{'VT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})                      ? "$rH_data->{$key}->{'DA'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})                ? "$rH_data->{$key}->{'dbSNP_rs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})               ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'})              ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})                     ? "$rH_data->{$key}->{'1KG'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})                    ? "$rH_data->{$key}->{'CG69'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})                 ? "$rH_data->{$key}->{'ESP5400'}\n"                  : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Ref | 4 - Var | 5 - source | 6 - score | 7 - strand | 8 - frame | 9 - observedAlleles | 10 - coverage | 11 - refAlleleCounts | 12 - refAlleleStarts | 13 - refAlleleMeanQV | 14 - novelAlleleCounts | 15 - novelAlleleStarts | 16 - novelAlleleMeanQV | 17 - diColor1 | 18 - diColor2 | 19 - het | 20 - flag | 21 - VariantClass | 22 - Variant_function_type | 23 - Gene | 24 - Variant_type | 25 - Detailed_annotation_main | 26 - dbSNP_rs | 27 - 1KG | 28 - CG69 |
        }

    } elsif ($VARIANT_CALLER eq "varscan") {
        my $resultsHeader = "Chrom\tPosition\tRef\tVar\tReads1\tReads2\tVarFreq\tStrands1\tStrands2\tQual1\tQual2\tPvalue\tMapQual1\tMapQual2\tReads1Plus\tReads1Minus\tReads2Plus\tReads2Minus\tVC\tVFT\tGene\tVT\tDA\tdbSNP_rs\t1KG\tCG69\tSIFT\tPPv2\tLRT\tMT\tPhyloP\tGERP++\tESP5400";
        print RESULTS "$resultsHeader\n";
        foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
            my $printLine;
            $printLine .= ($rH_data->{$key}->{'Chrom'})       ? "$rH_data->{$key}->{'Chrom'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Start'})       ? "$rH_data->{$key}->{'Start'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Ref'})         ? "$rH_data->{$key}->{'Ref'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Var'})         ? "$rH_data->{$key}->{'Var'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Reads1'})      ? "$rH_data->{$key}->{'Reads1'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'Reads2'})      ? "$rH_data->{$key}->{'Reads2'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'VarFreq'})     ? "$rH_data->{$key}->{'VarFreq'}\t"                  : "\t";
            $printLine .= ($rH_data->{$key}->{'Strands1'})    ? "$rH_data->{$key}->{'Strands1'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'Strands2'})    ? "$rH_data->{$key}->{'Strands2'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'Qual1'})       ? "$rH_data->{$key}->{'Qual1'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Qual2'})       ? "$rH_data->{$key}->{'Qual2'}\t"                    : "\t";
            $printLine .= ($rH_data->{$key}->{'Pvalue'})      ? "$rH_data->{$key}->{'Pvalue'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'MapQual1'})    ? "$rH_data->{$key}->{'MapQual1'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'MapQual2'})    ? "$rH_data->{$key}->{'MapQual2'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'Reads1Plus'})  ? "$rH_data->{$key}->{'Reads1Plus'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'Reads1Minus'}) ? "$rH_data->{$key}->{'Reads1Plus'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'Reads2Plus'})  ? "$rH_data->{$key}->{'Reads1Plus'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'Reads2Minus'}) ? "$rH_data->{$key}->{'Reads2Minus'}\t"              : "\t";
            $printLine .= ($rH_data->{$key}->{'VC'})          ? "$rH_data->{$key}->{'VC'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'VFT'})         ? "$rH_data->{$key}->{'VFT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'Gene'})        ? "$rH_data->{$key}->{'Gene'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'VT'})          ? "$rH_data->{$key}->{'VT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'DA'})          ? "$rH_data->{$key}->{'DA'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'dbSNP_rs'})    ? "$rH_data->{$key}->{'dbSNP_rs'}\t"                 : "\t";
            $printLine .= ($rH_data->{$key}->{'CdbSNP_rs'})   ? "$rH_data->{$key}->{'CdbSNP_rs'}\t"                : "\t";
            $printLine .= ($rH_data->{$key}->{'NCdbSNP_rs'})  ? "$rH_data->{$key}->{'NCdbSNP_rs'}\t"               : "\t";
            $printLine .= ($rH_data->{$key}->{'1KG'})         ? "$rH_data->{$key}->{'1KG'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'CG69'})        ? "$rH_data->{$key}->{'CG69'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'SIFT'})        ? "$rH_data->{$key}->{'SIFT'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'PPv2'})        ? "$rH_data->{$key}->{'PPv2'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'LRT'})         ? "$rH_data->{$key}->{'LRT'}\t"                      : "\t";
            $printLine .= ($rH_data->{$key}->{'MT'})          ? "$rH_data->{$key}->{'MT'}\t"                       : "\t";
            $printLine .= ($rH_data->{$key}->{'PhyloP'})      ? "$rH_data->{$key}->{'PhyloP'}\t"                   : "\t";
            $printLine .= ($rH_data->{$key}->{'GERP'})        ? "$rH_data->{$key}->{'GERP'}\t"                     : "\t";
            $printLine .= ($rH_data->{$key}->{'ESP5400'})     ? "$rH_data->{$key}->{'ESP5400'}\n"                  : "\n";
            print RESULTS $printLine;
#           | 1 - Chrom | 2 - Start | 3 - Ref | 4 - Var | 5 - Reads1 | 6 - Reads2 | 7 - VarFreq | 8 - Strands1 | 9 - Strands2 | 10 - Qual1 | 11 - Qual2 | 12 - Pvalue | 13 - MapQual1 | 14 - MapQual2 | 15 - Reads1Plus | 16 - Reads1Minus | 17 - Reads2Plus | 18 - Reads2Minus | 19 - VariantClass | 20 - Variant_function_type | 21 - Gene | 22 - Variant_type | 23 - Detailed_annotation_main | 24 - dbSNP_rs | 25 - 1KG | 26 - CG69 | 27 - SIFT | 28 - PPv2 | 29 - LRT | 30 - MT | 31 - PhyloP |
        }

    }
    close RESULTS;
}

sub printGff {
    my %arg = @_;
    my $rH_data = $arg{-data};

    # output results to a results file
    my $resultFile = $OUTFILE;
    (not $resultFile) and $resultFile = $OUTPATH . $ID . "_ANNOVAR_output." . $VARIANT_CALLER . "." . $BUILDVER . ".RESULTS." . $OUTFORMAT;
    open(RESULTS, ">$resultFile") || die "Could not open RESULTS file $resultFile :\n$!\n";

    print RESULTS "##gff-version 3\n";
    if ($HEADER) {
        printHeader( -fh => \*RESULTS );
    }
    print RESULTS $COMMON_HEADER;
    print RESULTS $DEFINITION_HEADER;
    print RESULTS $DIBAYES_HEADER;
    print RESULTS $DINDEL_HEADER;
    print RESULTS $GATK_HEADER;
    print RESULTS $MPILEUP_HEADER;
    print RESULTS $SMALLINDELS_HEADER;
    print RESULTS $VARSCAN_HEADER;

    my $resultsHeader = "#seq_id\tsource_tag\tprimary_tag\tstart\tend\tscore\tstrand\tframe\tattributes";
    print RESULTS "$resultsHeader\n";

    foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
#       | 1 - seq_id | 2 - source_tag | 3 - primary_tag | 4 - start | 5 - end | 6 - score | 7 - strand | 8 - frame | 8 - attributes (starting with ID=uniqueID;) |
        if ($BUILDVER =~ /v37|b37/) {
            print RESULTS  "$rH_data->{$key}->{Chrom}\t";
        } else {
            print RESULTS  "chr$rH_data->{$key}->{Chrom}\t";
        }
        print RESULTS  "$SOURCE\t$rH_data->{$key}->{VC}\t$rH_data->{$key}->{Start}\t$rH_data->{$key}->{Stop}\t$rH_data->{$key}->{score}\t$rH_data->{$key}->{strand}\t$rH_data->{$key}->{frame}\tID=$rH_data->{$key}->{Chrom}_$rH_data->{$key}->{Start}_$rH_data->{$key}->{Ref}/$rH_data->{$key}->{Var}";
        foreach my $head (sort {$a cmp $b} keys %{$rH_data->{$key}}) {
           next if ($head =~ /^(Chrom|source|VC|Start|Stop|score|strand|frame)$/);
            print RESULTS ";" . $head . "=" . $rH_data->{$key}->{$head} if ((defined $rH_data->{$key}->{$head}) && ($rH_data->{$key}->{$head} ne "") && ($rH_data->{$key}->{$head} !~ /^\s*$/));
        }
        print RESULTS "\n";
    }
    close RESULTS
}

sub printVcf {
    my %arg = @_;
    my $rH_data = $arg{-data};

    # output results to a results file
    my $resultFile = $OUTFILE;
    (not $resultFile) and $resultFile = $OUTPATH . $ID . "_ANNOVAR_output." . $VARIANT_CALLER . "." . $BUILDVER . ".RESULTS." . $OUTFORMAT;
    open(RESULTS, ">$resultFile") || die ("Could not open RESULTS file $resultFile :\n$!\n");

    print RESULTS "##fileformat=VCFv4.0\n";
    print RESULTS "##source=$SOURCE\n";
    if ($HEADER) {
        printHeader( -fh => \*RESULTS );
    }
    print RESULTS $COMMON_HEADER;
    print RESULTS $DEFINITION_HEADER;
    print RESULTS $DIBAYES_HEADER;
    print RESULTS $DINDEL_HEADER;
    print RESULTS $GATK_HEADER;
    print RESULTS $MPILEUP_HEADER;
    print RESULTS $SMALLINDELS_HEADER;
    print RESULTS $VARSCAN_HEADER;

    my $info_header;
    my $format_header;
    my $filter_header;
    print RESULTS $info_header;
    print RESULTS $format_header;
    print RESULTS $filter_header;

    my $resultsHeader = "#CHROMt\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
    print RESULTS "$resultsHeader\n";

    foreach my $key (sort {$a cmp $b} keys %{$rH_data}) {
#       | 1 - seq_id | 2 - source_tag | 3 - primary_tag | 4 - start | 5 - end | 6 - score | 7 - strand | 8 - frame | 8 - attributes (starting with ID=uniqueID;) |
        print RESULTS  "$rH_data->{$key}->{Chrom}\t";
        print RESULTS  "$rH_data->{$key}->{Start}\t";
        if ($BUILDVER =~ /v37|b37/) {
            print RESULTS  "$rH_data->{$key}->{Chrom}_$rH_data->{$key}->{Start}_$rH_data->{$key}->{Ref}/$rH_data->{$key}->{Var}\t";
        } else {
            print RESULTS  "chr$rH_data->{$key}->{Chrom}_$rH_data->{$key}->{Start}_$rH_data->{$key}->{Ref}/$rH_data->{$key}->{Var}\t";
        }
        print RESULTS  "$rH_data->{$key}->{Ref}\t";
        print RESULTS  "$rH_data->{$key}->{Var}\t";
        print RESULTS  "$rH_data->{$key}->{score}\t";
        print RESULTS  "$rH_data->{$key}->{Filter}\t";
        print RESULTS  "$rH_data->{$key}->{Info}\t";
        print RESULTS  "$rH_data->{$key}->{Format}\t";
        print RESULTS  "$rH_data->{$key}->{Sample}\n";
        foreach my $head (sort {$a cmp $b} keys %{$rH_data->{$key}}) {
           next if ($head =~ /^(Chrom|Start|Ref|Var|score|Filter|Info|Format|Sample)$/);
            print RESULTS ";" . $head . "=" . $rH_data->{$key}->{$head} if ((defined $rH_data->{$key}->{$head}) && ($rH_data->{$key}->{$head} ne "") && ($rH_data->{$key}->{$head} !~ /^\s*$/));
        }
        print RESULTS "\n";
    }
    close RESULTS
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

 aligner2Annovar.pl [arguments]

 Mandatory arguments to control input and output files
        -i,  --infile     <file_path>    path of the directory that contains the input file(s)
        -of, --outformat  <string>       output file format (csv or gff)
        -b,  --buildver   <string>       genome build version (hg18 (or 18), hg19 (or 19), b37 or v37)
        -vc, --vcaller    <string>       name of the software (algorithm) that did generate the input files (varscan, diBayes, smallIndel(s), mpileup, gatk, dindel, custom...)
        -id, --identifier <string>       specify the identificator (e.g. sample) on which the script will run (used to create proper working files...)

 Optional arguments:
        -o,  --outfile <file_path>       path of the directory that will contain the output file(s) (default: current user directory)
        -h,  --header  <file_path>       path of the header file that will be print at the beginning of the output result file (default: no header)
        -d,  --dbsnp   <string>          dbSNP version (snp130, snp131 or snp132)
        -e,  --exclude <string>          database filter(s) to exclude from analysis (e.g snp or snp,1kg,cg)
             --help                      print help message
        -m,  --man                       print complete documentation
        -v,  --verbose                   use verbose output

 Function: convert a csv (varscan) or gff (diBayes, smallIndels) input to an Annovar readable format,
 then execute Annovar for genome, dbSNP & 100genome annotations.

 Example:
       aligner2Annovar.pl -id=S00017 -b=hg18 -vc=diBayes -i=/home/NGS/S00184/BIOSCOPE.hg18_38Mbp/diBayes/S00184/S00184_SNP.gff -of=gff
       aligner2Annovar.pl -id=DxMol1 -b=hg19 -vc=smallindels -i=/home/NGS/DxMol1/smallindel/smallIndel-pe.gff -of=gff

 Version: $LastChangedDate: 2011-01-17 14:45:31 -0500 (Mon, 17 Jan 2011) $

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
untill now, seven variant callers are supported : varscan, diBayes, smallIndel(s), mpileup, gatk, dindel, custom

=item B<--outformat>

specifies the user desired output format.
supported choices are : gff (default) & csv.

=item B<--exclude>

specifies the database filters which the user wants to exclude from his analysis.
supported choices are : snp (exclude dbSNP filter), 1kg (exclude 1000 Genome filter), cg (exclude Complete Genomics filter), sift (exclude SIFT filter), pp2 (exclude PolyPhen.v2 filter), mt (exclude Mutation Taster filter), lrt (exclude LRT filter), phylop(exclude dbSNP filter)


=cut
