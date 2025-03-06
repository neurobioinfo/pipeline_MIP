#!/usr/bin/perl
######################################################
#written by Ian Byrell Stanaway bard@u.washington.edu#
######################################################
use strict;
my %arg;
&parseCommandLine;	
my $mips_file = $arg{-mips_file};
my $genome_sequence_dir = $arg{-genome_dir};
my $snp_file = $arg{-snp_file};
my @mips_file_contents;
&get_mips_file;

#restriction sites
my $ntalwl_seq_5 = "GGATC";
my $ntalwl_seq_3 = "GATCC";
my $nbbsrdi_seq_5 = "GCAATG";
my $nbbsrdi_seq_3 = "CATTGC";
my $ntalwl_seq_5_count = 0;
my $ntalwl_seq_3_count = 0;
my $nbbsrdi_seq_5_count = 0;
my $nbbsrdi_seq_3_count = 0;

my $bad_restriction_site_design_count = 0;
my $universal_left_end_mip_seq = "AGGACCGGATCAACT";#matches $ntalwl_seq_5 = "GGATC"; once
my $universal_right_end_mip_seq = "CATTGCGTGAACCGA";#matches $nbbsrdi_seq_3 = "CATTGC"; once
my $universal_middle_mip_seq = "CTTCAGCTTCCCGATATCCGACGGTAGTGT";
#example:
#AGGACCGGATCAACTacgcgtgccatctgccacccCTTCAGCTTCCCGATATCCGACGGTAGTGTcaacactcgttttgtgtcccCATTGCGTGAACCGA

my $rank;
my $chr;
my $ext_probe_start;
my $ext_probe_stop;
my $ext_probe_seq;
my $ext_copy_count;
my $lig_probe_start;
my $lig_probe_stop;
my $lig_probe_seq;
my $lig_copy_count;
my $scan_start;
my $scan_stop;
my $scan_seq;
my $feature_start_position;
my $feature_stop_position;
my $feature_mip_count;
my $mip_strand;
my $mip_info;
my $mip_redesign;

my $feature_mip_count = 0;
my $mip_count = 0;
my $snp_position;
my $total_snp_load_count = 0;
my $allele_1;
my $allele_2;
my $chr_fasta = "NA";
my $chr_sequence = "";
my %chr_snp_positions;
open (CHRSNPFILE, ">Chr_snps_size.txt") or die "can't create chr_snps_size\n";
my $outfile = $mips_file . ".fixed_snps";
open (OUTFILE, ">$outfile");

my $line_count = 0;
my $number_of_features = $#mips_file_contents;
foreach my $lines (@mips_file_contents)
{
	chomp $lines;
	if ($lines =~ m/^>|^#/)#header
	{
		print OUTFILE $lines . "\n";
	}
	elsif ($lines !~ m/snp/)
	{
		$mip_count++;
		$lines =~ s/^\d+/$mip_count/;
		print OUTFILE $lines . "\n";
	}
	else
	{
		$line_count++;
		#print $lines . "\n";
		my @line_contents = split(/\s+/, $lines);
		my $tmp_chr = $chr;#save the old one
		$rank = $line_contents[1];
		$chr = $line_contents[2];
		$ext_probe_start = $line_contents[3];
		$ext_probe_stop = $line_contents[4];
		$ext_probe_seq = $line_contents[5];
		$ext_copy_count = $line_contents[6];
		$lig_probe_start = $line_contents[7];
		$lig_probe_stop = $line_contents[8];
		$lig_probe_seq = $line_contents[9];
		$lig_copy_count = $line_contents[10];
		$scan_start = $line_contents[11];
		$scan_stop = $line_contents[12];
		$scan_seq = $line_contents[13];
		$feature_start_position = $line_contents[14];
		$feature_stop_position = $line_contents[15];
		$feature_mip_count = $line_contents[16];
		$mip_strand = $line_contents[17];
		$mip_info = $line_contents[18];

#		my $tmp_chr_fasta = "/nfs/home/bard/ucsc/hg18/chr$chr" . ".fa";
#		if ($tmp_chr_fasta ne $chr_fasta)
#		{
#			$chr_fasta = "/nfs/home/bard/ucsc/hg18/chr$chr" . ".fa";
#BJO
		my $tmp_chr_fasta = "$genome_sequence_dir/chr$chr" . ".fa";
		if ($tmp_chr_fasta ne $chr_fasta)
		{
			$chr_fasta = "$genome_sequence_dir/chr$chr" . ".fa";
			&get_chr_fasta_sequence;

			print "loading chr$chr snps\n";
			delete ($chr_snp_positions{$tmp_chr});
			&load_chr_snps;
		}
		&redesign_mip;
	}
}
close OUTFILE;
print "loaded total $total_snp_load_count snps\nbad_restriction_site_design_count\t$bad_restriction_site_design_count\n";
print CHRSNPFILE "loaded total $total_snp_load_count snps\n";
close CHRSNPFILE;

sub redesign_mip
{
	if ($mip_strand eq "+")
	{
		my $snp_count = 0;
		for my $int ($ext_probe_start .. $ext_probe_stop)
		{
			if (exists $chr_snp_positions{$chr}{$int})
			{
				$snp_count++;
				$snp_position = $int;
			}
		}
		for my $int ($lig_probe_start .. $lig_probe_stop)
		{
			if (exists $chr_snp_positions{$chr}{$int})
			{
				$snp_count++;
				$snp_position = $int;
			}
		}
		my $extension_arm_length = $ext_probe_stop - $ext_probe_start + 1;
		print $ext_probe_seq . "ext\n";
		$ext_probe_seq = uc(substr($chr_sequence,$ext_probe_start - 1,$extension_arm_length));
		print $ext_probe_seq . "ext\n";
		print $lig_probe_seq . "lig\n";
		my $ligation_arm_length = $lig_probe_stop - $lig_probe_start + 1;
		$lig_probe_seq = uc(substr($chr_sequence,$lig_probe_start - 1,$ligation_arm_length));
		print $lig_probe_seq . "lig\n";
		if ($snp_count == 0)
		{
			$mip_count++;
			print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t+\tmissing_snp_now\n";
			print "$line_count/$number_of_features\t$mip_count\n";
		}
		elsif ($snp_count == 1)
		{#make for allele_1
			#make two mips, on the same strand for each allele of the single snp found
			&get_alleles_for_position;#$snp_position
			if ($snp_position >= $ext_probe_start && $snp_position <= $ext_probe_stop)
			{
				my $allele_probe_seq_position = $snp_position - $ext_probe_start;#0 based for array
				my @arr_ext_probe_sequence =  split(//, $ext_probe_seq);
				#make allele 1
				$mip_count++;
				$arr_ext_probe_sequence[$allele_probe_seq_position] = $allele_1;
				$ext_probe_seq = "";
				foreach my $bases (@arr_ext_probe_sequence)
				{
					$ext_probe_seq = $ext_probe_seq . $bases;
				}
				&restriction_site_screen;		
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{			
					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $allele_probe_seq_position + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t+\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}
				#make allele 2
				$mip_count++;
				$arr_ext_probe_sequence[$allele_probe_seq_position] = $allele_2;
				$ext_probe_seq = "";
				foreach my $bases (@arr_ext_probe_sequence)
				{
					$ext_probe_seq = $ext_probe_seq . $bases;
				}

				&restriction_site_screen;		
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{
					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $allele_probe_seq_position + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t+\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}
			}
			elsif ($snp_position >= $lig_probe_start && $snp_position <= $lig_probe_stop)						
			{#must be between the ligerse target arms

				my $allele_probe_seq_position = $snp_position - $lig_probe_start;#0 based for array
				my @arr_lig_probe_sequence =  split(//, $lig_probe_seq);
				#make allele 1
				$mip_count++;
				$arr_lig_probe_sequence[$allele_probe_seq_position] = $allele_1;
				$lig_probe_seq = "";
				foreach my $bases (@arr_lig_probe_sequence)
				{
					$lig_probe_seq = $lig_probe_seq . $bases;
				}

				&restriction_site_screen;		
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{

					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $allele_probe_seq_position + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t+\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}

				#make allele 2
				$mip_count++;
				$arr_lig_probe_sequence[$allele_probe_seq_position] = $allele_2;
				$lig_probe_seq = "";
				foreach my $bases (@arr_lig_probe_sequence)
				{
					$lig_probe_seq = $lig_probe_seq . $bases;
				}
				&restriction_site_screen;
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{
					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $allele_probe_seq_position + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t+\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}
			}
		}
		else
		{
			#has extra snps
		}
	}
	else#is - strand
	{
		my $extension_arm_length = $ext_probe_stop - $ext_probe_start + 1;
		print $ext_probe_seq . "ext\n";
		$ext_probe_seq = uc(substr($chr_sequence,$ext_probe_start - 1,$extension_arm_length));
		print $ext_probe_seq . "ext\n";
		print $lig_probe_seq . "lig\n";
		my $ligation_arm_length = $lig_probe_stop - $lig_probe_start + 1;
		$lig_probe_seq = uc(substr($chr_sequence,$lig_probe_start - 1,$ligation_arm_length));
		print $lig_probe_seq . "lig\n";

		#do the rev/- strand snps
		my $snp_count = 0;
		for my $int ($ext_probe_start .. $ext_probe_stop)
		{
			if (exists $chr_snp_positions{$chr}{$int})
			{
				$snp_count++;
				$snp_position = $int;
			}
		}
		for my $int ($lig_probe_start .. $lig_probe_stop)
		{
			if (exists $chr_snp_positions{$chr}{$int})
			{
				$snp_count++;
				$snp_position = $int;
			}
		}
		if ($snp_count == 0)
		{
			#&reverse_complement_target;is already - strand from the first design round, were just fixing the change in spec for mip arm snps
			&reverse_complement_ext_arms;
			&reverse_complement_lig_arms;
			$mip_count++;
			print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t-\tmissing_snp_now\n";
			print "$line_count/$number_of_features\t$mip_count\n";
		}
		elsif ($snp_count == 1)
		{#make for allele_2
			#make two mips, one on each strand for each allele of the snp single snp found
			&get_alleles_for_position;#$snp_position
			if ($snp_position >= $ext_probe_start && $snp_position <= $ext_probe_stop)
			{
				my $allele_probe_seq_position = $snp_position - $ext_probe_start;#0 based for array
				my @arr_ext_probe_sequence =  split(//, $ext_probe_seq);
				#allele 1
				$mip_count++;
				$arr_ext_probe_sequence[$allele_probe_seq_position] = $allele_1;
				$ext_probe_seq = "";
				foreach my $bases (@arr_ext_probe_sequence)
				{
					$ext_probe_seq = $ext_probe_seq . $bases;
				}
				&reverse_complement_ext_arms;
				&reverse_complement_lig_arms;
				&restriction_site_screen;
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{

					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $extension_arm_length - ($allele_probe_seq_position + 1) + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t-\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}
				#allele 2
				$mip_count++;
				$arr_ext_probe_sequence[$allele_probe_seq_position] = $allele_2;
				$ext_probe_seq = "";
				foreach my $bases (@arr_ext_probe_sequence)
				{
					$ext_probe_seq = $ext_probe_seq . $bases;
				}
				&reverse_complement_ext_arms;
				#&reverse_complement_lig_arms;#lig arm is already the right direction
				&restriction_site_screen;
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{
					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $extension_arm_length - ($allele_probe_seq_position + 1) + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t-\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}
			}
			elsif ($snp_position >= $lig_probe_start && $snp_position <= $lig_probe_stop)						
			{#must be between the ligation target arms
				my $allele_probe_seq_position = $snp_position - $lig_probe_start;#0 based for array
				my @arr_lig_probe_sequence =  split(//, $lig_probe_seq);
				#allele 1
				$mip_count++;
				$arr_lig_probe_sequence[$allele_probe_seq_position] = $allele_1;
				$lig_probe_seq = "";
				foreach my $bases (@arr_lig_probe_sequence)
				{
					$lig_probe_seq = $lig_probe_seq . $bases;
				}
				&reverse_complement_ext_arms;
				&reverse_complement_lig_arms;
				&restriction_site_screen;
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{
					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $ligation_arm_length - ($allele_probe_seq_position + 1) + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.1\t-\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					#print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}

				#allele 2
				$mip_count++;
				$arr_lig_probe_sequence[$allele_probe_seq_position] = $allele_2;
				$lig_probe_seq = "";
				foreach my $bases (@arr_lig_probe_sequence)
				{
					$lig_probe_seq = $lig_probe_seq . $bases;
				}
				#&reverse_complement_ext_arms;arm is already the right direction
				&reverse_complement_lig_arms;
				&restriction_site_screen;
				if ($ntalwl_seq_5_count == 1 && $nbbsrdi_seq_3_count == 1 && $nbbsrdi_seq_5_count == 0 && $ntalwl_seq_3_count == 0)
				{
					$feature_mip_count =~ s/\.\d//g;
					my $allele_probe_position = $ligation_arm_length - ($allele_probe_seq_position + 1) + 1;
					print OUTFILE "$mip_count\t$rank\t$chr\t$ext_probe_start\t$ext_probe_stop\t$ext_probe_seq\t$ext_copy_count\t$lig_probe_start\t$lig_probe_stop\t$lig_probe_seq\t$lig_copy_count\t$scan_start\t$scan_stop\t$scan_seq\t$feature_start_position\t$feature_stop_position\t$feature_mip_count.2\t-\tsnp:$snp_position:$allele_probe_position:$allele_1/$allele_2\n";
					print "$line_count/$number_of_features\t$mip_count\n";
				}
				else
				{
					$mip_count++;
					print OUTFILE "$mip_count\tNA\t$chr\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$scan_stop\t$scan_seq\tNA\tNA\trestriction_site\n";
				}

			}
		}
	}
}

sub restriction_site_screen
{
	$ntalwl_seq_5_count = 0;
	$ntalwl_seq_3_count = 0;
	$nbbsrdi_seq_5_count = 0;
	$nbbsrdi_seq_3_count = 0;
	my $scan_target_sequence_with_probe_arms = $ext_probe_seq . $scan_seq . $lig_probe_seq;
	#print "$scan_target_sequence_with_probe_arms\n";
	#CCTCCCCTGGGCCCCCCACCGGCACCCTCCGCCGCCCCTTCTTGAACACACTCAGGAATGGCTGCTCCAGGTTTTCTCGCTGGTTCTCCAGGTCCAGGATGTCCTAGGAGGAGTAGAGCTCAGGGGAGGGGGCTTTTCCCAGGTTCCTCACA
	my $mip_seq = $universal_left_end_mip_seq . $lig_probe_seq . $universal_middle_mip_seq . $ext_probe_seq . $universal_right_end_mip_seq;
	#print "$mip_seq\n";
	#AGGACCGGATCAACTCTTTTCCCAGGTTCCTCACACTTCAGCTTCCCGATATCCGACGGTAGTGTCCTCCCCTGGGCCCCCCACCCATTGCGTGAACCGA
	#check the mip for restriction sites
	my $length_mip_seq = length($mip_seq) - 1;#0 based
	for my $int (0 .. ($length_mip_seq - 4))
	{
		my $sub_mip_seq_length_5 = substr($mip_seq, $int, 5);
		if ($sub_mip_seq_length_5 =~ m/$ntalwl_seq_5/gi)
		{
			$ntalwl_seq_5_count++;
		}
		if ($sub_mip_seq_length_5 =~ m/$ntalwl_seq_3/gi)
		{
			$ntalwl_seq_3_count++;
		}
	}
	for my $int (0 .. ($length_mip_seq - 5))
	{
		my $sub_mip_seq_length_6 = substr($mip_seq, $int, 6);
		if ($sub_mip_seq_length_6 =~ m/$nbbsrdi_seq_5/gi)
		{
			$nbbsrdi_seq_5_count++;
		}
		if ($sub_mip_seq_length_6 =~ m/$nbbsrdi_seq_3/gi)
		{
			$nbbsrdi_seq_3_count++;
		}
	}

	#check the target with probe arms
	my $length_target_sequence_with_probe_arms = length($scan_target_sequence_with_probe_arms) - 1;#0 based
	for my $int (0 .. ($length_target_sequence_with_probe_arms - 4))
	{
		my $sub_target_seq_length_5 = substr($length_target_sequence_with_probe_arms, $int, 5);
		if ($sub_target_seq_length_5 =~ m/$ntalwl_seq_5/gi)
		{
			$ntalwl_seq_5_count++;
		}
		if ($sub_target_seq_length_5 =~ m/$ntalwl_seq_3/gi)
		{
			$ntalwl_seq_3_count++;
		}
	}
	for my $int (0 .. ($length_target_sequence_with_probe_arms - 5))
	{
		my $sub_target_seq_length_6 = substr($length_target_sequence_with_probe_arms, $int, 6);
		if ($sub_target_seq_length_6 =~ m/$nbbsrdi_seq_5/gi)
		{
			$nbbsrdi_seq_5_count++;
		}
		if ($sub_target_seq_length_6 =~ m/$nbbsrdi_seq_3/gi)
		{
			$nbbsrdi_seq_3_count++;
		}
	}
}
sub reverse_complement_ext_arms
{
	my $flipped_seq;
	my $x;
	while ($ext_probe_seq ne "")
	{
		$x = chop($ext_probe_seq);
		$flipped_seq .= $x;
	}
	### complement sequence ###
	$flipped_seq =~ tr/ATCG/TAGC/;
	$ext_probe_seq = $flipped_seq;

}
sub reverse_complement_lig_arms
{
	my $flipped_seq = "";
	my $x = "";
	while ($lig_probe_seq ne "")
	{
		$x = chop($lig_probe_seq);
		$flipped_seq .= $x;
	}
	### complement sequence ###
	$flipped_seq =~ tr/ATCG/TAGC/;
	$lig_probe_seq = $flipped_seq;
}

sub reverse_complement_target
{
	my $flipped_seq;
	my $x;
	while ($scan_seq ne "")
	{
		$x = chop($scan_seq);
		$flipped_seq .= $x;
	}
	$flipped_seq =~ tr/ATCG/TAGC/;
	$scan_seq = $flipped_seq;
}

sub get_alleles_for_position
{
	my @alleles =  split(//, $chr_snp_positions{$chr}{$snp_position});
	$allele_1 = $alleles[0];
	$allele_2 = $alleles[1];
}
#replaced sub here BJO
sub load_chr_snps
{
	my $snp_load_count = 0;
	#my $snp_file = "/nfs/home/bard/dbSNP/build_129/xml_chr/gt_chr".$chr.".xml_trimmed_to_site_info.rs_orientation_info";
	#my $snp_file = "/nfs/home/bard/dbSNP/build_129/ds_xml_chr/ds_ch".$chr.".xml";
	
	############## COMMENTED OUT BY JOE ON 110613
	# my $snp_file = "/net/grc/vol1/references/human/snps/dbSNP/GATK_rod_hg18_build_130/dbsnp_130_hg18.rod";
	
	open (SNPFILE, $snp_file) or die "Can't open\n $snp_file\n";
	while (<SNPFILE>)
	{
		my $line = $_;
		chomp $line;
		my @line_contents = split(/\s+/, $line);
		my $chr = $line_contents[1];
		$chr =~ s/chr//;
		my $position = $line_contents[2]+1; #BJO
		my $snp_strand_orientation = $line_contents[6];
		my $alleles = $line_contents[9];
		my @alleles_content = split(/\//, $alleles);
		my $a1 = $alleles_content[0];
		my $a2 = $alleles_content[1];
		if ($snp_strand_orientation eq "-")
		{#flip them to be forward
			$a1 =~ tr/ATCG/TAGC/;
			$a2 =~ tr/ATCG/TAGC/;
			$chr_snp_positions{$chr}{$position} = $a1 . $a2;
			#print "$chr\t$position\t$a1\t$a2\tlig->ext\n";
			$snp_load_count++;
			$total_snp_load_count++;
			#print "$total_snp_load_count\n";
		}
		else#just load the forward strand
		{
			$chr_snp_positions{$chr}{$position} = $a1 . $a2;
			#print "$chr\t$position\t$a1\t$a2\text\n";
			$snp_load_count++;
			$total_snp_load_count++;
			#print "$total_snp_load_count\n";
		}
	}
	close SNPFILE;
}

sub parseCommandLine 
{
    my ( $useage ) = "

-mips_file

";

        for (my $i = 0; $i <= $#ARGV; $i++)
        {
                if ($ARGV[$i] =~ /^-/)
                {
                        $arg{$ARGV[$i]} = $ARGV[$i+1];
                }
        }
        die($useage) if (!($arg{-mips_file}));
}
sub get_mips_file
{
	open (FILE, $mips_file);
	@mips_file_contents = <FILE>;
	close FILE;
}
sub get_chr_fasta_sequence
{
	open (FILE, $chr_fasta) or die "Can't open $chr_fasta\n";
	my @contents = <FILE>;
	close FILE;
	$chr_sequence = "";
	foreach my $lines (@contents)
	{
		chomp $lines;
		if ($lines =~ m/>/g)
		{}
		else
		{
			$chr_sequence = $chr_sequence . uc($lines);
		}
	}
}
