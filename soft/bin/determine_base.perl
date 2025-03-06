#!/usr/bin/perl
use strict;

our $file = $ARGV[0];
my $fileFH;
if ( $file =~ /\.gz$/ ) {
  open ($fileFH, '-|', 'zcat', "$file") || die "Cannot Open File using gzip";
}
else {
  open ($fileFH, "$file") || die "Cannot Open File";
}

<$fileFH>;
<$fileFH>;
<$fileFH>;

while (my $line = <$fileFH>) {
  chomp $line;
#  print "\"" . $line . "\"\n";
  my @char_array = split(//,"$line");
  foreach my $char (@char_array){
#    print ord($char) . " ";
    if ( ord($char) < 64 ) { print "33\n"; close $fileFH; exit;}
    if ( ord($char) > 64+34 ) { print "64\n"; close $fileFH; exit;}
  }
#  print "\n";
  <$fileFH>;
  <$fileFH>;
  <$fileFH>;
}
print "Cannot determine base\n";
close($fileFH);
