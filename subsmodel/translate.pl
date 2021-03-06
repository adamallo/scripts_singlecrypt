use strict;
use warnings;

my ($infile,$dictfile,$out_file,$i_column,$o_column)=@ARGV;
scalar @ARGV != 5 and die "Usage script input_file dict_file out_file in_column out_column\n";
open (my $DICTFILE, $dictfile) or die "Incorrect input dictionary file $dictfile\n";
my @dict=<$DICTFILE>;
close $DICTFILE;

open(my $INFILE, $infile) or die "Incorrect input file $infile\n";
my $back=$\;
$/="";
my $content=<$INFILE>;
$/=$back;
close $INFILE;

shift @dict;
my %final_dict;
for (my $i=0; $i<scalar(@dict)-1; ++$i)
{
	my @columns=split(" ",$dict[$i]);
	my $key=$columns[$i_column-1];
	my $subc=$columns[$o_column-1];
	print("$key,$subc,$i\n");
	$content=~s/\b$key\b/$i/gme;
	$final_dict{$i}=$subc;
}
foreach my $key (keys %final_dict)
{
#	print "Debug: $key, $final_dict{$key}\n";
#	print "$content\n";
	$content=~s/\b$key\b/$final_dict{$key}/gme;
#	print "$content\n";

}

#foreach my $subc (sort {$final_dict{$a} <=> $final_dict{$b}} keys %final_dict)
#{
#		$content=~s/\b$final_dict{$subc}\b/$subc/gme;	
#}
open(my $OUTFILE,">$out_file") or die "Incorrect output file $out_file\n";
print $OUTFILE $content;
close $OUTFILE;
exit;
