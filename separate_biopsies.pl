#!/usr/bin/perl -w
use warnings;
use strict;

(scalar(@ARGV) >0 && scalar(@ARGV) < 2)  or die "Usage: script directory";

my $dir=$ARGV[0];

chdir($dir) or die "The input directory $dir is not accesible";

my @files=<*allA.txt>;
my $pool_query="whole_epi";
for my $filea (@files)
{
        ##File opening and reading
        my @dataA;
        my @dataB;
        my $fileb=$filea;
        my $id=$filea;
        $fileb=~s/(.*)allA.txt/$1allB.txt/;
        $id=~s/(.*)_phased_100_allA.txt/$1/;
        (-e $fileb) or die "There is no matching B allele file for $filea";
        my $FILEA;
        my $FILEB;
        open($FILEA,"<",$filea);
        open($FILEB,"<",$fileb);
        my @a_cont= <$FILEA>;
        my @b_cont= <$FILEB>;
        close($FILEA);
        close($FILEB);

        my @samples=(split("\t",$a_cont[0])); #Header
        my @samplesB=(split("\t",$b_cont[0]));
	chomp(@samples);
	chomp(@samplesB);
        ((join("",@samples) eq join("",@samplesB)) && (scalar(@samples)==scalar(@samplesB)) && (scalar(@a_cont) eq scalar(@b_cont))) or die "The parsed samples from files $filea and $fileb are not equivalent, check input files!\n"; #Checking if the samples are the same, the same number and the same number of loci 
	my @pool_columns=(0,1,2,3,4); #Fixed columns
	my @crypt_columns=(0,1,2,3,4);
	for (my $i=5; $i<scalar @samples; ++$i)
	{
		my $Apool=($samples[$i]=~/$pool_query/);
		my $Bpool=($samplesB[$i]=~/$pool_query/);
		$Apool != $Bpool and die "Input files for A and B are not compatible. Check the input files";
		if($Apool)
		{
			push(@pool_columns,$i);
		}
		else
		{
			push(@crypt_columns,$i);
		}
	}
	
	open(my $OUTAPOOL,">${id}_phased_100poolA.txt");
	open(my $OUTBPOOL,">${id}_phased_100poolB.txt");
	open(my $OUTACRYPT,">${id}_phased_100cryptsA.txt");
	open(my $OUTBCRYPT,">${id}_phased_100cryptsB.txt");
	print($OUTAPOOL join("\t",@samples[@pool_columns]),"\n");
	print($OUTBPOOL join("\t",@samplesB[@pool_columns]),"\n");
	print($OUTACRYPT join("\t",@samples[@crypt_columns]),"\n");
        print($OUTBCRYPT join("\t",@samplesB[@crypt_columns]),"\n");

	for (my $i=1;$i<scalar @a_cont;++$i)
	{
		@dataA=split("\t",$a_cont[$i]);
		chomp(@dataA);
		@dataB=split("\t",$b_cont[$i]);
		chomp(@dataB);
		print($OUTAPOOL join("\t",@dataA[@pool_columns]),"\n");
		print($OUTBPOOL join("\t",@dataB[@pool_columns]),"\n");
		print($OUTACRYPT join("\t",@dataA[@crypt_columns]),"\n");
                print($OUTBCRYPT join("\t",@dataB[@crypt_columns]),"\n");
	}

	close($OUTAPOOL);
	close($OUTBPOOL);
	close($OUTACRYPT);
	close($OUTBCRYPT);

}
		
