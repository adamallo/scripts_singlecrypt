#!/usr/bin/perl -w

use strict;
use warnings;
use Time::Piece;
use Time::Seconds;
use List::Util qw(sum);
use POSIX qw(floor ceil);
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

##Configuration variables
#########################
our $max_n_copies=6;##So far I am considering until 6 copies per allele. We could modify this to only consider the maximum actual number. Important if we would create the substitution model accordingly.
our $period_log=10000;
our $n_gen=100000000;

##ARGV variables
################
our $output_format="human";
our $input_dir="";
our $sample_prior;
our $seed=20;
our $fixed;
my $usage="\nUsage: $0 [options] -i input_directory -f format<nexus,xml,human> -s seed \n\nOptions:\n--------\n-p/--sample_prior: xml will make BEAST sample from the pior only (no data)\n\n";
my $help;

######################################################
##MAIN
######################################################

##Getopt
######################################################
(! GetOptions(
        'input_dir|i=s' => \$input_dir,
        'format|f=s' => \$output_format,
        'sample_prior|p' => \$sample_prior,
	'seed|s=i' => \$seed,
	'fixed_cnv' => \$fixed,
        'help|h' => \$help,
                )) or ((($input_dir eq "") || (! -d $input_dir)) || $help) and die $usage;

chdir($input_dir) or die "The input directory $input_dir is not accesible";

srand $seed;

my @files=<*A.txt>;
my @timestamp_files=<*_data_timestamps.txt>;

##Hash with patient: DOB timestamp
my $DOB_file="dobs.txt";
open(my $FILEDOB,$DOB_file);
my @line_dobs=<$FILEDOB>;
close($FILEDOB);
my %dobs;
my @temp;

##Stats file
my $stats_file="stats.csv";
open(my $STATS,">$stats_file");
print($STATS "file,n_leaves,delta_time,n_sites\n");

##DOBS parsing
foreach my $dobline (@line_dobs)
{
	@temp=split(" ",$dobline);
	$dobs{$temp[0]}=$temp[1];
	#print("Key: $temp[0], timestamp:$temp[1]\n"); #
}

##Main loop
for my $filea (@files)
{
	##File opening and reading
	my @dataA;
	my @dataB;
	my $fileb=$filea;
	my $fileTimestamps=$filea;
	my $id=$filea;
	$fileb=~s/(.*)A.txt/$1B.txt/;
	$id=~s/(.*)_phased.*.txt/$1/;
	my $dob=$dobs{$id};
	if (!defined($dob))
	{
		$dob=0;
		print("WARNING:There is no date of birth for the patient $id\n");
	}
	$fileTimestamps=~s/(.*)_phased.*.txt/$1_data_timestamps.txt/;
	(-e $fileb) or die "There is no matching B allele file for $filea";
	(-e $fileTimestamps) or die "There is no matching timestamps file $fileTimestamps for $filea";
	my $FILEA;
	my $FILEB;
	my $FILETIMES;
	open($FILEA,"<",$filea);
	open($FILEB,"<",$fileb);
	open($FILETIMES,$fileTimestamps);
	my @a_cont= <$FILEA>;
	my @b_cont= <$FILEB>;
	my @time_lines= <$FILETIMES>;
	close($FILEA);
	close($FILEB);
	close($FILETIMES);

	my @samples=split("\t",$a_cont[0]); #Header
	my @samplesB=split("\t",$b_cont[0]);
	chomp(@samples);
	chomp(@samplesB);
	#print("DEBUG: ",join(",",@samples));

	((join("",@samples) eq join("",@samplesB)) && (scalar(@samples)==scalar(@samplesB)) && scalar(@a_cont) eq scalar(@b_cont)) or die "The parsed samples from files $filea and $fileb are not equivalent, check input files!\n"; #Checking if the samples are the same, the same number and the same number of loci

	##Parsing the alleles into dataA and dataB arrays	
	splice(@samples,0,5);# Keeping only sample names
	splice(@samplesB,0,5);
	my @a_gens; #Genotypes
	my @b_gens;
	my $temp_a;
	my $temp_b;
	print("Working in $id, fileA $filea file B $fileb, n_loci= ",scalar(@a_cont)-1 ,", ",scalar(@b_cont)-1,"\n");
	
	for(my $i=0;$i<(scalar(@a_cont)-1);++$i) ###Locus. From 0 to size-2 instead of 0 to size-1 due to the header (using i+1 to read)
	{
		@a_gens=split("\t",$a_cont[$i+1]);
		@b_gens=split("\t",$b_cont[$i+1]);
		chomp(@a_gens);
		chomp(@b_gens);
		splice(@a_gens,0,5);#Keeping only the genotypes
		splice(@b_gens,0,5);

		#print("DEBUG: line $i ",scalar @a_gens," ",scalar @b_gens," ",scalar @samples," ", scalar @samplesB,"\n");
		for(my $j=0; $j<scalar(@samples); ++$j) ###Sample
		{
			$dataA[$j][$i]=$a_gens[$j];
			$dataB[$j][$i]=$b_gens[$j];
		}
	}

	correctBaseline(\@dataA,\@dataB);
	limitToModel(\@dataA,\@dataB);
	eliminateInvariableWildTypes(\@dataA,\@dataB);

	##Timestamp hash
	my %times;
	foreach my $timeline (@time_lines)
	{
        	@temp=split(" ",$timeline);
        	$times{$temp[0]}=$temp[1];
        	#print("Key: $temp[0], timestamp: $temp[1]\n"); #DEBUG
	}

	my $OUTFILE;
	my $outname=$filea;
	my $n_samples=scalar(@samples);
	my $n_char=scalar(@{$dataA[0]});
	my $max_date;
	my $min_date;

	##STATS output
	for (my $i=0;$i<$n_samples;++$i)
        {
               	$times{$samples[$i]} or die "There is not timestamp information for the sample $samples[$i] in the file $fileTimestamps\n";
		my $timestamp=$times{$samples[$i]}-$dob;
		my $t_object=Time::Seconds->new($timestamp); ##Times in years since DOB
		my $time=$t_object->years;
		unless ($max_date)
		{
			$max_date=$time;
		}
		unless ($min_date)
		{
			$min_date=$time;
		}
		if ($time > $max_date)
		{
			$max_date=$time;
		}
		if ($time < $min_date)
		{
			$min_date=$time;
		}
        }
	my $fileid=$outname;
	$fileid=~s/A\.txt//;
	my $deltaTime=$max_date-$min_date;
	print($STATS "$fileid,$n_samples,$deltaTime,$n_char\n");
	##

	#Data output depending on the selected format
	#############################################
	if ($output_format=~/nexus/i)
	{
		$outname=~s/A\.txt/\.nex/;
		open($OUTFILE,">$outname");
		print $OUTFILE "#NEXUS\n";
		print $OUTFILE "begin taxa;\n\tdimensions ntax=$n_samples;\n\ttaxlabels\n";
		for (my $i=0;$i<$n_samples;++$i)
		{
			print $OUTFILE "\t\t$samples[$i]\n";
		}
		print $OUTFILE ";\nend;\nbegin characters;\n\tdimensions nchar=$n_char;\n\tformat symbols=\"";
		
		for (my $a=0;$a<=$max_n_copies;++$a)
		{
			for (my $b=0;$b+$a<=$max_n_copies;++$b)
			{
				print $OUTFILE GenToAscii($a,$b);
			}
		}
		print $OUTFILE "\" missing=? gap=-;\nmatrix\n";
		for (my $i=0;$i<$n_samples;++$i)
                {
                        print $OUTFILE "\t\t$samples[$i]\t";
			for(my $j=0;$j<$n_char;++$j)
			{
				print $OUTFILE $sample_prior?"?":GenToAscii(($dataA[$i][$j]),($dataB[$i][$j]));
			}
			print $OUTFILE "\n";
                }
		print $OUTFILE ";\nend;";
	}
	elsif ($output_format=~/xml/i)
	{
		$outname=~s/A\.txt/\.xml/;
		my $beast_outname=$outname;
		$beast_outname=~s/\.xml//;
		open($OUTFILE,">$outname");
		print $OUTFILE "<?xml version=\"1.0\" standalone=\"yes\"?>\n<beast>\n<taxa id=\"taxa\">\n";
		my $timestamp;
		for (my $i=0;$i<$n_samples;++$i)
                {
                        print $OUTFILE "\t<taxon id=\"$samples[$i]\">\n";
			#print "DEBUG: taxon $samples[$i], timestamp $times{$samples[$i]}\n"; #DEBUG
			$times{$samples[$i]} or die "There is not timestamp information for the sample $samples[$i] in the file $fileTimestamps\n";
			$timestamp=$times{$samples[$i]}-$dob;
			my $t_object=Time::Seconds->new($timestamp); ##Times in years since DOB
			my $time=$t_object->years;
			unless ($max_date)
			{
				$max_date=$time;
			}
			unless ($min_date)
			{
				$min_date=$time;
			}
			if ($time > $max_date)
			{
				$max_date=$time;
			}
			if ($time < $min_date)
			{
				$min_date=$time;
			}
			#print "DEBUG: sample timestamp $times{$samples[$i]}, DOB $dob, epochtime $timestamp, epochtime(years) $time\n"; #Debug
			print $OUTFILE "\t\t<date value=\"$time\" direction=\"forwards\" units=\"years\"/>";
			print $OUTFILE "\t</taxon>\n";
                }
		print $OUTFILE "</taxa>\n\n<generalDataType id=\"cnv\">";
		my @states;

		my $state=0;
		for (my $a=0;$a<=$max_n_copies;++$a)
                {
                        for (my $b=0;$b+$a<=$max_n_copies;++$b)
                        {
				my $code=GenToAscii($a,$b);
				push(@states,$code);
                                print $OUTFILE "\n\t<state code=\"$code\"/> <!-- Genotype: $a,$b ; Beast State: $state -->";
				$state+=1;
                        }
		
                }
		
		print $OUTFILE "\n\t<ambiguity code=\"-\" states=\"",join("",@states),"\"/>";
		print $OUTFILE "\n\t<ambiguity code=\"?\" states=\"",join("",@states),"\"/>";
		print $OUTFILE "\n</generalDataType>";
		print $OUTFILE "\n\n<alignment id=\"alignment\">\n\t<dataType idref=\"cnv\"/>";

		for (my $i=0;$i<$n_samples;++$i)
                {
                        print $OUTFILE "\n\t<sequence>\n\t\t<taxon idref=\"$samples[$i]\"/>\n\t\t";
                        for(my $j=0;$j<$n_char;++$j)
                        {
                                print $OUTFILE $sample_prior?"?":GenToAscii(($dataA[$i][$j]),($dataB[$i][$j]));
                        }
                        print $OUTFILE "\n\t</sequence>";
                }
                print $OUTFILE "\n</alignment>\n";
		my $diff_date=$max_date-$min_date;
		my $cnv_operators=qq{<scaleOperator scaleFactor="0.25" weight="0.5">
			<parameter idref="cnv.loss"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.25" weight="0.5">
			<parameter idref="cnv.conversion"/>
		</scaleOperator>};
		my $cnv_priors=qq{<!-- Loss and conversion (relative to gain) rate priors-->
				<exponentialPrior mean="1.0" offset="0.0">
					<parameter idref="cnv.loss"/>
				</exponentialPrior>
				<exponentialPrior mean="1.0" offset="0.0">
					<parameter idref="cnv.conversion"/>
				</exponentialPrior>};
		if ($fixed)
		{
			$cnv_operators="";
			$cnv_priors="";
		}
		my $text = qq{
	
	<ascertainedCharacterPatterns id="patterns">
	        <alignment idref="alignment"/>
	        <state code='H'/>
	</ascertainedCharacterPatterns>	
	
	<!-- A prior assumption that the population size has remained constant       -->
	<!-- throughout the time spanned by the genealogy.                           -->
	<constantSize id="constant" units="years">
		<populationSize>
			<parameter id="constant.popSize" value="1" lower="0.0"/>
		</populationSize>
	</constantSize>

	<!-- Generate a random starting tree under the coalescent process            -->
	<coalescentSimulator id="startingTree">
		<taxa idref="taxa"/>
		<constantSize idref="constant"/>
	</coalescentSimulator>

	<!-- Generate a tree model                                  -->
	<treeModel id="treeModel">
		<coalescentTree idref="startingTree"/>
		<rootHeight>
			<parameter id="treeModel.rootHeight"/>
		</rootHeight>
		<nodeHeights internalNodes="true">
			<parameter id="treeModel.internalNodeHeights"/>
		</nodeHeights>
		<nodeHeights internalNodes="true" rootNode="true">
			<parameter id="treeModel.allInternalNodeHeights"/>
		</nodeHeights>
	</treeModel>

	<!-- Generate a coalescent likelihood                                        -->
	<coalescentLikelihood id="coalescent">
		<model>
			<constantSize idref="constant"/>
		</model>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</coalescentLikelihood>

	<!-- The strict clock (Uniform rates across branches)                        -->
 
	<strictClockCenancestorBranchRates id="branchRates">
		<rate>
			<parameter id="clock.rate" value="1"/>
		</rate>
	</strictClockCenancestorBranchRates>

	<frequencyModel id="frequencies">
		<dataType idref="cnv"/>
		<frequencies>
			<parameter id="cnv.frequencies" value="0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"/>
		</frequencies>
	</frequencyModel>
	
	<CNVModel id="cnv_subsmodel">
		<frequencies>
			<frequencyModel idref="frequencies"/>
        	</frequencies>
        	<gain_rate>
			<parameter id="cnv.gain" value="1" lower="0"/>
        	</gain_rate>
        	<loss_rate>
			<parameter id="cnv.loss" value="1" lower="0"/>
        	</loss_rate>
        	<conversion_rate>
			<parameter id="cnv.conversion" value="1" lower="0"/>
        	</conversion_rate>
	</CNVModel>

	<siteModel id="siteModel">
		<substitutionModel>
			<CNVModel idref="cnv_subsmodel"/>
		</substitutionModel>
	</siteModel>

	<cenancestorTreeLikelihood id="treeLikelihood" useAmbiguities="false">
		<patterns idref="patterns"/>
		<treeModel idref="treeModel"/>
		<siteModel idref="siteModel"/>
        	<cenancestorHeight>
        		<parameter id="luca_height" lower="$diff_date" upper="$max_date"/>
        	</cenancestorHeight>
		<cenancestorBranch>
			<parameter id="luca_branch" value="1" upper="$min_date" lower="0.0"/>
			<!-- Value 1 as a safe starting value -->
		</cenancestorBranch>
		<strictClockCenancestorBranchRates idref="branchRates"/>
	</cenancestorTreeLikelihood>
	
	<operators id="operators" optimizationSchedule="default">
		<scaleOperator scaleFactor="0.5" weight="10.0">
                        <parameter idref="clock.rate"/>
        	</scaleOperator>
		$cnv_operators
		<subtreeSlide size="2.5" gaussian="true" weight="15.0"> <!-- 2.5 years. They will be automatically optimized by BEAST though -->
			<treeModel idref="treeModel"/>
		</subtreeSlide>
		<narrowExchange weight="15.0">
			<treeModel idref="treeModel"/>
		</narrowExchange>
		<wideExchange weight="3.0">
			<treeModel idref="treeModel"/>
		</wideExchange>
		<wilsonBalding weight="3.0">
			<treeModel idref="treeModel"/>
		</wilsonBalding>
		<scaleOperator scaleFactor="0.75" weight="5.0">
			<parameter idref="treeModel.rootHeight"/>
		</scaleOperator>
		<uniformOperator weight="30.0">
			<parameter idref="treeModel.internalNodeHeights"/>
		</uniformOperator>
		
		<scaleOperator scaleFactor="0.2" weight="1.0"> <!-- We operate the branch since it is relative to the root. Operating luca_height is error prone, since it depends on the root -->
                        <parameter idref="luca_branch"/>
                </scaleOperator>

		<scaleOperator scaleFactor="0.5" weight="3.0">
			<parameter idref="constant.popSize"/>
		</scaleOperator>

                <upDownOperator scaleFactor="0.75" weight="5.0">
                        <up>
                                <parameter idref="clock.rate"/>
                        </up>
                        <down>
                                <parameter idref="treeModel.allInternalNodeHeights"/>
                        </down>
                </upDownOperator>

	</operators>

	<!-- Define MCMC                                                             -->
	<mcmc id="mcmc" chainLength="$n_gen" autoOptimize="true" operatorAnalysis="$beast_outname.ops">
		<posterior id="posterior">
			<prior id="prior">
                                <coalescentLikelihood idref="coalescent"/>
				<oneOnXPrior>
					<parameter idref="constant.popSize"/>
				</oneOnXPrior>
	                        
				<!-- Clock (gain) Rate Prior. -->
	                        <logNormalPrior mean="-4.0" stdev="2.5" offset="0.0" meanInRealSpace="false">
	                                <parameter idref="clock.rate"/>
	                        </logNormalPrior>
				$cnv_priors
				                                
				<!-- Cenancestor Prior on the height, since it is easier to have a meaningfull prior on it (time of the initial development of the BE fragment) -->
                                <uniformPrior lower="$diff_date" upper="$max_date">
                                	<parameter idref="luca_height"/>
                                </uniformPrior>
			</prior>
			<likelihood id="likelihood">
				<cenancestorTreeLikelihood idref="treeLikelihood"/>
			</likelihood>
		</posterior>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="$period_log">
			<column label="Posterior" dp="4" width="12">
				<posterior idref="posterior"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			<column label="rel_loss_rate" sf="6" width="12">
				<parameter idref="cnv.loss"/>
			</column>
			<column label="rel_conv_rate" sf="6" width="12">
				<parameter idref="cnv.conversion"/>
			</column>
			<column label="gain_rate" sf="6" width="12">
				<parameter idref="clock.rate"/>
			</column>

			<column label="rootHeight" sf="6" width="12">
				<parameter idref="treeModel.rootHeight"/>
			</column>
			
			<column label="luca_height" sf="6" width="12">
				<parameter idref="luca_height"/>
			</column>
			
			<column label="luca_branch" sf="6" width="12">
				<parameter idref="luca_branch"/>
			</column>
		</log>

		<!-- write log to file                                                       -->
		<log id="fileLog" logEvery="$period_log" fileName="$beast_outname.log" overwrite="false">
			<posterior idref="posterior"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="cnv.loss"/>
			<parameter idref="cnv.conversion"/>
			<parameter idref="treeModel.rootHeight"/>
			<parameter idref="luca_height"/>
			<parameter idref="luca_branch"/>
			<parameter idref="constant.popSize"/>
			<parameter idref="clock.rate"/>
			<coalescentLikelihood idref="coalescent"/>
		</log>

		<!-- write tree log to file                                                  -->
		<logTree id="treeFileLog" logEvery="$period_log" nexusFormat="true" fileName="$beast_outname.trees" sortTranslationTable="true">
			<treeModel idref="treeModel"/>
			<trait name="rate" tag="rate">
				<strictClockCenancestorBranchRates idref="branchRates"/>
			</trait>
			<posterior idref="posterior"/>
		</logTree>
	</mcmc>
	<report>
		<property name="timer">
			<mcmc idref="mcmc"/>
		</property>
	</report>
</beast>
};
		print $OUTFILE "$text\n";
	}
	else
	{
		$outname=~s/A\.txt/\.txt/;
		open($OUTFILE,">$outname");
		for(my $i=0;$i<scalar(@samples);++$i)
		{
			print $OUTFILE "$samples[$i] ";
			for(my $j=0;$j<$n_char;++$j)
			{
				print $OUTFILE $sample_prior?"??":"$dataA[$i][$j]$dataB[$i][$j],";
			}
			print $OUTFILE "\n";
		}
	}	
	close($OUTFILE);
}

close($STATS);
exit;

##FUNCTIONS
################

sub GenToAscii
{
	my ($A,$B)=@_;
	my $tot_a=0;
	for (my $i=1;$i<=$A;++$i)
	{
		$tot_a+=$max_n_copies-$i+2;
	}
	return chr($tot_a+$B+64);#Ascii
}

sub correctBaseline
{
	my ($A,$B)=@_;
	my $nLoci=scalar @{${$A}[0]};
	my $nSamples=scalar @{$A};
#	print("DEBUG: $nLoci, $nSamples\n");
	for (my $sample=0;$sample<$nSamples;++$sample)
	{
		my @temp=(@{${$A}[$sample]},@{${$B}[$sample]});
		my $estBaseline=sum(@temp)/@temp;
		my $intBaseline=int($estBaseline + 0.5); #Round 
#		print("DEBUG: $estBaseline, $intBaseline\n");
		if ($intBaseline > 1)
		{
			print("Sample $sample has estimated integer baseline $intBaseline\n");
		}
		for (my $locus=0;$locus<$nLoci;++$locus)
		{
#			print("DEBUG: sample $sample, locus $locus, ${${$A}[$sample]}[$locus], $intBaseline\n");
			if (${${$A}[$sample]}[$locus] % $intBaseline == 0) ##
			{
				${${$A}[$sample]}[$locus]/=$intBaseline;
			}
			elsif(rand() <0.5)
			{
				${${$A}[$sample]}[$locus]=floor(${${$A}[$sample]}[$locus]/$intBaseline);
			}
			else
			{
				${${$A}[$sample]}[$locus]=ceil(${${$A}[$sample]}[$locus]/$intBaseline);
			}
			if (${${$B}[$sample]}[$locus] % $intBaseline == 0) ##
                        {
                                ${${$B}[$sample]}[$locus]/=$intBaseline;
                        }
                        elsif(rand() <0.5)
                        {
                                ${${$B}[$sample]}[$locus]=floor(${${$B}[$sample]}[$locus]/$intBaseline);
                        }
                        else
                        {
                                ${${$B}[$sample]}[$locus]=ceil(${${$B}[$sample]}[$locus]/$intBaseline);
                        }

		}
	}
}

sub limitToModel
{
	my($A,$B)=@_;
	my $nLoci=scalar @{${$A}[0]};
        my $nSamples=scalar @{$A};
	my $temp_a=1;
	my $temp_b=1;
	for (my $sample=0;$sample<$nSamples;++$sample)
        {
		for (my $locus=0;$locus<$nLoci;++$locus)
		{
		        if(${${$A}[$sample]}[$locus]+${${$B}[$sample]}[$locus] > $max_n_copies) ##We need to truncate the values
	                {
	                        if(${${$A}[$sample]}[$locus] == ${${$B}[$sample]}[$locus]) ##Both are equal, we set the maximum value
	                        {
	                                $temp_a=$max_n_copies/2;
	                                $temp_b=$max_n_copies/2;
	                        }
	                        else #Not really needed. Here we truncate alleles if they are higher than the maximum by theirselves
	                        {
	                                $temp_a=(${${$A}[$sample]}[$locus] > $max_n_copies ? $max_n_copies : ${${$A}[$sample]}[$locus]);
	                                $temp_b=(${${$B}[$sample]}[$locus] > $max_n_copies ? $max_n_copies : ${${$B}[$sample]}[$locus]);
	                        }
	                        while($temp_a+$temp_b>$max_n_copies) ##Truncate the one that has more. If they have the same, the one that had originally less
	                        {
	                                if($temp_a>$temp_b)
	                                {
	                                        $temp_a-=1;
	                                }
	                                elsif($temp_b>$temp_a)
	                                {
	                                        $temp_b-=1;
	                                }
	                                else
	                                {
	                                        if(${${$A}[$sample]}[$locus]<${${$B}[$sample]}[$locus])
	                                        {
	                                                $temp_a-=1;
	                                        }
	                                        else
	                                        {
	                                                $temp_b-=1;
	                                        }
	                                }
	                        }
	                        print("WARNING: Segment with states ${${$A}[$sample]}[$locus],${${$B}[$sample]}[$locus], while the maximum state is $max_n_copies. The output state will be $temp_a,$temp_b\n");
	                        ${${$A}[$sample]}[$locus]=$temp_a;
	                        ${${$B}[$sample]}[$locus]=$temp_b;
	                }
		}
        }

}

sub eliminateInvariableWildTypes
{
	my($A,$B)=@_;
        my $nLoci=scalar @{${$A}[0]};
        my $nSamples=scalar @{$A};
	my @tokeep;
	for (my $locus=0;$locus<$nLoci;++$locus)
	{
		my $variable=0;
#		print("DEBUG: locus $locus: ");
        	for (my $sample=0;$sample<$nSamples;++$sample)
        	{
#			print("${${$A}[$sample]}[$locus],${${$B}[$sample]}[$locus] ");
			if(${${$A}[$sample]}[$locus] != 1 || ${${$B}[$sample]}[$locus] != 1)
			{
				$variable=1;
				last;
			}
			
		}
		if ($variable == 1)
		{
#			print(" --> variable\n");
			push(@tokeep,$locus);
		}
#		else
#		{
#			print("\n");
#		}
	}
	for (my $sample=0;$sample<$nSamples;++$sample)
	{
#		print("DEBUG: tokeep @tokeep\n");
		my @tempA=@{${$A}[$sample]}[@tokeep];
		${$A}[$sample]=\@tempA;
		my @tempB=@{${$B}[$sample]}[@tokeep];
		${$B}[$sample]=\@tempB;
	}
}


