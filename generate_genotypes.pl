#!/usr/bin/perl -w

use strict;
use warnings;

our $max_n_copies=6;##So far I am considering until 6 copies per allele. We could modify this to only consider the maximum actual number. Important if we would create the substitution model accordingly.
 

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

(scalar(@ARGV) >0 && scalar(@ARGV) < 3)  or die "Usage: script input_directory [format(nexus,xml,human)]";

my $input_dir=$ARGV[0];
my $output_format;
if (scalar @ARGV eq 2)
{
	$output_format=$ARGV[1];
}
else
{
	$output_format="human";
}

chdir($input_dir) or die "The input directory $input_dir is not accesible";

my @files=<*A.txt>;

for my $filea (@files)
{
	my @dataA;
	my @dataB;
	my $fileb=$filea;
	$fileb=~s/(.*)allA.txt/$1allB.txt/;
	(-e $fileb) or die "There is no matching B allele file for $filea";
	my $FILEA;
	my $FILEB;
	open($FILEA,"<",$filea);
	open($FILEB,"<",$fileb);
	my @a_cont= <$FILEA>;
	my @b_cont= <$FILEB>;
	close($FILEA);
	close($FILEB);

	my @samples=(split(" ",$a_cont[0]));
	my @samplesB=(split(" ",$b_cont[0]));
	
	((join("",@samples) eq join("",@samplesB)) && (scalar(@samples)==scalar(@samplesB)) && scalar(@a_cont) eq scalar(@b_cont)) or die "The parsed samples from files $filea and $fileb are not equivalent, check input files!\n";
	splice(@samples,0,5);
	splice(@samplesB,0,5);
	
	my @a_gens;
	my @b_gens;
	
	for(my $i=0;$i<scalar(@a_cont)-1;++$i)
	{
		@a_gens=split(" ",$a_cont[$i+1]);
		@b_gens=split(" ",$b_cont[$i+1]);
		splice(@a_gens,0,5);
		splice(@b_gens,0,5);
		for(my $j=0; $j<scalar(@samples); ++$j)
		{
			$dataA[$j][$i]=$a_gens[$j];
			$dataB[$j][$i]=$b_gens[$j];
		}
	}
	my $OUTFILE;
	my $outname=$filea;
	my $n_samples=scalar(@samples);
	my $n_char=scalar(@a_gens);
	
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
				print $OUTFILE GenToAscii(($dataA[$i][$j]),($dataB[$i][$j]));
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
		for (my $i=0;$i<$n_samples;++$i)
                {
                        print $OUTFILE "\t<taxon id=\"$samples[$i]\">\n";
			#print $OUTFILE "\t\t<date value=\"$date\" direction=\"forwards\" units=\"years\"/>"; DM if eventually I use this to generate the final xml I would need to put here the tip dates.
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
                                print $OUTFILE GenToAscii(($dataA[$i][$j]),($dataB[$i][$j]));
                        }
                        print $OUTFILE "\n\t</sequence>";
                }
                print $OUTFILE "\n</alignment>\n";
		my $text = qq{
	<patterns id="patterns" from="1" strip="false">
		<alignment idref="alignment"/>
	</patterns>

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
			<parameter id="treeModel.rootHeight" lower="30"/>
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
<!--
			<parameter id="clock.rate" value="2.3E-5" lower="0.0" upper="100.0"/>
-->
			<parameter id="clock.rate" value="1"/>
		</rate>
	</strictClockCenancestorBranchRates>

	<frequencyModel id="frequencies">
		<dataType idref="cnv"/>
		<alignment idref="alignment"/>
		<frequencies>
			<parameter id="cnv.frequencies" dimension="28"/>
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
        	<cenancestor>
        		<parameter id="luca" value="30.0" lower="30.0" upper="30.0"/>
        	</cenancestor>
		<strictClockCenancestorBranchRates idref="branchRates"/>
	</cenancestorTreeLikelihood>
	
	<operators id="operators" optimizationSchedule="default">
		<scaleOperator scaleFactor="0.75" weight="3">
                        <parameter idref="cnv.gain"/>
        	</scaleOperator>
        	<scaleOperator scaleFactor="0.75" weight="3">
                        <parameter idref="cnv.loss"/>
        	</scaleOperator>
        	<scaleOperator scaleFactor="0.75" weight="3">
                        <parameter idref="cnv.conversion"/>
        	</scaleOperator>
<!--		
		<scaleOperator scaleFactor="0.75" weight="3">
                        <parameter idref="clock.rate"/>
        	</scaleOperator>
-->
		<subtreeSlide size="0.06" gaussian="true" weight="15">
			<treeModel idref="treeModel"/>
		</subtreeSlide>
		<narrowExchange weight="15">
			<treeModel idref="treeModel"/>
		</narrowExchange>
		<wideExchange weight="3">
			<treeModel idref="treeModel"/>
		</wideExchange>
		<wilsonBalding weight="3">
			<treeModel idref="treeModel"/>
		</wilsonBalding>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="treeModel.rootHeight"/>
		</scaleOperator>
		<uniformOperator weight="30">
			<parameter idref="treeModel.internalNodeHeights"/>
		</uniformOperator>

		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="constant.popSize"/>
		</scaleOperator>

                <upDownOperator scaleFactor="0.75" weight="3">
                        <up>
<!--
                                <parameter idref="clock.rate"/>
-->
                        </up>
                        <down>
                                <parameter idref="treeModel.allInternalNodeHeights"/>
                        </down>
                </upDownOperator>

	</operators>

	<!-- Define MCMC                                                             -->
	<mcmc id="mcmc" chainLength="20000000" autoOptimize="true" operatorAnalysis="$beast_outname.ops">
		<posterior id="posterior">
			<prior id="prior">
                                <coalescentLikelihood idref="coalescent"/>
								<!-- rate priors???-->
                                	<!-- Constrains root height, makes prior negative inf if root height gets beyond 30 -->
                                <uniformPrior lower="0.0" upper="30.0">
                                	<parameter idref="treeModel.rootHeight"/>
                                </uniformPrior>
				<oneOnXPrior>
					<parameter idref="constant.popSize"/>
				</oneOnXPrior>
			</prior>
			<likelihood id="likelihood">
				<cenancestorTreeLikelihood idref="treeLikelihood"/>
			</likelihood>
		</posterior>
		<operators idref="operators"/>

		<!-- write log to screen                                                     -->
		<log id="screenLog" logEvery="1000">
			<column label="Posterior" dp="4" width="12">
				<posterior idref="posterior"/>
			</column>
			<column label="Prior" dp="4" width="12">
				<prior idref="prior"/>
			</column>
			<column label="Likelihood" dp="4" width="12">
				<likelihood idref="likelihood"/>
			</column>
			<column label="gain_rate" sf="6" width="12">
				<parameter idref="cnv.gain"/>
			</column>
			<column label="loss_rate" sf="6" width="12">
				<parameter idref="cnv.loss"/>
			</column>
			<column label="conversion_rate" sf="6" width="12">
				<parameter idref="cnv.conversion"/>
			</column>
<!--			<column label="clock.rate" sf="6" width="12">
				<parameter idref="clock.rate"/>
			</column>
-->
			<column label="rootHeight" sf="6" width="12">
				<parameter idref="treeModel.rootHeight"/>
			</column>
		</log>

		<!-- write log to file                                                       -->
		<log id="fileLog" logEvery="1000" fileName="$beast_outname.log" overwrite="false">
			<posterior idref="posterior"/>
			<prior idref="prior"/>
			<likelihood idref="likelihood"/>
			<parameter idref="cnv.gain"/>
			<parameter idref="cnv.loss"/>
			<parameter idref="cnv.conversion"/>
			<parameter idref="treeModel.rootHeight"/>
			<parameter idref="constant.popSize"/>
			<parameter idref="clock.rate"/>
			<cenancestorTreeLikelihood idref="treeLikelihood"/>
			<coalescentLikelihood idref="coalescent"/>
		</log>

		<!-- write tree log to file                                                  -->
		<logTree id="treeFileLog" logEvery="1000" nexusFormat="true" fileName="$beast_outname.trees" sortTranslationTable="true">
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
			for(my $j=0;$j<scalar(@a_gens)-1;++$j)
			{
				print $OUTFILE "$dataA[$i][$j]$dataB[$i][$j],";
			}
			print $OUTFILE "\n";
		}
	}	
	close($OUTFILE);
}

