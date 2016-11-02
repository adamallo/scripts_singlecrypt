#!/bin/bash
SCRIPTS_DIR=/Users/Diego/Desktop/singlecrypt/scripts/subsmodel
mkdir outlogs
mkdir trees
cd input_files 
../bin/CoalEvol7.3.5
perl -p0e 's/D.*?= //sg' Results/trees | perl -p0e 's/\n+/\n/sg' > ../trees/trees.tree ##Eliminates extra information printed by CoalEvol
mv Results ../outlogs
cd ..
mkdir xml
python $SCRIPTS_DIR/../add_cenan_years.py -i trees/trees.tree -o trees/rooted_scaled.trees -gt 365 -od 20 ##Adds the outgroup and scales the tree in time units
mkdir alignments
rates=(0.0001 0.001 0.01);gds=(0 0.001 "random");mods=("none" "baseline" "max2");i=0;seed=1;while read tree; do echo $tree > trees/tree${i}.tree;for gd in ${gds[*]}; do for rate in ${rates[*]}; do for mod in ${mods[*]}; do echo "Simulating tree${i}_r${rate}_gd${gd}_mod${mod}.xml"; if [[ $gd == "random" ]]; then python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -pgd 0 --randomGD 2 -ggd 0 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed --period 2000 --ngen 20000000; else python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -pgd $gd -ggd 0.31 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed --period 2000 --ngen 20000000;fi;done;seed=$(($seed+1));done;done;i=$(($i+1));done < trees/rooted_scaled.trees > outlogs/log_simulation.txt

#Estimation
for i in *; do sbatch beast.sbatch -seed 15 $i > $i.log_stdout;done

#Tree analysis
for i in *.trees; do name=$(echo $i | sed "s/\.xml.trees*//"); sbatch -p private --mem 10000 <(echo -e '#!/bin/bash'"\ntreeannotator -burnin 1000000 $i ${name}_MCC.trees");done

#Parameter analysis
echo -e "rep r gd mod cnv.loss cnv.conversion treeModel.rootHeight luca_height luca_branch constant.popSize clock.rate" > parameters.out;for i in tree*.xml.log;do rep=$(echo $i | sed "s/tree\(.*\)_r.*_gd.*_mod.*\.xml\.log/\1/");r=$(echo $i | sed "s/tree.*_r\(.*\)_gd.*_mod.*\.xml\.log/\1/");gd=$(echo $i | sed "s/tree.*_r.*_gd\(.*\)_mod.*\.xml\.log/\1/");mod=$(echo $i | sed "s/tree.*_r.*_gd.*_mod\(.*\)\.xml\.log/\1/");awk -v burning=500 -v rep=$rep -v r=$r -v gd=$gd -v mod=$mod 'BEGIN{FS="\t";a=b=c=d=e=f=g=0}{if(NR>3+burning)a+=$5;b+=$6;c+=$7;d+=$8;e+=$9;f+=$10;g+=$11}END{print rep,r,gd,mod,a/(NR-3),b/(NR-3),c/(NR-3),d/(NR-3),e/(NR-3),f/(NR-3),g/(NR-3)}' $i >> parameters.out;done

