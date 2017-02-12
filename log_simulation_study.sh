#!/bin/bash

###ATTENTION: This script has been put together as a log and it will not work as expected if run at once (it does not wait for cluster jobs to finish, for example)

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
rates=(0.0003 0.003 0.03);gds=(0 0.001 "random");mods=("none" "baseline" "max2");i=0;seed=1;while read tree; do echo $tree > trees/tree${i}.tree;for gd in ${gds[*]}; do for rate in ${rates[*]}; do for mod in ${mods[*]}; do echo "Simulating tree${i}_r${rate}_gd${gd}_mod${mod}.xml"; if [[ $gd == "random" ]]; then python $SCRIPTS_DIR/simulate.py -r $rate -c 1 -g 1 -d 1 -pgd 0 --randomGD 2 -ggd 0 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed --period 2000 --ngen 20000000; else python $SCRIPTS_DIR/simulate.py -r $rate -c 1 -g 1 -d 1 -pgd $gd -ggd 0.31 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed --period 2000 --ngen 20000000;fi;done;seed=$(($seed+1));done;done;i=$(($i+1));done < trees/rooted_scaled.trees > outlogs/log_simulation.txt
#rates=(0.0001 0.001 0.01);gds=(0 0.001 "random");mods=("none" "baseline" "max2");i=0;seed=1;while read tree; do echo $tree > trees/tree${i}.tree;for gd in ${gds[*]}; do for rate in ${rates[*]}; do for mod in ${mods[*]}; do echo "Simulating tree${i}_r${rate}_gd${gd}_mod${mod}.xml"; if [[ $gd == "random" ]]; then python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -pgd 0 --randomGD 2 -ggd 0 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed --period 2000 --ngen 20000000; else python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -pgd $gd -ggd 0.31 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed --period 2000 --ngen 20000000;fi;done;seed=$(($seed+1));done;done;i=$(($i+1));done < trees/rooted_scaled.trees > outlogs/log_simulation.txt

##
#xmls are used to be run three times in folders r1-r3
#####################################################

#Estimation
###########
cd r1
for i in *; do sbatch beast.sbatch -seed 25 $i;done
cd ../r2
for i in *; do sbatch beast.sbatch -seed 26 $i;done
cd ../r3
for i in *; do sbatch beast.sbatch -seed 27 $i;done
cd..

#Double check
cd r1
for i in *.xml; do name=$(echo $i | sed "s/.xml//"); nsamples=$(wc -l ${name}.xml.log | awk '{print $1}'); if [[ $nsamples -ne 10004 ]];then sbatch beast.sbatch -seed 25 -overwrite $i ;fi;done
cd ../r2
for i in *.xml; do name=$(echo $i | sed "s/.xml//"); nsamples=$(wc -l ${name}.xml.log | awk '{print $1}'); if [[ $nsamples -ne 10004 ]];then sbatch beast.sbatch -seed 26 -overwrite $i ;fi;done
cd ../r3
for i in *.xml; do name=$(echo $i | sed "s/.xml//"); nsamples=$(wc -l ${name}.xml.log | awk '{print $1}'); if [[ $nsamples -ne 10004 ]];then sbatch beast.sbatch -seed 27 -overwrite $i ;fi;done

#Tree analysis per replicate
############################
#for j in r*
#do
#    cd $j
#        for i in *.trees; do name=$(echo $i | sed "s/\.xml.trees*//"); sbatch -p private --mem 10000 <(echo -e '#!/bin/bash'"\ntreeannotator -burninTrees 500 $i ${name}_MCC.trees");done
#    cd ..
#done
#Double check
#for j in r*
#do 
#    cd $j
#    for i in *.xml; do name=$(echo $i | sed "s/.xml//"); if [[ ! -f ${name}_MCC.trees ]];then sbatch -p private --mem 10000 <(echo -e '#!/bin/bash'"\ntreeannotator -burninTrees 500 ${i}.trees ${name}_MCC.trees") ;fi;done
#    cd ..
#done


#Old parameter analysis by replicate
####################################
echo -e "rep r gd mod cnv.loss cnv.conversion treeModel.rootHeight luca_height luca_branch constant.popSize clock.rate" > parameters.out;for i in tree*.xml.log;do rep=$(echo $i | sed "s/tree\(.*\)_r.*_gd.*_mod.*\.xml\.log/\1/");r=$(echo $i | sed "s/tree.*_r\(.*\)_gd.*_mod.*\.xml\.log/\1/");gd=$(echo $i | sed "s/tree.*_r.*_gd\(.*\)_mod.*\.xml\.log/\1/");mod=$(echo $i | sed "s/tree.*_r.*_gd.*_mod\(.*\)\.xml\.log/\1/");awk -v burning=500 -v rep=$rep -v r=$r -v gd=$gd -v mod=$mod 'BEGIN{FS="\t";a=b=c=d=e=f=g=0}{if(NR>3+burning)a+=$5;b+=$6;c+=$7;d+=$8;e+=$9;f+=$10;g+=$11}END{print rep,r,gd,mod,a/(NR-3),b/(NR-3),c/(NR-3),d/(NR-3),e/(NR-3),f/(NR-3),g/(NR-3)}' $i >> parameters.out;done


#Reorganization for RWTY
########################
for i in r*;do r=$(echo $i|sed "s/r//");for j in $i/*;do mv $j $(echo $j|sed "s/\(.*\)\.\(...\)/\1.r$r.\2/");done;done
mkdir results;mv r*/* results/;mkdir results/MCC;mv results/*MCC* results/MCC
mkdir results/out; mv results/*.out results/out
#rmdir r*
for i in $(ls results/*.trees | xargs -n1 basename | sed "s/\.r.\.trees//" | uniq | sort);do mkdir results/$i; for j in results/$i*; do if [[ -f $j ]];then mv $j results/$i;fi;done;done
cd results
for i in *; do if [[ -d $i ]];then sbatch Rscript.sbatch ../../../../scripts_singlecrypt/rwty.R $i ;fi;done

#Tree analysis combined replicates
##################################
cd results
for i in tree*; do if [[ -d $i ]]; then name=$(echo $i|sed "s/.xml//g"); if [[ ! -f $i/${name}_MCC.trees ]]; then cd $i;sbatch -p private --mem 10000 <( echo -e '#!/bin/bash'"\nlogcombiner -trees -burnin 1000000 *xml.*.trees combined.trees.mcc\ntreeannotator combined.trees.mcc ${name}_MCC.trees" );cd..;fi;fi;done

##When the jobs finished
mkdir MCC; for i in tree*; do if [[ -d $i ]]; then rm -f $i/combined.trees.mcc;mv $i/*MCC* MCC;fi;done
sbatch -p private Rscript.sbatch ../../../../../scripts_singlecrypt/RF.R MCC ../trees/rooted_scaled.trees rf.csv

#Parameter analysis combined replicates
#######################################
burnin=500
tailn=$(( $burnin+4 ))
tailvar="+$tailn"
for i in tree*; do if [[ -d $i ]]; then name=$(echo $i|sed "s/.xml//g");rm -f $i/$name.combined.log;for j in $i/$name*r*.log; do if [[ ! -f $i/$name.combined.log ]];then head -n 3 $j > $i/$name.combined.log;fi; cat $j | tail -n $tailvar >>$i/$name.combined.log;done;fi;done
echo -e "rep r gd mod cnv.loss cnv.conversion treeModel.rootHeight luca_height luca_branch constant.popSize clock.rate" > parameters.out;for i in */tree*.combined.log;do rep=$(echo $i | sed "s/.*\/tree\(.*\)_r.*_gd.*_mod.*\.log/\1/");r=$(echo $i | sed "s/.*\/tree.*_r\(.*\)_gd.*_mod.*\.log/\1/");gd=$(echo $i | sed "s/.*\/tree.*_r.*_gd\(.*\)_mod.*\.log/\1/");mod=$(echo $i | sed "s/.*\/tree.*_r.*_gd.*_mod\(.*\)\.combined\.log/\1/");awk -v burnin=0 -v rep=$rep -v r=$r -v gd=$gd -v mod=$mod 'BEGIN{FS="\t";a=b=c=d=e=f=g=0}{if(NR>3+burnin)a+=$5;b+=$6;c+=$7;d+=$8;e+=$9;f+=$10;g+=$11}END{print rep,r,gd,mod,a/(NR-3),b/(NR-3),c/(NR-3),d/(NR-3),e/(NR-3),f/(NR-3),g/(NR-3)}' $i >> parameters.out;done


