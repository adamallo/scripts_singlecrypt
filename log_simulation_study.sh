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
rates=(0.0001 0.001 0.01);gds=(0 0.05 "random");mods=("none" "baseline" "max2");i=0;seed=1;while read tree; do echo $tree > trees/tree${i}.tree;for gd in ${gds[*]}; do for rate in ${rates[*]}; do for mod in ${mods[*]}; do echo "Simulating tree${i}_r${rate}_gd${gd}_mod${mod}.xml"; if [[ $gd == "random" ]]; then python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -gd 0 --randomGD 2 -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed; else python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -gd $gd -n 100 --xml -a 40 --mod $mod -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}_mod${mod}.xml --seed $seed;fi;done;seed=$(($seed+1));done;done;i=$(($i+1));done < trees/rooted_scaled.trees

