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
python $SCRIPTS_DIR/../add_cenan_years.py -i trees/trees.tree -o trees/rooted_scaled.trees -gt 365 -od 20 -ox xml/tree -a 40 ##Adds the outgroup and scales the tree in time units
mkdir alignments

rates=(0.0001 0.001 0.01);gds=(0 0.05);i=0;while read tree; do echo $tree > trees/tree${i}.tree;for gd in ${gds[*]}; do for rate in ${rates[*]}; do python $SCRIPTS_DIR/simulate.py -d $rate -c $rate -g $rate -gd $gd -n 100 --xml True -i trees/tree${i}.tree -o xml/tree${i}_r${rate}_gd${gd}.xml.temp; cat xml/tree${i}.xml > xml/tree${i}_r${rate}_gd${gd}.xml; cat xml/tree${i}_r${rate}_gd${gd}.xml.temp >> xml/tree${i}_r${rate}_gd${gd}.xml; rm -f tree${i}_r${rate}_gd${gd}.xml.temp;done;done;i=$(($i+1));done < trees/rooted_scaled.trees
