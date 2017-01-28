for i in *.xml; do sbatch beast.sbatch;done
for i in *.trees; do name=$(echo $i | sed "s/\.xml\.trees//"); sbatch -p private --mem 10000 <( echo -e '#!/bin/bash'"\ntreeannotator -burninTrees 1000 $i ${name}_MCC.trees" );done
for i in r*;do r=$(echo $i|sed "s/r//");for j in $i/*;do mv $j $(echo $j|sed "s/\(.*\)\.\(...\)/\1.r$r.\2/");done;done
mkdir results;mv r*/* results/;mkdir results/MCC;mv results/*MCC* results/MCC
mkdir results/out; mv results/*.out results/out
rmdir r*
for i in $(ls results/*.trees | xargs -n1 basename | sed "s/\.r.\.trees//" | uniq | sort);do mkdir results/$i; for j in results/$i*; do if [[ -f $j ]];then mv $j results/$i;fi;done;done
for i in *phased*; do sbatch -c 3 Rscript.sbatch ../../../../scripts_singlecrypt/rwty.R $i 5 3;done

##Incompleted runs with dump_states
for i in *.xml ; do name=$(head $i | tail -n 1 | sed 's/<.*\"\(.*\)\">/\1/g');file=$(grep -l $name * | tail -n1);echo sbatch beast.sbatch -seed 25 -dump_every 10000000 -load_dump $file $i;done

##Putting incompleted runs together
for i in 500m/r1*
do
    ./scripts_singlecrypt/merge_incomplete_beast_runs.sh $i
done
