perl ../scripts/generate_genotypes.pl -i data/ -f xml -s 20 > data_generation.log

for i in *.xml; do sbatch beast.sbatch;done
for i in r*;do r=$(echo $i|sed "s/r//");for j in $i/*;do mv $j $(echo $j|sed "s/\(.*\)\.\(...\)/\1.r$r.\2/");done;done
mkdir results/out; mv results/*.out results/out
rmdir r*
for i in $(ls results/*.trees | xargs -n1 basename | sed "s/\.r.\.trees//" | uniq | sort);do mkdir results/$i; for j in results/$i*; do if [[ -f $j ]];then mv $j results/$i;fi;done;done
for i in *phased*; do sbatch -c 3 Rscript.sbatch ../../../../scripts_singlecrypt/rwty.R $i 5 3;done

for i in *phased*; do name=$(echo $i|sed "s/\///g");cd $i;sbatch -p private --mem 10000 <( echo -e '#!/bin/bash'"\nlogcombiner -trees -burnin 25000000 *.r*.trees combined.trees.mcc\ntreeannotator combined.trees.mcc ${name}_MCC.trees" );cd..;done 

mkdir MCC; for i in *phased*; do rm -f $i/combined.trees.mcc;mv $i/*MCC* MCC;done

##Incompleted runs with dump_states
for i in *.xml ; do name=$(head $i | tail -n 1 | sed 's/<.*\"\(.*\)\">/\1/g');file=$(grep -l $name * | tail -n1);echo sbatch beast.sbatch -seed 25 -dump_every 10000000 -load_dump $file $i;done

##Putting incompleted runs together
for i in 500m/r1*
do
    ./scripts_singlecrypt/merge_incomplete_beast_runs.sh $i
done

##Generating ages for statistical analyses
cd data
echo "patient type age" > ages.tsv;for i in *.xml; do patient=$(echo $i | sed "s/_phased.*//g");data=$(echo $i|sed "s/.*_100\(.*\)\.xml/\\1/g");if [[ $data != "_all" ]];then age=$(grep "<date" $i | sed "s/.*value=\"\([^\"]*\)\".*/\\1/g" | sort -n -r | uniq | head -1); echo $patient $data $age >> ages.tsv; fi;done

##Minor stats for the paper
#for i in `find . -regex ".*/[^29].*pool\.xml"`; do sed -n "/<taxon[^\/]*>/p" $i | wc -l ;done | awk 'BEGIN{a=0} {a=a+$1} END{print a}'
for i in `find . -regex ".*/.*pool\.xml"`; do sed -n "/<taxon[^\/]*>/p" $i | wc -l ;done | awk 'BEGIN{a=0} {a=a+$1} END{print a}'
#for i in `find . -regex ".*/[^29].*crypts\.xml"`; do sed -n "/<taxon[^\/]*>/p" $i | wc -l ;done | awk 'BEGIN{a=0} {a=a+$1} END{print a}'
for i in `find . -regex ".*/.*crypts\.xml"`; do sed -n "/<taxon[^\/]*>/p" $i | wc -l ;done | awk 'BEGIN{a=0} {a=a+$1} END{print a}'
