for i in *.xml; do sbatch beast.sbatch;done
for i in *.trees; do name=$(echo $i | sed "s/\..*//"); treeannotator -burnin 500000 $i ${name}_MCC.trees;done
for i in run*;do run=$(echo $i|sed "s/run//");for j in $i/*;do mv $j $(echo $j|sed "s/\(.*\)\.\(...\)/\1.run$run.\2/");done;done
mkdir results;mv run*/* results/;mkdir results/MCC;mv results/*MCC* results/MCC
mkdir results/out; mv results/*.out results/out
rmdir run*
for i in $(ls results/*.trees | xargs -n1 basename | sed "s/\.run.\.trees//" | uniq | sort);do mkdir results/$i; for j in results/$i*; do if [[ -f $j ]];then mv $j results/$i;fi;done;done
