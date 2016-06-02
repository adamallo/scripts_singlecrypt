for i in *.xml; do sbatch beast.sbatch;done
for i in *.trees; do name=$(echo $i | sed "s/\..*//"); treeannotator -burnin 500000 $i ${name}_MCC.trees;done
