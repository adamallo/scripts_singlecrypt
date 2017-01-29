##Putting incompleted runs together

if [[ $# -ne 1 ]]
then
    echo "Usage: $0 directory"
    exit
fi

cd $1

eliminate_extra=500050000

for i in *.xml
do
	name=$(echo $i | sed "s/.xml//")
        unset states
        declare -a states
        istate=0
        for j in $name.log*
        do
            echo $j
            states[istate]=$(head -n 4 $j | tail -n 1 | awk '{print $1}')
            istate=$(( $istate + 1 ))
        done

        istate=0
        for j in $name.log*
        do
            stop=${states[$(($istate + 1 ))]}
            if [[ $eliminate_extra != "" ]] && [[ $stop -eq "" ]]
            then
                stop=$eliminate_extra
            fi
            
            if [[ $stop -eq "" ]]
            then
                stop_tree=""
            else
                stop_tree="STATE_$stop"
            fi

            echo $istate $j $stop

            id=$(echo $j | sed "s/$name.log_//")

            if [[ $istate -eq 0 ]]
            then
               awk -v stop=$stop '{if ($1 != stop ) {print $0} else {exit}}' $j > $name.merged.log
               awk -v stop=$stop_tree '{if ($2 != stop ) {print $0} else {exit}}' $name.trees_$id > $name.merged.trees
            else
                tail -n +4 $j | awk -v stop=$stop '{if ($1 != stop ) {print $0} else {exit}}' >> $name.merged.log
                awk -v stop=$stop_tree 'BEGIN{skip=1}{if(skip==1 && $1 == "tree"){skip=0} if (skip==0){if($2 != stop ) {print $0} else {exit}}}' $name.trees_$id >> $name.merged.trees
            fi
            istate=$(( $istate + 1 )) 
        done
	echo "End;" >> $name.merged.trees ###Last line. Otherwise the Nexus file does not match the standard and fails with some software
done
