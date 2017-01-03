##Putting incompleted runs together

if [[ ${#ARGV[*]} -ne 1 ]]
then
    echo "Usage: $0 directory"
    exit
fi

cd $1

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
            echo $istate $j $stop
            if [[ $istate -eq 0 ]]
            then
               awk -v stop=$stop '{if ($1 != stop ) {print $0} else {exit}}' $j > $name.test
            else
                tail -n +4 $j | awk -v stop=$stop '{if ($1 != stop ) {print $0} else {exit}}' >> $name.test
            fi
            istate=$(( $istate + 1 )) 
        done

done
