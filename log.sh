#CSV pre-processing
for i in *_data.csv;do name=$(basename -s .csv $i);sed -e 's/".*"//g' -e 's/ /_/g' -e 's/^M//' $i >${name}_final.csv;done #Remove comments from csv (they content commas).

#Generation of files with leaf dates
mv 852_data_final.csv 852_prepro.csv; awk 'BEGIN{FS=",";OFS=","}{$17=",";print $0}' 852_prepro.csv > 852_data_final.csv ## These files do not have the second column of ids in the same column. This adds columns to fix that #Padding 1
mv 391_data_final.csv 391_prepro.csv; awk 'BEGIN{FS=",";OFS=","}{$17=",,";print $0}' 391_prepro.csv > 391_data_final.csv #Padding 2
for i in *_final.csv; do name=$(basename -s _final.csv $i);gawk 'BEGIN{FS=","}{
if($16!="" && $4!="" && $4!="EndoDate")
{
	date=$4
	gsub(/\-/, FS,date)
	split(date,array)
	month= 1 + (index("JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC", toupper(array[2])) - 1) / 3
	if(array[3]>=16)
	{
		year=1900+array[3]
	}else{
		year=2000+array[3]
	}
	timestamp= mktime(sprintf("%d %d %d %d %d %d",year,month,array[1],0,0,0))
	print($16,timestamp)
}
if($24!="" && $4!="" && $4!="EndoDate")
{
        date=$4
        gsub(/\-/, FS,date)
        split(date,array)
        month= 1 + (index("JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC", toupper(array[2])) - 1) / 3
        if(array[3]>=16)
        {
                year=1900+array[3]
        }else{
                year=2000+array[3]
        }
        timestamp= mktime(sprintf("%d %d %d %d %d %d",year,month,array[1],0,0,0))
        print($24,timestamp)
}
	
}' $i > ${name}_timestamps.txt;done


