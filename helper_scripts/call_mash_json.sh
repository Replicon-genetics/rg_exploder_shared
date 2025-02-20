# Shell script to run mash_json.py from curation directory

if [ "$rootRG" == "" ]
then
    echo "rootRG is unset - exiting script"
    exit 
fi

rootapplicationdir="$rootRG"/

rootdatadir=$rootapplicationdir"data_sources/"

#input_seq37=$rootdatadir"exploder_input_37_1000"
#input_seq38=$rootdatadir"exploder_input_38_1000"

curation_seq37=$rootdatadir"GRCH37_sequences_1000_curation"
curation_seq38=$rootdatadir"GRCH38_sequences_1000_curation"

app="helper_python/mash_json.py"

stopit="stop"
input37="37"
input38="38"

if [ "$1" == $input37 ]
then
    inputdir=$input_seq37
    curationdir=$curation_seq37
elif [ "$1" == $input38 ]
    then
        inputdir=$input_seq38
        curationdir=$curation_seq38
    else
        inputdir=$stopit
fi

if [[ "$inputdir" != "$stopit" ]]
then
    #echo "valid directory $inputdir"
    #echo " python3 $rootapplicationdir$app"
    #cd $inputdir
    echo "writing loci.json to $curationdir"
    cd $curationdir
    python3 $rootapplicationdir$app

else
    echo "Parameter entered should be $input37 or $input38, not $1"
fi

