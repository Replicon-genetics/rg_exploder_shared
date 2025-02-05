# Shell script to run mash_json.py from curation directory
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"
curation="_curation"
app="helper_python/mash_json.py"

input37="37"
input38="38"
datadir37="GRCH37_sequences_1000"
datadir38="GRCH37_sequences_1000"


if [ "$1" == $input37 ]
then
    datadir=$input37
elif [ "$1" == $input38 ]
    then
      datadir=$input38
else
    datadir="stop"
fi

if [ $datadir != "stop" ]
then
    echo "valid parameter $1"
    #cd $rootdatadir$datadir$curation
    #python3 $rootapplicationdir$app
else
    echo "invalid parameter '$1': choose 37 or 38"
fi

