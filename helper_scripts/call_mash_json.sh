# Shell script to run mash_json.py from curation directory
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"
curation="_curation"
app="helper_python/mash_json.py"

GRCHver="stop"

input37="GRCH37"
input38="GRCH38"
datadir37="GRCH37_sequences_1000"
datadir38="GRCH38_sequences_1000"

appropriate37=$rootdatadir$datadir37$curation
appropriate38=$rootdatadir$datadir38$curation

if echo "$PWD" | grep -q "$input37"; then GRCHver=input37 ; fi 
if echo "$PWD" | grep -q "$input38"; then GRCHver=input38 ; fi 


if [ $GRCHver != "stop" ]
then
    echo "valid directory $PWD"
    ##cd $rootdatadir$datadir$curation
    python3 $rootapplicationdir$app
else
    echo "$PWD is not an appropriate directory to run this."
    echo " For $input37 cd to $appropriate37"
    echo " For $input38 cd to $appropriate38"
fi

