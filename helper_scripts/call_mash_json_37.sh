# Shell script to run mash_json.py from curation directory
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"
datadir="GRCH38_sequences_1000"
curation="_curation"
app="helper_python/mash_json.py"


cd $rootdatadir$datadir$curation

python3 $rootapplicationdir$app
