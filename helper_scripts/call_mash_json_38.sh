# Shell script to run mash_json.py from curation directory
dataroot="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/data_sources/"
pythonroot="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/helper_python"
targetdir="GRCH38_sequences_1000"
curation="_curation"

cd $dataroot$targetdir$curation
python3 $pythonroot/mash_json.py
