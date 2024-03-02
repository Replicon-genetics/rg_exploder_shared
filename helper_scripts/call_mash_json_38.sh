# Shell script to run mash_json.py from curation directory
datadir="GRCH38_sequences_1000"
thispath="/Users/caryodonnell/Desktop/Replicon/"
curation="_curation"

cd $thispath$datadir$curation
python3 /Users/caryodonnell/mytools/mash_json.py