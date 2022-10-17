# Shell script to create output directory heirarchies from input ones
outputdir="/Users/caryodonnell/mravn_rg_exploder/python/output"
inputdir="/Users/caryodonnell/mravn_rg_exploder/python/input"
cd $outputdir

/bin/ls $inputdir | grep -v json | while read locus
do
mkdir $locus
done