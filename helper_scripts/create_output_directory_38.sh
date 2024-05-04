# Shell script to create output directory heirarchies from input ones
rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootdatadir=$rootapplicationdir"data_sources/"

inputdir=$rootdatadir"exploder_input_38_1000"
outputdir=$rootdatadir"exploder_output_38_1000"

cd $outputdir

/bin/ls $inputdir | grep -v json | while read locus
do
mkdir $locus
done