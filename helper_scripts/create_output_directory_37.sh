# Shell script to create output directory heirarchies from input ones
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"

inputdir=$rootdatadir"exploder_input_37_1000"
outputdir=$rootdatadir"exploder_output_37_1000"

mkdir $outputdir
cd $outputdir

/bin/ls $inputdir | grep -v json | while read locus
do
mkdir $locus
done