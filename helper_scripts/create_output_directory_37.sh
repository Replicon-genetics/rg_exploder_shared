# Shell script to create output directory heirarchies from input ones
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"
enddatadir="_37_1000"

inputdir=$rootdatadir"exploder_input"$enddatadir
outputdir=$rootdatadir"exploder_output"$enddatadir

cd $outputdir

/bin/ls $inputdir | grep -v json | while read locus
do
mkdir $locus
done