# Shell script to create output directory heirarchies from input ones
outputdir="/Users/caryodonnell/Replicon/exploder_output_37_1000"
inputdir="/Users/caryodonnell/input_sequences_ex03"
cd $outputdir

/bin/ls $inputdir | grep -v json | while read locus
do
mkdir $locus
done