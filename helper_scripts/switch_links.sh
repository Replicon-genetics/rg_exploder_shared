# Shell script to switch links of input & output directories between 37 & 38 - takes 1 parameter, defaults to 38
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
#rootapplicationdir="/Users/caryodonnell/Desktop/Replicon/"

rootprogdir=$rootapplicationdir"exploder_python/"


data_sources="data_sources/"
inputdirhead="exploder_input_"
inputdirtail="_1000"
input37="37"
input38="38"
outputdirhead="exploder_output_"
outputdirtail="_1000"
cd $rootprogdir
rm output # Don't worry about 'not found' warnings
rm input # Don't worry about 'not found' warnings
echo $1
if [ "$1" == $input37 ]
then
    outputdir=$input37
    inputdir=$input37
else
    outputdir=$input38
    inputdir=$input38
fi

cd $rootprogdir
ln -s "../"$data_sources$outputdirhead$outputdir$outputdirtail output
ln -s "../"$data_sources$inputdirhead$inputdir$inputdirtail input
