# Shell script to create output directory heirarchies from input ones
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
#rootapplicationdir="/Users/caryodonnell/Desktop/Replicon/"
rootdatadir=$rootapplicationdir"data_sources/"

stopit="stop"
input37="37"
input38="38"

enddatadir37="_37_1000"
enddatadir38="_38_1000"

if [ "$1" == $input37 ]
then
    enddatadir=$enddatadir37
elif [ "$1" == $input38 ]
then
    enddatadir=$enddatadir38
else
    enddatadir=$stopit
fi


if [ $enddatadir == $stopit ]
then
    echo "Parameter entered should be $input37 or $input38, not $1"
else

    inputdir=$rootdatadir"exploder_input"$enddatadir
    outputdir=$rootdatadir"exploder_output"$enddatadir

 if [ -e "$outputdir" ]
    then 
        echo $outputdir "exists"
    else
        echo  "creating "$outputdir
        mkdir $outputdir
    fi

    cd $outputdir

    /bin/ls $inputdir | while read dirname
    do
    mkdir $dirname
    done
    ls -l $dirname
fi