# Shell script to create output directory heirarchies from input ones

if [ "$rootRG" == "" ]
then
    echo "rootRG is unset - exiting script"
    exit 
fi

rootapplicationdir="$rootRG"/
data_sources="data_sources/"
rootdatadir=$rootapplicationdir$data_sources

output_37="exploder_output_37_1000"
output_dir37=$rootdatadir$output_37
input_dir37=$rootdatadir"exploder_input_37_1000"

output_38="exploder_output_38_1000"
output_dir38=$rootdatadir$output_38
input_dir38=$rootdatadir"exploder_input_38_1000"

stopit="stop"
input37="37"
input38="38"

jsontxt="json"

if [ "$1" == $input37 ]
then
    output=$output_37
    inputdir=$input_dir37
    outputdir=$output_dir37

elif [ "$1" == $input38 ]
then
    output=$output_38
    inputdir=$input_dir38
    outputdir=$output_dir38
else
    inputdir=$stopit
fi

if [[ "$inputdir" == "$stopit" ]]
then
    echo "Parameter entered should be $input37 or $input38, not $1"
else

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
        isjson=`echo $dirname  | cut -f2 -d"."`
        if [[ "$isjson" != "$jsontxt" ]]
        then
            mkdir $dirname
        fi
    done
    ls -l $dirname

    cd "$rootapplicationdir"/exploder_python
    /bin/rm output
    linkdir="../$data_sources$output"
    ln -s $linkdir output
    cd $pw

    cd $pw
fi