# Shell script to create output directory heirarchies from input ones

if [ "$rootRG" == "" ]
then
    echo "rootRG is unset - exiting script"
    exit 
fi

rootapplicationdir="$rootRG"/

rootdatadir=$rootapplicationdir"data_sources/"

output_seq37=$rootdatadir"exploder_output_37_1000"
input_seq37=$rootdatadir"exploder_input_37_1000"

output_seq38=$rootdatadir"exploder_output_38_1000"
input_seq38=$rootdatadir"exploder_input_38_1000"

stopit="stop"
input37="37"
input38="38"

jsontxt="json"

if [ "$1" == $input37 ]
then
    inputdir=$input_seq37
    outputdir=$output_seq37

elif [ "$1" == $input38 ]
then
    inputdir=$input_seq38
    outputdir=$output_seq38
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
    ln -s $outputdir output

    cd $pw
fi