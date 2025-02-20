# Shell script to create input directory heirarchies from curation ones

if [ "$rootRG" == "" ]
then
    echo "rootRG is unset - exiting script"
    exit 
fi

rootapplicationdir="$rootRG"/
data_sources="data_sources/"
rootdatadir=$rootapplicationdir$data_sources

targetdir37="GRCH37_sequences_1000"
input_seq37="exploder_input_37_1000"

targetdir38="GRCH38_sequences_1000"
input_seq38="exploder_input_38_1000"

curation="_curation"
dirtxt="_dir"

hap0="_000000.gb"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"
jsonfile="_transcripts.json"

stopit="stop"
input37="37"
input38="38"

if [ "$1" == $input37 ]
then
    targetdir=$targetdir37
    inputseq=$input_seq37
    
elif [ "$1" == $input38 ]
then
    targetdir=$targetdir38
    inputseq=$input_seq38
else
    targetdir=$stopit
fi

# Do this for 37 & 38
if [ "$targetdir" != $stopit ] 
then
    inputdir=$rootdatadir$inputseq
    if [ -e "$inputdir" ]
    then 
        echo $inputdir "exists"
    else
        echo  "creating "$inputdir
        mkdir $inputdir
    fi
    
    cd $inputdir
    pw=$PWD
    curdir=$targetdir$curation

    /bin/ls ../$targetdir | while read dir
    do
    modlocus=`echo $dir | sed -e "s/$dirtxt//" `
    locus=`echo $dir  | cut -f1 -d"_"`
    echo "modlocus:" $modlocus  "locus:" $locus

    mkdir $modlocus
    cd $modlocus

    ln -s ../../$curdir/$modlocus$curation/ .
    ln -s $modlocus$curation/*.gb .
    ln -s $modlocus$curation/$modlocus$jsonfile .
    /bin/mv $modlocus$noseq $modlocus$hap0
    cd $pw
    done

    cd $pw
    /bin/ls | while read dir
    do
    locus=`echo $dir  | cut -f1 -d"_"`
    echo "$dir :"
    grep VERSION $dir/*.gb
    done

    cd "$rootapplicationdir"/exploder_python
    /bin/rm input

    linkdir="../$data_sources$inputseq"
    ln -s $linkdir input
    cd $pw

else
    echo "Parameter entered should be $input37 or $input38, not $1"
fi

