# Shell script to create curation heirarchies

if [ "$rootRG" == "" ]
then
    echo "rootRG is unset - exiting script"
    exit 
fi

rootapplicationdir="$rootRG"/

rootdatadir=$rootapplicationdir"data_sources/"

targetdir37="GRCH37_sequences_1000"
targetdir38="GRCH38_sequences_1000"

curation="_curation"
dirtxt="_dir"
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
elif [ "$1" == $input38 ]
then
    targetdir=$targetdir38
else
    targetdir=$stopit
fi

# Do this for 37 & 38
if [ "$targetdir" != $stopit ] 
then
    if [ -e "$rootdatadir$targetdir$curation" ]
    then 
        echo $rootdatadir$targetdir$curation "exists"
    else
        echo  "creating "$rootdatadir$targetdir$curation
        mkdir $rootdatadir$targetdir$curation
    fi

    cd $rootdatadir$targetdir$curation
    pw=$PWD

    /bin/ls ../$targetdir | while read dir
    do
    modlocus=`echo $dir | sed -e "s/$dirtxt//" `
    locus=`echo $dir  | cut -f1 -d"_"`
    echo "modlocus:" $modlocus  "locus:" $locus

    mkdir $modlocus$curation
    cd $modlocus$curation
    ln -s ../../$targetdir/$modlocus$dirtxt/ .
    ln -s $modlocus$dirtxt/$modlocus$noseq .
    ln -s $modlocus$dirtxt/$modlocus$locseq .
    ln -s $modlocus$dirtxt/$locus$ensembl $modlocus$ensembl
    /bin/cp -p $modlocus$dirtxt/$modlocus$jsonfile .
    cd $pw
    done

    cd $pw
    /bin/ls | while read dir
    do
    locus=`echo $dir  | cut -f1 -d"_"`
    echo "$dir :"
    grep VERSION $dir/$locus"_ensembl"  $dir/*.gb
    echo $dir/$locus"_transcripts.json" && grep Region $dir/*".json"
    done

else
    echo "Parameter entered should be $input37 or $input38, not $1"
fi

