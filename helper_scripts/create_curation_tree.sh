# Shell script to create curation heirarchies

make_curation()
{
    prelocus=$1
    modlocus=$2
    PW=$3
    if [[ "$modlocus" == "$dironly" ]]
        then
            newlocus=$prelocus
        else
            newlocus=$prelocus$undertxt$modlocus
    fi

    mkdir $newlocus$curation
    cd $newlocus$curation
    ln -s ../../$targetdir/$newlocus$dirtxt/ .
    ln -s $newlocus$dirtxt/$prelocus$noseq $newlocus$noseq
    ln -s $newlocus$dirtxt/$prelocus$locseq $newlocus$locseq 
    ln -s $newlocus$dirtxt/$prelocus$ensembl $newlocus$ensembl

    if [[ "$modlocus" == "$dironly" ]]
    then
        /bin/cp -p $newlocus$dirtxt/$prelocus$jsonfile barf
        sed -e "s/$prelocus/$newlocus/g" barf > $newlocus$jsonfile
        /bin/rm barf
    else
        /bin/cp -p $newlocus$dirtxt/$prelocus$jsonfile $newlocus$jsonfile
    fi
    cd $PW
}

#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
#rootapplicationdir="/Users/caryodonnell/Desktop/Replicon/"
rootdatadir=$rootapplicationdir"data_sources/"

targetdir37="GRCH37_sequences_1000"
targetdir38="GRCH38_sequences_1000"

curation="_curation"
dirtxt="_dir"
dironly="dir"
undertxt="_"
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

    /bin/ls ../$targetdir | while read dirname
    do
    first=`echo $dirname  | cut -f1 -d"_"`
    second=`echo $dirname  | cut -f2 -d"_"`
    make_curation $first $second $pw
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

