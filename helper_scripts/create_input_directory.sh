# Shell script to create input directory heirarchies from curation ones

fill_curation()
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

    mkdir $newlocus
    cd $newlocus

    ln -s ../../$curdir/$newlocus$curation/ .
    ln -s $newlocus$curation/*.gb .
    /bin/mv $newlocus$noseq $newlocus$hap0
    cd $PW
}

#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
#rootapplicationdir="/Users/caryodonnell/Desktop/Replicon/"
rootdatadir=$rootapplicationdir"data_sources/"

targetdir37="GRCH37_sequences_1000"
input_seq37=$rootdatadir"exploder_input_37_1000"

targetdir38="GRCH38_sequences_1000"
input_seq38=$rootdatadir"exploder_input_38_1000"

curation="_curation"
dirtxt="_dir"
dironly="dir"
undertxt="_"

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
    input_seq=$input_seq37
elif [ "$1" == $input38 ]
then
    targetdir=$targetdir38
    input_seq=$input_seq38

else
    targetdir=$stopit
fi

# Do this for 37 & 38
if [ "$targetdir" != $stopit ] 
then
    if [ -e "$input_seq" ]
    then 
        echo $input_seq "exists"
    else
        echo  "creating "$input_seq
        mkdir $input_seq
    fi

    cd $input_seq
    pw=$PWD
    curdir=$targetdir$curation

    /bin/ls ../$targetdir | while read dirname
    do
    first=`echo $dirname  | cut -f1 -d"_"`
    second=`echo $dirname  | cut -f2 -d"_"`
    fill_curation $first $second $pw
    done

    cd $pw
    /bin/ls | while read dir
    do
    locus=`echo $dir  | cut -f1 -d"_"`
    echo "$dir :"
    grep VERSION $dir/*.gb
    done

else
    echo "Parameter entered should be $input37 or $input38, not $1"
fi

