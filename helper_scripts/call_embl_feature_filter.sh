# Shell script to list the Ensemble downloaded GB zipped files and run the filter program 

if [ "$rootRG" == "" ]
then
    echo "rootRG is unset - exiting script"
    exit 
fi

rootapplicationdir="$rootRG"/

rootdatadir=$rootapplicationdir"data_sources/"
targetdir37="GRCH37_sequences_1000"
targetdir38="GRCH38_sequences_1000"

filter8="helper_python/embl_feature_filter8.py"
filter_revise="helper_python/embl_feature_filter_revise.py"

dirtxt="_dir"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"
jsonfile="_transcripts.json"
stopit="stop"

# There are two filter programs ...
python_filter1=$rootapplicationdir$filter8
# embl_feature_filter8.py is a difficult-to-follow line-by-line parsing of the input file
python_filter2=$rootapplicationdir$filter_revise
# embl_feature_filter_revise.py uses Biopython and the RG_exploder modules to parse the input file
# The output differs only in the genbank-format output files, whereby there are additional, but blank,
# fields in KEYWORDS, SOURCE, ORGANISM for the Biopython. There seems no way to remove these
# The COMMENT file output is also constitutive, but includes a reference to the processing in that COMMENT
# There is a problem that appeared after Biopython 1.80, first noticed 1.83 
# reported, by this author as https://github.com/biopython/biopython/issues/4876
# Despite the 1.84 update, the only way to stop Biopython tripping up on the -non-standard- variant definitions
# is to set the input sequence length to zero on the lOCUS line
# eg: from 
# LOCUS       1 75011 bp DNA HTG 19-AUG-2022
# to
# LOCUS       1 0 bp DNA HTG 19-AUG-2022
# which embl_feature_filter_revise.py does by creating a temporary input file
# Biopython actually corrects this on the output, but this will become superfluous if Biopython ever gets this one fixed.

# ... Pick one:
#python_filter=$python_filter1
python_filter=$python_filter2

# -j for shorter json file without mRNA and CDS join data
# jopt="-j"
jopt=""

sourcedata=$PWD
echo "sourcedata:" $sourcedata

if [[ "$sourcedata" == "$rootdatadir$targetdir37" ]]
then
    targetdir=$targetdir37
elif [[ "$sourcedata" == "$rootdatadir$targetdir38" ]]
then
    targetdir=$targetdir38
else
    targetdir=$stopit
fi

if [ "$targetdir" != $stopit ] 
then
    /bin/ls $1 | while read dir
    do
    echo " ************* "
    modlocus=`echo $dir | sed -e "s/$dirtxt//" `
    locus=`echo $dir  | cut -f1 -d"_"`
    echo "modlocus:" $modlocus  "locus:" $locus
    cd $sourcedata/$1/$dir
    mv ensembl.txt.gz $locus$ensembl.gz
    gunzip $locus$ensembl.gz
    python3 $python_filter -i $locus$ensembl -a $jopt
    if [[ "$modlocus" != "$locus" ]]
    then
        sed -i -- "s/$locus/$modlocus/g" $locus$jsonfile
        /bin/rm $locus$jsonfile--
        echo "renaming $locus$jsonfile to $modlocus$jsonfile"
        /bin/mv $locus$jsonfile $modlocus$jsonfile
        echo "renaming $locus$noseq to $modlocus$noseq"
        /bin/mv $locus$noseq $modlocus$noseq
        echo "renaming $locus$locseq to $modlocus$locseq"
        /bin/mv $locus$locseq $modlocus$locseq 
    fi
    done
    cd  $pw
else
    echo "This is not a good starting directory. Try $rootdatadir$targetdir37 or $rootdatadir$targetdir38"
fi
# Next step is: 'sh call_mash_json.sh'