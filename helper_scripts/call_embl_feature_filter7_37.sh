# Shell script to list the Ensemble downloaded GB zipped files and run the filter program 
rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootdatadir=$rootapplicationdir"data_sources/"

jsonfile="_transcripts.json"
ensembl="_ensembl"
curation="_curation"
dirtxt="_dir"
# -j for longer json file with mRNA and CDS join data
jopt="-j"
#jopt=""
# -gopt to correctly label release
gopt="37"

cd $rootdatadir$thisdata
pw=$PWD
/bin/ls $1 | while read dir
do
modlocus=`echo $dir | sed -e "s/$dirtxt//" `
locus=`echo $dir  | cut -f1 -d"_"`
echo $modlocus  $locus
cd $pw/$1/$dir
mv ensembl.txt.gz $locus$ensembl.gz
gunzip $locus$ensembl.gz
python3 /Users/caryodonnell/mytools/embl_feature_filter7.py -i $locus$ensembl -a $jopt -g $gopt
cp -p $locus$jsonfile ../../$thisdata$curation/$modlocus$curation/.
done
cd  $pw


locus="KRAS"
modlocus="KRAS_minus"
cd $modlocus$dirtxt
sed -i -- "s/$locus/$modlocus/g" $locus$jsonfile
cp -p $locus$jsonfile ../../$thisdata$curation/$modlocus$curation/$modlocus$jsonfile
/bin/rm ../../$thisdata$curation/$modlocus$curation/$locus$jsonfile
/bin/rm $locus$jsonfile--
cd  $pw

# Next step is: sh call_mash_json_37.sh