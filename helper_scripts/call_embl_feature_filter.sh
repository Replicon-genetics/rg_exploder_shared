# Shell script to list the Ensemble downloaded GB zipped files and run the filter program 
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"

jsonfile="_transcripts.json"
ensembl="_ensembl"
curation="_curation"
dirtxt="_dir"

#Choose one of two filter programs
python_filter1="/Users/caryodonnell/mytools/embl_feature_filter8.py"
# embl_feature_filter8.py is a difficult-to-follow line-by-line parsing of the input file
python_filter2="/Users/caryodonnell/mytools/embl_feature_filter_revise.py"
# embl_feature_filter_revise.py uses Biopython and the RG_exploder modules to parse the input file
# The output differs only in the genbank-format output files, whereby there are additional 
# (blank) fields in KEYWORDS, SOURCE & ORGANISM for the Biopython. There seems no way to remove these
# The COMMENT file output is also constitutive, but includes a reference to the processing in that COMMENT
# There's a problem that appeared after Biopython 1.80, first noticed 1.83 
# reported, by this author as https://github.com/biopython/biopython/issues/4876
# Despite the 1.84 update, the only way to stop Biopython tripping up on the "non-standard" variant definitions
# is to set the input sequence length to zero on the lOCUS line
# eg: from 
#LOCUS       1 75011 bp DNA HTG 19-AUG-2022
#to
#LOCUS       1 0 bp DNA HTG 19-AUG-2022
# which embl_feature_filter_revise.py does by creating a temporary input file
# Biopython actually corrects this on the output, but this will become superfluous if Biopython ever gets this one fixed.

# Pick one!
# python_filter=$python_filter1
python_filter=$python_filter2

# -j for shorter json file with mRNA and CDS join data
#jopt="-j"
jopt=""

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
python3 $python_filter -i $locus$ensembl -a $jopt
cp -p $locus$jsonfile ../../$thisdata$curation/$modlocus$curation/.
done
cd  $pw

locus="PTEN"
modlocus="PTEN_a"
cd $modlocus$dirtxt
sed -i -- "s/$locus/$modlocus/g" $locus$jsonfile
cp -p $locus$jsonfile ../../$thisdata$curation/$modlocus$curation/$modlocus$jsonfile
/bin/rm ../../$thisdata$curation/$modlocus$curation/$locus$jsonfile
/bin/rm $locus$jsonfile--
cd  $pw

locus="KRAS"
modlocus="KRAS_minus"
cd $modlocus$dirtxt
sed -i -- "s/$locus/$modlocus/g" $locus$jsonfile
cp -p $locus$jsonfile ../../$thisdata$curation/$modlocus$curation/$modlocus$jsonfile
/bin/rm ../../$thisdata$curation/$modlocus$curation/$locus$jsonfile
/bin/rm $locus$jsonfile--
cd  $pw

locus="BRCA2"
modlocus="BRCA2_minus"
cd $modlocus$dirtxt
sed -i -- "s/$locus/$modlocus/g" $locus$jsonfile
cp -p $locus$jsonfile ../../$thisdata$curation/$modlocus$curation/$modlocus$jsonfile
/bin/rm ../../$thisdata$curation/$modlocus$curation/$locus$jsonfile
/bin/rm $locus$jsonfile--
cd  $pw
# Next step is: 'sh call_mash_json.sh'