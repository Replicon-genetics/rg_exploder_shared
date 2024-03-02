# Shell script to list the Ensemble downloaded GB zipped files and run the filter program 
datadir="GRCH37_sequences_1000"
thispath="/Users/caryodonnell/Desktop/Replicon/"
jsonfile="_transcripts.json"
ensembl="_ensembl"
curation="_curation"
dirtxt="_dir"
# -j for longer json file with mRNA and CDS join data
jopt="-j"
#jopt=""
# -g to correctly label release
gopt="37"

cd $thispath$datadir
pw=$PWD
/bin/ls $1 | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
cd $pw/$1/$dir
mv ensembl.txt.gz $locus$ensembl.gz
gunzip $locus$ensembl.gz
python3 /Users/caryodonnell/mytools/embl_feature_filter7.py -i $locus$ensembl -a $jopt -g $gopt
cp -p $locus$jsonfile ../../$datadir$curation/$locus$curation/.
done
# Now a fix for the KRAS / KRAS_minus setup
cd $thispath$datadir
cp -p KRAS_dir/KRAS_transcripts.json ../$datadir$curation/KRAS_curation/.
cp -p KRAS_minus_dir/KRAS_transcripts.json ../$datadir$curation/KRAS_minus_curation/.
cd ../$datadir$curation/KRAS_minus_curation/
sed -e 's/KRAS/KRAS_minus/g' KRAS_transcripts.json > KRAS_minus_transcripts.json
/bin/rm KRAS_transcripts.json
cd  $pw
# Next step is: sh call_mash_json_37.sh
