# Shell script to list the Ensemble downloaded GB zipped files and run the filter program 
datadir="GRCH38_sequences_1000"
thispath="/Users/caryodonnell/Desktop/Replicon/"
jsonfile="_transcripts.json"
ensembl="_ensembl"
curation="_curation"
dirtxt="_dir"
# -j for longer json file with mRNA and CDS join data
#jopt="-j"
jopt=""
# -gopt to correctly label release
gopt="38"

cd $thispath$datadir
pw=$PWD
/bin/ls $1 | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
cd $pw/$1/$dir
mv ensembl.txt.gz $locus$ensembl.gz
gunzip $locus$ensembl.gz
#python3 /Users/caryodonnell/mytools/embl_feature_filter7.py -i $locus$ensembl -a

python3 /Users/caryodonnell/mytools/embl_feature_filter7.py -i $locus$ensembl -a $jopt -g $gopt
cp -p $locus$jsonfile ../../$datadir$curation/$locus$curation/.
done

prelocus="PTEN"
locus="PTEN_a"
cd $locus$dirtxt
python3 /Users/caryodonnell/mytools/embl_feature_filter7.py -i $prelocus$ensembl -a $jopt -g $gopt
sed -i -- 's/PTEN/PTEN_a/g' $prelocus$jsonfile
cp -p $prelocus$jsonfile ../../$datadir$curation/$locus$curation/$locus$jsonfile
/bin/rm $prelocus$jsonfile--
cd ..

prelocus="KRAS"
locus="KRAS_minus"
cd $locus$dirtxt
python3 /Users/caryodonnell/mytools/embl_feature_filter7.py -i $prelocus$ensembl -a $jopt -g $gopt
sed -i -- 's/KRAS/KRAS_minus/g' $prelocus$jsonfile
cp -p $prelocus$jsonfile ../../$datadir$curation/$locus$curation/$locus$jsonfile
/bin/rm $prelocus$jsonfile--
cd ..

cd  $pw
