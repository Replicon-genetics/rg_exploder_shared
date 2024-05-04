# Shell script to create curation heirarchies
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"

targetdir="GRCH38_sequences_1000"

curation="_curation"
dirtxt="_dir"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"
cd $rootdatadir$targetdir$curation
pw=$PWD

jsonfile="_transcripts.json"
/bin/ls ../$targetdir | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
mkdir $locus$curation
cd $locus$curation
ln -s ../../$targetdir/$locus$dirtxt/ .
ln -s $locus$dirtxt/$locus$noseq .
ln -s $locus$dirtxt/$locus$locseq .
ln -s $locus$dirtxt/$locus$ensembl .
/bin/cp -p $locus$dirtxt/$locus$jsonfile .
cd ..
done
prelocus="PTEN"
locus="PTEN_a"
mkdir $locus$curation
cd $locus$curation
ln -s ../../$targetdir/$locus$dirtxt/ .
ln -s $locus$dirtxt/$prelocus$noseq $locus$noseq
ln -s $locus$dirtxt/$prelocus$locseq $locus$locseq 
ln -s $locus$dirtxt/$prelocus$ensembl $locus$ensembl
/bin/cp -p $locus$dirtxt/$prelocus$jsonfile barf
sed -e "s/$prelocus/$locus/g" barf > $locus$jsonfile
/bin/rm barf
cd ..
prelocus="KRAS"
locus="KRAS_minus"
mkdir $locus$curation
cd $locus$curation
ln -s ../../$targetdir/$locus$dirtxt/ .
ln -s $locus$dirtxt/$prelocus$noseq $locus$noseq
ln -s $locus$dirtxt/$prelocus$locseq $locus$locseq 
ln -s $locus$dirtxt/$prelocus$ensembl $locus$ensembl
/bin/cp -p $locus$dirtxt/$prelocus$jsonfile barf
sed -e "s/$prelocus/$locus/g" barf > $locus$jsonfile
/bin/rm barf
cd ..

prelocus="BRCA2"
locus="BRCA2_minus"
mkdir $locus$curation
cd $locus$curation
ln -s ../../$targetdir/$locus$dirtxt/ .
ln -s $locus$dirtxt/$prelocus$noseq $locus$noseq
ln -s $locus$dirtxt/$prelocus$locseq $locus$locseq 
ln -s $locus$dirtxt/$prelocus$ensembl $locus$ensembl
/bin/cp -p $locus$dirtxt/$prelocus$jsonfile barf
sed -e "s/$prelocus/$locus/g" barf > $locus$jsonfile
/bin/rm barf
cd ..

cd $pw