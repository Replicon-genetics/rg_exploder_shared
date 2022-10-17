# Shell script to create curation heirarchies for GRCH37

root="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/data_sources/"
datadir="GRCH37_sequences_1000"
curation="_curation"
dirtxt="_dir"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"
mkdir $root$datadir$curation
cd $root$datadir$curation
pw=$PWD

jsonfile="_transcripts.json"
/bin/ls ../$datadir | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
mkdir $locus$curation
cd $locus$curation
ln -s ../../$datadir/$locus$dirtxt/ .
ln -s $locus$dirtxt/$locus$noseq .
ln -s $locus$dirtxt/$locus$locseq .
ln -s $locus$dirtxt/$locus$ensembl .
/bin/cp -p $locus$dirtxt/$locus$jsonfile .

done
prelocus="KRAS"
locus="KRAS_minus"
mkdir $locus$curation
cd $locus$curation
ln -s ../../$datadir/$locus$dirtxt/ .
ln -s $locus$dirtxt/$prelocus$noseq $locus$noseq
ln -s $locus$dirtxt/$prelocus$locseq $locus$locseq 
ln -s $locus$dirtxt/$prelocus$ensembl $locus$ensembl
/bin/cp -p $locus$dirtxt/$prelocus$jsonfile barf # ?
sed -e "s/$prelocus/$locus/g" barf > $locus$jsonfile # ?
/bin/rm barf # ?

cd $pw

