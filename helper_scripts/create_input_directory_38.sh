# Shell script to create input directory heirarchies from curation ones
#baseindir="GRCH37_sequences_1000"
#base_input_seq="/Users/caryodonnell/input_sequences_ex03"
baseindir="GRCH38_sequences_1000"
base_input_seq="/Users/caryodonnell/exploder_input_38_1000"
cd $base_input_seq
pw=$PWD
thispath="/Users/caryodonnell/Desktop/Replicon/"
curation="_curation"
targetdir=$baseindir$curation
dirtxt="_dir"
hap0="_00000.gb"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"

jsonfile="_transcripts.json"
/bin/ls ../$targetdir | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
mkdir $locus
cd $locus
ln -s ../../$targetdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus$noseq $locus$hap0
ln -s $locus$dirtxt/$locus$locseq .
cd ..
done

rm -r loci.json

locus="PTEN_a"
mkdir $locus
cd $locus
ln -s ../../$targetdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus$noseq $locus$hap0
ln -s $locus$dirtxt/$locus$locseq .
cd ..

locus="KRAS_minus"
mkdir $locus
cd $locus
ln -s ../../$targetdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus$noseq $locus$hap0
ln -s $locus$dirtxt/$locus$locseq .
cd ..
cd $pw

locus="BRCA2_minus"
mkdir $locus
cd $locus
ln -s ../../$targetdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus$noseq $locus$hap0
ln -s $locus$dirtxt/$locus$locseq .
cd ..
cd $pw

/bin/ls | while read dir
do
echo "$dir :"
grep VERSION $dir/*$hap0
done


