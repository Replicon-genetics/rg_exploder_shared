# Shell script to create input directory heirarchies from curation ones
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootdatadir=$rootapplicationdir"data_sources/"

baseindir="GRCH37_sequences_1000"
base_input_seq=$rootdatadir"exploder_input_37_1000"

cd $base_input_seq
pw=$PWD
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

locus="KRAS_minus"
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


