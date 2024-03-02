# Shell script to create input directory heirarchies from curation ones
baseindir="GRCH37_sequences_1000"
base_input_seq="/Users/caryodonnell/exploder_input_37_1000"
cd $base_input_seq
pw=$PWD
thispath="/Users/caryodonnell/Desktop/Replicon/"
curation="_curation"
curationdir=$baseindir$curation
dirtxt="_dir"
hap0="_000000.gb"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"

/bin/ls ../$curationdir | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
mkdir $locus
cd $locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0
rm *modified*.*
cd ..
done

/bin/rm -r loci.json/*
/bin/rmdir loci.json

locus="KRAS_minus"
mkdir $locus
cd $locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0
rm *modified*.*
cd ..
cd $pw

/bin/ls | while read dir
do
echo "$dir :"
grep VERSION $dir/*$hap0
done

/bin/ln -s ../$curationdir/loci.json .
