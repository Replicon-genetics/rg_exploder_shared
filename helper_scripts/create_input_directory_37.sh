# Shell script to create input AND output directory trees from curation directories

root="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/data_sources/"
datadir="GRCH37_sequences_1000"
base_input_seq=$root"exploder_input_37_1000"
base_output_seq=$root"exploder_output_37_1000"
mkdir $base_input_seq
mkdir $base_output_seq
pw=$PWD
cd $base_input_seq

curation="_curation"
curationdir=$datadir$curation
dirtxt="_dir"
hap0="_000000.gb"
noseq="_noseq.gb"
locseq="_locseq.gb"
ensembl="_ensembl"

/bin/ls $root$curationdir | while read dir
do
locus=`echo $dir  | cut -f1 -d"_"`
mkdir $base_output_seq/$locus
mkdir $base_input_seq/$locus
cd $base_input_seq/$locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0
done

cd $base_input_seq
locus="KRAS_minus"
mkdir $base_output_seq/$locus
mkdir $base_input_seq/$locus
cd $base_input_seq/$locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0

# Custom changes
locus="AK2"
cd $base_input_seq/$locus
rm AK2_locseq_bak.gb
rm AK2_locseq_modified.gb
# end of Custom changes

cd $base_output_seq
/bin/rmdir loci.json

cd $base_input_seq
/bin/rm -r loci.json/*
/bin/rmdir loci.json

/bin/ls | while read dir
do
echo "$dir :"
grep VERSION $dir/*$hap0
done

/bin/ln -s ../$curationdir/loci.json .

cd $pw
# Make sure you have run mash_json.py to create loci.json
# Now run RG_exploder_globals_make
