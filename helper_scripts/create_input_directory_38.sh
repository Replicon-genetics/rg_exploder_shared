# Shell script to create input AND output directory trees from curation directories

root="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/data_sources/"
datadir="GRCH38_sequences_1000"
base_input_seq=$root"exploder_input_38_1000"
base_output_seq=$root"exploder_output_38_1000"
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

locus="KRAS_minus"
mkdir $base_output_seq/$locus
mkdir $base_input_seq/$locus
cd $base_input_seq/$locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0

locus="BRCA2_minus"
mkdir $base_output_seq/$locus
mkdir $base_input_seq/$locus
cd $base_input_seq/$locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0

locus="PTEN_a"
mkdir $base_output_seq/$locus
mkdir $base_input_seq/$locus
cd $base_input_seq/$locus
ln -s ../../$curationdir/$locus$curation $locus$dirtxt
ln -s $locus$dirtxt/$locus*.gb .
mv $locus$noseq $locus$hap0

# Custom changes
locus="ATM"
cd $base_input_seq/$locus
mv ATM_parent1.gb ATM_hap1.gb
mv ATM_parent2.gb ATM_hap2.gb
mv ATM_parent1_som1.gb ATM_som1.gb
mv ATM_parent2_som2.gb ATM_som2.gb
mv ATM_var_set1.gb ATM_som3.gb

locus="EGFR"
cd $base_input_seq/$locus
mv EGFR_parent1.gb EGFR_hap1.gb
mv EGFR_parent2.gb EGFR_hap2.gb
mv EGFR_parent1_somatic_var1.gb EGFR_som1.gb
mv EGFR_parent2_somatic_var14.gb EGFR_som14.gb
mv EGFR_parent1_somatic_var15.gb EGFR_som15.gb
mv EGFR_parent2_somatic_var16.gb EGFR_som16.gb
mv EGFR_parent1_somatic_var32.gb EGFR_som32.gb
mv EGFR_parent2_somatic_var39.gb EGFR_som39.gb

locus="KRAS"
cd $base_input_seq/$locus
mv KRAS_var2.gb KRAS_hap1.gb
mv KRAS_var3.gb KRAS_hap2.gb
rm KRAS_plus_testsetA_var.gb

locus="KRAS_minus"
cd $base_input_seq/$locus
mv KRAS_minus_var2.gb KRAS_minus_hap1.gb
mv KRAS_minus_var3.gb KRAS_minus_hap2.gb
#rm KRAS_plus_testsetA_var.gb
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
