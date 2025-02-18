# Copy latest versions across to shared python area
rootdevdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
devpythondir=$rootdevdir"exploder_python/"
devhelp_pythondir=$rootdevdir"helper_python/"

rootsharedir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
sharepythondir=$rootsharedir"exploder_python/"
sharehelp_pythondir=$rootsharedir"helper_python/"

same="/."
cd $devpythondir
cp -p RG_exploder_globals.py $sharepythondir$same
cp -p RG_exploder_globals_make.py $sharepythondir$same
cp -p RG_exploder_gui.py $sharepythondir$same
cp -p RG_exploder_io.py $sharepythondir$same
cp -p RG_exploder_main.py $sharepythondir$same
cp -p RG_exploder_process.py $sharepythondir$same
cp -p RG_exploder_builder.py $sharepythondir$same
cp -p Biopython_fix_metro.py $sharepythondir/Biopython_fix.py

cd $devhelp_pythondir
cp -p embl_feature_filter8.py $sharehelp_pythondir$same
cp -p embl_feature_filter_revise.py $sharehelp_pythondir$same
cp -p mash_json.py $sharehelp_pythondir$same
cp -p RG_exploder_io.py $sharehelp_pythondir$same
cp -p gaussian_test.py $sharehelp_pythondir$same
cp -p rayleigh_test.py $sharehelp_pythondir$same
cp -p List_LRGs_GRCh38.txt $sharehelp_pythondir$same
cp -p List_LRGs_GRCh37.txt $sharehelp_pythondir$same
cp -p config_setup.json $sharehelp_pythondir/config.json

cd $sharehelp_pythondir
ln -s ../exploder_python/RG_exploder_globals.py .
ln -s ../exploder_python/RG_exploder_io.py .
ln -s ../exploder_python/RG_exploder_process.py . 
ln -s ../exploder_python/Biopython_fix.py .








