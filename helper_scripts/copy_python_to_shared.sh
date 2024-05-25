# Copy latest versions across to shared python area
rootdevdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootdevpythondir=$rootdevdir"exploder_python/"
rootdevhelp_pythondir=$rootdevdir"helper_python/"
rootdevhelp_scriptdir=$rootdevdir"helper_scripts/"

rootsharedir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootsharepythondir=$rootsharedir"exploder_python/"
rootsharehelp_pythondir=$rootsharedir"helper_python/"
rootsharehelp_scriptdir=$rootsharedir"helper_scripts/"

same="/."
cd $rootdevpythondir
cp -p RG_exploder_globals.py $rootsharepythondir$same
cp -p RG_exploder_globals_make.py $rootsharepythondir$same
cp -p RG_exploder_gui.py $rootsharepythondir$same
cp -p RG_exploder_gui_addon.py $rootsharepythondir$same
cp -p RG_exploder_io.py $rootsharepythondir$same
cp -p RG_exploder_main.py $rootsharepythondir$same
cp -p RG_exploder_process.py $rootsharepythondir$same
cp -p Biopython_fix_metro.py $rootsharepythondir/Biopython_fix.py

cd $rootdevhelp_pythondir
cp -p embl_feature_filter7.py $rootsharehelp_pythondir$same
cp -p mash_json.py $rootsharehelp_pythondir$same
cp -p RG_exploder_io.py $rootsharehelp_pythondir$same
cp -p gaussian_test.py $rootsharehelp_pythondir$same
cp -p rayleigh_test.py $rootsharehelp_pythondir$same
cp -p List_LRGs_GRCh38.txt $rootsharehelp_pythondir$same
cp -p List_LRGs_GRCh37.txt $rootsharehelp_pythondir$same

cd $rootdevhelp_scriptdir
cp -p call_embl_feature_filter7_37.sh $rootsharehelp_scriptdir$same
cp -p call_embl_feature_filter7_38.sh $rootsharehelp_scriptdir$same
cp -p call_mash_json_37.sh $rootsharehelp_scriptdir$same
cp -p call_mash_json_38.sh $rootsharehelp_scriptdir$same

cp -p create_curation_tree_37.sh $rootsharehelp_scriptdir$same
cp -p create_curation_tree_38.sh $rootsharehelp_scriptdir$same
cp -p create_input_directory_37.sh $rootsharehelp_scriptdir$same
cp -p create_input_directory_38.sh $rootsharehelp_scriptdir$same
cp -p create_output_directory_37.sh $rootsharehelp_scriptdir$same
cp -p create_output_directory_38.sh $rootsharehelp_scriptdir$same

cp -p switch_links.sh $rootsharehelp_scriptdir$same
cp -p switch_Biopython_fix.sh $rootsharehelp_scriptdir$same

cp -p copy_python_to_shared.sh $rootsharehelp_scriptdir$same








