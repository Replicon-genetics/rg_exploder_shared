# Copy latest versions across to shared python area
rootdevdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootdevpythondir=$rootdevdir"exploder_python/"
rootdevhelp_scriptdir=$rootdevdir"helper_scripts/"

rootsharedir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
rootsharepythondir=$rootsharedir"exploder_python/"
rootsharehelp_scriptdir=$rootsharedir"helper_scripts/"

same="/."

cd $rootdevhelp_scriptdir
cp -p call_embl_feature_filter.sh $rootsharehelp_scriptdir$same
cp -p call_mash_json.sh $rootsharehelp_scriptdir$same

cp -p create_curation_tree_37.sh $rootsharehelp_scriptdir$same
cp -p create_curation_tree_38.sh $rootsharehelp_scriptdir$same
cp -p create_input_directory_37.sh $rootsharehelp_scriptdir$same
cp -p create_input_directory_38.sh $rootsharehelp_scriptdir$same
cp -p create_output_directory_37.sh $rootsharehelp_scriptdir$same
cp -p create_output_directory_38.sh $rootsharehelp_scriptdir$same

cp -p switch_links.sh $rootsharehelp_scriptdir$same
cp -p switch_Biopython_fix.sh $rootsharehelp_scriptdir$same

cp -p copy_python_to_shared.sh $rootsharehelp_scriptdir$same
cp -p copy_scripts_to_shared.sh  $rootsharehelp_scriptdir$same







