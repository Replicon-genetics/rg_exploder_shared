# Copy latest versions across to shared scripts area

rootdevdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
devhelp_scriptdir=$rootdevdir"helper_scripts/"
devdocsdir=$rootdevdir"documents/"

devAWSdocdir=$devdocsdir"AWS_data_management/"
devnotesdir=$devdocsdir"development_journals/"

devmaintdir=$devdocsdir"maintenance_guides/"
devgenedatacurdir=$devmaintdir"Gene_data_curation/"
devpresentationdir=$devdocsdir"presentations/"

rootproddir="/Users/caryodonnell/Documents/repositories/mravn/rg_exploder/"
prodpyodidedir=$rootproddir"pyodide/"
prodwebdistdir=$rootproddir"webdist/"

rootsharedir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
sharehelp_scriptdir=$rootsharedir"helper_scripts/"
sharedocsdir=$rootsharedir"documents/"

shareAWSdocdir=$sharedocsdir"AWS_data_management/"
sharenotesdir=$sharedocsdir"development_notes/"

sharemaintdir=$sharedocsdir"maintenance_guides/"
sharegenedatacurdir=$sharemaintdir"Gene_data_curation/"

sharepyodidedir=$rootsharedir"pyodide/"
sharewebdistdir=$rootsharedir"webdist/"

same="."

cd $devhelp_scriptdir
cp -p call_embl_feature_filter.sh $sharehelp_scriptdir$same
cp -p call_mash_json.sh $sharehelp_scriptdir$same

cp -p create_curation_tree.sh $sharehelp_scriptdir$same
cp -p create_input_directory.sh $sharehelp_scriptdir$same
cp -p create_output_directory.sh $sharehelp_scriptdir$same

cp -p switch_links.sh $sharehelp_scriptdir$same
cp -p switch_Biopython_fix.sh $sharehelp_scriptdir$same

cp -p copy_python_to_shared.sh $sharehelp_scriptdir$same
cp -p copy_scripts_to_shared.sh  $sharehelp_scriptdir$same

cp -p make_rootRG.sh $sharehelp_scriptdir$same
cp -p test_rootRG.sh $sharehelp_scriptdir$same


cp -p $devAWSdocdir"AWS_deployments_register.ods"  $shareAWSdocdir$same
cp -p $devnotesdir"exploder_issues_log.ods" $sharenotesdir$same
cp -p $devgenedatacurdir"Replicon-specified_variations.ods" $sharegenedatacurdir$same
cp -p $devmaintdir"Synthetic_Reads_Generator_webApp_updating."* $sharemaintdir$same
cp -p $devmaintdir"Data_management."* $sharemaintdir$same

cp -p $prodpyodidedir"copy_build_to_webdist.sh" $sharepyodidedir$same
cp -p $prodpyodidedir"delete_old_build.sh" $sharepyodidedir$same
cp -p $prodpyodidedir"tidy_dist.sh" $sharepyodidedir$same

cp -p $prodpyodidedir"build/webworker.js" $sharepyodidedir"build/."
cp -p $prodpyodidedir"packages/rg_exploder/rg_exploder/RG_exploder_io.py" $sharepyodidedir"packages/rg_exploder/rg_exploder/."
cp -p $prodpyodidedir"packages/rg_exploder/rg_exploder/RG_exploder_webmain.py" $sharepyodidedir"packages/rg_exploder/rg_exploder/."
cp -p $prodpyodidedir"packages/rg_exploder/rg_exploder/setup.py" $sharepyodidedir"packages/rg_exploder/rg_exploder/."


cp -p $prodwebdistdir"public/webworker.js" $sharewebdistdir"public/."
cp -p $prodwebdistdir"src/App.vue" $sharewebdistdir"src/."
cp -p $prodwebdistdir"src/components/"* $sharewebdistdir"src/components/."




