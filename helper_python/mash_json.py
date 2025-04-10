#!/usr/local/bin/python3
Progver="mash_json.py"
ProgverDate="13-Feb-2025"

'''
This merges the content of multuple {locus}_transcrupts.json files as part of the pre-processing required to run the RG_exploder_*.py suite.
Consider using embl_feature_filter_revise.py instead

Copyright © 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025 ; Cary O'Donnell 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see the original repository at
    https://github.com/snowlizardz/rg_exploder_shared, or the licences at <https://www.gnu.org/licenses/>.


    Contact: syrgenreads@gmail.com

Created using python3 and BioPython

Longer description:
This merges the content of all the individual {locus}_transcrupts.json files in directories below the current directory
into a single file loci.json. It is a source for RG_exploder_globals_make.py to create config.json
which in turn is a source of variable values required to run RG_exploder_main.py
with a bit more development work, the need for loci.json can be eliminated,
at the cost of more reads from the input directory. But as this was set up for running from AWS,
it was envisaged that this would keep reads down.
'''
import json
from datetime import date
import RG_exploder_io as RG_io
from collections import OrderedDict

global config_file,json_transcripts,merged_transcriptconstants
config_file="config.json"
json_transcripts="_transcripts.json"
curation="_curation"
#curation="" # Was originally curation="_curation"
merged_transcriptconstants=OrderedDict()
#merged_transcriptconstants={}

today = date.today()
# Textual month, day and year	
r_date = today.strftime("%Y %B %d")
#print("r_date =%s"%r_date)

def config_file_in(in_file):
    config_exists = RG_io.is_file(in_file)
    config_in_data=OrderedDict()
    if config_exists:
        #print("config_exists %s as %s"%(config_exists,in_file))
        with RG_io.open_read(in_file) as json_data_file:
            config_in_data = json.load(json_data_file)
    else:
        print("%s does not exists"%in_file)
    return config_exists,config_in_data
# end of config_file_in()

def get_reference_seqs_configs():
    Reference_sequences=config_in_data["Reference_sequences"]
    return Reference_sequences
# end of set_reference_seqs_configs()

def get_locus_configs():
    global Reference_sequences
    for item in Reference_sequences:
        this_locus=config_in_data["Reference_sequences".item]

def save_merged_json():
    global merged_transcriptconstants
    sorted_merged=dict(sorted(merged_transcriptconstants.items()))
    out_data={
        "Reference_sequences":sorted_merged
        }
    config_file="loci.json"
    with RG_io.open_write(config_file) as write_file:
        json.dump(out_data, write_file,indent=4)    

directories=RG_io.list_files('.','*%s'%curation)
#directories=directories.sort()
print("directories %s"%directories)
for folder in directories:
    infolder=str(folder)
    print("Folder %s"%infolder)
    locus=infolder.split(curation)[0]
    #locus=infolder
    #print("Locus %s"%locus)
    #target_file="./%s/%s%s"%(infolder,locus,json_transcripts)
    target_file="./%s/%s%s"%(infolder,locus,json_transcripts)
    exists,data=config_file_in(target_file)
    if exists:
        locusentry=data["Reference_sequences"][locus]
        #print("found %s"%locusentry)
        #merged_transcriptconstants.update({locus:locusentry})
        merged_transcriptconstants[locus]=locusentry
    else:
        print("fail %s"%target_file)
save_merged_json()
