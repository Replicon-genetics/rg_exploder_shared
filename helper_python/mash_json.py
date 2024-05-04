#!/usr/local/bin/python3
Progver="mash_json.py"
ProgverDate="02-Mar-2022"
# This merges the content of all the individual {locus}_transcrupts.json files in the 'curation' directory
# into a single file loci.json. It is a source for RG_exploder_globals_make.py to create config.json
# which in turn is a source of variable values required to run RG_exploder_main.py
# with a bit more development work, the need for loci.json can be eliminated,
# at the cost of more reads from the input directory. But as this was set up for running from AWS,
# it was envisaged that this would keep reads down.
import json
from datetime import date
import RG_exploder_io as RG_io
from collections import OrderedDict

global config_file,json_transcripts,merged_transcriptconstants
config_file="config.json"
json_transcripts="_transcripts.json"
curation="_curation"
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
        print("config_exists %s as %s"%(config_exists,in_file))
        with RG_io.open_read(in_file) as json_data_file:
            config_in_data = json.load(json_data_file)
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
    #print("Folder %s"%infolder)
    locus=infolder.split(curation)[0]
    #print("Locus %s"%locus)
    target_file="./%s/%s%s"%(infolder,locus,json_transcripts)
    exists,data=config_file_in(target_file)
    if exists:
        locusentry=data["Reference_sequences"][locus]
        #print(locusentry)
        #merged_transcriptconstants.update({locus:locusentry})
        merged_transcriptconstants[locus]=locusentry

save_merged_json()
