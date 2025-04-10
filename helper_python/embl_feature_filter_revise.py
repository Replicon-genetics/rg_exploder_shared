#!/usr/local/bin/python3
Progver="embl_feature_filter_revise.py"
ProgverDate="09-Mar-2025"
'''

This reads in Genbank format files as part of the pre-processing required to run the RG_exploder_*.py suite

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
This processes the {locus}_Ensembl_download.gz file to eliminate unwanted items from the feature table
creating, optionally
a) A file with the target locus, mRNA features and sequence. {locus}_locseq.gb - this is the reference sequence for main application 
b) A file with the target locus, variant features, omits sequence. {locus}_var.gb
c) A file with the target locus, mRNA features, variant features and sequence. {locus}_filtered.gb

And constitutively:
d) A file with the target locus, no variants, no sequence {locus}_noseq.gb - this is the {locus}_000000.gb file for the main application
e) A file with the transcript ids in json format   {locus}_transcripts.json - contributor to config.json file

Unwanted only in that these features (primarily the database cross-reference
definitions, but also CDS, misc_RNA and exons) are unnecessary for the variants-processing applications.
The prime interests are:
   the mRNA join features to identify specific transcripts.
   the "variations" features, with dbsnp ids

This uses BioPython to produce (almost) the same output as embl_feature_filter8.py

####### NB: the command-line option -c supresses CDS inclusion in the {locus}_transcripts.json file ######


'''
import sys, getopt
import binascii  #needed by is_gz_file
import gzip
from functools import partial
import json
from datetime import date
from io import StringIO #python 3
import RG_exploder_io as RG_io
import RG_exploder_globals as RG_globals
import RG_exploder_process as RG_process
import copy
import os
import inspect

from RG_exploder_globals import process_file_configs_GUI

#BioPython imports

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord # Used in writing protein sequence

def read_config():
    #Initialise config
    #print ("%s"%inspect.getmodule(process_file_configs_GUI))
    splits=str(inspect.getmodule(process_file_configs_GUI)).split('/')
    path_to="/"
    for i in range(1,len(splits)-1):
        path_to+="%s%s"%(splits[i],"/")
    #print("%s"%pathonly)
    if RG_globals.process_file_configs_GUI("%sconfig.json"%path_to):
        print("config read")
        #RG_globals.is_print_diagnostics=True # Switch on diasgnostic prints from RG_process
        pass
    else:
        print("*** WARNING: config NOT read ***")
        pass

global empty_string,inputfile,outfeatfile,outvarfile
global outtempfile,outreffile,outtmphandle,r_date,mRNA_transcript_list,outfeatfile
global mRNA_join_list,last_generange
global target_locus,is_append_tmp_file,Locus

target_locus=""
last_generange="undefined"
empty_string=''
is_append_tmp_file=False
Locus="LOCUS       "

mRNA_transcript_list=dict()
mRNA_join_list=dict()
CDS_join_list=dict()

today = date.today()
# Textual month, day and year	
r_date = today.strftime("%Y %B %d")
this_date = today.strftime("%Y %B %d")
#print("r_date =%s"%r_date)

MANE_Select_dict={
          "AK2":"672715",
          "ATM":"675843",
          "BARD1":"260947",
          "BRCA1":"357654",
          "BRCA2":"380152",
          "BRIP1":"259008",
          "CDH1":"261769",
          "CDK12":"447079",
          "CHEK1":"438015",
          "CHEK2":"404276",
          "CIITA":"324288",
          "EGFR":"275493",
          "EPCAM":"263735",
          "FANCL":"233741",
          "KRAS":"311936",
          "KRAS_minus":"311936",
          "MLH1":"231790",
          "MSH2":"233146",
          "MSH6":"234420",
          "NBN":"265433",
          "NCF1":"289473",
          "NF1":"358273",
          "PALB2":"261584",
          "PPP2R2A":"380737",
          "PMS2":"265849",
          "PTEN":"380152",
          "PTEN_a":"undefined",
          "RAD51B":"471583",
          "RAD51C":"337432",
          "RAD51D":"345365",
          "RAD54L":"371975",
          "STK11":"326873",
          "TP53":"269305"
          }

LRG_id_dict={
          "AK2":"133",
          "ATM":"135",
          "BARD1":"297",
          "BRCA1":"292",
          "BRCA2":"293",
          "BRIP1":"300",
          "CDH1":"301",
          "CDK12":"1413",
          "CHEK1":"0000",
          "CHEK2":"302",
          "CIITA":"49",
          "EGFR":"304",
          "EPCAM":"215",
          "FANCL":"501",
          "KRAS":"344",
          "KRAS_minus":"344",
          "MLH1":"216",
          "MSH2":"218",
          "MSH6":"219",
          "NBN":"158",
          "NCF1":"87",
          "NF1":"214",
          "PALB2":"308",
          "PMS2":"161",
          "PPP2R2A":"0000",
          "PTEN":"311",
          "PTEN_a":"0000",
          "RAD51B":"0000",
          "RAD51C":"314",
          "RAD51D":"516",
          "RAD54L":"1414",
          "STK11":"319",
          "TP53":"321"
          }

def modify_mane_select_key(mane_select_found,mane_select_id,mane_select_id_key,sorted_list):
    if mane_select_found:
        key_found=False
        for key in sorted_list:
            #print("key %s; value %s"%(key,sorted_mRNA_transcript_list[key]))
            mod_key=mane_select_id+"m"
            if key==mane_select_id or key==mod_key :
                key_found=True
                new_key=mane_select_id_key
                sorted_list[new_key] = sorted_list.pop(key)
                break
        if key_found:
            new_list=dict()
            new_list.setdefault(new_key,sorted_list[new_key])
            for key in sorted_list:
                #print("key %s; value %s"%(key,sorted_list[key]))
                if key != new_key:
                    new_list.setdefault(key,sorted_list[key])
            sorted_list=new_list
    return sorted_list
            
def make_transcriptconstants():
    global transcript_ids,ensembl_geneid,source_mapping,last_generange,LRG_id_dict,MANE_Select_dict
    global outj2,grchver

    
    if grchver=="GRCh37":
        Release="Ensembl Release 75 (Variation data updated April 2021)"
        
    elif grchver=="GRCh38":
        Release= "Ensembl Release 105 (Dec 2021)"
        
    else:
        Release= "Undefined"

    try:
        LRGid=LRG_id_dict[target_locus]
    except:
        LRGid="00000"

    sorted_mRNA_transcript_list=dict(sorted(mRNA_transcript_list.items()))
 
    sorted_mRNA_join_list=dict(sorted(mRNA_join_list.items()))
    sorted_CDS_join_list=dict(sorted(CDS_join_list.items()))
    
    try:
        mane_select=MANE_Select_dict[target_locus]
        mane_select_id="%s-%s"%(target_locus,mane_select)
        mane_select_id_key="%s(MANE_Select)"%mane_select_id
        mane_select_found=True
        version="v0.95"
    except:
        mane_select_found=False
        mane_select_id="undefined"
        mane_select_id_key="undefined"
        version="null"
        
    sorted_mRNA_transcript_list=modify_mane_select_key(mane_select_found,mane_select_id,mane_select_id_key,sorted_mRNA_transcript_list)
    
    if not outj2:

        sorted_mRNA_join_list=modify_mane_select_key(mane_select_found,mane_select_id,mane_select_id_key,sorted_mRNA_join_list)
        sorted_CDS_join_list=modify_mane_select_key(mane_select_found,mane_select_id,mane_select_id_key,sorted_CDS_join_list)

        transcriptconstants = {
            target_locus: {
                "Release":Release,
                "Retrieval_date":r_date,
                "Region":source_mapping,
                "Locus_range":last_generange,
                "is_join_complement":is_join_complement,
                "LRG_id":LRGid,
                "Ensembl_id":ensembl_geneid,
                "mRNA":sorted_mRNA_transcript_list,
                "mRNA_join":sorted_mRNA_join_list,
                "CDS_join":sorted_CDS_join_list,
            "MANE_Select":{
                "version":version,
                "mRNA":mane_select_id_key
                }
             }
            }
    else:

        transcriptconstants = {
            target_locus: {
                "Release":Release,
                "Retrieval_date":r_date,
                "Region":source_mapping,
                "Locus_range":last_generange,
                "is_join_complement":is_join_complement,
                "LRG_id":LRGid,
                "Ensembl_id":ensembl_geneid,
                "mRNA":sorted_mRNA_transcript_list,
            "MANE_Select":{
                "version":version,
                "mRNA":mane_select_id_key
                }
             }
            }
        
    config_out_data={
        "Reference_sequences":transcriptconstants
        }
    write_transcript_configs(target_locus,config_out_data)


def write_transcript_configs(locus,out_data):
    global last_generange
    config_file="%s_transcripts.json"%locus
    with RG_io.open_write(config_file) as write_file:
        json.dump(out_data, write_file,indent=4)
    print ('Output transcript ids file is: %s '%config_file)
    if last_generange =="undefined" or not found_locus:
        print('*** Looks like target locus %s was NOT FOUND ***'%target_locus)

def main(argv):
    global empty_string, inputfile,outtempfile,outreffile,outfeatfile,outvarfile,outnoseqfile,target_locus,ensembl_geneid,exclude_CDS
    global outj2,grchver
    global target_locus
    inputfile = empty_string
    outtempfile = empty_string
    outreffile =empty_string
    outfeatfile = empty_string
    outvarfile= empty_string
    outj2 = False

    outtempfile_postname="_temp.gb"
    outreffile_postname="_filtered.gb"
    outfeatfile_postname="_locseq.gb"
    outvarfile_postname="_filtervar.gb"
    outnoseqfile_postname="_noseq.gb"
    exclude_CDS=False # Default for -c command line switch for the exclusion or inclusion of the CDS feature

    opts, args = getopt.getopt(argv,"hacji:f:l:v:g:",["ifile=","ffile=","lfile=","vfile=","grch="])
    '''
    try:
        opts, args = getopt.getopt(argv,"hacij:f:l:v:",["ifile=","ffile=","lfile=","vfile="])
    except getopt.GetoptError:
        err_and_quit()
    '''
    for opt, arg in opts:
        if opt == '-h':
            print ('usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile> -c {CDS-suppress} -j {alternative json}'%Progver)
            #print ('usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile>'%Progver)
            sys.exit()
        elif opt == '-c':
            exclude_CDS=True #  Use this switch for the exclusion or inclusion of CDS
        elif opt == '-j':
            outj2= True
        elif opt in ("-i", "--ifile"):
            inputfile = arg
            splits=inputfile.split("_")
            target_locus=splits[0]
        elif opt == '-a':
            outreffile =target_locus+outreffile_postname
            outfeatfile =target_locus+outfeatfile_postname
            outvarfile =target_locus+outvarfile_postname
        elif opt in ("-f", "--ffile"):
            if arg:
                outreffile = arg
            else:
                outreffile =target_locus+outreffile_postname
        elif opt in ("-l", "--lfile"):
            if arg:
                outfeatfile = arg
            else:
                outfeatfile =target_locus+outfeatfile_postname
        elif opt in ("-v", "--vfile"):
            if not arg:
                outvarfile =target_locus+outvarfile_postname
            else:
                outvarfile = arg
                
    if inputfile == empty_string:
        err_and_quit()
    if (outreffile == empty_string) and (outfeatfile == empty_string) and (outvarfile == empty_string):
        outreffile =target_locus+outreffile_postname

    outtempfile=target_locus+outtempfile_postname

    outnoseqfile=target_locus+outnoseqfile_postname
            
    print ('Running: %s, version date %s'%(Progver,ProgverDate))
    print ('Target locus is %s'%target_locus)
    print ('exclude_CDS %s'%exclude_CDS)

    print ('Input file is: %s '%inputfile)
    print ('Output filtered file (features, variants and seq) file is: %s '%outreffile)
    print ('Output locus file is: %s '%outfeatfile)
    print ('Output variants file is: %s '%outvarfile)
    print ('Output noseq file is: %s '%outnoseqfile)

def err_and_quit():
    print ('error: usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile> -c {CDS-suppress} -j {alternative json}'%Progver)
    sys.exit(2)
    return

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

# ================================ 
#           Above here - old flat-file scan routines             
# ================================

# =================================================
# Sequence-input routines
# ================================================
global is_show_infilepath
is_show_infilepath=False

def clean_up():
    close_temp_file()
    delete_tempfile()
    is_append_tmp_file=False

def close_temp_file():
    global is_append_tmp_file,outtmphandle
    if is_append_tmp_file:
        outtmphandle.close()
    return

def delete_tempfile():
    global is_append_tmp_file,outtmphandle
    if is_append_tmp_file:
        outtmphandle.close()
    print(" Deleting temporary file %s - you may have to authorise"%outtempfile)
    os.system("/bin/rm %s"%outtempfile) 
    return

def update_temp_file(instring):
    global outtempfile,is_append_tmp_file,outtmphandle
    if outtempfile !=empty_string:
        if not is_append_tmp_file:
            outtmphandle = open(outtempfile,"w")
            is_append_tmp_file=True
            outtmphandle.write(instring)
        else:
            outtmphandle.write(instring)
    return

def read_as_text():
    with open(inputfile) as f:
        for line in f:
            if Locus in line:
                #print("Line: %s"%line)
                splits=line.split()
                splits[0]=Locus
                splits[2]='0' #Set sequence length to zero to overcome Biopython bug with variants that are replace="-/N"
                line=' '.join(splits)+"\n"
                #print("Line: %s"%line)
            update_temp_file(line)
    warning="/Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/Bio/GenBank/__init__.py:841: BiopythonParserWarning: Expected sequence length 0,"
    print("Expect to get warning from Biopython like:\n***%s***"%warning)
    print("This is because have deliberately created a temporary input file with sequence length set to zero")
    close_temp_file()

def update_journal(instring):
    print("%s"%instring)

# Use this instead of the one below to allow for reading of gzip files when using Python
# This is **not** suitable for webapp repository because there is no support for gzip files
def read_single_record_inputgz(infile,seqform):
    seq_record,success=RG_io.read_single_record_input_gz_support(infile,seqform)
    return seq_record,success

# This is safe to use in Python or for webapp, but assumes no input files are gz
def read_single_record_input(infile,seqform):
    with RG_io.open_read(infile) as f:
        seq_record = SeqIO.read(f,seqform)
        seq_record,success = RG_process.add_extra_objects(seq_record)
    return seq_record,success
# end of read_single_record_input(infile,seqform)

def read_seqrecord(infile,embl_or_genbank):
    SEQ_record=""
    exists = RG_io.is_file(infile)
    if exists:
        update_journal("\nReading file %s"%infile)
        # second returned parameter "exists" is really a "success" for the file-read, but want to overwrite "exists" to return an overall success/fail condition here
        SEQ_record,exists = read_single_record_inputgz(infile,embl_or_genbank)
        if exists:
            MaxVarPos=len(SEQ_record.seq)
            SEQ_record.howmany="all"
            ''' By definition Refseq CIGAR is all Ms'''
            SEQ_record.mutlabel=""
            SEQ_record.cigar=str(MaxVarPos)+"M"
            SEQ_record.cigarbox=[0, int(MaxVarPos), 'M', int(MaxVarPos)]
            SEQ_record.mutbox=RG_process.haveamutbox(SEQ_record.cigarbox)
            # Copying in a modified version of RG_globals.bio_parameters["target_build_variant"]["exonplus_lookup"] which will be inherited by each variant entry
            # This is later modified in each variant haplotype to create a lookup table for creating paired-end exomic reads 
            SEQ_record.exonplus_lookup=[]
            #print("SEQ_record.exonplus_lookup[:10]%s\n[-10:]%s "%(SEQ_record.exonplus_lookup[:10],SEQ_record.exonplus_lookup[-10:]))
            SEQ_record.Headclip=0
            SEQ_record.Tailclip=0
            SEQ_record.splicecount=0
            SEQ_record.endclipcount=0
            SEQ_record.is_ref_spliced=False
            SEQ_record.unclipped_length=MaxVarPos
            SEQ_record.clipped_length=MaxVarPos
            SEQ_record.spliced_length=MaxVarPos
            SEQ_record.endlocus=len(SEQ_record.seq)
            SEQ_record.rev_endlocus=0
            SEQ_record.Generated_Fragcount=0
            SEQ_record.Saved_Fragcount=0
            SEQ_record.Saved_Unpaired_Fragcount=0
            SEQ_record.Saved_Paired_Fragcount=0
            
            update_journal("Success reading %s"%infile)
        else:
            update_journal("Failure reading ACCESSION line for %s"%infile)
    else:
        update_journal(" \n*** %s file not found ***"%infile)
        exists=False
    return SEQ_record,exists
# end of read_seqrecord(infilename,embl_or_genbank,out_head)


# ================================ 
#           Sequence file output           
# ================================
    from RG_exploder_globals import Seq_Format,Seq_IO_file_ext # CONSTANTS: defines Genbank or Embl as file format; path to refseq file
    global Seq_Format,Seq_IO_file_ext




def write_gb_features2(SeqRec,varfeats,out_file,readmetxt,style):
    # Derived from RG_exploder_main.py.write_gb_features 
    Outfilepath=""
    print("out_file %s; style %s"%(out_file,style))

    # Several fields in SeqRec.annotations seem to be ignored by SeqIO.write, so 'date' always comes out as '01-JAN-1980'
    # I have tried!
    # Any COMMENT fields, likewise.
    # "full" and "ref" contain the sequence; "ref" omits variants

    if style=="ref": #  Without variants
        CopySeqRec=copy.copy(SeqRec)
    elif style == "full": # Full sequence record including variants
        CopySeqRec=copy.copy(SeqRec)
        for item in varfeats:
            CopySeqRec.features.append(item) 
    elif style == "noseq" or style == "vars":
        #Take clipped sequence and only the base annotations into a new sequence record
        CopySeqRec=RG_process.modify_seq_in_record(SeqRec.seq[0:0],SeqRec)
        CopySeqRec.features=RG_process.get_source_gene_features(SeqRec,target_locus) # Overwrites all features with basic
        if style == "vars": #  With variants
            for item in varfeats:
                CopySeqRec.features.append(item)
    else:
        print("Invalid parameter %s for write_gb_features2"%style)
        exit()
             
    out_gbfile="%s%s"%(Outfilepath,out_file)
    update_journal(" Writing %s"%out_file)
    gbout = RG_io.open_write(out_gbfile)
    out_handle=StringIO()
    SeqIO.write(CopySeqRec,out_handle,"genbank")
    out_data = out_handle.getvalue()
    gbout.write(out_data)
    out_handle.close()
    gbout.close()
    return
# end of write_gb_features2(SeqRec,out_file,readmetxt,style)

# ================================ 
#           This is it             
# ================================

if __name__ == "__main__":
   main(sys.argv[1:])

def get_joinstring(in_feature_location):
    location_string=str(in_feature_location).strip("join").strip("{}")       
    #print("location_string %s"%(location_string))
    splits=location_string.split(",")
    #print("splits %s"%splits)
    out_joinstring=""
    for item in splits:
        duff = "".join([i for i in item if i not in "()+-[]"])
        begin,end=duff.split(":")
        ibegin=int(begin)+1
        #print("duff%s"%duff)
        out_joinstring+=str(ibegin)+":"+end+","
    return out_joinstring[:-1]

def modify_annotations(annotations):
    unwanted_annotations=['keywords','taxonomy']
    new_annotations={}
    for (key) in annotations:
        #print("key %s, CopySeqRec.annotations.get(key) %s"%(key,CopySeqRec.annotations.get(key)))
        if any(x in key for x in unwanted_annotations):
            pass
        elif key=='comment':
            comments=str(annotations.get(key))
            newcomment=comments[:comments.find(")")+1]+", then filtered using "+Progver+" on "+ this_date +"."
            new_annotations.update({key:newcomment})
        else:
            #print("passing: key %s, CopySeqRec.annotations.get(key) %s"%(key,CopySeqRec.annotations.get(key)))
            new_annotations.update({key:annotations.get(key)})         
    return new_annotations


def embl_filter_record(SeqRec,locus_tag):
    CopySeqRec=copy.deepcopy(SeqRec)
    CopySeqRec.features=filter_all_features(SeqRec,locus_tag,True)# Overwrite the features; True retains variants
    CopySeqRec.annotations =modify_annotations(CopySeqRec.annotations)         
    return CopySeqRec

def embl_filter_record_novars(SeqRec,locus_tag):
    CopySeqRec=copy.deepcopy(SeqRec)
    CopySeqRec.features=filter_all_features(SeqRec,locus_tag,False)# Overwrite the features; False omits variants
    CopySeqRec.annotations =modify_annotations(CopySeqRec.annotations)  
    return CopySeqRec


def filter_all_features(seq_record,target_loc,with_var):
    wanted_features=['source','gene','mRNA','CDS']
    wanted_CDS_db_xref=['CCDS:']
    '''
    a) For each feature, check it's one we want to keep
    b) For each item in the variants table, list positions of wanted ones
    '''
    #print("At filter_varfeature_index")
    #print("seq_record.id %s"%seq_record.id)
    #print("seq_record.features %s"%seq_record.features)

    found_locus=False
    features=[]
    
    out_index_count=0
    #target_locus_match="'"+RG_globals.target_locus+"'"
    target_locus_match="'"+target_loc+"'"
    
    for (index, feature) in enumerate(seq_record.features):
        feature=seq_record.features[index]
        if any(x in feature.type for x in wanted_features): 
            if feature.type == 'source':
                features.append(feature); out_index_count+=1
            elif feature.type == 'gene':
                for key in feature.qualifiers:
                    if key=='locus_tag': # Matching the feature groupings to the desired locus, otherwise
                        locus_txt=str(feature.qualifiers.get(key))
                        if target_locus_match in locus_txt:
                            #print("feature.type %s, key %s, value %s; target_locus_match %s"%(feature.type,key,locus_txt,target_locus_match))
                            found_locus=True
                        else:
                            found_locus=False
                if found_locus:
                    #print("feature.type %s: feature %s"%(feature.type,feature))
                    features.append(feature); out_index_count+=1
                # print("seq_record.features[index]%s "%seq_record.features[index])
            elif found_locus:
                if feature.type == 'mRNA':
                    features.append(feature); out_index_count+=1
                    #print("feature.type %s; feature %s"%(feature.type,feature))
                else: #feature.type == 'CDS':
                    new_feature=copy.copy(feature)
                    new_CDS_qualifiers= {}   
                    qualifiers=str(feature.qualifiers)
                    #print("feature.type: %s; qualifiers %s"%(feature.type,qualifiers))                    
                    for key in feature.qualifiers:
                        if key !='translation':
                            if key== 'db_xref':
                                new_xrefs=[]
                                for item in feature.qualifiers.get(key):                                
                                    if any(x in item for x in wanted_CDS_db_xref):
                                        new_xrefs.append(item)
                                if len(new_xrefs) > 0:
                                    new_CDS_qualifiers.update({key:new_xrefs}) 
                            else:
                                new_CDS_qualifiers.update({key:feature.qualifiers.get(key)})  
                    new_feature.qualifiers=new_CDS_qualifiers
                    features.append(new_feature); out_index_count+=1 
    if with_var:
        varindex=RG_process.get_varfeature_index(seq_record)# Comes back sorted: highest position first
        for (index) in varindex:
            features.append(seq_record.features[index])
    #return features,source_feature_index,gene_feature_index,mRNA_feature_index,CDS_feature_index
    return features
# end of filter_all_features(seq_record,with_var)


def make_files():
    global mRNA_transcript_list,mRNA_join_list,CDS_join_list,grchver,source_mapping,last_generange,is_join_complement,ensembl_geneid,r_date,empty_string,found_locus
    global outtempfile
    read_as_text() # Creates outtempfile
    #SEQ_record,exists = read_seqrecord(inputfile,"genbank")
    SEQ_record,exists = read_seqrecord(outtempfile,"genbank")
    found_locus=False
    is_join_complement=False
    ensembl_geneid="undefined"
    match_id="no"
    if exists:
        wanted_features=['source','gene','mRNA','CDS']
        #wanted_transcript_qualifiers=['gene','note','standard_name']
        target_locus_match="'"+target_locus+"'"
        print("target_locus:%s"%target_locus)
        NewRec=embl_filter_record_novars(SEQ_record,target_locus) # Clean it up to this target_locus
        varfeatures=RG_process.get_varfeatures(SEQ_record) # Features in featstart highest to lowest order because calls order_featindex(varfeat_index,seq_record,True) - reverses
        varfeatures=RG_process.order_featstart(varfeatures,False)  # Puts var features in featstart order (lowest to highest)
        
        #print("NewRec: %s"%(NewRec))
        #print("NewRec.id %s"%(NewRec.id))
        ##print("NewRec.Description: %s"%(NewRec.description))
        #print("NewRec.numFeatures: %s"%(NewRec.number_of_features))
        #print("NewRec.annotations: %s"%(NewRec.annotations))
        #print("NewRec.id %s"%(NewRec.id))
        #print("NewRec.annotations.get('date'): %s"%(NewRec.annotations.get('date')))

        r_date=NewRec.annotations.get('date')

        splits=NewRec.id.split(":")
    
        NewRec.GRChver=splits[1]
        NewRec.chrom=splits[2]
        NewRec.abs_start=splits[3]
        NewRec.abs_end=splits[4]
        NewRec.polarity=splits[5]

        grchver=splits[1]
        source_mapping="%s:%s:%s:%s:%s"%(splits[1],splits[2],splits[3],splits[4],splits[5])
        #print("Region: %s"%source_mapping)
        #print("grchver:%s"%grchver)
                    
        for (index, feature) in enumerate(NewRec.features):
            feature=NewRec.features[index]
            if any(x in feature.type for x in wanted_features): 
                if feature.type == 'source':
                    pass
                    #print("feature %s"%feature)
                    #ensembl_geneid 
                elif feature.type == 'gene':
                    for key in feature.qualifiers:
                        if key=='gene':
                            match_id=str(feature.qualifiers.get(key)).strip("[]'")
                        if key=='locus_tag': # Matching the feature groupings to the desired locus, otherwise
                            locus_txt=str(feature.qualifiers.get(key))
                            print("locus_txt %s"%locus_txt)
                            if target_locus_match in locus_txt:
                                found_locus=True
                                ensembl_geneid=match_id
                                print("ensembl_geneid %s"%ensembl_geneid)
                                # Now add location
                                locus_begin,locus_end=str(feature.location).strip("[]").split(":")
                                locus_end,pol=locus_end.split("]")
                                begin=int(locus_begin)+1
                                if pol =="(-)":
                                    is_join_complement=True
                                else:
                                    is_join_complement=False
                                last_generange=str(begin)+":"+locus_end
                                #print("locus_begin %s,locus_end %s, pol %s,complement %s"%(locus_begin,locus_end,pol,complement))
                            else:
                                found_locus=False
                elif found_locus :
                    if feature.type == 'mRNA':
                        #pass
                        #print("feature: %s"%(feature))
                        #print("feature.location: %s"%(feature.location))
                        #print("feature.type: %s"%(feature.type))
                        #print("gene: %s"%feature.qualifiers.get('gene'))
                        mRNA_gene_name=str(feature.qualifiers.get('gene')).strip("[]'")
                        #print("mRNA_gene_name: %s"%mRNA_gene_name)
                        #print("standard_name: %s"%feature.qualifiers.get('standard_name'))
                        if mRNA_gene_name==match_id:
                            mRNA_transcript_name=str(feature.qualifiers.get('standard_name')).strip("[]'")
                            #print("standard_name: %s"%mRNA_transcript_name)
                            #mRNA_transcript_list.append(mRNA_transcript_name)
                            tidname=mRNA_transcript_name.split('.')[0]
                            truncated_transcript_id=tidname.split('ENST')[1].lstrip("0")

                            mRNA_key=target_locus+"-"+truncated_transcript_id
                            mRNA_transcript_list.update({mRNA_key:mRNA_transcript_name})
                
                            #print("feature.location: %s"%(feature.location))

                            mRNA_joinstring=get_joinstring(feature.location)
                            mRNA_join_list.update({mRNA_key:mRNA_joinstring})
 
                    elif feature.type == 'CDS' and not exclude_CDS: 
                        #print("feature.type: %s"%(feature.type))
                        #print("gene: %s"%feature.qualifiers.get('gene'))
                        CDS_gene_name=str(feature.qualifiers.get('gene')).strip("[]'")
                        #print("CDS_gene_name: %s"%CDS_gene_name)
                        #print("note: %s"%feature.qualifiers.get('note'))
                        if CDS_gene_name==match_id:
                            CDS_note=str(feature.qualifiers.get('note')).strip("[]'")
                            tidtxt,CDS_transcript_name=CDS_note.split("=")
                            truncated_CDS_transcript_name=CDS_transcript_name.split('.')[0].split('ENST')[1].lstrip("0")
                            #print("truncated_transcript_id %s"%truncated_transcript_id)
                            #print("truncated_CDS_transcript_name: %s"%truncated_CDS_transcript_name)
                            if truncated_transcript_id == truncated_CDS_transcript_name:
                                CDS_joinstring=get_joinstring(feature.location)
                                CDS_join_list.update({mRNA_key:CDS_joinstring})
            elif feature.type !="variation":
                print("Unexpected feature %s"%feature)
        
        #print("mRNA_transcript_list:%s"%mRNA_transcript_list)
        #print("mRNA_join_list:%s"%mRNA_join_list)
        #if not exclude_CDS:
        #    print("CDS_join_list:%s"%mRNA_join_list)
        make_transcriptconstants()
        if found_locus:
            if outfeatfile != empty_string: # "_locseq.gb" 
                write_gb_features2(NewRec,varfeatures,outfeatfile,"","ref") # Don't understand why, if this goes after "full", it includes variants
            if outreffile != empty_string: # "_filtered.gb"
                write_gb_features2(NewRec,varfeatures,outreffile,"","full")
            if outnoseqfile != empty_string:
                write_gb_features2(NewRec,varfeatures,outnoseqfile,"","noseq")
            if outvarfile != empty_string: # "_filtervar.gb"
                write_gb_features2(NewRec,varfeatures,outvarfile,"","vars")

    else:
        exit()

read_config()
make_files()
clean_up()
#read_as_text()



    

