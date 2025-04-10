#!/usr/local/bin/python3
Progver="embl_feature_filter8.py"
ProgverDate="21-Feb-2024"
'''

This reads in Genbank format files as part of the pre-processing required to run the RG_exploder_*.py suite.
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
This processes the {locus}_Ensembl_download.gz file to eliminate unwanted items from the feature table
creating, optionally
a) A file with the target locus, mRNA features and sequence. {locus}_locseq.gb - this is the reference sequence for main application 
b) A file with the target locus, variant features, omits sequence. {locus}_var.gb
c) A file with the target locus, mRNA features, variant features and sequence. {locus}_filtered.gb

And constitutively:
d) A file with the target locus, no variants, no sequence {locus}_noseq.gb - this is the hap0 file for the main application
e) A file with the transcript ids in json format   {locus}_transcripts.json - contributor to config.json file

Unwanted only in that these features (primarily the database cross-reference
definitions, but also CDS, misc_RNA and exons) are unnecessary for the Replicon
variants-processing applications.
The prime interests are:
   the mRNA join features to identify specific transcripts.
   the "variations" features, with dbsnp ids


####### NB: the command-line option -c to  flip between retaining mRNA join(n..m) ... and mRNA join(complement(n..m) is deprecated ######
            because it is automatically selected by the setjoin() function
####### NB: the command-line option -c now supresses CDS inclusion######


'''
import sys, getopt
import binascii  #needed by is_gz_file
import gzip
from functools import partial
import json
from datetime import date
import RG_exploder_io as RG_io

global empty_string,inputfile,outfeatfile,is_append_featfile,outvarfile,is_append_varfile, is_append_noseqfile,nocomment
global outreffile,is_append_ref_file,r_date,mRNA_transcript_list,CDS_transcript_list,CDS_protein_list,first_transcript,is_first_ts
global mrna_store,is_mrna_store,mRNA_join_list,mrna_feat,mrna_feat_line,cds_feat,cds_feat_line,codon_start,cds_store,is_cds_store,last_generange
is_append_featfile=False
is_append_varfile=False
is_append_ref_file=False
is_append_noseqfile=False

last_generange="undefined"
mrna_store=""
is_mrna_store=False
empty_string=''
first_transcript="undefined"
is_first_ts=False
mrna_feat="mRNA   "
mrna_feat_line=""
mRNA_transcript_list=dict()
mRNA_join_list=dict()
CDS_transcript_list=dict()
CDS_protein_list=dict()
CDS_join_list=dict()
cds_store=""
is_cds_store=False
cds_feat="CDS   "
cds_feat_line=""
codon_start="codon_start"
dots=".."
carat="^"

today = date.today()
# Textual month, day and year	
r_date = today.strftime("%Y %B %d")
#print("r_date =%s"%r_date)

'''
unwanted=['xref="CCDS','xref="RefSeq','xref="Uni','xref="UCSC','xref="HPA','db_xref="HGNC',
          'db_xref="EMBL','xref="GO','db_xref="protein','xref="KEGG','xref="ChEMBL','xref="PDB',
          'xref="Reactome','xref="Vega','xref="goslim','xref="OTTT']
'''

unwanted=['xref="RefSeq_mRNA_predicted:','xref="RefSeq_peptide_predicted:','xref="Uni','xref="UCSC','xref="HPA','db_xref="HGNC',
          'db_xref="EMBL','xref="GO','db_xref="protein','xref="KEGG','xref="ChEMBL','xref="PDB',
          'xref="Reactome','xref="Vega','xref="goslim','xref="OTTT','xref="BioGRID','xref="RNAcentral',
          'db_xref="RefSeq_peptide:','db_xref="RefSeq_ncRNA:','db_xref="RefSeq_ncRNA_predicted:'
          ]           #'xref="RefSeq_mRNA:',

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
def save_generange():
    global generange
    generange=line.split(geneline)
    generange=generange[1]

def check_mrnastorelines():
    global mrna_store,is_mrna_store,mrna_feat,mrna_feat_line
    if is_mrna_store:
        #print("mrna_store %s"%(mrna_store))
        splitline=mrna_store.split(mrna_feat)
        mrna_feat_line=splitline[1]
        mrna_feat_line=mrna_feat_line.replace(' ', '')
        mrna_store=""
    else:
        mrna_feat_line="wot?"
    is_mrna_store=False

def check_cdsstorelines():
    global cds_store,is_cds_store,cds_feat,cds_feat_line
    if is_cds_store:
        #print("cds_store %s"%(cds_store))
        splitline=cds_store.split(cds_feat)
        cds_feat_line=splitline[1]
        cds_feat_line=cds_feat_line.replace(' ', '')
        cds_store=""
    else:
        cds_feat_line="wot?"
    is_cds_store=False


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

def modify_joinlist(sorted_list):
    # Edit out all the join information and create list of lists
    for key in sorted_list:
        edited_item=''.join( c for c in sorted_list[key] if  c not in 'joincomplet()' )
        edited_item=edited_item.replace("..", ":")
        sorted_list.update({key:edited_item})
    return sorted_list
            
def make_transcriptconstants():
    global transcript_ids,ensembl_geneid,source_mapping,first_transcript,last_generange,LRG_id_dict,MANE_Select_dict
    global outj2,grchver

    splits=source_mapping.split(":")
    grchver=splits[0]

    print("grchver %s"%grchver)
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

    if "complement" in last_generange:
        is_join_complement=True
        edited_generange=''.join( c for c in last_generange if  c not in 'complent()' )
        last_generange=edited_generange
    else:
        is_join_complement=False
    
    last_generange = last_generange.replace("..", ":")
    
    if not outj2:

        sorted_mRNA_join_list=modify_mane_select_key(mane_select_found,mane_select_id,mane_select_id_key,sorted_mRNA_join_list)
        sorted_CDS_join_list=modify_mane_select_key(mane_select_found,mane_select_id,mane_select_id_key,sorted_CDS_join_list)

        sorted_mRNA_join_list=modify_joinlist(sorted_mRNA_join_list)
        #print("sorted_mRNA_join_list %s"%(sorted_mRNA_join_list))
        sorted_CDS_join_list=modify_joinlist(sorted_CDS_join_list)


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
    if last_generange =="undefined":
        print('*** Looks like target locus %s was NOT FOUND ***'%target_locus)


def main(argv):
    global empty_string, inputfile,outreffile,outfeatfile,outvarfile,outnoseqfile,target_locus,ensembl_geneid,exclude_complement,exclude_CDS
    global outj2,grchver
    inputfile = empty_string
    outreffile =empty_string
    outfeatfile = empty_string
    outvarfile= empty_string
    outj2 = False
    grchver = empty_string

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
            usage_message("")
            #print ('usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile> -g <grchver:37/38> -c {CDS-suppress} -j {alternative json}'%Progver)
            #print ('usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile>'%Progver)
            sys.exit()
        # Former command line switch for complement, now automatically selected by setjoin()
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
        '''
        elif opt in ("-g", "--grch"):
            if not arg:
                grchver ="undef"
            else:
                grchver = arg
        '''

                
    if inputfile == empty_string:
        err_and_quit()
    if (outreffile == empty_string) and (outfeatfile == empty_string) and (outvarfile == empty_string):
        outreffile =target_locus+outreffile_postname

    outnoseqfile=target_locus+outnoseqfile_postname
            
    print ('Running: %s, version date %s'%(Progver,ProgverDate))
    print ('Target locus is %s'%target_locus)
    print ('exclude_CDS %s'%exclude_CDS)

    print ('Input file is: %s '%inputfile)
    print ('Output filtered file (features, variants and seq) file is: %s '%outreffile)
    print ('Output locus file is: %s '%outfeatfile)
    print ('Output variants file is: %s '%outvarfile)
    print ('Output noseq file is: %s '%outnoseqfile)


def usage_message(addtxt):
    print ('%susage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile> -c {CDS-suppress} -j {alternative json}'%(addtxt,Progver))
    
def err_and_quit():
    usage_message("error: ")
    #print ('error: usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile> -c {CDS-suppress} -j {alternative json}'%Progver)
    #print ('usage: %s -i <inputfile> -a {all} -f <filterfile> -l <locusfile> -v <outvarfile> -g <grchver:37/38> -c {CDS-suppress} -j {alternative json}'%Progver)

    sys.exit(2)
    return

def store_featfile(instring):
    '''
    '''
    global empty_string,outfeatfile,is_append_featfile,outfeathandle
    if outfeatfile !=empty_string:
        if not is_append_featfile:
            outfeathandle = open(outfeatfile,"w")
            is_append_featfile=True
        outfeathandle.write(instring)
    update_ref_file(instring)
    return

def update_ref_file(instring):
    '''
    '''
    global outreffile,is_append_ref_file,outrefhandle
    if outreffile !=empty_string:
        if not is_append_ref_file:
            outrefhandle = open(outreffile,"w")
            is_append_ref_file=True
        outrefhandle.write(instring)
    return

def update_featfile(instring):
    '''
    '''
    global empty_string,outfeatfile,is_append_featfile,outfeathandle,mrna_store,is_mrna_store,cds_store,is_cds_store
    instring=test_for_comment(instring)
    if outfeatfile !=empty_string:
        if not is_append_featfile:
            outfeathandle = open(outfeatfile,"w")
            is_append_featfile=True
        is_keep_line=True
        for (item) in unwanted:
            if (item in instring):
                is_keep_line=False
                break
        if is_keep_line:
            outfeathandle.write(instring)
            update_varfile(instring)
    else:
        update_varfile(instring)
        
    if is_mrna_store:
                mrna_store+=line[:-1]
    if is_cds_store:
                cds_store+=line[:-1]
    return

def update_varfile1(instring):
    '''
    '''
    global empty_string,outvarfile,is_append_varfile,outvarhandle,ignore_var,nocomment
    instring1=instring
    if not ignore_var and (outvarfile != empty_string):
        if not is_append_varfile:
            outvarhandle = open(outvarfile,"w")
            is_append_varfile=True
        if Locus in instring:
            instring=get_newlocline(instring)
        elif Origin in instring:
            nocomment=False
            ignore_var=True
            instring=instring+"//\n"
        outvarhandle.write(instring)
        update_noseqfile(instring)
    update_ref_file(instring1) 
    return

def test_for_comment(instring):
    global is_comment,source_mapping
    if comment in instring or organism in instring or keywords in instring:
        is_comment=True
        instring=""        
    elif features in instring:
        is_comment=False
    elif is_comment:
        ''' This switches off sending SOURCE & organism content in instring  '''
        instring=""
    elif Accession in instring:
        colonpos=instring.find(':')
        source_mapping=instring[colonpos+1:]
        source_mapping=source_mapping.rstrip()
    return(instring)
  
def update_varfile(instring):
    '''
    '''
    global empty_string,outvarfile,is_append_varfile,outvarhandle,ignore_var,nocomment
    instring1=instring
    if not ignore_var and (outvarfile != empty_string):
        if not is_append_varfile:
            outvarhandle = open(outvarfile,"w")
            is_append_varfile=True
        if Locus in instring:
            instring=get_newlocline(instring)
        elif Origin in instring:
            nocomment=False
            ignore_var=True
            instring=instring+"//\n"
        else:
            instring=test_for_comment(instring)
        outvarhandle.write(instring)
        update_noseqfile(instring)
    update_ref_file(instring1) 
    return


def update_noseqfile(instring):
    '''
    This is just a cut-down version of update_varfile
    '''
    global empty_string,outnoseqfile,is_append_noseqfile,outnoseqhandle
    if outnoseqfile != empty_string:
        if not is_append_noseqfile:
            outnoseqhandle = open(outnoseqfile,"w")
            is_append_noseqfile=True
        if not is_variants or (Origin in instring):
            outnoseqhandle.write(instring) 
    return

def get_newlocline(instring):
    global Locus,r_date
    '''
    instring looks like
    LOCUS       12 46148 bp DNA HTG 9-OCT-2018

    want to set the bp value to 0 by returning the line

    LOCUS       12 0 bp DNA HTG 9-OCT-2018
    '''
    listpos=0
    slist=instring.split()
    r_date=slist[-1]
    for string in slist:
        if string=='bp':
            slist[listpos-1]='0'
            break
        listpos+=1
    newstring=Locus
    slist.pop(0)
    for string in slist:
        newstring=newstring+string+" "
    newstring=newstring[:-1]
    newstring+="\n"
    return newstring

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def close_outfiles():
    global is_append_ref_file,outrefhandle,is_append_featfile,outfeathandle,is_append_varfile,outvarhandle
    make_transcriptconstants()
    if is_append_ref_file:
        outrefhandle.close()
    if is_append_varfile:
        outvarhandle.close()
    if is_append_featfile:
        outfeathandle.close()
    return


def test_complement():
    global exclude_complement,joincomplement,complement,line,ignore,exclude_joinAC,joinAC,ignorest,joincomplementAC,complementAC
    global reading_join
    
    if exclude_complement:         # exclude_complement=True
        #if joincomplement in line:
        if complement in line:
            ignore=True            # Screens out join(complement and join(complement(AC
        elif (exclude_joinAC and joinAC in line): 
            ignore=True            # Screens out join(AC if exclude_joinAC ==True
            reading_join=False
        else: 
            ignore=False           # Effectively allows join( or join(AC if exclude_joinAC==False
            reading_join=True
            if not ignorest:
                update_featfile(line)
    else:                          # exclude_complement=False
        #if joincomplement in line: # We are looking to keep join(complement
        if complement in line: # We are looking to keep complement
            #if (exclude_joinAC and joincomplementAC in line):
            if (exclude_joinAC and complementAC in line):
                ignore=True            # Screens out join(complement(AC if exclude_joinAC ==True
                reading_join=False
            else:
                ignore=False           # Effectively allows join(complement or join(complement(AC if exclude_joinAC==False
                reading_join=True
                if not ignorest:
                    update_featfile(line)
        else:
            ignore=True
    return


def setjoin():
    global exclude_complement
    # exclude_complement is a switch for the exclusion or inclusion of any "join(complement" lines'''
    # False means join(complement is retained
    # True means the reverse: join(complement is rejected
    # see also exclude_joinAC
    
    if complement in storelines:
        exclude_complement=False
    else:
        exclude_complement=True
    print("exclude_complement %s"%exclude_complement)

def test_variants():
    global line,is_CDS,ignore_var,is_variants,finish,variation,dots,carat
    if variation in line:
        is_CDS=False
        ignore_var=False
        is_variants=True
        # modify line if start is higher than end - this will be an insert
        # change format in line from eg: "variation       166..165" to "variation       165^165"
        # This appears to be BioPython regarding this as 'non-standard' Genbank format, and trying to
        # trap inconsistencies in linear vs circular DNA. Somehow the high..low format gets caight in this trap
        # Not an issue in Biopython 1.80, but manifested in Biopython 1.83. This is intended for future protection.
        if dots in line:
            splitline=line.split()
            values=splitline[1].split(dots)
            #print("splitline %s | values %s"%(splitline,values))
            if values[0] > values[1]:
                newline="%s%s%s%s\n"%(variation,values[1],carat,values[0])
                update_varfile(newline)
            else:
                update_varfile(line)  
        else:
            update_varfile(line)
    elif Origin in line:
        is_CDS=False
        finish=True
        update_featfile(line)
    elif is_variants:
        update_varfile(line)

# ================================ 
#           This is it             
# ================================

if __name__ == "__main__":
   main(sys.argv[1:])

start=False
firstgene=False
geneid_found=False
foundloc=False
store=False
ignore=False
ignorest=False
finish=False
skip_to_variants=False
ignore_var=False
nocomment=False
global is_comment
is_comment=False
is_variants=False
is_CDS=False
is_exon=False
in_translate=False
exclude_translation=True # ''' Use this switch for the exclusion or inclusion of the translation in CDS feature'''
exclude_joinAC=True # ''' Use this switch for the exclusion or inclusion of any "join(AC" or join(complement(AC lines'''
#exclude_joinAC=False # ''' This kills '''
exclude_complement=True # Reset elsewhere

reading_join=False

ensembl_geneid="**** EMPTY ****" 

startline='/db_xref="taxon'
geneline="gene   "
generange="0..0"
genetag="/gene="
gene_equal_quote='/gene="'
GRChid="unset"
locus_tag="/locus_tag="
storelines=""
variation="     variation       "
translation="/translation="
db_xref="/db_xref="
misc_rna="misc_RNA   "
exon="exon   "
keywords="KEYWORDS"
features="FEATURES   "
comment="COMMENT     "
organism="  ORGANISM"
Locus="LOCUS       "
Origin="ORIGIN"
Accession="ACCESSION"
complement="complement"
joincomplement="join(complement"
joinAC="join(AC"
joincomplementAC="join(complement(AC"
complementAC="complement(AC"

standard_name='/standard_name='
note_transcript_id='/note="transcript_id='


'''splits=inputfile.split("_")
target_locus=splits[0]'''


# Read in matching cvs if present
enstfile="./%s_ENSTID.csv"%target_locus
enstids= {}
is_enstids=False
'''
if RG_io.is_file(enstfile):
     with RG_io.open_read(enstfile) as idfile:
         for line in idfile:
             splitline=line.split(',')
             enstids[splitline[1]]=splitline[0]
     is_enstids=True
'''
is_gz = is_gz_file(inputfile)
if not is_gz:
    _open = open
    _openmode='r'
elif is_gz:
    _open = partial(gzip.open, mode='rt')
    _openmode='rb'
else:
    program_exit('Unknown file encoding at %s'%inputfile)

with _open(inputfile) as f:
    
    for line in f:
        if finish:
            update_featfile(line)
        elif foundloc:
            if locus_tag in line: # It's a non-matching one after matching previously found, so skip to variants-processing
                #print(line)
                ignore=True
                ignorest=True
                is_CDS=False
                in_translate=False
                skip_to_variants=True
            elif skip_to_variants:
                test_variants()
            elif geneline in line:
                is_CDS=False
                ignore=True
                in_translate=False
            elif cds_feat in line:
                ignore_var=True
                in_translate=False
                if not exclude_CDS:
                    is_CDS=True
                    is_cds_store=True
                    test_complement()
                    if ignore:
                        is_CDS=False
                else:
                    ignore=True
            elif exon in line:
                is_CDS=False
                ignore=True
                in_translate=False
            elif misc_rna in line:
                is_CDS=False
                ignore=True
                in_translate=False
            elif mrna_feat in line:
                is_mrna_store=True
                is_CDS=False
                ignore_var=True
                in_translate=False
                test_complement()
            elif codon_start in line and is_CDS: # Rare occurrence where join(complement ... followed by /codon_start gets included
                line=''
            elif gene_equal_quote in line and not ignore : # Have to ignore gene id in misc_rna entries
                # Catch the gene_id in a mRNA or CDS definition
                check_mrnastorelines()
                check_cdsstorelines()
                update_featfile(line)
     
            elif (standard_name in line) and not ignore: # Have to ignore in misc_rna entries
                # We have a transcript id in a mRNA definition
                splitline=line.split(standard_name)
                #print("splitline %s"%splitline)
                splitline=splitline[1].split('"')
                #print("splitline %s"%splitline)
                transcript_id=splitline[1]
                tidname="unset"
                if is_enstids:
                    for key,value in enstids.items():
                        if key==transcript_id:
                            tidname=value
                            break
                else:
                    tidname=transcript_id.split('.')[0]
                    #print("line %s; standard_name %s; transcript_id %s; tidname %s"%(line,standard_name,transcript_id,tidname))
                    truncated_transcript_id=tidname.split('ENST')[1]
                    tidname=truncated_transcript_id.lstrip("0")
                
                tkey2="%s-%s"%(target_locus,tidname)
                #print("tkey %s,transcript_id %s"%(tkey,transcript_id))
                if not is_first_ts:
                    first_transcript=tkey2
                mRNA_transcript_list.setdefault(tkey2,transcript_id)
                mRNA_join_list.setdefault(tkey2,mrna_feat_line)
                update_featfile(line)

            elif (note_transcript_id in line) and not ignore: # Have to ignore in misc_rna entries
                # We have a transcript id in a CDS definition
                splitline=line.split(note_transcript_id)
                #print("splitline %s"%splitline)
                #splitline=splitline[1].split('"')
                transcript_id=splitline[1]
                tidname="unset"
                if is_enstids:
                    for key,value in enstids.items():
                        if key==transcript_id:
                            tidname=value
                            break
                else:
                    tidname=transcript_id.split('.')[0]
                    #print("line %s; note_transcript_id %s; transcript_id %s; tidname %s"%(line,note_transcript_id,transcript_id,tidname))
                    truncated_transcript_id=tidname.split('ENST')[1]
                    tidname=truncated_transcript_id.lstrip("0")
                #print("note_transcript_id %s; transcript_id %s: tidname %s"%(note_transcript_id,transcript_id,tidname))
                tkey2="%s-%s"%(target_locus,tidname)
                #print("tkey %s,transcript_id %s"%(tkey,transcript_id))
                CDS_transcript_list.setdefault(tkey2,transcript_id)
                CDS_join_list.setdefault(tkey2,cds_feat_line)
                update_featfile(line)
            
            elif variation in line:
                test_variants()
            elif Origin in line:
                is_CDS=False
                finish=True
                update_featfile(line)
            elif is_variants:
                update_varfile(line)
            elif is_CDS:
                '''
                if db_xref in line:
                    pass
                '''
                if in_translate:
                    ignore=True
                elif translation in line:
                    in_translate=True
                    if exclude_translation:
                        ignore=True
                    else:
                        update_featfile(line)
                        '''update_varfile(line)'''  
                else:
                    update_featfile(line)
                    ''' update_varfile(line) '''
            elif not ignorest and not ignore:
                update_featfile(line)
        # end of elif foundloc       
        elif not start:
            if startline in line:
                start =True
            update_featfile(line)
        elif not firstgene:
            if geneline in line:
                save_generange() # Might be unnecessary here
                firstgene =True
                store=True
                storelines+=line
            else:
                update_featfile(line)
        elif locus_tag in line:
            # To distinguish, eg: /locus_tag="BRCA1" from /locus_tag="BRCA1P1"
            #if target_locus in line:
            splitters=line.split('"')
            if target_locus == splitters[1] :
                setjoin()
                foundloc=True
                last_generange=generange[:-1].strip()
                '''
                Use this instead if you don't want the target locus in the feats file 
                store_featfile(storelines)
                '''
                if genetag in storelines:
                    splitline=storelines.split('=')
                    if len(splitline)>1:
                        ensembl_geneid=splitline[1].rstrip()
                        geneid_found=True
                    else:
                        ensembl_geneid="**** WOT ******"
                    print("ensembl_geneid for json file = %s"%ensembl_geneid)
                storelines+=line
                update_featfile(storelines)
            store=False
            ''' This way we actually skip the non-target loci. Is this really what we want? '''
            storelines=""
        elif geneline in line:
            save_generange() # Definitely necessary here
            storelines=""
            store=True
            storelines+=line
        elif store:
            storelines+=line

close_outfiles()

'''
# Shell script to list the Ensemble downloaded GB zipped files and run the filter program 
pw=$PWD
/bin/ls $1 | while read dir
do
cd $pw/$1/$dir
ls *.gz | while read file
do
python3 /Users/caryodonnell/Documents/Replicon/python_scripts/embl_feature_filter3.py -i $file 
done
done
cd $pw

'''


    

