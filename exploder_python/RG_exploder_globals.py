#!/usr/local/bin/python3
#Prg_ver="RG_exploder_globals_11"
#Prg_verDate="07-Mar-2024"
# © author: Cary O'Donnell for Replicon Genetics 2018, 2019, 2020, 2021, 2022, 2023, 2024

import time  # Used in getime()
import json
import copy
# import os #needed to verify file presence using os.path # mravn change
import RG_exploder_io as RG_io # File input/output
global config_file
config_file="input/config.json"

def get_strand_extra_label():
    def get_strand_label(inval):
        label="_%s"%inval["label"].split("-")[0].lower()
        return label

    if is_frg_paired_end:
        strandlabel=get_strand_label(bio_parameters["is_frg_paired_end"])
        extralabel=strandlabel
    else:
        if is_duplex:
            strandlabel=get_strand_label(bio_parameters["is_duplex"])
            extralabel=strandlabel
        else:
            strandlabel=get_strand_label(bio_parameters["is_simplex"])
            extralabel=strandlabel
        if is_onefrag_out:
            extralabel+="_each"     
    if is_muts_only:
        extralabel+="_varonly"
    return extralabel,strandlabel

def get_locus_transcript():
# What's the current value of Reference Sequence?

    locus_transcript=target_transcript_name    
    if '(' in locus_transcript:
        splitters=locus_transcript.split('(')
        locus_transcript=splitters[0]
    if locus_transcript == empty_transcript_name or is_frg_paired_end:# Forcing when is_frg_paired_end
        locus_transcript="%s-%s"%(target_locus,empty_transcript_name.lower())
    elif is_CDS:
        locus_transcript="%s-%s"%(locus_transcript,CDS_trigger)
    else:
        locus_transcript="%s-%s"%(locus_transcript,mRNA_trigger)
        
    return locus_transcript

def save_user_configs():
    global GUI_ConfigText,config_json_end,config_json_head
    extralabel,strandlabel=get_strand_extra_label()
    #user_config_file="%s%s/%s%s%s"%(outfilepathroot,target_locus,config_json_head,get_locus_transcript(),config_json_end)
    user_config_file="%s%s/%s%s%s%s"%(outfilepathroot,target_locus,get_locus_transcript(),config_json_head,extralabel,config_json_end)
    GUI_ConfigText="Configuration at %s"%(getime())
    config_out_data=make_out_user_data()
    with RG_io.open_write(user_config_file) as write_file:
        json.dump(config_out_data, write_file,indent=4)
    write_file.close()
# end of save_user_configs()

def config_file_in():
    global config_in_data,config_json_head,config_file
    config_exists = RG_io.is_file(config_file)
    if config_exists:
        #print("config_exists %s"%config_exists)
        with RG_io.open_read(config_file) as json_data_file:
            config_in_data = json.load(json_data_file)
        set_configs()
    return config_exists
# end of config_file_in()

def make_out_user_data():
    make_custom_stringconstants()
    bio_parameters2=make_bio_parameters_configs3()
    Reference_sequences2=make_Reference_sequences2()
    config_out_user_data={
        "custom_stringconstants":custom_stringconstants,
        "bio_parameters":bio_parameters2,
        "Reference_sequences":Reference_sequences2
        }
    return config_out_user_data
# end of make_out_user_data()

def make_custom_stringconstants():
    global custom_stringconstants
    custom_stringconstants = {
        "DatasetIDText":DatasetIDText,
        "CustomerIDText":CustomerIDText,
        "GUI_ConfigText":GUI_ConfigText
         }
    #print("custom_stringconstants %s"%(custom_stringconstants))
# end of make_custom_stringconstants

def make_bio_parameters_configs3():
    # Update all the latest values worth saving for end-use examination, that haven't already been put in bio_parameters 
    global bio_parameters,mutlabels,mutfreqs
    
    # Changes to other objects in bio_parameters["target_build_variant"] should have been saved already
    # These originals were converted into distinct globals, and may have changed, so save explicitly
    
    bio_parameters["target_locus"]["value"]=target_locus
    bio_parameters["target_transcript_name"]["value"]=target_transcript_name
    bio_parameters["target_transcript_id"]["value"]=target_transcript_id
    bio_parameters["is_CDS"]["value"]=is_CDS
    bio_parameters["Fraglen"]["value"]=Fraglen
    bio_parameters["Fragdepth"]["value"]=Fragdepth
    bio_parameters["Exome_extend"]["value"]=Exome_extend
    bio_parameters["is_flip_strand"]["value"]=is_flip_strand
    bio_parameters["is_frg_paired_end"]["value"]=is_frg_paired_end
    bio_parameters["is_duplex"]["value"]=is_duplex
    #bio_parameters["is_simplex"]["value"] does not change. Only used as a label
    bio_parameters["is_fasta_out"]["value"]=is_fasta_out
    bio_parameters["is_onefrag_out"]["value"]=is_onefrag_out
    bio_parameters["is_muts_only"]["value"]=is_muts_only
    bio_parameters["is_frg_label"]["value"]=is_frg_label
    bio_parameters["is_use_absolute"]["value"]=is_use_absolute
    bio_parameters["is_fastacigar_out"]["value"]=is_fastacigar_out
    bio_parameters["is_vars_to_lower"]["value"]=is_vars_to_lower
    bio_parameters["is_journal_subs"]["value"]=is_journal_subs
    bio_parameters["is_fastq_out"]["value"]=is_fastq_out
    bio_parameters["Qualmin"]["value"]=Qualmin
    bio_parameters["Qualmax"]["value"]=Qualmax
    bio_parameters["is_write_ref_fasta"]["value"]=is_write_ref_fasta
    bio_parameters["is_mut_out"]["value"]=is_mut_out
    bio_parameters["is_write_ref_ingb"]["value"]=is_write_ref_ingb
    bio_parameters["is_sam_out"]["value"]=is_sam_out
    bio_parameters["gauss_mean"]["value"]=gauss_mean
    bio_parameters["gauss_SD"]["value"]=gauss_SD

    bio_parameters_tmp=copy.deepcopy(bio_parameters)
    bio_parameters_tmp["target_build_variant"]["mrnapos_lookup"]="hidden" # mrnapos_lookup is not typically end-use informative
 
    # Set varfreqs as a single object instead of two: mutlabels & mutfreqs (legacy stuff)
    varfreqs=dict()
    for index in range(len(mutlabels)):
                       varfreqs.setdefault(str(mutlabels[index]),mutfreqs[index])
    bio_parameters_tmp["mutlabelfreqs"]=varfreqs
    # Not using 'for item in mutlabels' here because pyodide (when running in Vue.js) treats the item "00000" as a Boolean and barfs
    # Likewise  'for item in mutfreqs', the item could be 0, and a Boolean
    
    # Or use the two. Vue.js does not barf
    #bio_parameters_tmp["mutlabels"]=mutlabels
    #bio_parameters_tmp["mutfreqs"]=mutfreqs
    return bio_parameters_tmp
# end of make_bio_parameters_configs3()
    

def make_Reference_sequences2():
    global Reference_sequences

    try: # Are the joins defined in globals via the json file?
        barf=Reference_sequences[target_locus]["CDS_join"][target_transcript_name] # tests if it exists
        CDS_join=Reference_sequences[target_locus]["CDS_join"]
        mRNA_join=Reference_sequences[target_locus]["mRNA_join"]
    except:
        CDS_join=bio_parameters["target_build_variant"]["CDS_join"]              
        mRNA_join=bio_parameters["target_build_variant"]["mRNA_join"]

    Reference_sequences_tmp={
        target_locus : {
            "Release":Reference_sequences[target_locus]["Release"],
            "Retrieval_date": Reference_sequences[target_locus]["Retrieval_date"],
            "Region":Reference_sequences[target_locus]["Region"],
            "Locus_range":Reference_sequences[target_locus]["Locus_range"],
            "is_join_complement":Reference_sequences[target_locus]["is_join_complement"],
            "LRG_id":Reference_sequences[target_locus]["LRG_id"],
            "Ensembl_id":  Reference_sequences[target_locus]["Ensembl_id"],
            "mRNA":Reference_sequences[target_locus]["mRNA"],
            "mRNA_join":mRNA_join,
            "CDS_join":CDS_join,
            "MANE_Select":Reference_sequences[target_locus]["MANE_Select"]
            }
        }
    return Reference_sequences_tmp

def process_file_configs_GUI(): # COD read config.json only for GUI read-only version
    config_file_in()
# end of process_file_configs_webUI()

def set_configs():
    set_reference_seqs_configs()
    set_string_configs()
    set_io_configs()
    set_diagnostic_configs()
    set_constants_configs()
    set_bio_parameters_configs()
    set_dependant_configs()
# end of set_configs

def set_reference_seqs_configs():
    global Reference_sequences
    Reference_sequences=config_in_data["Reference_sequences"]
# end of set_reference_seqs_configs()

def set_dependant_configs():
    # =================================================================
    # Initialising menu lists - used by GUI and web UI
    # =================================================================
    global GeneList,LRG_GeneList,Reference_sequences
    GeneList=[]
    LRG_GeneList=[]
    for locus in Reference_sequences:
        GeneList.append(locus)
        if Reference_sequences[locus]["LRG_id"]=="LRG_000":
            LRG_GeneList.append(locus)
        else:
            LRG_GeneList.append("%s (%s)"%(locus,Reference_sequences[locus]["LRG_id"]))
    
    global target_locus,ensembl_geneid
    target_locus=GeneList[0] # set the target locus variable initial value
    ensembl_geneid=Reference_sequences[target_locus]['Ensembl_id'] # set initial value
    
# end of set_dependant_configs

# ==========================
# Setting global constants
# ==========================

def set_string_configs():
        global config_in_data
        global CustomerIDText,DatasetIDText,CopyrightText
        global seq_polarity_plus,seq_polarity_minus,seq_polarity_none
        global unwanted_replace, mut_types,mRNA_trigger,CDS_trigger,empty_transcript_name,empty_transcript_id,skip_trigger
        global title_label,reference_gene,reference_haplotype,variants_header,variants_label,frequency_label,options_label,reads_list_label,reads_type_label,help_label,about_label,more_label
        global GRCh38_txt,GRCh37_txt,LRG_txt,ensembl38_gene_url,ensembl38_transcript_url,ensembl38_view_url,ensembl37_view_url,lrg_view_url,help_url
        global read_annotation,frag_annotation,spaceman,rev_label,fwd_label,IUPAC_codes
        global ensembl37_gene_url,ensembl37_transcript_url
        global ReadsList

        ReadsList=config_in_data["ReadsList"]
        stringconstants=config_in_data["stringconstants"]
        DatasetIDText=stringconstants["DatasetIDText"]
        CustomerIDText=stringconstants["CustomerIDText"]
        CopyrightText=stringconstants["CopyrightText"]
        seq_polarity_plus=stringconstants["seq_polarity_plus"]
        seq_polarity_minus=stringconstants["seq_polarity_minus"]
        seq_polarity_none=stringconstants["seq_polarity_none"]

        unwanted_replace=stringconstants["unwanted_replace"]
        mut_types=stringconstants["mut_types"]
        mRNA_trigger=stringconstants["mRNA_trigger"]
        CDS_trigger=stringconstants["CDS_trigger"]
        empty_transcript_name=stringconstants["empty_transcript_name"]
        empty_transcript_id=stringconstants["empty_transcript_id"]
        skip_trigger=stringconstants["skip_trigger"]
        title_label=stringconstants["title_label"]
        reference_gene=stringconstants["reference_gene"]
        reference_haplotype=stringconstants["reference_haplotype"]
        variants_header=stringconstants["variants_header"]
        variants_label=stringconstants["variants_label"]
        frequency_label=stringconstants["frequency_label"]
        options_label=stringconstants["options_label"]
        reads_list_label=stringconstants["reads_list_label"]
        reads_type_label=stringconstants["reads_type_label"]
        help_label=stringconstants["help_label"]
        about_label=stringconstants["about_label"]
        more_label=stringconstants["more_label"]
        
        GRCh38_txt=stringconstants["GRCh38_txt"]
        GRCh37_txt=stringconstants["GRCh37_txt"]
        LRG_txt=stringconstants["LRG_txt"]
        ensembl38_gene_url=stringconstants["ensembl38_gene_url"]
        ensembl38_transcript_url=stringconstants["ensembl38_transcript_url"]
        ensembl38_view_url=stringconstants["ensembl38_view_url"]
        ensembl37_gene_url=stringconstants["ensembl37_gene_url"]
        ensembl37_transcript_url=stringconstants["ensembl37_transcript_url"]
        ensembl37_view_url=stringconstants["ensembl37_view_url"]
        lrg_view_url=stringconstants["lrg_view_url"]
        help_url=stringconstants["help_url"]
        read_annotation=stringconstants["read_annotation"]
        frag_annotation=stringconstants["frag_annotation"]
        spaceman=stringconstants["spaceman"]
        rev_label=stringconstants["rev_label"]
        fwd_label=stringconstants["fwd_label"]
        IUPAC_codes=stringconstants["IUPAC_codes"]
# end of set_string_configs()

def set_io_configs():
    global config_in_data
    global infilepathroot,outfilepathroot,pygui_outfilepathroot
    global journhead,journend,readmehead,readmeend,config_json_head,config_json_end
    global Seq_Format,Seq_IO_file_ext,Ref_file_name
    global mutlabels

    IOconstants=config_in_data["IOconstants"]
    infilepathroot=IOconstants["infilepathroot"]
    outfilepathroot=IOconstants["outfilepathroot"]
    pygui_outfilepathroot=IOconstants["pygui_outfilepathroot"]
    journhead=IOconstants["journhead"]
    journend=IOconstants["journend"]
    readmehead=IOconstants["readmehead"]
    readmeend=IOconstants["readmeend"]
    config_json_head=IOconstants["config_json_head"]
    config_json_end=IOconstants["config_json_end"]
    Seq_Format=IOconstants["Seq_Format"]
    Seq_IO_file_ext=IOconstants["Seq_IO_file_ext"]
    Ref_file_name=IOconstants["Ref_file_name"]
# end of set_io_configs()

def set_diagnostic_configs():
    global is_include_indels, is_dels_to_dots, is_dels_to_dots_override
    global is_print_diagnostics
    global is_mutate_ref,is_make_exome,is_trim_to_gene
    global is_fastq_random
    global is_show_infilepath
    global MaxVarPos,is_htm_journal,is_make_reference_files1
    global config_in_data
    
    diagnostic_vars=config_in_data["diagnostic_vars"]
    is_include_indels=diagnostic_vars["is_include_indels"]
    is_dels_to_dots=diagnostic_vars["is_dels_to_dots"]
    is_dels_to_dots_override=diagnostic_vars["is_dels_to_dots_override"]
    is_print_diagnostics=diagnostic_vars["is_print_diagnostics"]
    is_mutate_ref=diagnostic_vars["is_mutate_ref"]
    is_make_exome=diagnostic_vars["is_make_exome"]
    is_trim_to_gene=diagnostic_vars["is_trim_to_gene"]
    is_fastq_random=diagnostic_vars["is_fastq_random"]
    is_show_infilepath=diagnostic_vars["is_show_infilepath"]
    MaxVarPos=diagnostic_vars["MaxVarPos"]
    is_htm_journal=diagnostic_vars["is_htm_journal"]
    is_make_reference_files1=diagnostic_vars["is_make_reference_files1"]

# end of set_diagnostic_configs()

def set_constants_configs():
    global config_in_data,NormaliseUpper,Mapq_Min,Mapq_Max,mutlabels,mutfreqs,DOC_precision,Paired_insert_Min,Paired_insert_Max
    admin_constants=config_in_data["admin_constants"]
    NormaliseUpper=admin_constants["NormaliseUpper"]
    Mapq_Min=admin_constants["Mapq_Min"]
    Mapq_Max=admin_constants["Mapq_Max"]
    mutlabels=admin_constants["mutlabels"]
    mutfreqs=admin_constants["mutfreqs"]
    DOC_precision=admin_constants["DOC_precision"]
    Paired_insert_Min=admin_constants["Paired_insert_Min"]
    Paired_insert_Max=admin_constants["Paired_insert_Max"]
    
# end of set_constants_configs

def set_bio_parameters_configs():
    global config_in_data,bio_parameters
    global target_locus,is_fastacigar_out,is_onefrag_out,is_muts_only,is_frg_paired_end,is_flip_strand,is_duplex,is_frg_label,is_journal_subs,is_make_exome,is_trim_to_gene
    global is_vars_to_lower, is_fasta_out,is_fastq_out,is_sam_out,is_mut_out,is_write_ref_fasta,is_write_ref_ingb,is_use_absolute
    global gauss_mean,gauss_SD,gauss_mean_Max,gauss_mean_Min,gauss_SD_Max,gauss_SD_Min
    #global mutfreqs,mutlabels
    global target_transcript_name,target_transcript_id,is_CDS

    global Exome_extend,Exome_extend_Min,Exome_extend_Max,Fraglen,FraglenMin,FraglenMax
    global Fragdepth,FragdepthMin,FragdepthMax,Qualmin,QualMIN,Qualmax,QualMAX
    
    bio_parameters=config_in_data["bio_parameters"]
    
    target_locus=bio_parameters["target_locus"]["value"]
    target_transcript_name=bio_parameters["target_transcript_name"]["value"]
    target_transcript_id=bio_parameters["target_transcript_id"]["value"]
    
    is_CDS=bio_parameters["is_CDS"]["value"]
    is_fastacigar_out=bio_parameters["is_fastacigar_out"]["value"]
    is_onefrag_out=bio_parameters["is_onefrag_out"]["value"]
    is_muts_only=bio_parameters["is_muts_only"]["value"]
    is_frg_paired_end=bio_parameters["is_frg_paired_end"]["value"]
    is_flip_strand=bio_parameters["is_flip_strand"]["value"]
    is_duplex=bio_parameters["is_duplex"]["value"]
    is_frg_label=bio_parameters["is_frg_label"]["value"]
    is_journal_subs=bio_parameters["is_journal_subs"]["value"]
    is_vars_to_lower=bio_parameters["is_vars_to_lower"]["value"]
    is_fasta_out=bio_parameters["is_fasta_out"]["value"]
    is_fastq_out=bio_parameters["is_fastq_out"]["value"]
    is_sam_out=bio_parameters["is_sam_out"]["value"]
    is_mut_out=bio_parameters["is_mut_out"]["value"]
    is_write_ref_fasta=bio_parameters["is_write_ref_fasta"]["value"]
    is_write_ref_ingb=bio_parameters["is_write_ref_ingb"]["value"]   
    is_use_absolute=bio_parameters["is_use_absolute"]["value"]
    gauss_mean=bio_parameters["gauss_mean"]["value"]
    gauss_SD=bio_parameters["gauss_SD"]["value"]
    gauss_mean_Max=bio_parameters["gauss_mean"]["max"]
    gauss_mean_Min=bio_parameters["gauss_mean"]["min"]
    gauss_SD_Max=bio_parameters["gauss_SD"]["max"]
    gauss_SD_Min=bio_parameters["gauss_SD"]["min"]
    
    #mutfreqs=bio_parameters["mutfreqs"]["value"]
    #mutlabels=bio_parameters["mutfreqs"]["label"]
    
    Exome_extend=bio_parameters["Exome_extend"]["value"]
    Exome_extend_Min=bio_parameters["Exome_extend"]["min"]
    Exome_extend_Max=bio_parameters["Exome_extend"]["max"]
    
    Fraglen=bio_parameters["Fraglen"]["value"]
    FraglenMin=bio_parameters["Fraglen"]["min"]
    FraglenMax=bio_parameters["Fraglen"]["max"]
    
    Fragdepth=bio_parameters["Fragdepth"]["value"]
    FragdepthMin=bio_parameters["Fragdepth"]["min"]
    FragdepthMax=bio_parameters["Fragdepth"]["max"]
    
    
    Qualmin=bio_parameters["Qualmin"]["value"]
    QualMIN=bio_parameters["Qualmin"]["min"]
    
    Qualmax=bio_parameters["Qualmax"]["value"]
    QualMAX=bio_parameters["Qualmax"]["max"]    
# end of set_bio_parameters_configs

def getime():
    localtime = time.asctime(time.localtime(time.time()))
    return localtime


# =================================================================
# Main
# =================================================================
process_file_configs_GUI()
