#!/usr/local/bin/python3
#Prg_ver="RG_exploder_globals_make
#Prg_verDate="10-Apr-2025"
'''
This module creates the config.json file from all the contributing input directories and is part of the data-preparation. 
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
'''
import sys
import json
import subprocess
import os
'''
RG module imports
'''
import RG_exploder_io as RG_io # File input/output

######################  SET THE GLOBALS IN set_config_consts() BEFORE RUNNING ################

def set_config_consts():
    global exploder_root
    global is_use_TKinter,CustomerIDText,PublicIDText,MAINFILE #   These are most likely to need resetting between runs

    global is_pair_monitor
    ######## Revisit these each time a data set is renewed. Now down to 2 as MAINVER, DATEVER, GRCH_dataset derived ########
    #CustomerIDText="EBaker"
    PublicIDText="Public"
    CustomerIDText=PublicIDText # Set to this to override a CustomerIDText
    #CustomerIDText="C Roos" # See also: is_exome_paired_end ; is_pair_monitor
    vuedir="mravn" # Solely to define how is_use_TKinter is set
    sharedir="shared" # To define how 
    ######## End of: Revisit these each time a data set is renewed ########

    #### Should not need changing unless modify throughout other code ####

    global config_file_output,config_file_reference_seqs
    global rootdatadir,rootapplicationdir,roothelpscriptdir

    set_exploder_root="../" # Do it relative
    exploder_root="%s/"%os.path.abspath(set_exploder_root) # Ensure this is defined as full path, required to support the hyperlinks to results in the Python GUI

    root_sources="data_sources"
    rootdatadir="%sdata_sources"%exploder_root
    rootapplicationdir="%sexploder_python"%exploder_root # Not used for App.vue setup            
    roothelpscriptdir="%shelper_scripts"%exploder_root

    config_file_output="input/config.json"
    config_file_reference_seqs="input/loci.json"
    MAINFILE="RG_exploder_main.py"
    #MAINFILE="RG_exploder_main_roos.py"

    if vuedir in exploder_root:
        is_use_TKinter=False  # When setting up to use the Vue.js GUI
        is_pair_monitor=False
    else:
        is_use_TKinter=True  # When setting up to use the Python GUI : RG_exploder_gui.py
        is_pair_monitor=True # Sets to False when is_use_TKinter=True

    if sharedir in exploder_root:
        is_pair_monitor=False # True to save monitor file with paired ends; False - no monitor file 
        
    set_defaults()

def config_file_out():
    #set_defaults()
    out_data=make_out_data()
    with RG_io.open_write(config_file_output) as write_file:
        json.dump(out_data, write_file,indent=4)
# end of config_file_out()

def config_file_in(in_file):
    config_exists = RG_io.is_file(in_file)
    config_in_data=""
    if config_exists:
        #print("config_exists %s"%config_exists)
        with RG_io.open_read(in_file) as json_data_file:
            config_in_data = json.load(json_data_file)
    return config_exists,config_in_data
# end of config_file_in()


def main_file_read(in_file):
    main_exists = RG_io.is_file(in_file)
    line2=""; line3=""
    if main_exists:
        #print("main_exists %s"%main_exists)
        with RG_io.open_read(in_file) as main_file:
            main_file.readline()
            line2=main_file.readline()
            line3=main_file.readline()
        main_file.close()
    return main_exists,line2,line3
# end of main_file_read()

def process_file_configs_python(): # Only create config-file if file not already present (manual deletion necessary)
    global infilepathroot
    config_exists = RG_io.is_file(config_file_output)
    if config_exists:
        subprocess.run(["rm","%s%s"%(infilepathroot,config_file_output)])
        print("Config file %s%s RE-write"%(infilepathroot,config_file_output))
    else:
        print("Config file %s%s NEW write"%(infilepathroot,config_file_output)) 
    config_file_out()
# end of process_file_configs_python()

def set_defaults():
    global CustomerIDText,PublicIDText
    set_string_defaults()
    set_io_defaults()
    set_diagnostic_defaults()
    set_constants_defaults()
    set_gui_vars_defaults()
    set_vars_seq_limits_defaults()
    if CustomerIDText==PublicIDText:
        set_Reference_sequences_configs_public3(PublicIDText)
    else:
        print("No set up for CustomerIDText==%s"%CustomerIDText)
        set_Reference_sequences_configs_public3(CustomerIDText)
# end of set_defaults

def set_defaults1():
    global GRCH_dataset, CustomerIDText
    set_string_defaults()
    set_io_defaults()
    set_diagnostic_defaults()
    set_constants_defaults()
    set_gui_vars_defaults()
    set_vars_seq_limits_defaults()
    if GRCH_dataset=="GRCh37":
        if CustomerIDText=="Public":
            #set_Reference_sequences_configs_GRCh37_1000_public()
            set_Reference_sequences_configs_GRCh37_1000_public2()
        else:
            set_Reference_sequences_configs_GRCh37_1000()
    else:
        if CustomerIDText=="Public":
            #set_Reference_sequences_configs_GRCh38_1000_public()
            set_Reference_sequences_configs_GRCh38_1000_public2()
        else:
            set_Reference_sequences_configs_GRCh38_1000()
    set_dynamic_vars_defaults()
# end of set_defaults1

def make_out_data():
    global ReadsList
    make_stringconstants()
    make_IOconstants()
    make_diagnostic_vars()
    make_admin_constants()
    make_bio_parameters_configs()
    config_out_data={
        "stringconstants":stringconstants,
        "IOconstants":IOconstants,
        "diagnostic_vars":diagnostic_vars,
        "admin_constants":admin_constants,
        "bio_parameters":bio_parameters,
        "ReadsList":ReadsList,
        "Reference_sequences":Reference_sequences,
        }
    return config_out_data
# end of make_out_data()

# ==========================
# Setting global constants
# ==========================
def set_string_defaults():
    global CopyrightText
    CopyrightText="Copyright © Replicon Genetics & Cary O'Donnell 2021-2025. Available under GNU Affero General Public License version 3"
    #There are several ways of declaring polarity so we define the known ones and use as needed
    global seq_polarity_plus,seq_polarity_minus,seq_polarity_none
    seq_polarity_plus="1"; seq_polarity_minus="-1";seq_polarity_none="0"

    global unwanted_replace # ''' list of undesired variant substitutions that must be filtered out. May want to expand this
    unwanted_replace=["HGMD_MUTATION","dbSNP_variant"]

    global mut_types,IUPAC_codes
    # mut_types is a pattern object for the mutation types in CIGAR that selects fragments containing them when is_muts_only = True  '''
    #           Used by is_mut_cigar'''
    mut_types = "DIX" # Mutation types declared in CIGAR: Deletion, Insertion, X - variant.
                      # N is a CIGAR item, but is ignored as a *mutation* type

    IUPAC_codes="[ACGTRYSWKMBDHVN]"

    global mRNA_trigger,skip_trigger,CDS_trigger,empty_transcript_name,empty_transcript_id
    mRNA_trigger="mRNA"       # ''' CONSTANT: External dependency. This word in place of join in the refseq feature table will initiate a splice: the join regions become exons'''
    skip_trigger="gap"          # ''' CONSTANT: Used internally to create an entry in the feature table of a "variant" sequence, in the position of the "splice"
                                # '''           and so skip this region when considering how to fragment it in MutateVarSeq
                                # '''           The clever bit comes when assessing whether another variant feature coincides, or overlaps with, this skip and how to handle it: in get_mutref_olap

    CDS_trigger="CDS"           # ''' CONSTANT: External dependency. Use in combination with is_CDS
    #empty_transcript_name="Locus"
    empty_transcript_name="Gene"
    empty_transcript_id=""

    
    global GRCh38_txt,GRCh37_txt,LRG_txt
    GRCh38_txt="GRCh38"  # CONSTANT - text id for genome build version
    GRCh37_txt="GRCh37"   #  CONSTANTS  - text id for genome build version
    LRG_txt="LRG" #  CONSTANTS - label for URL only
    
    global ensembl38_gene_url,ensembl38_transcript_url,ensembl38_view_url,ensembl37_view_url,lrg_view_url,help_url,about_url,more_url
    global ensembl37_gene_url,ensembl37_transcript_url
    ensembl38_gene_url="https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;"  # ''' CONSTANT for Ensembl hyperlink to gene
    #eg: https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000133703;r=12:25205246-25250936
    ensembl38_transcript_url="https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;"  # ''' CONSTANT for Ensembl hyperlink to transcript
    
    ensembl37_gene_url="https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;"  # ''' CONSTANT for Ensembl 37 hyperlink to gene
    ensembl37_transcript_url="https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;"  # ''' CONSTANT for Ensembl 37 hyperlink to transcript


    # Alternatively:
    ensembl38_view_url="https://www.ensembl.org/Homo_sapiens/Location/View?r=" # ''' CONSTANT for Ensembl hyperlink to gene
      #eg: https://www.ensembl.org/Homo_sapiens/Location/View?r=12:25205246-25250936
    ensembl37_view_url="https://grch37.ensembl.org/Homo_sapiens/Location/View?r="
    #eg:"https://grch37.ensembl.org/Homo_sapiens/Location/View?"
    #lrg_view_url="http://ftp.ebi.ac.uk/pub/databases/lrgex/"
    lrg_view_url="https://www.ensembl.org/Homo_sapiens/LRG/Summary?lrg=LRG_"

    help_url="https://repliconevaluation.wordpress.com/replicon-genetics" # The help page
    about_url="https://repliconevaluation.wordpress.com/about" # The about page
    more_url="https://repliconevaluation.wordpress.com/more" # The more page

    global spaceman
    #spaceman="\t"   # Space separator in Fastaheader is a tab
    spaceman=" "    # Space separator in Fastaheader is a single blank space

    global fwd_label,rev_label,read_annotation,frag_annotation
    fwd_label="f" # In case someone wants it to be blank or fwd or something else
    fwd_label="%sfwd"%spaceman
    fwd_label="/1" # David Orr wants this
    rev_label="r" # In case someone wants it to be blank or rev or something else
    rev_label="%srev"%spaceman #
    rev_label="/2" # David Orr wants this

    read_annotation="read" # Do we call the output sequences "read" or "fragment" or what?
    #read_annotation="fragment"
    frag_annotation="frg" # The prefix to the fragment-read number in fasta description eg: >frg_1 

    global hapname,title_label,template_label,reference_gene,reference_haplotype,variants_header,variants_label,frequency_label,options_label,reads_list_label,reads_type_label,help_label,about_label,more_label
    hapname="Haplotype"
    #title_label="NGS Read Simulator"
    title_label="Synthetic %ss Generator"%read_annotation.title() 
    #reference_gene="Reference Sequence"
    template_label="Template"
    #reference_gene="Reference Source"
    reference_gene="Reference Gene"
    reference_haplotype="Reference %s"%hapname
    #variants_header="Variant Sequences"
    #variants_header="%s Definition"%hapname
    variants_header="%s %s"%(hapname,template_label)
    #variants_label="Variant Source"
    #variants_label="Variations Source"
    #variants_label="%s Source"%hapname
    variants_label="%s Variant"%hapname
    #frequency_label="Source Frequency"
    frequency_label="Source Ratio"
    options_label="Output Options"
    reads_list_label="%ss Configuration"%read_annotation.title()
    reads_type_label="%ss type"%read_annotation.title()
    help_label="Need help..?"
    about_label="About..."
    more_label="More..."
# end of set_string_defaults()

def make_stringconstants():
    global stringconstants
    stringconstants = {
        "DatasetIDText":DatasetIDText,
        "CustomerIDText":CustomerIDText,
        "CopyrightText":CopyrightText,
        "seq_polarity_plus":seq_polarity_plus,"seq_polarity_minus":seq_polarity_minus,
        "seq_polarity_none":seq_polarity_none,
        "unwanted_replace":unwanted_replace,
        "mut_types":mut_types,
        "mRNA_trigger":mRNA_trigger,
        "skip_trigger":skip_trigger,
        "CDS_trigger":CDS_trigger,
        "empty_transcript_name":empty_transcript_name,"empty_transcript_id":empty_transcript_id,
        "GRCh38_txt":GRCh38_txt,"GRCh37_txt":GRCh37_txt,"LRG_txt":LRG_txt,
        "title_label":title_label,
        "reference_gene":reference_gene,
        "reference_haplotype":reference_haplotype,
        "variants_header":variants_header,
        "variants_label":variants_label,
        "frequency_label":frequency_label,
        "options_label":options_label,
        "reads_list_label":reads_list_label,
        "reads_type_label":reads_type_label,
        "help_label":help_label,
        "about_label":about_label,
        "more_label":more_label,
        "ensembl38_gene_url":ensembl38_gene_url,
        "ensembl38_transcript_url":ensembl38_transcript_url,
        "ensembl38_view_url":ensembl38_view_url,
        "ensembl37_gene_url":ensembl37_gene_url,
        "ensembl37_transcript_url":ensembl37_transcript_url,
        "ensembl37_view_url":ensembl37_view_url,
        "lrg_view_url":lrg_view_url,
        "help_url":help_url,
        "about_url":about_url,
        "more_url":more_url,
        "read_annotation":read_annotation,
        "frag_annotation":frag_annotation,
        "spaceman":spaceman,
        "rev_label":rev_label,"fwd_label":fwd_label,
        "IUPAC_codes":IUPAC_codes
         }
    #print("stringconstants %s"%(stringconstants))
# end of make_stringconstants

def set_io_defaults():
    # =================================================================
    # Setting global constants for IO
    # =================================================================
    global infilepathroot,outfilepathroot,pygui_outfilepathroot # CONSTANTS: set filepaths for location of input and output files. Obviously better as command-line parameters or config-file setting

    global journhead,journend,readmehead,readmeend,config_json_head,config_json_end
    journhead="_001"; journend="_journal"  # CONSTANTS for output journal file
    readmehead="_000"; readmeend="_readme" # CONSTANTS for output readme file
    config_json_head="_002"                # CONSTANT for output config file
    config_json_end="_config.txt"          # CONSTANT for output config file

    global Seq_Format,Seq_IO_file_ext,Ref_file_name
    Seq_Format ="genbank" #  CONSTANT: For sequence_format file IO '''
    Seq_IO_file_ext=".gb" #  CONSTANT: File extension name to coincide with Seq_Format '''
    #Ref_file_name="ref"   #  CONSTANT: File name extension xxx as in EGFR_xxx.gb
    #Ref_file_name="loc"  # switching back to "ref" and will use "hap0" for "variant" that doesn't contain any variants ... April 2021
    Ref_file_name="locseq"   #  CONSTANT: File name extension xxx as in EGFR_xxx.gb 

    global mutlabels # These are the root names for the expected input files.
                 # When combined with infilepathroot and target_locus deliver the full filename path
    mutlabels=["hap0","hap1","hap2"]
    for count in range(1,10):
        mutlabels.append("som%s"%(str(count)))

    #make_IOconstants()
# end of set_io_defaults()

def make_IOconstants():
    global IOconstants
    IOconstants = {
        "infilepathroot":"input/",
        "outfilepathroot":"output/",
        "pygui_outfilepathroot":pygui_outfilepathroot,
        "journhead":journhead,"journend":journend,"readmehead":readmehead,"readmeend":readmeend,
        "config_json_head":config_json_head,
        "config_json_end":config_json_end,
        "Seq_Format":Seq_Format,"Seq_IO_file_ext":Seq_IO_file_ext,"Ref_file_name":Ref_file_name,
        }
    #print("IOconstants %s"%IOconstants)
# end of make_IOconstants()

def set_diagnostic_defaults():
    # ===========================================================================
    # Setting admin-configurable diagnostic switches for development-testing-only
    # Ensure each is set to "User-config default" or expect trouble
    # ===========================================================================
    global is_include_indels, is_dels_to_dots, is_dels_to_dots_override,is_exome_paired_end,is_pair_monitor,is_use_TKinter

    is_exome_paired_end=True # False to force paired-end for mRNA & CDS; True to enable exomic paired-end
    if not is_use_TKinter: # Way too slow for high parameter values giving big files in Vue.js
        is_pair_monitor=False

    is_include_indels=True      #  Set is_include_indels False to ignore indel definitions in variant feature table.
                                #  User-config default should be True
    is_dels_to_dots=False       # Deleted bases are converted to dots in MutateVarSeq when is_dels_to_dots flag is True .
                                #  User-config default should be False
    is_dels_to_dots_override=False   # is_dels_to_dots_override is only used in get_mutrecords(REF_record,embl_or_genbank) to temporarily
                                     #  switch is_dels_to_dots to is_dels_to_dots_override to prevent the dots appearing in output files
                                     # But with is_dels_to_dots_override=True, it will override and the dots will appear
                                     #  in {target}mut.fasta as well as the fragments. Useful for diagnostics, but not for real data output
                                     #  16-Aug-2020: wondering if is_dels_to_dots_override is now vestigial and epistatic
                                     #  User-config default should be False

    global is_print_diagnostics      # Switch to print more verbose runtime messages for diagnostic purposes.'''
    is_print_diagnostics=False       #Previously also had is_print_progress messages which were converted to journal output '''
                                     #  User-config default should be False

    global is_make_exome
    is_make_exome =True   # is_make_exome =True the Reference sequence will be processed for exons, (but only if the reference file has "splice" features)
                          # is_make_exome = False Reference is not processed for exons
                          # This switch is solely here as an override to permit processing of the same data & parameters set with- and without- splicing
                          # It is deprecated as a user-selectable option

    global is_trim_to_gene
    is_trim_to_gene =True # is_trim_to_gene =True the source for the Reference sequence will be trimmed to the gene boundary, in the same way as used for exons)
                          # is_trim_to_gene = False the source for the Reference sequence is not trimmed 
                          # This switch is solely here as an override to permit processing of the same data & parameters set with- and without- trimming
                          # It is not available as a user-selectable option
    
    global is_mutate_ref
    is_mutate_ref= False  #  To include the variants defined in the Refseq file included within the Reference sequence set: is_mutate_ref= True
                          #  Variants in Refseq file should NOT include inserts or deletions where any haplotype variants are to be introduced
                          #  as all position-definitions resulting will be incorrect.
                          #  User-config default should therefore be False
                          #  This has potential for re-use as a feature - but more development and testing work needed
                          #  Jan-2022 - is_mutate_ref not used as originally intended with make_ref_varseq() deprecated

    global is_fastq_random  # When set to True, this implements an all-positions random in FASTQ output
    is_fastq_random= True   # With Qualmin==Qualmax, random is NOT used despite this setting
    #is_fastq_random= False # When False it implements a faster "slice" routine than a random: using random strings cut at a random position, but this is still slower than Qualmin==Qualmax

    global is_show_infilepath
    is_show_infilepath=False # The Journal file will show full path name of input files when True. Use True for debugging.
                            # Only the stub (without infilepathroot) is shown when set to False.
                            #  User-config default should be False

    global is_htm_journal # Save an htm version of the journal file when set to True 
    is_htm_journal =True

    global MaxVarPos # The maximum length of Refseq used. It is set to the length of the reference sequence in read_refseqrecord, unless set here > 0''' 
    MaxVarPos=0      #  Set it here for testing purposes eg: MaxVarPos=250 to truncate the sequence that is read in'''
                     #  User-config default should be zero
                     # NB: Unresolved bug when MaxVarPos > 0
                     #     File ... RG_exploder_process.py", line 756, in get_absolute_position
                     #     abs_pos=seq_record.offset+seq_record.strand_mod*seqpos
                     #     AttributeError: 'SeqRecord' object has no attribute 'offset'


    global is_reads_upper # FASTA reads output case
    is_reads_upper = True    # True:  FASTA sequence in upper case; inserts & SNP are lower case
    #is_reads_upper = False  # False: FASTA sequence in lower case; inserts & SNP are upper case    


# end of set_diagnostic_defaults()

def make_diagnostic_vars():
    global diagnostic_vars
    diagnostic_vars ={
        "is_exome_paired_end":is_exome_paired_end,
        "is_pair_monitor":is_pair_monitor,
        "is_include_indels":is_include_indels,
        "is_dels_to_dots":is_dels_to_dots,
        "is_dels_to_dots_override":is_dels_to_dots_override,
        "is_print_diagnostics":is_print_diagnostics,
        "is_mutate_ref":is_mutate_ref,
        "is_make_exome":is_make_exome,
        "is_trim_to_gene":is_trim_to_gene,
        "is_fastq_random":is_fastq_random,
        "is_show_infilepath":is_show_infilepath,
        "is_htm_journal":is_htm_journal,
        "MaxVarPos":MaxVarPos,
        "is_reads_upper":is_reads_upper
        }
#print("diagnostic_vars %s"%diagnostic_vars)
#end of make_diagnostic_vars()

def set_constants_defaults():
    # ====================================================================
    # Setting admin-configurable constants required across modules
    # ====================================================================

    global QualMIN,QualMAX  # Absolute range values for Sanger FASTQ quality Qualmin Qualmax variables
    QualMIN=0
    QualMAX=93

    global Exome_extend_Max, Exome_extend_Min # Absolute range values for Exome_extend
    Exome_extend_Max=50
    Exome_extend_Min=0

    global FraglenMin,FraglenMax,FragdepthMin,FragdepthMax # Absolute range values for Fraglen and Fragdepth
    #FraglenMin=8 - the usual minimum
    FraglenMin=4  # Smaller numbers for test purposes
    FraglenMax=2000
    FragdepthMin=1
    #FragdepthMax=500 # original
    FragdepthMax=800 # Christophe Roos

    global NormaliseUpper # Figure for normalising mutfreqs values
    NormaliseUpper=1000

    global Mapq_Min,Mapq_Max # Used in write_frag_samout range mapq=randint(Mapq_Min,Mapq_Max)
    Mapq_Min=15
    Mapq_Max=30

    global DOC_precision # Number of decimal places to report DOC
    DOC_precision=2
    
# end of set_constants_defaults

def make_admin_constants():
    global admin_constants,NormaliseUpper,Mapq_Min,Mapq_Max,mutlabels,mutfreqs,DOC_precision,Paired_insert_Min,Paired_insert_Max
    admin_constants={
        "NormaliseUpper":NormaliseUpper,"Mapq_Min":Mapq_Min,"Mapq_Max":Mapq_Max,
        "mutlabels":mutlabels,"mutfreqs":mutfreqs,"DOC_precision":DOC_precision,
        "Paired_insert_Min":Paired_insert_Min,"Paired_insert_Max":Paired_insert_Max
        }
#end of make_admin_constants

def set_gui_vars_defaults():
    # =======================================================================
    # Initialising gui-configurable global flags for sequence output features
    # These are all user-changeable variables
    # =======================================================================
    ''' def initialise_flags_output_configs()'''

    global is_CDS
    is_CDS = False      # is_CDS = True uses the join configuration for CDS instead of mRNA when target_transcript_name is not == empty_transcript_name ie: not genomic
                        # is_CDS = False uses the mRNA join configuration when not genomic

    global is_fastacigar_out,is_onefrag_out,is_muts_only,is_frg_paired_end,is_flip_strand, is_duplex,is_frg_label,is_journal_subs,is_use_absolute,is_vars_to_lower
    is_fastacigar_out = True #  is_fastacigar_out = True appends the CIGAR string to the fragmented fasta, fastq and {target}_mut.fasta sequences.
                             #  is_fastacigar_out = False to prevent this
    is_onefrag_out = False   #  is_onefrag_out = False to create random start positions, in random order to meet a fragdepth target.
                             #  is_onefrag_out = True : sequential start positions, just once
    is_muts_only = False  # is_muts_only = True : suppresses output of non-variant-containing fragments (ie: suppress where CIGAR contains only M).
                          # is_muts_only = False: allows all fragments to be produced.
    is_frg_paired_end= True  # "None" will hide the option in GUI menu.
    is_flip_strand = False # Testing
    is_duplex=False           # is_duplex = True  - 50% of created fragments are delivered in reverse-complement, so there is a mix of forward and reverse strands
                              # is_duplex = False - all fragments are in forward orientation
    is_frg_label = True   # is_frg_label = True labels the fragment source  (hapN, somN)
                          # is_frg_label = False for no labelling

    is_journal_subs = False  # is_journal_subs = True to write status of each subsituted variant to journal whether match or mismatch to definition
                          # is_journal_subs = False - suppress status of **matched** subsituted variant to journal
                          
    is_use_absolute =True # is_use_absolute =True - show variants in absolute genomic coordinates rather than relative, as defined in the ref file
                          # is_use_absolute =False - use relative positions only

    is_vars_to_lower=False  # is_vars_to_lower=True: Variant nucleotides (SNVs and Inserts) are converted to lower case, in contrast to the rest of the upper-case seqeunce
                            # is_vars_to_lower=False: all fragment-nucleotides are shown in upper case
                            # Note interaction with is_reads_upper which switches these around when is_reads_upper=False

    global is_fasta_out,is_fastq_out,is_sam_out,is_mut_out
    is_fasta_out=True # is_fasta_out = True creates a FASTA output of the fragmented mutated sequences.
                      # is_fasta_out = False suppresses output of FASTA file
    is_fastq_out=False # is_fastq_out = True creates a FASTQ output of the fragmented mutated sequences.
                      # is_fastq_out = False suppresses output of FASTQ file
    is_sam_out=False   # is_sam_out = True creates a SAM output of the fragmented mutated sequences.
                      # is_sam_out = False suppresses output of SAM file
    is_mut_out=False   # is_mut_out = True creates a Fasta output of the **un-fragmented** mutated sequences
                      # is_mut_out = False suppresses this

    global is_write_ref_fasta, is_write_ref_ingb
    is_write_ref_fasta = False  # is_write_ref_fasta = True writes out the input reference file as a fasta format file
                               # is_write_ref_fasta = False suppresses this ... is this beung used differently?

    is_write_ref_ingb = False   # is_write_ref_ingb = True writes out the input variant files in genbank-format, feature table (no sequence).
                               #  is_write_ref_ingb = False suppresses this.

# end of set_gui_vars_defaults()

def make_bio_parameters_configs():
    global bio_parameters,ReadsList,read_annotation,variants_header,empty_transcript_name,is_reads_upper
    if is_reads_upper:
        case_subs="lower"
    else:
        case_subs="upper"
    bio_parameters={
        "target_locus":{
            "label": empty_transcript_name,
            "value": target_locus
            },
        "target_transcript_name":{
            "label":template_label,
            "value":target_transcript_name
            },
        "target_transcript_id":{
            "value":target_transcript_id
            },
        "target_build_variant": {
            "is_get_ref": False,
            "is_save_var": False,
            "is_get_muttranscripts":False,
            "mRNA_join":"",
            "CDS_join":"",
            "headclip":0,
            "tailclip":0,
            "exonplus_lookup":[0],
            "transcript_view":"",
            "abs_offset":0,
            "ref_strand":1,
            "local_begin":0,
            "local_end":0,
            "joinlist":"",
            "refseq":"ATGC",  #  Test value. Should be ""

            "abs_Begin":{
                "label":"Absolute Begin",
                "value":0
                },
            
            "abs_End":{
                "label":"Absolute End",
                "value":0, 
                },
            
            "trans_Begin":{
                "label":"Locus",
                "value":1,
                "min":1, 
                "max":15
                },

            "trans_End":{
                "label":"Locus",
                "value":1,
                "min":1,
                "max":15
                },
            
            "trans_Begin_ext": {
                "label":"Begin extension",
                "value":0,
                "min":-50,
                "max":50
                },
            
            "trans_End_ext": {
                "label":"End extension",
                "value":0, 
                "min":-50,
                "max":50
                },
            
            "ref_subseq": {
                "label":"Reference Sequence", 
                "value":"", 
                "length":0,
                "viewstring":""
                },

            "var_subseq": {
                "label":"Variant Sequence",
                "length":0,
                "value":""
                },

            "hap_name":{
                "label":"Source Name",
                "value":"hapx"
                },
 
             "var_name":{
                "label":"Variant Name",
                "value":""
                },
            # "variant_view":"", # Might be useful - see RG_exploder_main.save_add_var()
            "AddVars": []
            },
        "is_CDS": {
            "label": "CDS only",
            "value":is_CDS
            },
        "mutfreqs": {
            "label":[],
            "value":[]
            },
        "Fraglen": {
            "label": "%s length"%read_annotation.title(),
            "value": Fraglen,
            "min": FraglenMin,
            "max": FraglenMax
            },
        "Fragdepth": {
            "label": "Depth of cover",
            "value": Fragdepth,
            "min": FragdepthMin,
            "max": FragdepthMax
            },
        "Exome_extend": {
            "label": "(Exome) Extension",
            "value": Exome_extend,
            "min": Exome_extend_Min,
            "max": Exome_extend_Max
            },
        "is_flip_strand": {
            "label": "Flip polarity", 
            "value": is_flip_strand
            },
        "is_frg_paired_end": {
            "label": "Paired-end", # must be two words separated by - ; see def get_strand_label in main
            "value": is_frg_paired_end
            },
        "is_duplex": {
            "label": "Dual-strand", # must be two words separated by - ; see def get_strand_label in main
            "value": is_duplex
            },
        "is_simplex": {
            "label": "Single-strand", # must be two words separated by - ; see def get_strand_label in main
            "value": None # Never changed, it's only here for the label
            },
        "is_fasta_out": {
            "label": "%ss in FASTA format"%read_annotation.title(),
            "value": is_fasta_out
            },
        "is_onefrag_out": {
            "label": "- Each possible %s"%read_annotation,
            "value": is_onefrag_out
            },
        "is_muts_only": {
            "label": "- Variant %ss only"%read_annotation,
            "value": is_muts_only
            },
        "is_frg_label": {
            "label": "- Add source coordinates...",
            "value": is_frg_label
            },
        "is_use_absolute": {
            "label": "- ... plus genome coordinates",
            "value": is_use_absolute
            },
        "is_fastacigar_out": {
            "label": "- CIGAR annotation",
            "value": is_fastacigar_out
            },
        "is_vars_to_lower": {
            "label": "- Substitutions in %s case"%case_subs,
            "value": is_vars_to_lower
            },
        "is_journal_subs": {
            "label": "- Journal the substitutions",
            "value": is_journal_subs
            },
        "is_fastq_out": {
            "label": "%ss in FASTQ format"%read_annotation.title(),
            "value": is_fastq_out
            },
        "Qualmin": {
            "label": "- FASTQ quality min",
            "value": Qualmin,
            "min": QualMIN,
            "max": QualMAX
            },
        "Qualmax": {
            "label": "- FASTQ quality max",
            "value": Qualmax,
            "min": QualMIN,
            "max": QualMAX
            },
        "is_write_ref_fasta": {
            "label": "Save %s"%(reference_haplotype),
            "value": is_write_ref_fasta
            },
        "is_mut_out": {
            "label": "Save %ss"%(variants_label),
            "value": is_mut_out
            },
        "is_write_ref_ingb": {
            "label": "Save Source Features",
            "value": is_write_ref_ingb
            },
        "is_sam_out": {
            "label": "%ss in SAM format"%read_annotation.title(),
            "value": is_sam_out
            },
        "gauss_mean": {
            "label":"Mean insert size",
            "value":gauss_mean,
            "min":gauss_mean_Min,
            "max":gauss_mean_Max
            },
        "gauss_SD": {
            "label":"SD insert size",
            "value":gauss_SD,
            "min":gauss_SD_Min,
            "max":gauss_SD_Max
            }
        }

    ReadsList=[]
    ReadsList.append(bio_parameters["is_frg_paired_end"]["label"])
    ReadsList.append(bio_parameters["is_duplex"]["label"])
    ReadsList.append(bio_parameters["is_simplex"]["label"])
# end of make_bio_parameters_configs()


def set_Reference_sequences_configs_public3(custtext):
    global CustomerIDText,DatasetIDText,Reference_sequences,MAINVER,DATEVER
    global GeneList,GRCH_dataset,MAINFILE
    exists,stuff=config_file_in(config_file_reference_seqs)
    print("Setting up for CustomerIDText==%s"%CustomerIDText)
    if exists:
        Reference_sequences=stuff["Reference_sequences"]
        set_dynamic_vars_defaults()
        # Get the GRCH_dataset version from loci.json file
        region=Reference_sequences[GeneList[0]]["Region"]
        GRCH_dataset=region.split(":")[0]
        DatasetIDText="Demo %s"%GRCH_dataset
        main_exists,line2,line3=main_file_read(MAINFILE)
        # Get the MAINVER and DATEVER directly from the main python code 
        if main_exists:
            date=line3.split('"')[1]
            DATEVER=date.split("-")[1]+" "+date.split("-")[2]
            
            main=line2.split('"')[1]
            main2=main.split('.')[0]
            MAINVER=main2.split('_')[3]+'_'+main2.split('_')[4]
            CustomerIDText="%s v.%s, %s"%(CustomerIDText,MAINVER,DATEVER)
        else:
            print("Fail: %s does not exist"%MAINFILE)
            sys.exit() 
    else:
        print("Fail: %s does not exist"%config_file_reference_seqs)
        sys.exit()

def set_vars_seq_limits_defaults():
    # ====================================================================
    # Initialising gui-configurable variables for sequence output features
    # These are all user-changeable variables
    # ====================================================================
    ''' def initialise_seq_limits_configs '''

    global Qualmin,Qualmax  # ''' Initialise the range of FASTQ quality values that are set at random in get_quality_list '''
                            # NB: QualMIN,QualMAX are the absolute limits for when choosing
    #Qualmin=35
    #Qualmax=50

    Qualmin=42 # Default set to same value as Qualmax to make default faster than 
    Qualmax=42

    global Exome_extend  # If is_make_exome, this value enables extension of the join region by this value.
    Exome_extend=0 # This is intended to have eg +3 nt and -3 nt each side of an intron/exon boundary within the "exome" sequence

    global Fraglen,Fragdepth
    Fraglen=200 # Set a value for the length of a fragment eg: Fraglen=50
    Fragdepth=3 # Set a value for the depth-of-cover eg: Fragdepth=50 # NB: related variable Fragcount

     # For a quick runtime to minimise output size for testing purposes, use these:
     #Fraglen=8
     #Fragdepth=1

    #''' Long catch for out-of-range initialised values
    global FraglenMin,FraglenMax
    if Fraglen > FraglenMax:
        Fraglen=FraglenMax
    if Fraglen < FraglenMin:
        Fraglen = FraglenMin

    if Fragdepth > FragdepthMax:
       Fragdepth = FragdepthMax
    if Fragdepth < FragdepthMin:
       Fragdepth = FragdepthMin

    global mutfreqs  # Initialising the frequencies of fragments by source in mutlabels
    mutfreqs=[0,50]

    # To prevent out-of-range errors, auto-extend the mutfreqs array to the same length as mutlabels by the value in last position
    # COD 6th Mar-2021: looks vestigial as this test is also done in main.mutfreqs_extend.But, placed here, it supports pre-setting of different frequencies in config.json
    while len(mutfreqs) < len(mutlabels):
        mutfreqs.append(mutfreqs[-1])

    global Paired_insert_Min,Paired_insert_Max
    Paired_insert_Min=40
    Paired_insert_Max=850
    
    global gauss_mean,gauss_SD, gauss_mean_Max,gauss_mean_Min,gauss_SD_Max,gauss_SD_Min
    gauss_mean=200
    gauss_SD=2
    gauss_mean_Min=Paired_insert_Min
    gauss_mean_Max=Paired_insert_Max

    if gauss_mean > Paired_insert_Max:
        gauss_mean=Paired_insert_Max
    elif gauss_mean < Paired_insert_Min:
        gauss_mean=Paired_insert_Min
        
    gauss_SD_Min=0
    gauss_SD_Max=20

    
#end of set_vars_seq_limits_defaults

def set_dynamic_vars_defaults():
    # =================================================================
    # Initialising dynamic variables
    # =================================================================
    global GeneList,Reference_sequences,is_CDS,CDS_trigger,mRNA_trigger,splice_trigger,splice_trigger
    GeneList=[]
    for locus in Reference_sequences:
        GeneList.append(locus)
    
    global target_locus,target_transcript_name,target_transcript_id,empty_transcript_name,empty_transcript_id
    target_locus=GeneList[0] # set the target locus variable initial value
    target_transcript_name=empty_transcript_name
    target_transcript_id=empty_transcript_id

    if is_CDS:
        splice_trigger=CDS_trigger
    else:
        splice_trigger=mRNA_trigger

    # Python testing outside Python GUI
    #target_locus="ATM"
    #target_transcript_name="ATM-201"
    #target_transcript_id="ENST00000278616"

    #target_locus="KRAS"
    #target_transcript_name="KRAS-201"
    #target_transcript_id="ENST00000256078"

    # end of set_dynamic_vars_defaults    

# end of set_dynamic_vars_defaults


# =================================================================
# Main
# =================================================================
if __name__ == "__main__":

    global is_use_TKinter,GRCH_dataset,infilepathroot,outfilepathroot,pygui_outfilepathroot
    global exploder_root,rootdatadir,rootapplicationdir,roothelpscriptdir
    global MAINVER,DATEVER
    set_config_consts()

    if GRCH_dataset == "GRCh38":
        infilepathroot = "%s/exploder_input_38_1000/"%rootdatadir
        outfilepathroot = "%s/exploder_output_38_1000/"%rootdatadir
    elif GRCH_dataset == "GRCh37":
        infilepathroot = "%s/exploder_input_37_1000/"%rootdatadir
        outfilepathroot = "%s/exploder_output_37_1000/"%rootdatadir
    else:
        print(" GRCH_dataset is unset or incorrect - retry!")
        exit()
    
    pygui_outfilepathroot=outfilepathroot  
    print("\nCreating %s for genome build version %s\n"%(config_file_output,GRCH_dataset))
    process_file_configs_python()
    print(" \nDid you set the following correctly before running?:\n 'exploder_root':%s\n 'is_use_TKinter':%s\n 'GRCh_dataset':%s\n 'CustomerIDText':%s\n'DatasetIDText': %s\n 'MAINVER':%s\n 'DATEVER':%s\n 'is_exome_paired_end':%s\n 'is_pair_monitor':%s\n 'is_reads_upper':%s"
          %(exploder_root,is_use_TKinter,GRCH_dataset,CustomerIDText,DatasetIDText,MAINVER,DATEVER,is_exome_paired_end,is_pair_monitor,is_reads_upper))
