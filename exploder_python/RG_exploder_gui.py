#!/usr/bin/python3
Progver="RG_exploder_gui_19.py"
ProgverDate="25-Aug-2024"
'''
© author: Cary O'Donnell for Replicon Genetics 2020, 2021, 2022, 2023, 2024

    The code in the section "Could be in another module" were once part of this GUI code,
    but moved to RG_main when these sections starting calling RG_main objects.
    
    When it became clear that much of these blocks of code added no value to RG_main as part of the App.vue / Vue.js GUI version, they
    were moved out of RG_main back here. This to lighten the code burden on RG_main

    The 'build' function is the original developer's name for the GUI function allowing the addition of further haplotypes,
    created after the initial 'explode' function ie: the fragmentation of sequence into reads

'''
import re
import RG_exploder_io as RG_io  # File input/output
import RG_exploder_globals as RG_globals  # Configurable constants and preset defaults for user-configurable variables
import RG_exploder_main as RG_main
import RG_exploder_process as RG_process
import RG_exploder_builder as RG_builder

import tkinter as tk
from tkinter import font
from tkinter import ttk
from PIL import ImageTk, Image
import webbrowser

def initialise_stuff():
   #initialise_global_vars() #  Later dev explicitly setting variables to RG_globals made practically all of the definitions in this function redundant
    initialise_modulevalues()
    return
# End of initialise_stuff()

# =================   GUI module setups =========================
def initialise_modulevalues():
    # Globals limited to this module
    global Progver,ProgverDate
    global local_begin,local_end,trans_Begin,trans_End,abs_Begin,abs_End
    global Chrom_txt
    local_begin=RG_globals.bio_parameters["target_build_variant"]["local_begin"]
    local_end=RG_globals.bio_parameters["target_build_variant"]["local_end"]
    abs_Begin=RG_globals.bio_parameters["target_build_variant"]["abs_Begin"]["value"]
    abs_End=RG_globals.bio_parameters["target_build_variant"]["abs_End"]["value"]
    trans_Begin=RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["value"]
    trans_End=RG_globals.bio_parameters["target_build_variant"]["trans_End"]["value"]
    Chrom_txt="X"

    global run_count
    run_count=0 # Count how many runs
    
    global IUPAC_sub
    IUPAC_sub=RG_globals.IUPAC_codes[:1]+"^"+RG_globals.IUPAC_codes[1:] # Insert a "^" to get a re.sub string for parsing pasted text. Bizarre hack
    global previous_target_locus, previous_target_transcript_name
    previous_target_locus=""
    previous_target_transcript_name=RG_globals.empty_transcript_name


    global ENS_target_url,ENS_target_url_txt,LRG_target_url,LRG_target_url_txt,ENS_ts_target_url,ENS_ts_target_url_txt
    ENS_target_url="" #  variable assigned on reading Refseq in python GUI
    ENS_target_url_txt="" #  variable assigned on reading Refseq in python GUI
    LRG_target_url="" #  variable assigned on reading Refseq in python GUI
    LRG_target_url_txt="" #  variable assigned on reading Refseq in python GUI
    ENS_ts_target_url=""
    ENS_ts_target_url_txt=""

    global icondir
    icondir="./icons/"

    global results_links,results_idx
    results_links=['']
    results_idx=[]

    global flags_output_labels,flags_output_bools
    global config_limit_labels,config_limit_vals,config_limit_mins,config_limit_maxs

    ''' Use below for testing set_tog_bools(test_labels,test_booleans)
    test_labels=["One","Two","Three"]
    test_booleans=[False,True,False]
    '''
    set_flags_output_labels()
    set_flags_output_bools()
    set_flags_dict()
    set_tog_bools(flags_output_labels,flags_output_bools)
    set_tog_bools2(flags_output_labels2,flags_output_bools2)

    set_config_limits()
    set_spin_vals(config_limit_vals,config_limit_labels,config_limit_mins,config_limit_maxs)
    
    set_src_slider_values_build([50,50],["Begin","End"],[0,0]) # New for builder
    init_set_src_slider_values(18)# Set number of frequency sliders to a value that doesn't obscure messages panel # Hidden for builder

    global max_seqlength,min_seqlength,max_ext_add,min_ext_add
    max_seqlength=RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["max"]
    min_seqlength=RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["min"]
    max_ext_add=RG_globals.bio_parameters["target_build_variant"]["trans_Begin_ext"]["max"]
    min_ext_add=RG_globals.bio_parameters["target_build_variant"]["trans_Begin_ext"]["min"]
    
    return
# End of def initialise_modulevalues()    

def set_flags_dict():
    return # Available for testing
    global flags_output_labels,flags_output_bools # externally pre-defined
    global flags_dict
    flags_dict={}
    i=0
    for label in flags_output_labels:
        flags_dict[label]=flags_output_bools[i]
        i+=1
    # print("flags_dict: %s"%flags_dict)

def set_flags_output_labels():
    global flags_output_labels,flag_widget_order,flags_output_labels2,flag_widget_order2

   # Set values in flag_widget_order to set the order-positions in flags_output_labels and flags_output_bools
   # this becomes the order in which "Output options" is presented for the boolean widgete
    flag_widget_order={#"is_duplex":1,            # Duplex reads
                       #"is_frg_paired_end":0,    #  Paired end reads
                       "is_flip_strand":0,
                       "is_fasta_out":1,         #  Reads in FASTA format
                       "is_onefrag_out":2,       #  Each possible read
                       "is_muts_only":3,         # Variant reads only
                       "is_frg_label":4,         # Annotate source positions ...
                       "is_use_absolute":5,      # ...plus absolute position
                       "is_fastacigar_out":6,    # CIGAR annotation
                       "is_vars_to_lower":7,     # Substitutions in lower case
                       "is_journal_subs":8,      # Journal the substitutions - Verbose documentation of variant substitutions
                       "is_fastq_out":9,         # Reads in FASTQ format
                       "is_write_ref_fasta":10,  # Save Template Sequence
                       "is_mut_out":11,          # Save Haplotype Sequences
                       "is_write_ref_ingb":12,   # Save Source features
                       "is_sam_out":13          #  Reads in SAM format
                       }
   
    flags_output_labels=[""] * len(flag_widget_order)

    for item in flag_widget_order:
        flags_output_labels[flag_widget_order[item]]=RG_globals.bio_parameters[item]["label"]

    flag_widget_order2={"is_CDS":0,            # 
                       }

    flags_output_labels2=[""] * len(flag_widget_order2)

    for item in flag_widget_order2:
        flags_output_labels2[flag_widget_order2[item]]=RG_globals.bio_parameters[item]["label"]


def set_flags_output_bools():
    global flags_output_bools,flag_widget_order,flags_output_bools2,flag_widget_order2    
    flags_output_bools=[""] * len(flag_widget_order)

    for item in flag_widget_order:
        flags_output_bools[flag_widget_order[item]]=RG_globals.bio_parameters[item]["value"]

    flags_output_bools2=[""]* len(flag_widget_order2)

    for item in flag_widget_order2:
        flags_output_bools2[flag_widget_order2[item]]=RG_globals.bio_parameters[item]["value"]

    #print("flags_output_bools2 %s \n"%flags_output_bools2)


def set_config_limits():
    # Not necessary to set all these as globals if incorporate set_spin_vals in this routine.
    # Done here for consistency with the set_flags routines
    global config_limit_labels,config_limit_vals,config_limit_mins,config_limit_maxs
    config_limit_labels=[RG_globals.bio_parameters["Exome_extend"]["label"],
                       RG_globals.bio_parameters["Fraglen"]["label"],
                       RG_globals.bio_parameters["Fragdepth"]["label"],
                       RG_globals.bio_parameters["Qualmin"]["label"],
                       RG_globals.bio_parameters["Qualmax"]["label"],
                       RG_globals.bio_parameters["gauss_mean"]["label"],
                       RG_globals.bio_parameters["gauss_SD"]["label"]
                         ]             
    config_limit_vals=[RG_globals.Exome_extend,RG_globals.Fraglen,RG_globals.Fragdepth,RG_globals.Qualmin,RG_globals.Qualmax,RG_globals.gauss_mean,RG_globals.gauss_SD]
    config_limit_mins=[RG_globals.Exome_extend_Min,RG_globals.FraglenMin,RG_globals.FragdepthMin,RG_globals.QualMIN,RG_globals.QualMIN,RG_globals.gauss_mean_Min,RG_globals.gauss_SD_Min]
    config_limit_maxs=[RG_globals.Exome_extend_Max,RG_globals.FraglenMax,RG_globals.FragdepthMax,RG_globals.QualMAX,RG_globals.QualMAX,RG_globals.gauss_mean_Max,RG_globals.gauss_SD_Max]

def set_tog_bools(in_labels,in_booleans):
    global  tog_labels,tog_booleans
    tog_labels=in_labels
    tog_booleans=[]
    for val in in_booleans:
        tog_booleans.append(val)
    #print("A0:tog_booleans %s"%tog_booleans)


def set_tog_bools2(in_labels,in_booleans):
    global  tog_labels2,tog_booleans2
    tog_labels2=in_labels
    tog_booleans2=[]
    for val in in_booleans:
        tog_booleans2.append(val)
    #print("A0:tog_booleans2 %s"%tog_booleans2)


def set_src_slider_values_build(invalues,inlabels,invalues2):
    global src_sliders_vals_build,src_sliders_strvals_build,src_sliders_labels_build,src_sliders_straddvals_build
    src_sliders_vals_build=[]
    src_sliders_strvals_build=[]
    src_sliders_labels_build=[]

    src_sliders_straddvals_build=[]

    count=0
    for v in invalues:
        src_sliders_vals_build.append(tk.IntVar(value=v))
        src_sliders_strvals_build.append(tk.StringVar(value=str(v)))
        src_sliders_labels_build.append(inlabels[count])
        count+=1
        
    for v in invalues2:
        src_sliders_straddvals_build.append(tk.StringVar(value=str(v)))


def init_set_src_slider_values(length): # Not used in builder
    global src_sliders_init_vals
    global src_sliders_vals,src_sliders_strvals,src_sliders_labels
    src_sliders_vals=[]
    src_sliders_init_vals=[]
    src_sliders_strvals=[]
    src_sliders_labels=[]

    preval=-1
    for v in range(length):
        src_sliders_vals.append(tk.IntVar(value=preval))
        src_sliders_init_vals.append(tk.IntVar(value=preval))
        src_sliders_strvals.append(tk.StringVar(value=str(preval)))
        src_sliders_labels.append("hap%s"%v)

def set_spin_vals(in_values,in_labels,mins,maxes):
    global spinvals,spin_labels,spin_mins,spin_maxs
    #global spin_safe_vals
    spinvals=[]
    #spin_safe_vals=[]
    spin_labels=[]
    spin_mins=mins
    spin_maxs=maxes
    count=0
    for v in in_labels:
        spinvals.append(tk.StringVar(value=in_values[count]))
        #spin_safe_vals.append(tk.StringVar(value=in_values[count]))
        spin_labels.append(v)
        count+=1

def to_bool(inflag):
    # If these are text instead of Boolean, resets to Boolean
    outflag={'True': True, 'False': False, True: True, False: False}.get(inflag)
    #print("flag_to_bool: inflag %s , outflag %s"%(inflag,outflag))
    return outflag

def flags_to_bools(inlist):
    i=0
    for item in inlist:
        inlist[i]=to_bool(item)
        i+=1
    return inlist

def get_flags_output_bools():
    global flags_output_bools,tog_booleans,flag_widget_order,tog_booleans2,flag_widget_order2
    #print("A:flags_output_bools %s"%flags_output_bools)
    #print("A1:tog_booleans %s"%tog_booleans)

    #print("GUI flags1: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)
    tog_booleans=flags_to_bools(tog_booleans)
    tog_booleans2=flags_to_bools(tog_booleans2)
    flags_output_bools=tog_booleans
    #print("1:flags_output_bools %s"%flags_output_bools)
    #print("1:tog_booleans       %s"%tog_booleans)

    RG_globals.is_flip_strand=tog_booleans[flag_widget_order["is_flip_strand"]]
    RG_globals.is_fasta_out=tog_booleans[flag_widget_order["is_fasta_out"]]
    RG_globals.is_fastq_out=tog_booleans[flag_widget_order["is_fastq_out"]]
    RG_globals.is_sam_out=tog_booleans[flag_widget_order["is_sam_out"]]

    RG_globals.is_write_ref_fasta=tog_booleans[flag_widget_order["is_write_ref_fasta"]]
    RG_globals.is_mut_out=tog_booleans[flag_widget_order["is_mut_out"]]
    #print("GUI flags2: RG_globals.is_mut_out %s; flags_output_bools[4]%s"%(RG_globals.is_mut_out,flags_output_bools[4]))
    RG_globals.is_write_ref_ingb=tog_booleans[flag_widget_order["is_write_ref_ingb"]]
    RG_globals.is_onefrag_out=tog_booleans[flag_widget_order["is_onefrag_out"]]
    RG_globals.is_muts_only=tog_booleans[flag_widget_order["is_muts_only"]]
    ##RG_globals.is_frg_paired_end=tog_booleans[flag_widget_order["is_frg_paired_end"]]
    ##RG_globals.is_duplex=tog_booleans[flag_widget_order["is_duplex"]]
    RG_globals.is_frg_label=tog_booleans[flag_widget_order["is_frg_label"]]
    RG_globals.is_use_absolute=tog_booleans[flag_widget_order["is_use_absolute"]]
    RG_globals.is_fastacigar_out=tog_booleans[flag_widget_order["is_fastacigar_out"]]
    RG_globals.is_vars_to_lower=tog_booleans[flag_widget_order["is_vars_to_lower"]]
    RG_globals.is_journal_subs=tog_booleans[flag_widget_order["is_journal_subs"]]
    RG_globals.is_CDS=tog_booleans2[flag_widget_order2["is_CDS"]]

    set_flags_dict()

def get_config_limits():
    # AKA: get_spin_vals
    global spinvals,config_limit_vals
    for i in range(len(spinvals)):
        #print("1: label %s, value %s"%(config_limit_labels[i],config_limit_vals[i]))
        #print("2: label %s, value %s"%(spin_labels[i],spinvals[i].get()))
        config_limit_vals[i]=int(float(spinvals[i].get()))


    RG_globals.Exome_extend=config_limit_vals[0]
    RG_globals.Fraglen=config_limit_vals[1]
    RG_globals.Fragdepth=config_limit_vals[2]
    RG_globals.Qualmin=config_limit_vals[3]
    RG_globals.Qualmax=config_limit_vals[4]
    RG_globals.gauss_mean=config_limit_vals[5]
    RG_globals.gauss_SD=config_limit_vals[6]

    # Checks if the min & max require inversion
    if RG_globals.Qualmax < RG_globals.Qualmin:
        RG_globals.Qualmin,RG_globals.Qualmax=RG_globals.Qualmax,RG_globals.Qualmin
        config_limit_vals[3],config_limit_vals[4]=config_limit_vals[4],config_limit_vals[3]
        correct_spin_vals() #see comments in object
        #refresh_gui() # Has no effect
    return

def correct_spin_vals():
    global previous_target_locus,previous_target_transcript_name
    # Sledgehammer to re-set everything where the Qualmin/max are swapped
    # Called conditionally only after a swap because the frame-refresh caused a screen-flicker when build absent
    set_config_limits()
    set_spin_vals(config_limit_vals,config_limit_labels,config_limit_mins,config_limit_maxs)
    # Now update the spinboxes to show the correction by refreshing frame they are in
    # However: since addition of the build frame the layout gets rearranged. Tried fixing ... no luck
    #refresh_main_options_list() # Need to change because it messes up the build layout.

def get_mutfreqs():
    # Set from src_sliders_init_vals (integers)
    #       or src_sliders_vals (tk.IntVar)
    #       or src_sliders_strvals (tk.StringVar)
    # src_sliders_init_vals would be user-saved, but the others user-set
    # could call save_src_sliders_vals() to force a save to src_sliders_init_vals, but would need a version without user-prompt
    # which to force?
    global src_sliders_init_vals,src_sliders_vals,src_sliders_strvals

    count=0
    for i in RG_globals.mutfreqs:
         # For user-saved
         #RG_globals.mutfreqs[count]=src_sliders_init_vals[count]
         # For user-current-selection
         #RG_globals.mutfreqs[count]=int(src_sliders_vals[count].get())
         # or
         RG_globals.mutfreqs[count]=int(src_sliders_strvals[count].get())
         count += 1
    #print("        src_sliders_init_vals: %s"%  src_sliders_init_vals)
    #print("RG_globals.mutfreqs: %s"%RG_globals.mutfreqs)
    return

def save_and_leave():
    #print("GUI save_and_leave1: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)
    get_mutfreqs()
    #print("GUI save_and_leave2: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)
    get_flags_output_bools()
    #print("GUI save_and_leave3: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)
    get_config_limits()
    return

def increment_run():
    global run_count,results_links, results_idx
    run_count+=1
    thisfilepath=RG_globals.pygui_outfilepathroot+RG_globals.target_locus+"/" # pygui_outfilepathroot as "output/" not suitable for this, must be full path
    #results_links.append("file://%s"%RG_globals.outfilepathroot)
    i=0 ; is_new = True
    for item in results_links:
        if item == thisfilepath:
            results_idx.append(i)
            is_new = False
            break
        i+=1

    if is_new:
        results_links.append(thisfilepath)
        results_idx.append(len(results_links)-1)

def set_refseq_target_url():
    # Call (eg: from GUI) to set URL without having read the reference file, by using config.json contents
    RG_globals.ensembl_geneid=RG_globals.Reference_sequences[RG_globals.target_locus]["Ensembl_id"]
    geneid=RG_globals.ensembl_geneid
    LRG_id=RG_globals.Reference_sequences[RG_globals.target_locus]["LRG_id"]
    
    global ENS_target_url,ENS_target_url_txt,LRG_target_url,LRG_target_url_txt,ENS_ts_target_url,ENS_ts_target_url_txt
    global Chrom_txt
        
    ENS_target_url=""
    ENS_target_url_txt=""
    LRG_target_url=""
    LRG_target_url_txt=""
    ENS_ts_target_url=""
    ENS_ts_target_url_txt=""

    GRChver_txt=RG_globals.Reference_sequences[RG_globals.target_locus]["Region"].split(":")[0]
    Chrom_txt=RG_globals.Reference_sequences[RG_globals.target_locus]["Region"].split(":")[1]
    # Special catch for non-standard chromosome identifiers such as PTEN_a
    if len(Chrom_txt)>2:
        Chrom_txt=Chrom_txt[:3]
    
    if GRChver_txt== RG_globals.GRCh37_txt:
        target="37"
    elif GRChver_txt== RG_globals.GRCh38_txt:
        target="38"
    else:
        target=""

    if LRG_id != "":
        #LRG_target_url=RG_globals.lrg_view_url+LRG_id+".xml"
        LRG_target_url=RG_globals.lrg_view_url+LRG_id
        #LRG_target_url_txt=RG_globals.bio_parameters["target_locus"]["label"]+" "+RG_globals.target_locus+": "+LRG_id
        LRG_target_url_txt="LRG "+RG_globals.target_locus+": "+LRG_id
        #print("LRG_target_url: %s"%LRG_target_url)
            
    if "." in geneid:
        geneid,version=RG_globals.ensembl_geneid.split(".")
        if target=="37":
            ENS_target_url=RG_globals.ensembl37_gene_url+"g="+geneid
            #ENS_target_url_txt="Ensembl "+RG_globals.target_locus+": "+RG_globals.GRCh37_txt+":"+geneid
            ENS_target_url_txt=RG_globals.GRCh37_txt+":"+geneid
            
        elif target=="38":
            ENS_target_url=RG_globals.ensembl38_gene_url+"g="+geneid
            #ENS_target_url_txt="Ensembl "+RG_globals.target_locus+": "+RG_globals.GRCh38_txt+":"+geneid
            ENS_target_url_txt=RG_globals.GRCh38_txt+":"+geneid
        else:
            ENS_target_url=""
            #ENS_target_url_txt="Ensembl "+RG_globals.target_locus+": "+"no build "+":"+geneid
            ENS_target_url_txt="GRCh??:"+geneid


    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
        ENS_ts_target_url=ENS_target_url
        ENS_ts_target_url_txt=ENS_target_url_txt
    else:
        tsid,version=RG_globals.Reference_sequences[RG_globals.target_locus]["mRNA"][RG_globals.target_transcript_name].split(".")
        if target=="37":
            # Special cases where transcript is not defined in build 37, is present in 38, and has been edited into the 37 source data
            # eg: in config.json AK2-672715m(MANE_Select) and in genbank file "ENST00000672715m.1 , with 'm' as the triggering identifier
            if tsid[-1:]=="m": # identify the modified id
                tsid=tsid[:-1] # chop it for the hyperlink
                target="38" #  Switch to 38
        if target=="37":
            ENS_ts_target_url=RG_globals.ensembl37_transcript_url+"g="+geneid+";t="+tsid
        elif target=="38":
            ENS_ts_target_url=RG_globals.ensembl38_transcript_url+"g="+geneid+";t="+tsid
        else:
            ENS_ts_target_url=""
        ENS_ts_target_url_txt="GRCh"+target+":"+tsid
    #print("ENS_ts_target_url %s"%ENS_ts_target_url)

def results_hyperLink(event):
    global results_links
    #print("event.widget.tag_names(tk.CURRENT)[1]: %s"%event.widget.tag_names(tk.CURRENT)[1])
    '''
    Intermittently get the following message written to the Python shell when using the above after clicking on link
       - or possibly the error is clicking just off the link, but it's not easily reproducible:

       Exception in Tkinter callback
       Traceback (most recent call last):
       File "/Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/tkinter/__init__.py", line 1883, in __call__
       return self.func(*args)
       File "/Users/caryodonnell/Desktop/Replicon/python_exploder_scripts/RG_exploder_gui_06.py", line 392, in results_hyperLink
       idx= int(event.widget.tag_names(tk.CURRENT)[1])
       ValueError: invalid literal for int() with base 10: 'link'

    Substitututed the try: except: construct to evade this error report.
    '''
    try:
        idx= int(event.widget.tag_names(tk.CURRENT)[1])
        #print ("idx %s, results_idx[idx-1]%s, results_links[results_idx[idx-1]] %s"%(idx,results_idx[idx-1],results_links[results_idx[idx-1]]))
        webbrowser.open("file:///%s"%results_links[results_idx[idx-1]])
    except:
        # Yup - goes here instead: a silent fail without diagnostic print
        # print("results_hyperLink(event) fail")
        pass

def save_and_go():
    global run_count,icondir
    global go_img,gobutton
    wait_img=ImageTk.PhotoImage(Image.open(icondir+"wait45.png"))
    gobutton.configure(image=wait_img)
    #gobutton.configure(text="WAIT")
    gobutton.configure(state="disable")
    save_and_leave()
    #print("GUI save_and_go: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)

    increment_run()
    write_GUI_text("Running Exploder on %s- run number %s\n"%(RG_globals.target_locus,run_count))
    ##RG_globals.bio_parameters["target_build_variant"]["exonpos_lookup"]=[0] # testing that setting exonpos_lookup to empty will force a re-run in main
    is_successful=RG_main.call_exploder_main()
    time_stamp=RG_globals.getime()
    #print("is_successful %s"%is_successful)
    if is_successful:
        write_GUI_text("Run %s Complete at %s\n"%(run_count,time_stamp))
    else:
        write_GUI_text("Run %s Failed at %s\n"%(run_count,time_stamp))

    extralabel,strandlabel=RG_globals.get_strand_extra_label()
    write_GUI_results_link("%s_%s run results. Check time stamp!"%(RG_globals.get_locus_transcript(),extralabel))
        
    gobutton.configure(state="normal")
    gobutton.configure(image=go_img)
    #gobutton.configure(text="GO")

# =================  END GUI module setups =========================


# =================  GUI module setups =========================


def set_pygui_defaults010():
    global pygui_frame_labels,pygui_button_labels # Used as tk object labels
    #"Save Variant to %s"%(RG_globals.bio_parameters["target_build_variant"]["hap_name"]["label"])]
    pygui_frame_labels=["Messages Panel"]
    pygui_button_labels=["Clear messages",
                         "Retrieve Reference Sequence",
                         "Save"]
# end of set_pygui_defaults010()

def tk_background_config(scope):
    # Here in case it's needed
    return
    scope.configure(background="white")
    
def ttk_background_config(scope):
    # Here in case it's needed
    style_4 = {'fg': 'white', 'bg': 'coral2', 'activebackground':
               'gray71', 'activeforeground': 'gray71'}
    scope.configure(**style_4)


def initialise_tk():
    set_pygui_defaults010()
    global pygui_frame_labels,pygui_button_labels
    
    global exploder_win,topbar,main,bottom_left,bottom_right
    global exploder_win,topbar,main,main_left,main_vars_list,main_build,main_options_list
    exploder_win= tk.Tk() # tk only, not ttk
    tk_background_config(exploder_win)
    exploder_win.title("Synthetic Reads Generator - Python tkinter GUI (v15)") 
    # Laying out grid
    topbar=ttk.Frame(exploder_win)
    main=ttk.Frame(exploder_win)
    bottom_left=ttk.Frame(exploder_win)
    bottom_centre=ttk.Frame(exploder_win)
    bottom_right=ttk.Frame(exploder_win)
    # Trying to set background all same
    tk_background_config(main)
    tk_background_config(topbar)
    tk_background_config(main)
    tk_background_config(bottom_left)
    tk_background_config(bottom_right)

    topbar.pack(side="top", fill="x")
    main.pack(side="top", fill="both", expand=True)
    bottom_left.pack(side="left", fill="both", expand=False)
    bottom_centre.pack(side="left", fill="both", expand=False)
    bottom_right.pack(side="right", fill="both", expand=False)

    initialise_main_left(main)
    initialise_main_vars_list(main)
    initialise_main_options_list(main)
    initialise_main_build(main)

    main_left.pack(side="left", fill="both", expand=True)
    main_vars_list.pack(side="left", fill="both", expand=True)
    main_options_list.pack(side="left", fill="both", expand=True) 
    main_build.pack(side="left", fill="both", expand=True)

    #topbar.label=tk.Label(topbar, width=30,text="Topbar Panel",foreground="blue")
    #topbar.header=tk.Label(topbar,text="%s | %s"%(RG_globals.CustomerIDText,RG_globals.DatasetIDText))
    topbar.header=ttk.Label(topbar,text="%s | %s"%(RG_globals.CustomerIDText,RG_globals.DatasetIDText))
    tk_background_config(topbar.header)
    topbar.header.pack()
    #topbar.help=tk.Label(topbar,text="%s"%(RG_globals.help_label))
    topbar.help=ttk.Label(topbar,text="%s"%(RG_globals.help_label))
    topbar.help.pack()
    set_topbar_help()
    ### Next lines add an Ensembl and LRG hyperlink to top of page, but these duplicated in titles of selection boxes
    #topbar.ENS_label=tk.Label(topbar,text="Ensembl:%s"%RG_globals.target_locus)
    #topbar.ENS_label.pack()
    #topbar.LRG_label=tk.Label(topbar,text="LRG:%s"%RG_globals.target_locus)
    #topbar.LRG_label.pack()
    ### If reinstate these, need to reinstate lines in refresh_genelabels as well

    #bottom_left.label=tk.Label(bottom_left, width=30,text="Messages Panel",)
    #bottom_left.label=tk.Label(bottom_left, width=30,text=pygui_frame_labels[0])
    bottom_left.label=ttk.Label(bottom_left, width=30,text=pygui_frame_labels[0])
    bottom_left.label.grid()

    bottom_centre.label=ttk.Label(bottom_centre, width=30,text="")
    bottom_centre.label.grid()

    #bottom_right.label=tk.Label(bottom_right, width=30,text="Gene structure")
    bottom_right.label=ttk.Label(bottom_right, width=30,text="Gene structure")
    bottom_right.label.grid()

    #main_vars_list.label=tk.Label(main_vars_list, width=10,text="Centre Panel",)
    #main_vars_list.label.pack(side="top", fill="x")
    # END Laying out grid
    return
# END initialise_tk()

def set_generic_label(framelabel,hyperlink,url,urltxt):
    if urltxt != "":
        framelabel['text']=urltxt
    f2 = font.Font(framelabel, framelabel.cget("font"))
    if url != "":
        f2.configure(underline = True)
        framelabel.configure(font=f2)
        framelabel['foreground']="blue"
        framelabel['cursor']="hand2"
        framelabel.bind('<Button-1>', hyperlink)
    else:
        #NB: Underline does not disappear
        f2.configure(underline = False)
        framelabel['foreground']="black"
        framelabel['cursor']="arrow"
        framelabel.unbind('<Button-1>')

def ensembl_hyperLink(event):
    global ENS_target_url
    if ENS_target_url != "":
        webbrowser.open("%s"%ENS_target_url)
 
def set_topbar_ENS_label_link():
    set_generic_label(topbar.ENS_label,ensembl_hyperLink,ENS_target_url,ENS_target_url_txt)

def lrg_hyperLink(event):
    global LRG_target_url
    #print("LRG_target_url %s"%LRG_target_url)
    if LRG_target_url != "":
        webbrowser.open("%s"%LRG_target_url)

def set_topbar_LRG_label_link():
    set_generic_label(topbar.LRG_label,lrg_hyperLink,LRG_target_url,LRG_target_url_txt)

def set_listbox_locus_label():
    #set_generic_label(main_left.genelistlabel,lrg_hyperLink,LRG_target_url,"") # LRG now
    set_generic_label(main_left.genelistlabel,ensembl_hyperLink,ENS_target_url,ENS_target_url_txt)    

def ensembl_ts_hyperLink(event):
    global ENS_ts_target_url
    if ENS_ts_target_url != "":
        webbrowser.open("%s"%ENS_ts_target_url)

def set_listbox_transcript_label():
    #set_generic_label(main_left.tslistlabel,ensembl_ts_hyperLink,ENS_ts_target_url,"")
    set_generic_label(main_left.tslistlabel,ensembl_ts_hyperLink,ENS_ts_target_url,ENS_ts_target_url_txt)

def help_hyperLink(event):
    if RG_globals.help_url != "":
        webbrowser.open("%s"%RG_globals.help_url)
        
def set_topbar_help():
    set_generic_label(topbar.help,help_hyperLink,RG_globals.help_url,"")

def initialise_main_left(frame):
    global main_left
    #main_left=tk.LabelFrame(main,text="Main_Left",padx=10,pady=10)
    main_left=ttk.LabelFrame(main,text="Main_Left")
    #tk_background_config(main_left)
    set_gene_list(main_left)

def initialise_main_options_list(frame):
    global main_options_list
    #main_options_list=tk.LabelFrame(frame,text="   Output options",padx=10,pady=10)
    #main_options_list=tk.LabelFrame(frame,text=RG_globals.options_label,padx=10,pady=10)
    main_options_list=ttk.LabelFrame(frame,text=RG_globals.options_label)
    #tk_background_config(main_options_list)

def initialise_main_vars_list(frame):
    global main_vars_list
    #main_vars_list=tk.LabelFrame(main,padx=10,pady=10) # Compatible with font setting in class source_sliders
    main_vars_list=ttk.LabelFrame(main) # Not compatible with font setting in class source_sliders
    #tk_background_config(main_vars_list)
    set_slider_gene_label()
    

def initialise_main_build(frame):
    global main_build
    #main_build=tk.LabelFrame(main,padx=10,pady=10)
    #main_build=tk.LabelFrame(main)
    main_build=ttk.LabelFrame(main)
    #tk_background_config(main_build)
    set_slider_build_label()
    #panel_build = tk.Frame(main_build)
    panel_build = ttk.Frame(main_build)
    #tk_background_config(panel_build)
    panel_build.grid()
    #panel_build.label=tk.Label(panel_build,text="    Source             Local_pos       Global_pos     Extension")
    panel_build.label=ttk.Label(panel_build,text=" Source               Local_coord  Extension  Genome_coord")
    panel_build.label.grid(row=0,column=0,sticky=tk.W)

def set_slider_gene_label():
    set_gene_label(main_vars_list)

def set_slider_build_label():
    set_builder_label(main_build)

def set_gene_label(frame):
    #frame['text']="  %s     %s"%(RG_globals.variants_header,RG_globals.frequency_label)
    frame['text']=" %ss"%RG_globals.variants_label
    frame_label=ttk.Label(frame,text=" %s             %s"%(RG_globals.bio_parameters["target_build_variant"]["hap_name"]["label"],RG_globals.frequency_label))
    frame_label.grid()
    #frame['text']=" Source Name             %s"%(RG_globals.frequency_label)
    
def set_builder_label(frame):
    #frame['text']="               Label         Position                 Extension         Global_mapping"  # Addition in builder
    frame['text']=" Create a new %s %s"%(RG_globals.target_locus,RG_globals.variants_label)  # Addition in builder

def set_gene_list(frame):
    #frame['text']="Gene List:%s"%RG_globals.target_locus
    #frame['text']=" Reference Gene"
    #frame['text']=pygui_frame_labels[1]
    frame['text']=RG_globals.reference_gene

# Source label/slider definitions
class source_sliders:
    def __init__(self, master1,v_index,**kwargs):
        global src_sliders_vals,src_sliders_strvals,src_sliders_labels
        #self.panel2 = tk.Frame(master1)
        self.panel2 = ttk.Frame(master1)
        self.panel2.grid()
        self.label=ttk.Label(self.panel2, width=15,text="%s_%s"%(RG_globals.target_locus,src_sliders_labels[v_index]))
        self.label.grid(row=0,column=0,sticky=tk.W)
        #f2 = font.Font(master1, master1.cget("font")) # Looks superfluous here
        #self.label.configure(font=f2)# Not supported ttk.LabelFrame, so looks superfluous
        ####### Extra to make this a hyperlink to the source file
        # However unless this is made READ only, implementing this allows the user to edit the file before a run.
        # This is Dangerous because it will change the source file without verification of constraint to correct number range
        # Best way to manage this is to use the 'Create a New Variations Source' feature aka 'Builder' 
        #f2.configure(underline = True)
        #self.label['fg']="blue"
        #self.label['cursor']="hand2"
        #self.label.bind('<Button-1>', lambda e: self.haplotype_hyperLink(self.label['text']))
        #######

        self.entryvar=src_sliders_strvals[v_index]
        #self.entry = tk.Entry(self.panel2, textvariable=self.entryvar, width=3)
        self.entry = ttk.Entry(self.panel2, textvariable=self.entryvar, width=3)
        self.entry.grid(row=0,column=1,sticky=tk.W)

        self.blank=''
        self.old_value = src_sliders_strvals[v_index].get()
        self.zero='0'
        self.entryvar.trace('w', self.checkentry)
        self.get, self.set = self.entryvar.get, self.entryvar.set

        self.scalevar=src_sliders_vals[v_index]
        #self.scale2 = tk.Scale(self.panel2,variable = self.scalevar,orient=tk.HORIZONTAL,sliderlength=25,**kwargs)
        self.scale2 = ttk.Scale(self.panel2,variable = self.scalevar,orient=tk.HORIZONTAL,**kwargs)
        self.scale2.grid(row=0,column=2,sticky=tk.W)
        self.scalevar.trace('w', self.setentryfromscale)

    def checkentry(self, *args):
       if self.get().isdigit():
          # the current value is only digits; allow this
          self.old_value = self.get()
          self.scalevar.set(int(self.get()))
          self.entryvar.set(str(self.scalevar.get()))
       elif self.get() == self.blank:
         # the current value is blank - set to zero
          self.old_value = self.zero
          self.scalevar.set(0)
          # self.entryvar.set(str(self.scalevar.get()))
       else:
          # there's non-digit characters in the input; reject this
          self.set(self.old_value)

    def setentryfromscale(self, *args):
       self.entryvar.set(str(self.scalevar.get()))


    def haplotype_hyperLink(event,w_txt):
        try:
            #idx= int(event.widget.tag_names(tk.CURRENT)[1])
            #print ("idx %s, results_idx[idx-1]%s, results_links[results_idx[idx-1]] %s"%(idx,results_idx[idx-1],results_links[results_idx[idx-1]]))
            in_file="%s%s/%s%s"%(RG_globals.infilepathroot,RG_globals.target_locus,w_txt,RG_globals.Seq_IO_file_ext)
            file_exists = RG_io.is_file(in_file)
            if file_exists:
                webbrowser.open("file:///%s"%in_file)
            else:
                clear_GUI_text()
                write_GUI_text("File does not exist - this will occur after creating a new %s"%RG_globals.variants_label)
        except:
            # Yup - goes here instead: a silent fail without diagnostic print
            # print("results_hyperLink(event) fail")
            pass

class source_sliders_builder:# The Locus_Begin / Locus_End entry panels
    def __init__(self, master1,v_index,**kwargs):
        global src_sliders_vals_build,src_sliders_strvals_build,src_sliders_labels_build
        global src_sliders_straddvals_build,max_seqlength,min_seqlength,max_ext_add,min_ext_add
        global Chrom_txt
        #self.panel2 = tk.Frame(master1)
        self.panel2 = ttk.Frame(master1)
        self.panel2.grid()
        self.label=ttk.Label(self.panel2, width=11,text="%s"%(src_sliders_labels_build[v_index]))
        self.label.grid(row=0,column=0,sticky=tk.W)
        
        self.entryvar=src_sliders_strvals_build[v_index]
        #self.entry = tk.Spinbox(self.panel2, textvariable=self.entryvar, width=6,from_=min_seqlength,to=max_seqlength)
        self.entry = ttk.Spinbox(self.panel2, textvariable=self.entryvar, width=6,from_=min_seqlength,to=max_seqlength)
        self.entry.grid(row=0,column=1,sticky=tk.W)

        self.abspos=0
        self.localpos=0
        
        self.blank=''
        self.old_value = src_sliders_strvals_build[v_index].get()
        self.zero='1'
        self.entryvar.trace('w', self.checkentry)
        self.get, self.set = self.entryvar.get, self.entryvar.set

        self.scalevar=src_sliders_vals_build[v_index]
        self.scalevar.trace('w', self.checkentry2)
        
        self.entry2var=src_sliders_straddvals_build[v_index]
        #self.entry2 = tk.Spinbox(self.panel2, textvariable=self.entry2var, width=4,from_=min_ext_add,to=max_ext_add,state='readonly')
        self.entry2 = ttk.Spinbox(self.panel2, textvariable=self.entry2var, width=4,from_=min_ext_add,to=max_ext_add,state='readonly')
        self.entry2.grid(row=0,column=2,sticky=tk.W)

        self.entry2var.trace('w', self.checkentry2)# Adds any change here to the absolute / genome coord

        #self.complement_check(self)
        self.complement_check() # Extra in builder
        
        #self.abslabel=tk.Label(self.panel2, width=10,text=str(self.abspos))
        self.abslabel=ttk.Label(self.panel2, width=11,text=str(self.abspos))
        self.abslabel.grid(row=0,column=3,sticky=tk.W)

    def checkentry(self, *args):
       if self.get().isdigit():
          # the current value is only digits; allow this
          self.old_value = self.get()
          self.scalevar.set(int(self.get()))
          self.entryvar.set(str(self.scalevar.get()))
          if int(self.entry.get()) > max_seqlength:
                 self.entryvar.set(str(max_seqlength))
                 self.set(str(max_seqlength))
       elif self.get() == self.blank:
         # the current value is blank - set to zero (or minimum)
          self.old_value = self.zero
          self.scalevar.set(1)
          self.set(1)
          #self.entryvar.set(str(self.scalevar.get()))
       else:
          # there's non-digit characters in the input; reject this
          self.set(self.old_value)
        
    def checkentry2(self, *args):
        global max_seqlength,min_seqlength
        self.entry.config(from_=min_seqlength, to=max_seqlength)
        self.entryvar.set(str(self.scalevar.get()))
        #self.complement_check(self)
        self.complement_check() # Extra in builder
        self.abslabel.config(text=Chrom_txt +":"+str(self.abspos))
        self.entry.config(from_=min_seqlength, to=max_seqlength)

    # Extra in builder
    def complement_check(self,*args):
        abs_offset=RG_globals.bio_parameters["target_build_variant"]["abs_offset"]
        #exonpos_lookup=RG_globals.bio_parameters["target_build_variant"]["exonpos_lookup"]
        self.extension=int(self.entry2var.get())
        if self.entry.get() == self.blank:
            self.entry.delete(0,"end")
            self.entry.insert(0,1)
        self.transpos=int(self.entry.get())
        if self.transpos > max_seqlength:
            self.transpos = max_seqlength
        if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
            headclip=RG_globals.bio_parameters["target_build_variant"]["headclip"] # This should be == the REFSEQ.Headclip value
            #headclip=int(RG_globals.bio_parameters[RG_globals.target_locus]["Locus_range"].split(":")[0])-1 # Do not need the previous!
            self.localpos=self.transpos+self.extension+headclip
        else:
            if RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]:
                modpos=-self.extension # alternative to headclip, but does same job here
            else:
                modpos=self.extension
            self.localpos=int(RG_globals.exonpos_lookup[self.transpos])+int(modpos)
        self.abspos=abs(abs_offset+self.localpos)
    # Extra in builder
        
# END class source_sliders_builder

def src_sliders_build_clear_seqs():
    global refseqtxt,varseqtxt,varnameseqtxt,out_seq,varnametxt,hapnameseqtxt,hapnametxt
    refseqtxt.configure(state="normal")
    refseqtxt.delete('1.0', tk.END)
    refseqtxt.configure(state="disable")

    varseqtxt.configure(state="normal")
    varseqtxt.delete('1.0', tk.END)
    varseqtxt.configure(state="disable")
    out_seq=""

    hapnameseqtxt.configure(state="normal")
    hapnameseqtxt.delete('1.0', tk.END)
    hapnameseqtxt.configure(state="disable")
    hapnametxt="hap01"

    varnameseqtxt.configure(state="normal")
    varnameseqtxt.delete('1.0', tk.END)
    varnameseqtxt.configure(state="disable")
    varnametxt=""
    

def src_sliders_instantiate(window):
    global src_slider_widgets,src_sliders_vals,src_slider_widgets_active_count,max_src_slider_widgets
    src_slider_widgets=[]
    src_slider_widgets_active_count=0
    # Instantiating the label/slider display
    for r in range(len(src_sliders_vals)):
        # src_slider_widgets.append(source_sliders(window,r,from_=0, to=100,showvalue=0)) # Use showvalue to hide(0), or show(1) the slider values above the slider
        src_slider_widgets.append(source_sliders(window,r,from_=0, to=100)) # showvalue not supported in ttk, but defaults to not showing above slider
    src_slider_widgets_active_count=len(src_slider_widgets)
    max_src_slider_widgets=len(src_slider_widgets)
    # End Instantiating the label/slider display
    return

def src_sliders_instantiate_build(window):
    global src_slider_widgets_build,src_sliders_vals_build,max_seqlength,src_slider_widgets,max_src_slider_widgets
    global local_begin,local_end,abs_Begin,abs_End,trans_Begin,trans_End
    global refseqtxt,varseqtxt,out_seq,varnameseqtxt,hapnameseqtxt,hapnametxt
    
    seq_region=""
    seq_region_success=False
    out_seq=""
    varnametxt=""
    hapnametxt="hap01"
    abs_maptxt=""
    
    src_slider_widgets_build=[]
    # Instantiating the label/slider display
    for r in range(len(src_sliders_vals_build)):
        src_slider_widgets_build.append(source_sliders_builder(window,r, to=max_seqlength,showvalue=0))

    def get_ref_range():
        global local_begin, local_end
        '''
        if local_begin > local_end:
           local_begin,local_end=local_end,local_begin #swap
        '''
        #print("local_end %s; local_begin %s; local_end-local_begin %s"%(local_end,local_begin,local_end-local_begin))

        success=False

        RG_globals.bio_parameters["target_build_variant"]["local_begin"]=local_begin
        RG_globals.bio_parameters["target_build_variant"]["local_end"]=local_end

        RG_globals.bio_parameters["target_build_variant"]["abs_Begin"]["value"]=abs_Begin
        RG_globals.bio_parameters["target_build_variant"]["abs_End"]["value"]=abs_End

        RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["value"]=trans_Begin
        RG_globals.bio_parameters["target_build_variant"]["trans_End"]["value"]=trans_End
        
        RG_globals.bio_parameters["target_build_variant"]["trans_Begin_ext"]["value"]=trans_Begin_ext
        RG_globals.bio_parameters["target_build_variant"]["trans_End_ext"]["value"]=trans_End_ext

        subseq,run_success=RG_builder.get_ref_subseq3(local_begin,local_end)
        if run_success:
            RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]=str(subseq)
        else:
            RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]="FAIL"
        
        ref_range=RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]
        
        ref_string=RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["viewstring"]
        #print("ref_range %s"%ref_range)
        if ref_range !="FAIL":
            success= True
        return ref_range,ref_string,success
    
    def calc_seq():
        global varnametxt
        global trans_text
        nonlocal seq_region,seq_region_success
        get_slider_vals()
        #print("trans_Begin %s; local_begin %s| trans_End %s;local_end %s"%(trans_Begin,local_begin,trans_End,local_end))
        seq_region,ref_string,seq_region_success=get_ref_range()        

        if seq_region_success:
            clear_GUI_text2()
            #varnametxt="%s_%s_to_%s"%(RG_globals.target_locus,abs_Begin,abs_End)
            #varnametxt="%s_to_%s"%(abs_Begin,abs_End)
            varnametxt=RG_globals.bio_parameters["target_build_variant"]["var_name"]["value"]

            #write_GUI_text2("[trans_Begin %s;local_begin %s, trans_End %s; local_end %s Strln: %s Calclen %s\n"%(trans_Begin,local_begin,trans_End,local_end,len(seq_region),abs(local_begin-local_end)+1))
            write_GUI_text2("[%s_Begin %s;local_begin %s, %s_End %s; local_end %s Strln: %s Calclen %s\n"%(trans_text,trans_Begin,local_begin,trans_text,trans_End,local_end,len(seq_region),abs(local_begin-local_end)+1))

            write_GUI_text2("seq_region\n%s\n"%ref_string)
            
            refseqtxt.configure(state="normal")
            refseqtxt.delete('1.0', tk.END)
            refseqtxt.insert(tk.END,ref_string)
            refseqtxt.configure(state="disable")

            varseqtxt.configure(state="normal")
            varseqtxt.delete('1.0', tk.END)
            varseqtxt.insert(tk.END,seq_region)

            hapnameseqtxt.configure(state="normal")
            hapnameseqtxt.delete('1.0', tk.END)
            hapnameseqtxt.insert(tk.END,hapnametxt)

            varnameseqtxt.configure(state="normal")
            varnameseqtxt.delete('1.0', tk.END)
            varnameseqtxt.insert(tk.END,varnametxt)
        else:
            src_sliders_build_clear_seqs()
            clear_GUI_text2()
            write_GUI_text2("Selected region has failed. Try again")        
        
    def save_seq():
        nonlocal seq_region,seq_region_success
        global local_begin,local_end
        global src_slider_widgets, src_slider_widgets_active_count,src_slider_widgets_old_active_count
        global src_sliders_vals,src_sliders_strvals
        global IUPAC_sub

        ref_seq=seq_region
        clear_GUI_text2()
        if ref_seq == "\n" or ref_seq == "" or not seq_region_success : # This should not happen
            write_GUI_text2("Empty reference sequence or error\n")
        elif max_src_slider_widgets == len(RG_globals.mutlabels):
            write_GUI_text2("Full! No more additional variants this time.\n")
        else:
            var_seq=varseqtxt.get('0.3',tk.END)
            var_seq = var_seq[:-1]
        
            if var_seq=="":
                var_seq="-"

            var_name=varnameseqtxt.get('0.3',tk.END)
            var_name = var_name[:-1]

            hap_name=hapnameseqtxt.get('0.3',tk.END)
            hap_name = hap_name[:-1]

            RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]=str(ref_seq)
            RG_globals.bio_parameters["target_build_variant"]["var_subseq"]["value"]=str(var_seq)
            RG_globals.bio_parameters["target_build_variant"]["var_name"]["value"]=var_name
            RG_globals.bio_parameters["target_build_variant"]["hap_name"]["value"]=hap_name
            
            ## Delete?: RG_globals.bio_parameters["target_build_variant"]["is_save_var"]=True
            ## Delete?:is_success=RG_gui_addon.call_buildstuff() # Success only if the ref_seq and var_seq differ
            is_success=save_add_var()
            if is_success:
                if abs_End > abs_Begin:
                    abs_maptxt="%s..%s"%(abs_Begin,abs_End)
                else:
                    abs_maptxt="%s..%s"%(abs_End,abs_Begin)
                spacer="                     /"
                ref_save,var_save,var_name=get_add_saves()
                if var_name=="Quack":
                    local_begin=1; local_end=1
                msgtxt= 'Haplotype Name: %s\n     variation       %s..%s\n%sreplace="%s/%s"\n%sabsolute_map="%s"\n%scomment="%s"\n'\
                 %(hap_name,local_begin,local_end,spacer,ref_save,var_save,spacer,abs_maptxt,spacer,var_name)
                write_GUI_text2(msgtxt)
                
                #RG_globals.AddVars.append(addvar)
                old_mutlabs_tot=len(RG_globals.mutlabels)
                RG_globals.mutlabels=get_mutlabs()
                new_mutlabs_tot=len(RG_globals.mutlabels)
                if  new_mutlabs_tot > old_mutlabs_tot:
                    #refresh_src_sliders() # Need better, because refresh_src_sliders() resets all defaults,only need to set a default for new mutlab
                    # This seems to work ...
                    setfreq=50
                    src_sliders_vals[new_mutlabs_tot-1]=tk.IntVar(value=setfreq)
                    src_sliders_strvals[new_mutlabs_tot-1]=tk.StringVar(value=str(setfreq))

                    RG_globals.mutfreqs[new_mutlabs_tot-1]=setfreq

                    src_count=len(RG_globals.mutlabels)
                    src_slider_widgets_old_active_count=src_slider_widgets_active_count
                    for i in range(src_count):
                        w_text="%s_%s"%(RG_globals.target_locus,RG_globals.mutlabels[i])
                        refresh_source_slider(i,w_text)
                    src_slider_widgets_active_count=src_count
            else:
                write_GUI_text("Either: Reference is unmodified; duplicate name; empty name(s); no Haplotype saved")
                               
            #print("RG_globals.AddVars: %s"%RG_globals.AddVars)
    # End Instantiating the label/slider display

    ############## Calc button  ##############

    calcseq=ttk.Button(window,text=pygui_button_labels[1], command=calc_seq)
    calcseq.grid(row=5,column=0,sticky=tk.W)

    ############## refseqtxt  ##############
    reflabel=ttk.Label(window,text=RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["label"]).grid(row=6,column=0) # Redundant?
    refseqtxt=tk.Text(window,height=4,width=50,background="grey")
    refseqtxt.bind("<1>", lambda event: refseqtxt.focus_set())
    refseqtxt.grid(row=8,column=0,sticky=tk.W)
    refseqtxt.insert(tk.END,seq_region)
    refseqtxt.configure(state="disable")

    ############## varseqtxt  ##############
    def handle_DNA_input(s):
        varseqtxt.configure(state = tk.NORMAL)
        
        if  s.keysym.lower() in { "delete"}:
            varseqtxt.delete(1.0, tk.END)    
        
        elif re.match(RG_globals.IUPAC_codes,s.char.upper()):
            var_seq=varseqtxt.get('0.3',tk.END)
            if var_seq[:-1]=="-":
                varseqtxt.delete(1.0, tk.END)
            varseqtxt.insert(
                varseqtxt.index(tk.INSERT),s.char.upper()
            )

        elif s.keysym.lower() in {"backspace"}:
            varseqtxt.delete(
                varseqtxt.index(tk.INSERT)
                + "-1c" * (s.keysym.lower() == "backspace")
            )
        elif  s.char.lower() in { "-"}:
            varseqtxt.delete(1.0, tk.END)
            varseqtxt.insert(tk.END,"-")

        elif s.state == 4 and s.keysym == 'c': # Detect CTRL-C for copy
            content = varseqtxt.selection_get()
            window.clipboard_clear()
            window.clipboard_append(content)
        
        elif s.state == 4 and s.keysym == 'v': # Detect CTRL-V for paste
            content=window.selection_get(selection='CLIPBOARD')
            varseqtxt.insert('end',re.sub(IUPAC_sub,'',content.upper())) # Filter for IUPAC-allowed chars
        
        else:
            pass
        
        var_seq=varseqtxt.get('0.3',tk.END)
        if var_seq[:-1]=="":
            varseqtxt.delete(1.0, tk.END)
            varseqtxt.insert(tk.END,"-")
    
        varseqtxt.configure(state = tk.DISABLED)

    varlabel=ttk.Label(window,text=RG_globals.bio_parameters["target_build_variant"]["var_subseq"]["label"]).grid(row=9,column=0)
    varseqtxt=tk.Text(window,height=4,width=50,background="white",foreground="black")
    varseqtxt.bind("<1>", lambda event: varseqtxt.focus_set())
    varseqtxt.bind("<Key>", lambda e: handle_DNA_input(e))
    varseqtxt.grid(row=10,column=0,sticky=tk.W)

    ##varseqscroll = tk.Scrollbar(window,bd=1,orient="vertical")
    #varseqscroll = ttk.Scrollbar(window,orient="vertical")
    #varseqscroll.grid(row=10,column=1,sticky=tk.E)

    #varseqtxt.config(yscrollcommand=varseqscroll.set)
    #varseqscroll.config(command=varseqtxt.yview)
    varseqtxt.insert(tk.END,seq_region)

    ############## varnameseqtxt  ##############
    def handle_varname_input(s):
        varnameseqtxt.configure(state = tk.NORMAL)
        if  s.keysym.lower() in { "delete"}:
            varnameseqtxt.delete(1.0, tk.END)
            
        elif re.search('[a-zA-Z0-9_.()>]',s.char):
            varnameseqtxt.insert(
                varnameseqtxt.index(tk.INSERT),s.char
            )
        elif s.keysym.lower() in {"backspace"}:
            varnameseqtxt.delete(
                varnameseqtxt.index(tk.INSERT)
                + "-1c" * (s.keysym.lower() == "backspace")
            )     
        else:
            pass
        varnameseqtxt.configure(state = tk.DISABLED)
    
    varnameseqlabel=ttk.Label(window,text=RG_globals.bio_parameters["target_build_variant"]["var_name"]["label"]).grid(row=11,column=0)
    varnameseqtxt=tk.Text(window,height=2,width=20,background="white",foreground="black")
    varnameseqtxt.bind("<1>", lambda event: varnameseqtxt.focus_set())
    varnameseqtxt.bind("<Key>", handle_varname_input)
    varnameseqtxt.grid(row=12,column=0)
    varnameseqtxt.insert(tk.END,varnametxt)

    ############## hapnameseqtxt  ##############

    def handle_hapname_input(s):
        hapnameseqtxt.configure(state = tk.NORMAL)
        if  s.keysym.lower() in { "delete"}: # s.char.lower() behaves the same as s.keysym.lower() here
            hapnameseqtxt.delete(1.0, tk.END)
            
        elif re.match('[a-zA-Z0-9_]',s.char):
            hapnameseqtxt.insert(
                hapnameseqtxt.index(tk.INSERT),s.char
            )
        elif s.keysym.lower() in {"backspace"}:
            hapnameseqtxt.delete(
                hapnameseqtxt.index(tk.INSERT)
                + "-1c" * (s.keysym.lower() == "backspace")
            )
        else:
            pass
        hapnameseqtxt.configure(state = tk.DISABLED)

    hapnameseqlabel=ttk.Label(window,text=RG_globals.bio_parameters["target_build_variant"]["hap_name"]["label"]).grid(row=13,column=0)
    hapnameseqtxt=tk.Text(window,height=2,width=20,background="white",foreground="black")
    #hapnameseqtxt.configure(bg="#1C6C0B", insertbackground='white')
    hapnameseqtxt.bind("<1>", lambda event: hapnameseqtxt.focus_set())
    hapnameseqtxt.bind("<Key>", handle_hapname_input)
    hapnameseqtxt.grid(row=14,column=0)
    hapnameseqtxt.insert(tk.END,hapnametxt)


    ############## Save button  ##############
    #save=tk.Button(buttons,text='Save',bg="blue", command=save_src_sliders_vals)
    #saveseq=tk.Button(window,text=pygui_button_labels[2],bg="blue", command=save_seq)
    saveseq=ttk.Button(window,text=pygui_button_labels[2], command=save_seq)
    saveseq.grid(row=15,column=0,sticky=tk.W)    

    # END Instantiating the save buttons to label/slider
    return

###### Start of: Could be in another module ####
# Save the additional variants
def is_make_addvar_complement():
    # Note that any transformation enacted in get_ref_subseq3 by calling is_make_ref_complement will have already been made
    # There's a different determination to decide whether it's flipped back again before saving to Addvar
    # If it's not transcript sequence, it's already on +, or flipped to +, so do not change
    # If the strand is + and the transcript is on the opposite strand, we previously flipped the default to -; so we flip the edited one back to +
    # If the strand is - and the transcript is on the same strand, we left it - ; so we flip the edited to +
    # We leave the others unchanged as they are already +
    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name: # Locus
        if RG_globals.bio_parameters["target_build_variant"]["ref_strand"]==1:  # Strand-check
            is_complement=False
        else:
            is_complement=None # Not really False: this is a special case. Changes will have incorrect specification!!
    elif RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]:  # Transcripts that are is_join_complement
        if RG_globals.bio_parameters["target_build_variant"]["ref_strand"]==1:  # Strand-check
            is_complement=True
        else:
            is_complement=False
    elif RG_globals.bio_parameters["target_build_variant"]["ref_strand"]==1:   # Transcripts Not is_join_complement
        is_complement=False
    else:
        is_complement=True
    return is_complement

def save_add_var():
    is_success=False
    if RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"] != RG_globals.bio_parameters["target_build_variant"]["var_subseq"]["value"]:
        is_success=True
    if is_success and (RG_globals.bio_parameters["target_build_variant"]["hap_name"]["value"] !="") and (RG_globals.bio_parameters["target_build_variant"]["var_name"]["value"] !=""):
        is_success=save_add_var2()
            
        # Not essential, but this may be needed for App.js message panel report
        '''
        spacer="                     /"
        if  int(addvar["abs_End"]) > int(addvar["abs_Begin"]):
            abs_maptxt="%s..%s"%(addvar["abs_Begin"],addvar["abs_End"])
        else:
            abs_maptxt="%s..%s"%(addvar["abs_End"],addvar["abs_Begin"])
        RG_globals.bio_parameters["target_build_variant"]["variant_view"]=\
                'Haplotype Name: %s\n     variation       %s..%s\n%sreplace="%s/%s"\n%sabsolute_map="%s"\n%scomment="%s"\n'\
                 %(addvar["hapname"],addvar["local_begin"],addvar["local_end"],spacer,addvar["ref_seq"],addvar["var_seq"],spacer,abs_maptxt,spacer,addvar["varname"])

        RG_globals.bio_parameters["target_build_variant"]["variant_view"]="Reference is unmodified; no Haplotype generated" # Set this when is_success=False

        '''
        #  End of: Not essential ...
    else:
        is_success=False        
    return is_success

def get_add_saves():
    is_do_complement=is_make_addvar_complement()
    varname_save=RG_globals.bio_parameters["target_build_variant"]["var_name"]["value"],
    if is_do_complement== None:
        ref_save="-"
        var_save="-"
        varname_save="Quack"
    elif is_do_complement:
        ref_save=RG_process.get_revcomp(RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"])
        var_save=RG_process.get_revcomp(RG_globals.bio_parameters["target_build_variant"]["var_subseq"]["value"])
    else:
        ref_save=RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]
        var_save=RG_globals.bio_parameters["target_build_variant"]["var_subseq"]["value"]
    return (ref_save,var_save,varname_save)

def save_add_var2():
    match = False
    for label in RG_globals.mutlabels:
        if label == RG_globals.bio_parameters["target_build_variant"]["hap_name"]["value"]:
            match = True
            break
    if not match:

        ref_save,var_save,varname_save=get_add_saves()
        addvar={"locus":RG_globals.target_locus,
                "hapname":RG_globals.bio_parameters["target_build_variant"]["hap_name"]["value"],
                "varname":varname_save,
                "local_begin":RG_globals.bio_parameters["target_build_variant"]["local_begin"],
                "local_end":RG_globals.bio_parameters["target_build_variant"]["local_end"],
                "ref_seq":ref_save,
                "var_seq":var_save,
                "abs_Begin":RG_globals.bio_parameters["target_build_variant"]["abs_Begin"]["value"],
                "abs_End":RG_globals.bio_parameters["target_build_variant"]["abs_End"]["value"]
                }
        RG_globals.bio_parameters["target_build_variant"]["AddVars"].append(addvar)
    return (not match)

###### End of: Could be in another module ####

# Defining the output flags toggle options
class flags_toggler:
    def __init__(self, master2,v_index,**kwargs):
        global tog_labels,tog_booleans

        #self.panel2 = tk.Frame(master2)
        self.panel2 = ttk.Frame(master2)
        self.panel2.grid()
        self.vindex=v_index

        if tog_booleans[self.vindex] != None:
            # True or False values - display label and show default box ticked or not
            self.label=ttk.Label(self.panel2,width=21,anchor="w",text=tog_labels[v_index])
            self.label.grid(row=0,column=0,sticky=tk.W)
            self.bool = tk.BooleanVar()
            self.bool.set(tog_booleans[self.vindex])

            #self.btn=tk.Checkbutton(self.panel2, width=5,selectcolor="blue", highlightcolor="blue", var=self.bool,command=self.set_tog)
            self.btn=ttk.Checkbutton(self.panel2, width=5, var=self.bool,command=self.set_tog)

            self.btn.grid(row=0,column=1,sticky=tk.W)
        else:
            # Null settings hide the label and box but preserve layout
            
            #self.label=tk.Label(self.panel2)# retain this or refresh_genelabels kills
            self.label=ttk.Label(self.panel2)# retain this or refresh_genelabels kills
            #self.label=tk.Label(self.panel2, width=25,text="", bg="light grey")
            #self.label.grid(row=0,column=0,sticky=tk.W)
            #self.label.grid(row=0,column=0)

    def set_tog(self):
        tog_booleans[self.vindex]=format(self.bool.get())
        #print("clicked %s %s"%(tog_labels[self.vindex],format(self.bool.get())))

# END Defining the output flags toggle options


class flags_toggler2_builder:
    def __init__(self, master2,v_index,**kwargs):
        global tog_labels2,tog_booleans2

        #self.panel2 = tk.Frame(master2)
        self.panel2 = ttk.Frame(master2)
        self.panel2.grid()
        self.vindex=v_index

        if tog_booleans2[self.vindex] != None:
            # True or False values - display label and show default box ticked or not
            self.label=ttk.Label(self.panel2,width=22,anchor="w",text=tog_labels2[v_index])
            self.label.grid(row=0,column=0,sticky=tk.W)
            self.bool = tk.BooleanVar()
            self.bool.set(tog_booleans2[self.vindex])

            #self.btn=tk.Checkbutton(self.panel2, width=5,selectcolor="blue", highlightcolor="blue", var=self.bool,command=self.set_tog)
            self.btn=ttk.Checkbutton(self.panel2, width=5,var=self.bool,command=self.set_tog)
            self.btn.grid(row=0,column=1,sticky=tk.W)
        else:
            # Null settings hide the label and box but preserve layout
            
            #self.label=tk.Label(self.panel2)# retain this or refresh_genelabels kills
            self.label=ttk.Label(self.panel2)# retain this or refresh_genelabels kills
            #self.label=tk.Label(self.panel2, width=25,text="", bg="light grey")
            #self.label.grid(row=0,column=0,sticky=tk.W)
            #self.label.grid(row=0,column=0)

    def set_tog(self):
        tog_booleans2[self.vindex]=format(self.bool.get())
        # This sets the global is_CDS boolean immediately, rather than wait until "GO"
        tobools=flags_to_bools(tog_booleans2)
        RG_globals.is_CDS=tobools[flag_widget_order2["is_CDS"]]
        # As this is a toggle, call the update to resplice
        if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
            #Only toggled CDS, but still need to clear
            src_sliders_build_clear_seqs()
        else: 
            get_muttranscripts(False)       
        #print("clicked %s %s"%(tog_labels[self.vindex],format(self.bool.get())))
# END Defining the output flags toggler2 builder options


# Defining the initialising values spinbox options
class spinlimits:
    def __init__(self, master2,v_index,**kwargs):
        global spinvals,spin_labels,spin_mins,spin_maxs
        #self.panel2 = tk.Frame(master2)
        self.panel2 = ttk.Frame(master2)
        self.panel2.grid()

        self.vindex=v_index
        self.inMin=spin_mins[self.vindex]
        self.inMax=spin_maxs[self.vindex]
        self.labelplus=" (%s - %s)"%(self.inMin,self.inMax)
        self.label=ttk.Label(self.panel2, width=20,anchor="w",text=spin_labels[self.vindex]+self.labelplus)
        self.label.grid(row=0,column=0,sticky=tk.W)
        self.putMin=self.inMin
        self.putMax=self.inMax

        self.spinvar=spinvals[self.vindex]
        #self.spin=tk.Spinbox(self.panel2, width=4, textvariable=self.spinvar,from_=self.inMin,to=self.inMax)
        self.spin=ttk.Spinbox(self.panel2, width=4, textvariable=self.spinvar,from_=self.inMin,to=self.inMax)
        #self.spin=tk.Spinbox(self.panel2, width=4, textvariable=self.spinvar,from_=1,to=999)
        self.spin.grid(row=0,column=1,sticky=tk.W)
        self.blank=''
        # self.old_value=spinvals[self.vindex]
        self.old_value=self.spinvar.get()
        self.last_safe_value=self.spinvar.get()
        self.spinvar.trace('w', self.checkspin2a)
        self.get, self.set = self.spinvar.get, self.spinvar.set

    def checkspin2a(self, *args):
       # A range-check done before processing by get_config_limits()
       if self.get().isdigit():
           # the current value is only digits; allow this
           if int(self.putMin) <= int(self.get()) <= int(self.putMax):
               # only allow anything within the wider constraint
               self.old_value = self.get()
               self.spinvar.set(int(self.get()))
               spinvals[self.vindex]=self.spinvar

           # Now reject as it's outside the wider constraint

           elif int(self.get()) < int(self.putMin): # If it's less than wider constraint, set it to that
                                                    # this is permitted to allow building up string eg: from 1 to 199
                   self.old_value = self.putMin
                   self.spinvar.set(int(self.putMin))
           else:
               self.set(self.old_value)

       elif self.get() == self.blank:
         # the current value is blank - set to min
          self.old_value = self.inMin
          self.spinvar.set(int(self.inMin))
          spinvals[self.vindex]=self.spinvar
       else:
          #there's non-digit characters in the input; reject this
          self.set(self.old_value)

       # Would typically do a range-check here, but constraining between 8-200, it's impossible to
       # enter any 2 digit numbers between 21 and 79, and awkward for many others
       # Here's one way to save the last in-range value, and use these, or just blast through as in get_config_limits()
       #global spin_safe_vals
       #if int(self.inMin) <= int(self.old_value) <= int(self.inMax):
       #    self.last_safe_value=self.old_value
       #    spin_safe_vals[self.vindex].set(int(self.last_safe_value))
       #print("Last Safe value %s self.old_value %s "%(self.last_safe_value,self.old_value))

# END Defining the initialising values spinbox options

def tog_instantiate(window):
    global tog_labels,toggler_widgets
    # Instatiating the output flags toggle options
    r=0
    toggler_widgets=[]
    for v in tog_labels:
        toggler_widgets.append(flags_toggler(window,r))
        r+=1

def tog_instantiate2_builder(window):
    global tog_labels2,toggler_widgets2
    # Instatiating the output flags toggle options
    r=0
    toggler_widgets2=[]
    for v in tog_labels2:
        toggler_widgets2.append(flags_toggler2_builder(window,r))
        r+=1

def spin_instantiate(window):
    global spin_widgets,spinvals
    # Instatiating the initialising values spinbox options
    r=0
    spin_widgets=[]
    for v in spinvals:
        spin_widgets.append(spinlimits(window,r))
        r+=1

def refresh_genelabels_builder():
    global toggler_widgets,flag_widget_order,previous_target_locus,previous_target_transcript_name
    set_flags_output_labels()
    redo=False
    
    for i in range(3):
        toggler_widgets[i].label['text']=flags_output_labels[i]

    #print("previous_target_locus %s, RG_globals.target_locus %s, previous_target_transcript %s, RG_globals.target_transcript_name %s"
    #      %(previous_target_locus,RG_globals.target_locus,previous_target_transcript,RG_globals.target_transcript_name))

    if previous_target_transcript_name != RG_globals.target_transcript_name:
        previous_target_transcript_name =RG_globals.target_transcript_name
        redo=True
        redo_locus=False

    if previous_target_locus !=RG_globals.target_locus:              
        previous_target_locus = RG_globals.target_locus
        RG_globals.target_transcript_name = RG_globals.empty_transcript_name
        redo=True
        redo_locus=True
 
    if redo:
        get_muttranscripts(redo_locus)
        set_refseq_target_url()
        # Restore the next two if reinstating the links in function initialise_tk()
        #set_topbar_ENS_label_link()
        #set_topbar_LRG_label_link()
        set_slider_build_label()
        set_listbox_locus_label()
        set_listbox_transcript_label()
        src_sliders_build_clear_seqs()

def add_mutlabs(mutlabs):
    if len(RG_globals.bio_parameters["target_build_variant"]["AddVars"])>0:
        for AddVar in RG_globals.bio_parameters["target_build_variant"]["AddVars"]:
            #print("RG_globals.AddVars item: %s"%item["locus"])
            if AddVar["locus"]==RG_globals.target_locus:
                match=False
                for hapname in mutlabs:
                    if AddVar["hapname"] == hapname:
                        match=True
                if not match:       
                    mutlabs.append(AddVar["hapname"])
    return mutlabs

def initialise_mutfiles():
    global Infilepath,gbfiles
    #RG_main.general_initiate()
    Infilepath=RG_globals.infilepathroot+RG_globals.target_locus+"/"
    gbfiles=RG_io.list_files(Infilepath,'*%s'%RG_globals.Seq_IO_file_ext)
    RG_globals.mutlabels=get_mutlabs()    

def get_mutlabs(): # Get list of non-ref files in target-input folder
    global Infilepath,gbfiles
    mutlabs=[]
    for filename in gbfiles:
        #print("filename: %s"%filename)
        mutlabel=re.sub('%s|%s|%s_'%(Infilepath,RG_globals.Seq_IO_file_ext,RG_globals.target_locus),'',str(filename))
        if mutlabel != RG_globals.Ref_file_name:
            #print("mutlabel: %s"%mutlabel)
            mutlabs.append(mutlabel)
        mutlabs.sort()
    mutlabs=add_mutlabs(mutlabs)
    RG_process.mutfreqs_extend(mutlabs) # Extend mutfreqs if necessary
    #print("mutlabs %s"%mutlabs)
    return mutlabs


def get_muttranscripts(redo_locus): # Derive the transcript tables
    # NB: redo_locus isn't used!
    # Very importantly the globals set here, but calculated in main (simply to allow use of Vue.js GUI version) :
    # exonpos_lookup,abs_offset,ref_strand,max_seqlength,min_seqlength
    global max_seqlength,min_seqlength
    # Do the call to return the above global values, plus a journaling text 'transcript_view'.
    success,RG_globals.exonpos_lookup=RG_builder.get_muttranscripts2(False) # False is *not* to extend exome for exonpos_lookup
    # Setting those global values in this module from RG_globals calculated in get_muttranscripts
    max_seqlength=RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["max"]
    min_seqlength=RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["min"]
    transcript_view=RG_globals.bio_parameters["target_build_variant"]["transcript_view"] # Journaling/ annotation only
    clear_GUI_text2()
    write_GUI_text2("%s"%transcript_view)
    if max_seqlength < 1:
        write_GUI_text2("Error in reading sequence - max_seqlength is %s, setting to 1\n"%max_seqlength)
        max_seqlength=1
    else:
        refresh_src_sliders_builder()
    return

def get_slider_vals():
    global local_begin,local_end,src_slider_widgets,abs_Begin,abs_End
    global trans_Begin,trans_End,trans_Begin_ext,trans_End_ext

    trans_Begin=src_slider_widgets_build[0].transpos
    trans_End=src_slider_widgets_build[1].transpos

    trans_Begin_ext=src_slider_widgets_build[0].extension
    trans_End_ext=src_slider_widgets_build[1].extension

    local_begin=src_slider_widgets_build[0].localpos
    local_end=src_slider_widgets_build[1].localpos
    
    if local_begin > local_end:
        local_begin,local_end = local_end,local_begin# Swap
        #NB: This isn't enough- need to reset all these values, probably in refresh_source_slider_builder

    if trans_Begin > trans_End:
        trans_Begin,trans_End = trans_End,trans_Begin# Swap

    abs_Begin=src_slider_widgets_build[0].abspos
    abs_End=src_slider_widgets_build[1].abspos

def refresh_gui():
    refresh_genelabels_builder()
    initialise_mutfiles()
    refresh_src_sliders()
    refresh_src_sliders_builder()
    
def refresh_src_sliders():
    # This is the Haplotype variants list sliders
    #print("refresh_src_sliders")
    global src_slider_widgets, src_slider_widgets_active_count,src_slider_widgets_old_active_count
    global src_sliders_vals,src_sliders_strvals
    #print("RG_globals.mutlabels %s"%RG_globals.mutlabels)
    src_count=len(RG_globals.mutlabels)
    src_slider_widgets_old_active_count=src_slider_widgets_active_count
    def_freqlow=30
    def_freqstd=50
    setfreq=def_freqstd
    #print("src_count:%s, len(RG_globals.mutfreqs %s)"%(src_count,len(RG_globals.mutfreqs)))
    for i in range(src_count):
        endchars=RG_globals.mutlabels[i][-2:]
        #print("endchars %s"%endchars)
        if endchars =="00" :
            setfreq=0
        elif "som" in RG_globals.mutlabels[i]:
            setfreq=def_freqlow
        else:
            setfreq=def_freqstd
            #print("setfreq %s"%setfreq)
        src_sliders_vals[i]=tk.IntVar(value=setfreq)
        src_sliders_strvals[i]=tk.StringVar(value=str(setfreq))
        #print("i in range(src_count):%s"%i)
        if i <  len(RG_globals.mutfreqs):
            RG_globals.mutfreqs[i]=setfreq
        else:
            RG_globals.mutfreqs.append(setfreq)
        #print("label %s; freq %s"%(RG_globals.mutlabels[i],RG_globals.mutfreqs[i]))

    for i in range(src_count):
        w_text="%s_%s"%(RG_globals.target_locus,RG_globals.mutlabels[i])
        refresh_source_slider(i,w_text)
    
    for i in range(src_count,src_slider_widgets_old_active_count):
        forget_source_slider(i)
    src_slider_widgets_active_count=src_count
    #print("mutlabels: %s , mutfreqs: %s"%(RG_globals.mutlabels,RG_globals.mutfreqs))

def refresh_source_slider(v_index,w_text):
    # This is a single Haplotype variant slider - called for each one by refresh_src_sliders()
    #print("at refresh_src_slider - singular")
    self=src_slider_widgets[v_index]
    self.panel2.grid()
    self.label['text']=w_text
    self.label.grid(row=0,column=0,sticky=tk.W)
    #self.hypertextlink="%s%s/%s%s"%(RG_globals.infilepathroot,RG_globals.target_locus,w_text,RG_globals.Seq_IO_file_ext)
    #print("self.hypertextlink %s"%self.hypertextlink)
    self.entryvar=src_sliders_strvals[v_index]
    self.entry['textvariable']=self.entryvar
    self.entry.grid(row=0,column=1,sticky=tk.W)
    self.get, self.set = self.entryvar.get, self.entryvar.set
    self.entryvar.trace('w', self.checkentry)
    self.scalevar=src_sliders_vals[v_index]
    self.scale2['variable']= self.scalevar
    self.scale2.grid(row=0,column=2,sticky=tk.W)
    self.scalevar.trace('w', self.setentryfromscale)
          
def refresh_src_sliders_builder():
    # The Local_pos and Extension section
    # Attempts to restrict ranges by current values of other range ... too hard in current state of code 27-Feb-2024
    #print("at refresh_src_sliders_builder")
    #print("ref_strand %s"%RG_globals.bio_parameters["target_build_variant"]["ref_strand"])
    global src_slider_widgets_build
    global src_sliders_vals_build,src_sliders_strvals_build,src_sliders_labels_build,max_seqlength
    global trans_text

    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
        #w_text=RG_globals.target_locus
        trans_text="Locus"
        if RG_globals.bio_parameters["target_build_variant"]["ref_strand"]==-1:  # What was this here for?
            trans_text="Quack"
        
    elif RG_globals.is_CDS:
        trans_text="CDS"
    else:
        trans_text="mRNA"
    src_sliders_build_clear_seqs()
    
    src_sliders_vals_build[0]=tk.IntVar(value=1)
    src_sliders_strvals_build[0]=tk.StringVar(value=str(1))
    src_sliders_straddvals_build[0]=tk.StringVar(value=str(0))
    #Begin
    refresh_source_slider_builder(0,trans_text) #CDA/mRNA etc
    src_sliders_vals_build[1]=tk.IntVar(value=max_seqlength) # vals
    src_sliders_strvals_build[1]=tk.StringVar(value=str(max_seqlength)) # strvals
    
    src_sliders_straddvals_build[1]=tk.StringVar(value=str(0))
    #End
    refresh_source_slider_builder(1,trans_text)

def refresh_source_slider_builder(v_index,w_text):
    global max_seqlength,min_seqlength,Chrom_txt
    #print("At refresh_source_slider_builder")
    self=src_slider_widgets_build[v_index]
    self.panel2.grid()
    self.label['text']="%s_%s"%(w_text,src_sliders_labels_build[v_index])# This is the Source  x_Begin / x_End label 
    self.label.grid(row=0,column=0,sticky=tk.W)

    # Localpos
    self.entryvar=src_sliders_strvals_build[v_index] # strvals
    self.entry['textvariable']=self.entryvar
    self.entry.config(from_=min_seqlength, to=max_seqlength)# Reset new limits on spinbox
    self.entry.grid(row=0,column=1,sticky=tk.W)
    self.get, self.set = self.entryvar.get, self.entryvar.set
    self.entryvar.trace('w', self.checkentry)
    self.scalevar=src_sliders_vals_build[v_index] # vals
    self.scalevar.trace('w', self.checkentry2)

    # Extension
    self.entry2var=src_sliders_straddvals_build[v_index]
    self.entry2['textvariable']=self.entry2var
    self.entry2.grid(row=0,column=2,sticky=tk.W)
    self.entry2var.trace('w', self.checkentry2)# Adds any change here to the absolute 
    self.complement_check()
        
    self.abslabel.config(text=Chrom_txt +":"+str(self.abspos))
    #self.entry.config(from_=min_seqlength, to=max_seqlength)
    return

def forget_source_slider(v_index):
    self=src_slider_widgets[v_index]
    self.label.grid_forget()
    self.entry.grid_forget()
    self.scale2.grid_forget()

def refresh_main_options_list():
    # Deprecated because it messes up the build layout. Tried fixing different ways ... no luck
    global main, main_options_list
    global main_build
    main_options_list.destroy()
    main_build.destroy()
    initialise_main_options_list(main)
    main_options_list.pack(side="right", fill="both", expand=True)
    tog_instantiate(main_options_list)
    spin_instantiate(main_options_list)

    initialise_main_build(main)
    main_build.pack(side="left", fill="both", expand=True)
    src_sliders_instantiate(main_vars_list)
    tog_instantiate2_builder(main_left)
    src_sliders_instantiate_build(main_build)
    
    #refresh_gui()
    
def instantiate_widgets():
    global is_guitext_on,is_guitext2_on
    is_guitext_on=False
    is_guitext2_on=False
    src_sliders_instantiate(main_vars_list)
    src_sliders_instantiate_build(main_build)
    setup_reads_listbox(main_options_list)
    tog_instantiate(main_options_list) 
    spin_instantiate(main_options_list)
    fill_transcript_list()
    setup_listbox(main_left)
    tog_instantiate2_builder(main_left)
    #setup_reads_listbox(main_left)
    setup_gobutton(main_left)
    setup_GUI_text(bottom_left)
    setup_GUI_text2(bottom_right)
    is_guitext_on=True
    is_guitext2_on=True

def listbox_select(evt):
    global listbox,transcript_list
    selected=str(listbox.get())
    count=0
    locus="empty"
    #for lrg_locus in RG_globals.LRG_GeneList:
    for lrg_locus in RG_globals.GeneList:
        if lrg_locus == selected:
            locus=RG_globals.GeneList[count]
            break
        count+=1
    RG_globals.target_locus=locus
    
    refresh_gui()
    
    fill_transcript_list()
    ts_listbox.config(value=transcript_list)
    ts_listbox.current(0)

def ts_listbox_select(evt):
    selected=str(ts_listbox.get())
    RG_globals.target_transcript_name=selected
    if RG_globals.target_transcript_name != RG_globals.empty_transcript_name:
        RG_globals.target_transcript_id=RG_globals.Reference_sequences[RG_globals.target_locus]["mRNA"][RG_globals.target_transcript_name]
    else:
        RG_globals.target_transcript_id=RG_globals.empty_transcript_id
    #print("RG_globals.target_transcript_name %s, RG_globals.target_transcript_id %s, RG_globals.is_make_exome %s"%(RG_globals.target_transcript_name,RG_globals.target_transcript_id,RG_globals.is_make_exome))
    refresh_genelabels_builder()

def setup_listbox(inframe):
    global listbox,ts_listbox
    #frame=ttk.Frame(inframe,relief=tk.SUNKEN)
    frame=ttk.Frame(inframe)
    frame.grid(row=0,column=0,sticky='nsew')
    '''
    frame.columnconfigure(0,weight=1)
    frame.rowconfigure(0,weight=1)
    '''
    inframe.genelistlabel = ttk.Label(frame,text=RG_globals.bio_parameters["target_locus"]["label"])
    inframe.genelistlabel.grid()
    listbox=ttk.Combobox(frame,width=22,height=10,value=RG_globals.GeneList)
    #listbox=ttk.Combobox(frame,width=22,height=1,value=RG_globals.LRG_GeneList)
    #listbox.set(RG_globals.target_locus)
    #listbox.set(RG_globals.LRG_GeneList[0])
    listbox.set(RG_globals.GeneList[0])
    listbox.bind('<<ComboboxSelected>>',listbox_select)
    listbox['state']='readonly'
    listbox.focus()
    listbox.grid()

    inframe.haptemplatelabel=ttk.Label(frame,text=RG_globals.reference_haplotype)
    inframe.haptemplatelabel.grid()
                                       
    
    inframe.tslistlabel = ttk.Label(frame,text=RG_globals.bio_parameters["target_transcript_name"]["label"])
    inframe.tslistlabel.grid()

    ts_listbox=ttk.Combobox(frame,width=22,height=10,value=transcript_list)
    ts_listbox.set(RG_globals.target_transcript_name)
    ts_listbox.bind('<<ComboboxSelected>>',ts_listbox_select)
    ts_listbox['state']='readonly'
    ts_listbox.focus()
    ts_listbox.grid()
    
def fill_transcript_list():
    global transcript_list
    transcript_list=[]
    transcript_list.append(RG_globals.empty_transcript_name)
    for transcript in RG_globals.Reference_sequences[RG_globals.target_locus]["mRNA"]:
        transcript_list.append(transcript)

def setup_reads_listbox(inframe):
    global reads_listbox
    #frame=ttk.Frame(inframe,relief=tk.SUNKEN)
    frame=ttk.Frame(inframe)
    frame.grid(row=0,column=0,sticky='nsew')
    '''
    frame.columnconfigure(0,weight=1)
    frame.rowconfigure(0,weight=1)
    '''
    #inframe.readslistlabel = ttk.Label(frame,text=RG_globals.reads_type_label)
    #inframe.readslistlabel.grid()
    reads_listbox=ttk.Combobox(frame,width=22,height=3,value=RG_globals.ReadsList)
    reads_listbox.set(RG_globals.ReadsList[0])
    reads_listbox.bind('<<ComboboxSelected>>',reads_listbox_select)
    reads_listbox['state']='readonly'
    reads_listbox.focus()
    reads_listbox.grid()

def reads_listbox_select(evt):
    global reads_listbox
    selected=str(reads_listbox.get())
    if selected ==RG_globals.bio_parameters["is_frg_paired_end"]["label"]:
        RG_globals.is_frg_paired_end=True
        RG_globals.is_duplex = False
    elif selected == RG_globals.bio_parameters["is_duplex"]["label"]:
        RG_globals.is_frg_paired_end= False
        RG_globals.is_duplex = True
        
    elif selected == RG_globals.bio_parameters["is_simplex"]["label"]:
        RG_globals.is_frg_paired_end= False
        RG_globals.is_duplex = False
        
def setup_gobutton(inframe):
    global go_img,gobutton,icondir,Style
    # image-objects have to have global scope or the tkinter objects fail
    go_img = ImageTk.PhotoImage(Image.open(icondir+"./go.png"))
    
    Style = ttk.Style()
    # style.theme_use options: ('clam', 'alt', 'default', 'classic')
    Style.theme_use('alt')
    Style.configure('TButton', width = 20, borderwidth=1, focusthickness=3, focuscolor='none')

    gobutton=ttk.Button(inframe,command=save_and_go,image=go_img)
    #gobutton=ttk.Button(inframe,text="GO",command=save_and_go)

    #ttk_background_config(gobutton)
    #gobutton.grid(row=30,column=1,sticky=tk.W)# When using go_img
    gobutton.grid(row=30,column=0,sticky=tk.E)# When using text
    #gobutton.grid()

    #quitbutton.grid(row=30,column=1,sticky=tk.W)
    #quitbutton.grid()

def setup_GUI_text2(inframe):
    global guitextbox2

    #clearguitxt2=tk.Button(inframe,text="Clear Structure", bg='blue',command=clear_GUI_text2)
    clearguitxt2=ttk.Button(inframe,text="Clear Structure",command=clear_GUI_text2)
    clearguitxt2.grid()

    #guitextframe2=tk.Frame(inframe,bd=1,relief=tk.SUNKEN)
    guitextframe2=ttk.Frame(inframe)
    guitextbox2=tk.Text(guitextframe2)
    guitextbox2.bind("<1>", lambda event: guitextbox2.focus_set())
    guitextbox2.configure(state="disable",background="grey")
    guitextbox2.grid(row=0,column=0,sticky=tk.W)
    #guiscrollbar2=tk.Scrollbar(guitextframe2,bd=1,orient="vertical")
    guiscrollbar2=ttk.Scrollbar(guitextframe2,orient="vertical")
    guiscrollbar2.grid(row=0,column=1,sticky=tk.E)
    guitextbox2.config(yscrollcommand=guiscrollbar2.set)
    guiscrollbar2.config(command=guitextbox2.yview)
    guitextframe2.grid()
    guitextbox2.update()

def setup_GUI_text(inframe):
    global guitextbox
    #clearguitxt=tk.Button(bottom_left,text='Clear messages', bg='blue',command=clear_GUI_text)
    #clearguitxt=tk.Button(inframe,text=pygui_button_labels[0], bg='blue',command=clear_GUI_text)
    clearguitxt=ttk.Button(inframe,text=pygui_button_labels[0],command=clear_GUI_text)
    clearguitxt.grid()

    #guitextframe=tk.Frame(inframe,bd=1,relief=tk.SUNKEN)
    guitextframe=ttk.Frame(inframe)
    guitextbox=tk.Text(guitextframe)
    guitextbox.bind("<1>", lambda event: guitextbox.focus_set())
    guitextbox.configure(state="disable",background="grey")
    #guitextbox.bind("<Key>", lambda e: "break") # Disables text box entry and copy
    guitextbox.grid(row=0,column=0,sticky=tk.W)
    #guiscrollbar=tk.Scrollbar(guitextframe,bd=1,orient="vertical")
    guiscrollbar=ttk.Scrollbar(guitextframe,orient="vertical")
    guiscrollbar.grid(row=0,column=1,sticky=tk.E)
    guitextbox.config(yscrollcommand=guiscrollbar.set)
    guiscrollbar.config(command=guitextbox.yview)
    guitextbox.tag_config('link', foreground="blue",underline=1)
    guitextbox.tag_bind('link', '<Button-1>', results_hyperLink)
    guitextframe.grid()
    guitextbox.update()

def write_GUI_text(instring):
    global guitextbox
    # Suppress first call during instatiation - vestigial?
    global is_guitext_on
    if is_guitext_on:
        guitextbox.configure(state="normal")
        guitextbox.insert(tk.END, instring)
        guitextbox.update()
        guitextbox.configure(state="disable")
    else:
        is_guitext_on=True

def write_GUI_text2(instring):
    global guitextbox2
    # Suppress first call during instatiation - vestigial?
    global is_guitext2_on
    if is_guitext2_on:
        guitextbox2.configure(state="normal")
        guitextbox2.insert(tk.END, instring)
        guitextbox2.update()
        guitextbox2.configure(state="disable")
    else:
        is_guitext2_on=True

def write_GUI_results_link(instring):
    global guitextbox,run_count
    global is_guitext_on
    if is_guitext_on:
        guitextbox.configure(state="normal")
        #guitextbox['cursor']="hand2"
        guitextbox.insert(tk.END, instring,('link', str(run_count)))
        guitextbox.insert(tk.END, "\n")
        guitextbox.update()
        #guitextbox['cursor']="arrow"
        #guitextbox.update()
        guitextbox.configure(state="disable")

def clear_GUI_text():
    global guitextbox
    guitextbox.configure(state="normal")
    guitextbox.delete('1.0', tk.END)
    guitextbox.configure(state="disable")

def clear_GUI_text2():
    global guitextbox2
    guitextbox2.configure(state="normal")
    guitextbox2.delete('1.0', tk.END)
    guitextbox2.configure(state="disable")

#initialise_stuff()
initialise_tk()
initialise_stuff()

#print("GUI 1: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)
instantiate_widgets()
#print("GUI 3: RG_globals.is_mut_out %s"%RG_globals.is_mut_out)
refresh_gui() 


# Keeping the widgets live
exploder_win.mainloop()
