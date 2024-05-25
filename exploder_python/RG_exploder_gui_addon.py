#!/usr/local/bin/python3
Progver="RG_exploder_gui_addon.py"
ProgverDate="25-May-2024"
'''
© author: Cary O'Donnell for Replicon Genetics 2018, 2019, 2020, 2021, 2022, 2023, 2024
Calls from RG_exploder_gui.py that were formally housed by RG_exploder_main.py

© author: Cary O'Donnell for Replicon Genetics 2018, 2019, 2020, 2021, 2022, 2023, 2024
'''
# =================================================
# python_imports
# =================================================
'''
Python imports. These define global constants, functions, data structures.
As such they must be defined first, being impractical to be called by subroutine and declaring all the globals
'''

import copy
#end of python_imports
'''
BioPython imports
'''
from Bio.SeqRecord import SeqRecord # Used in writing protein sequence

# =================================================
# End of python_imports
# =================================================

'''
RG module imports
'''
import RG_exploder_globals as RG_globals  # Configurable constants
import RG_exploder_main as RG_main
import RG_exploder_process as RG_process # Seqrecord Processing

global is_joindata_in_json
is_joindata_in_json= None

# ===============================================================
# End of initialising routines for global variables and constants
# ===============================================================

############################################################################
# Used only in gotcha fix in get_ref_subseq3 and get_transcript_data_process
############################################################################
    
def flip_is_frg_paired_end():
    # A hack for the 'builder' functions: get_ref_subseq3 & get_transcript_data_process
    # Calling splice_refseq when is_frg_paired end==True leaves the sequence unspliced, apart from end-trims
    # This is NOT what we want for the 'builder' functions, and this is the fix to repurpose splice_refseq
    global pre_paired_end
    pre_paired_end=copy.copy(RG_globals.is_frg_paired_end)
    RG_globals.is_frg_paired_end=False
    
def reset_is_frg_paired_end():
    # Re-instates original value of is_frg_paired_end after use of flip_is_frg_paired_end()
    global pre_paired_end
    RG_globals.is_frg_paired_end=pre_paired_end

#########################################
def get_transcript_data_json():
    # The joinlist is present in json file, use it
    global is_joindata_in_json
    is_joindata_in_json=True
    splice_joinlist_txt=[]
    if  RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
        splice_joinlist_txt=RG_globals.Reference_sequences[RG_globals.target_locus]["Locus_range"]
    else:
        if RG_globals.is_CDS:
            splice_joinlist_txt=RG_globals.Reference_sequences[RG_globals.target_locus]["CDS_join"][RG_globals.target_transcript_name]
        else:
            splice_joinlist_txt=RG_globals.Reference_sequences[RG_globals.target_locus]["mRNA_join"][RG_globals.target_transcript_name]
    #print("data_json: splice_joinlist_txt: %s"%(splice_joinlist_txt))
    splice_joinlist_txt=splice_joinlist_txt.split(",")
    return True,splice_joinlist_txt

def get_transcript_data_process(redo_locus):
    global REF_record,is_joindata_in_json
    is_joindata_in_json=False
    # The joinlist is not present in json file, so create it from the feature table
    # Fronts splice_a_sequence2 
    # The only interest in calling splice_a_sequence2 is to get the splice_joinlist_txt
    def make_join(in_joinlist):
        out_joinlist=""
        for item in in_joinlist:
            out_joinlist+="%s,"%item

        x_join={
        RG_globals.target_transcript_name:
            out_joinlist[:-1]
        }
        return x_join

    flip_is_frg_paired_end()
    RG_main.general_initiate(False)
    exists=True
    if redo_locus:
        REF_record,exists=RG_main.read_refseqrecord(RG_globals.Seq_Format)
    if not exists:
        print("Failed to read Refseq at RG_main.get_transcript_data_process")
    else:
        
        is_spliced,splicetotal,spliceseq,extra_record_features,target_gene,message,transcript_id,last_first,last_second,this_ref_offset,this_ref_endclip,splice_joinlist_txt_CDS=RG_main.splice_a_sequence2(REF_record,RG_globals.target_locus,RG_globals.CDS_trigger)
        RG_globals.bio_parameters["target_build_variant"]["CDS_join"]=make_join(splice_joinlist_txt_CDS) # Superfluous?

        is_spliced,splicetotal,spliceseq,extra_record_features,target_gene,message,transcript_id,last_first,last_second,this_ref_offset,this_ref_endclip,splice_joinlist_txt_mRNA=RG_main.splice_a_sequence2(REF_record,RG_globals.target_locus,RG_globals.mRNA_trigger)
        RG_globals.bio_parameters["target_build_variant"]["mRNA_join"]=make_join(splice_joinlist_txt_mRNA) # Superfluous?

        if RG_globals.is_CDS:
            splice_joinlist_txt=splice_joinlist_txt_CDS
        else:
            splice_joinlist_txt=splice_joinlist_txt_mRNA
    # This reverse is a fix solely to reinstate consistency between reading from config.json and the section in splice_a_sequence2 that already flips the join for is_join_complement=True
    # The other way to do it is to change both config.json AND the mrnapos_lookup assignment
    if RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]:
        splice_joinlist_txt.reverse()

    reset_is_frg_paired_end()
    #print("data_process splice_joinlist_txt: %s"%(splice_joinlist_txt))
    return exists,splice_joinlist_txt

#########################################
# MAIN
# This section once only had call_exploder_main() for use by a GUI coded in Vue.js,
# which batch-called this object after allowing the setting of user-modifiable parameters and pressing "Run"
#
# The Vue.js GUI version is available at https://repliconevaluation.wordpress.com/replicon-genetics/ngs-emulator-sets
# 
# The first version of RG_exploder_gui.py did some pre-reading of the input data, but "Reference_sequences" 
# was coded into config.json to replace this because values could not be passed back to the pyodide-compiled version fronted by Vue.js
# 
# Some of this pre-reading has been re-introduced (for the mRNA and CDS join features), eliminating the need to include it within config.json
# The Python version tests ('try') whether the mRNA and CDS join features are in config.json. If not, it pre-reads.
# 'Try' cannot be done in the pyodide version
# 
# There may be some other features that would be improved in the same way and eliminate the need for more of the config.json file
#
#    *The 'build' function is the original developer's name for the GUI function allowing the addition of further haplotypes,
#     This was developed after the initial 'explode' function; the fragmentation of sequence into reads
#########################################

def get_muttranscripts(redo_locus):
    # Derive the transcript tables
    # It might be obvious to use Refseq.abs_start, Refseq.abs_end, Refseq.polarity
    # But at this point we have not necessarily read in the sequence file yet. Data are from config.json

    # NB: redo_locus isn't used!
    # Very importantly the values set here to pass back to the calling GUI are:
    # mrnapos_lookup,abs_offset,ref_strand,max_seqlength,GRChver_txt,headclip,tailclip
    # plus feature_titles (non-critical)
    
    is_join_complement=RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]
    Region=RG_globals.Reference_sequences[RG_globals.target_locus]["Region"]
    abs_start=int(Region.split(":")[2])  # REFSEQ.abs_start - available only when sequence file read. This comes from config.json
    abs_end=int(Region.split(":")[3])    # REFSEQ.abs_end - ditto
    ref_strand=int(Region.split(":")[4]) # REFSEQ.polarity - ditto
    maxreflen=abs(abs_end-abs_start)+1
    if ref_strand == 1:
        offset=abs_start-1
    else:
        offset=abs_end+1

    abs_offset=offset*ref_strand
        
    #print("abs_start: %s; abs_end: %s; abs_offset: %s"%(abs_start,abs_end,abs_offset))

    locus_begin,locus_end=RG_globals.Reference_sequences[RG_globals.target_locus]["Locus_range"].split(":")
    locus_begin=int(locus_begin)
    locus_end=int(locus_end)
    
    # Need to store an offset for case of generating variants from *clipped* genomic sequence, so the offset is the Headclip length
    headclip=locus_begin-1
    tailclip=maxreflen-locus_end
    
    RG_globals.bio_parameters["target_build_variant"]["headclip"]=headclip
    RG_globals.bio_parameters["target_build_variant"]["tailclip"]=tailclip

    exon_total=0; intron_total=0; feature_titles=["Pre-Locus"]; starts=[]; ends=[]
    abs_starts=[]; abs_ends=[]; mrnapos_lookup=[0]; exon_length=[0];max_seqlength=0;
    template_length=0;template_cumulative_length=["-"]

    run_success,splice_joinlist_txt=get_joinlist()

    # set Pre-Locus range
    if is_join_complement:
        #headclip,tailclip=tailclip,headclip
        starts.append(maxreflen)
        ends.append(locus_end+1)
        exon_length.append(tailclip)
    else:
        starts.append(1)
        ends.append(locus_begin-1)
        exon_length.append(headclip)

    # set Locus ranges
    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name: # Locus
        feature_titles.append(RG_globals.empty_transcript_name)
        if is_join_complement:
            starts.append(locus_end)
            ends.append(locus_begin)
        else:
            starts.append(locus_begin)
            ends.append(locus_end)
        locus_length=locus_end-locus_begin+1
        exon_length.append(locus_length)
        max_seqlength=locus_length
        template_length+=locus_length
        template_cumulative_length.append(template_length)
        
    # or set introns & exons ranges
    else:
        feature_titles.append("Upstream")
        #mrnapos_lookup.append(0)
        exon_length.append(0)
        template_cumulative_length.append("-")
        if is_join_complement:
            starts.append(locus_end)
        else:
            starts.append(locus_begin)
                                                                        
        if RG_globals.is_CDS:
            exon_text="CDS "
            intron_text="CDS_Intron "
        else:
            exon_text="Exon "
            intron_text="Intron "

        #print("splice_joinlist_txt %s"%splice_joinlist_txt)

        for count in range(len(splice_joinlist_txt)): # was for item in splice_joinlist_txt
            this_pair=splice_joinlist_txt[count]
            begin,end=this_pair.split(":")
            begin=int(begin)
            end=int(end)

            #print("begin %s \n end %s"%(begin,end)) # Validation check
            #Last intron / upstream
            if is_join_complement:
                if end == locus_end:
                    ends.append(end)
                else:
                    ends.append(end+1)
            elif begin == locus_begin:
                    ends.append(begin)
            else:
                    ends.append(begin-1)
            exon_len=abs(end-begin+1)
            max_seqlength+=exon_len
            
            #print("begin:%s; end: %s; read_strand: %s; exon_len: %s"%(begin,end,read_strand,exon_len))
        
            exon_total+=1
            feature_titles.append(exon_text+str(exon_total))
            #print("lookup_begin: %s,lookup_end+1 %s"%(lookup_begin,lookup_end+1))
            this_exon=[]
            for n in range(begin,end+1):
                this_exon.append(n)
            if is_join_complement:
                this_exon.reverse()
            #print("exon %s: %s"%(exon_total,this_exon))
            for item2 in this_exon:
                mrnapos_lookup.append(item2)
            
            #print("Locus: %s mrnapos_lookup %s"%(RG_globals.target_locus,mrnapos_lookup))
            #print(" locus %s ; transcript %s; begin %s ; end %s ; exon_len %s "%(RG_globals.target_locus,RG_globals.target_transcript_name,begin,end,exon_len))

            if is_join_complement:
                starts.append(end)
                ends.append(begin)
                starts.append(begin-1)
            else:
                starts.append(begin)
                ends.append(end)
                starts.append(end+1)
            exon_length.append(exon_len)
            template_length+=exon_len
            template_cumulative_length.append(template_length)
            #Next intron / downstream start
            intron_total+=1
            feature_titles.append(intron_text+str(intron_total)+"-"+str(intron_total+1))
            template_cumulative_length.append(template_length)
            
            #print("%s"%mrnapos_lookup[(len(mrnapos_lookup)-10):-1])
        ## End of: for item in splice_joinlist_txt
        ## End of: for count in range(len(splice_joinlist_txt))
        feature_titles[-1]="Downstream"
        intron_total-=1
        exon_length.append(exon_len)
        template_cumulative_length[-1]="-"

        if is_join_complement:
            ends.append(locus_begin)
        else:
            ends.append(locus_end)

    ## End of: if RG_globals.target_transcript_name == RG_globals.empty_transcript_name: # Locus
    # Finishing off    
    feature_titles.append("Post-Locus")
    if is_join_complement:
        starts.append(locus_begin-1)
        ends.append(1)
        exon_length.append(headclip)
    else:
        starts.append(locus_end+1)
        ends.append(maxreflen)
        exon_length.append(tailclip)
    template_cumulative_length.append("-")

    if starts[-2]==starts[-1]:# This to catch where there's no downstream because exon boundary coincides with locus boundary
        ends[-2]=starts[-2]
        feature_titles[-2]="No Downstream"

    if starts[1]==ends[1]:# This to catch where there's no upstream because exon boundary coincides with locus boundary
        feature_titles[1]="No Upstream"
         
    #print("feature_titles %s \n starts %s \n ends %s"%(feature_titles,starts,ends)) # Validation check

    # This is all that is essential to send back to GUI
    RG_globals.bio_parameters["target_build_variant"]["mrnapos_lookup"]=mrnapos_lookup
    RG_globals.bio_parameters["target_build_variant"]["abs_offset"]=abs_offset
    RG_globals.bio_parameters["target_build_variant"]["ref_strand"]=ref_strand
    RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["max"]=max_seqlength
    RG_globals.bio_parameters["target_build_variant"]["trans_End"]["max"]=max_seqlength

    # The following is only to build transcript_view for user-notification via RG_globals.bio_parameters["target_build_variant"]["transcript_view"]
    # It recreates the intron/exon table as seen in Ensembl transcript view. Ideally put in journal somehow
    for count in range(len(ends)):
        abs_starts.append(abs(starts[count]+abs_offset))
        abs_ends.append(abs(ends[count]+abs_offset))

    if RG_globals.is_CDS and (RG_globals.target_transcript_name != RG_globals.empty_transcript_name):
        add="_CDS"
    else:
        add=""

    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
        trans_text="DNA"
    elif RG_globals.is_CDS:
        trans_text="CDS"
    else:
        trans_text="mRNA"

    transcript_view="Locus: %s, %s: %s%s,\tLength: %s\n Source %s\n"%(RG_globals.target_locus,trans_text,RG_globals.target_transcript_name,add,max_seqlength,
                                                                             RG_globals.Reference_sequences[RG_globals.target_locus]["Region"])
    
    transcript_view+=" Label\t\t Local\t\t| Global\t\t| Length| Template (Cumulative)\n"

    for count in range(len(feature_titles)):
        loc_start=starts[count]
        loc_end=ends[count]
        glob_start=abs_starts[count]
        glob_end=abs_ends[count]
        this_len=abs(ends[count]-starts[count])+1
        if this_len==1: this_len = 0
        #print("%s; start %s; end %s; this_len %s; templ %s"%(feature_titles[count],starts[count],ends[count],this_len,template_cumulative_length[count]))
            
        transcript_view+=" %s:\t %s - %s\t| %s - %s\t| %s\t | %s\n"%(feature_titles[count],loc_start,loc_end,glob_start,glob_end,
                                                                     this_len,template_cumulative_length[count])
    transcript_view+="\n"

    RG_globals.bio_parameters["target_build_variant"]["transcript_view"]=transcript_view #Not essential - only used for GUI information
    # End of user-notification section. 
    return run_success


# Diversion from main 'exploder' function to extract a subsequence that supports user-defined variants (aka the 'builder' function)
def get_joinlist():
    global is_joindata_in_json
    success=False
    is_joindata_in_json=None
    try: # Are the joins defined in globals via the json file?
        barf=RG_globals.Reference_sequences[RG_globals.target_locus]["CDS_join"] # tests if it exists
        success,splice_joinlist_txt=get_transcript_data_json()
    except: # Failed, so read the Reference file # Works for python, but always define i json for pyodide because pyodide doesn't like a fail  
        success,splice_joinlist_txt=get_transcript_data_process(True)
    if success:
        #RG_globals.bio_parameters["target_build_variant"]["splice_joinlist_txt"]=splice_joinlist_txt
        if splice_joinlist_txt == "":
            success=False
    return success,splice_joinlist_txt

def is_make_ref_complement():
    # Was this, but obviously wrong because of AK2
    # is_complement=(RG_globals.target_transcript_name != RG_globals.empty_transcript_name) and RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]
    # Create the default variant 'reference' sequence;
    # If it's not transcript sequence, then determine the strand, and present as forward (so flip if -1)
    # If it is a transcript then flip if the transcript is on the opposite strand (is_join_complment)
    # We are sending back the instruction whether to flip (True) or not (False)
    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name: # Locus
        if RG_globals.bio_parameters["target_build_variant"]["ref_strand"]==1:  # Strand-check
            is_complement=False
        else:
            is_complement=None # Not really True: this is a special case where should leave it or flip polarity entirely. What to do?
            
    elif RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]:# Transcripts
        is_complement=True
    else:
        is_complement=False                                    
    return is_complement

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

def get_ref_subseq3():
    flip_is_frg_paired_end()
    RG_main.general_initiate(False)
    REF_record,exists=RG_main.read_refseqrecord(RG_globals.Seq_Format)
    REF_record=RG_main.splice_refseq(REF_record)
    reset_is_frg_paired_end()

    local_begin=RG_globals.bio_parameters["target_build_variant"]["local_begin"]
    local_end=RG_globals.bio_parameters["target_build_variant"]["local_end"]
    is_do_complement=is_make_ref_complement()

    reflen=local_end-local_begin+1
   
    if reflen > 100:
        first=REF_record.seq[local_begin-1:local_begin+2]
        last=REF_record.seq[local_end-3:local_end]
        if is_do_complement:
            first=RG_process.get_revcomp(first)
            last=RG_process.get_revcomp(last)
        viewstring="%s...%s"%(first,last)
        subref="NNN"
        success=True
    else:
        subref=REF_record.seq[local_begin-1:local_end]
        if is_do_complement:
            subref=RG_process.get_revcomp(subref)
        if len(subref)>0:
            success=True
        else:
            success=False
        viewstring="%s"%str(subref)

    if is_do_complement==None:
        subref="XXX"
        viewstring="Quack"
        
    RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]=str(subref)
    RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["viewstring"]="%s bases:%s"%(reflen,str(viewstring))
    
    # Building a default naming string for the variant

    trans_Begin=RG_globals.bio_parameters["target_build_variant"]["trans_Begin"]["value"]
    trans_End=RG_globals.bio_parameters["target_build_variant"]["trans_End"]["value"]

    trans_Begin_ext=RG_globals.bio_parameters["target_build_variant"]["trans_Begin_ext"]["value"]
    trans_End_ext=RG_globals.bio_parameters["target_build_variant"]["trans_End_ext"]["value"]
    if trans_Begin_ext==0:
        tbe=""
    else:            
        tbe=str(trans_Begin_ext)
        if trans_Begin_ext >0:
            tbe="+%s"%tbe
            
    if trans_End_ext==0:
        tee=""
    else:
        tee=str(trans_End_ext)
        if trans_End_ext >0:
            tee="+%s"%tee

    if trans_Begin==trans_End and (trans_Begin_ext==trans_End_ext):
        if trans_Begin_ext==0:
            varname="%s"%(trans_Begin)
        else:
            varname="%s%s"%(trans_Begin,tbe)
    else:
        varname="%s%s_%s%s"%(trans_Begin,tbe,trans_End,tee)
        
    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
        varfront="g"
    elif RG_globals.is_CDS:
        varfront="c"
    else:
        varfront="t"

    RG_globals.bio_parameters["target_build_variant"]["var_name"]["value"]="%s.%s"%(varfront,varname)
    return subref,success

# Save the additional variants
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
 
    
#########################################
# When called from the python GUI module

def showmutfreqs():
    # Debug on these values
    print("RG_globals.mutfreqs %s"%RG_globals.mutfreqs)
    print("RG_globals.mutlabels %s"%RG_globals.mutlabels)
    
def call_buildstuff():
    #print(" ***** New call_exploder_main %s *****"%RG_globals.getime())
    #showmutfreqs()
    if RG_globals.bio_parameters["target_build_variant"]["is_get_ref"]:
        subseq,run_success=get_ref_subseq3()
        if run_success:
            RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]=str(subseq)
        else:
            RG_globals.bio_parameters["target_build_variant"]["ref_subseq"]["value"]="FAIL"
        RG_globals.bio_parameters["target_build_variant"]["is_get_ref"]=False

    elif RG_globals.bio_parameters["target_build_variant"]["is_get_muttranscripts"]:
        run_success=get_muttranscripts(True)
        RG_globals.bio_parameters["target_build_variant"]["is_get_muttranscripts"]=False

    elif RG_globals.bio_parameters["target_build_variant"]["is_save_var"]:
        run_success=save_add_var()
        RG_globals.bio_parameters["target_build_variant"]["is_save_var"]=False
    else:
        print("Called incorrectly")
        
    return run_success

