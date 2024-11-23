#!/usr/local/bin/python3
Progver="RG_exploder_main_30_6.py"
ProgverDate="07-Nov-2024"
'''
© author: Cary O'Donnell for Replicon Genetics 2018, 2019, 2020, 2021, 2022, 2023, 2024
This module reads in Genbank format files and uses any variant feature definitions to create those variants from the reference sequence.
These variants are then split into mutiple shorter fragments.
Output can be made to SAM, fasta & fastq files
Created using python3 and BioPython
======================================

Changes that could be added are:
g) Add function to read vcf file instead of Genbank ... will probably need to write a genbank to VCF converter first, as they don't appear to exist
h) Work out how to detect GC-runs, AT-runs and the type of mutation to introduce & add into cigar string
j) mutrecs is the biggest memory-user. Is there another way without creating them all at the same time?

=====================================

Major modification Feb 2023 - attempting to keep one code base to maintain different Biopython versions.
a) Deprecation of 'generic_dna' in MutableSeq and SeqRecord
b) sequence.reverse_complement() change to sequence.reverse_complement(inplace=True)
c) Uses extra module import Biopython_fix to modulate

© author: Cary O'Donnell for Replicon Genetics 2018, 2019, 2020, 2021, 2022, 2023, 2024
'''
# =================================================
# python_imports
# =================================================
'''
Python imports. These define global constants, functions, data structures.
As such they must be defined first, being impractical to be called by subroutine and declaring all the globals
'''
from io import StringIO #python 3
#from StringIO import StringIO #python 2
import numpy as np

import time  # Used in a few functions
import sys # Used in program_exit for sys.exit
from random import randint # Used in write_samout; generate_multisource_random_frags
from random import random
from random import choice
import copy
#end of python_imports
'''
BioPython imports
'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord # Used in writing protein sequence
from Bio.SeqFeature import SeqFeature, FeatureLocation  # used in splice_a_sequence
import Biopython_fix  # From Feb 2023 to fix MutableSeq, Seq, Bio.Alphabet deprecated from Biopython 1.78 (September 2020)

# =================================================
# End of python_imports
# =================================================

'''
RG module imports
'''
import RG_exploder_io as RG_io
# File input/output needed for verifying file presence using pathlib.Path eg: read_seqrecord. A different RG_io is used for pyodide when building webapp
# Used in input:read_single_record_input,read_seqrecord
# Used in output: initialise_journal,initialise_journal_htm,initialise_readme,write_gb_features,write_fastaout,write_fastqout,writemut,writeprotmut,write_refseq,write_samheader

import RG_exploder_globals as RG_globals  # Configurable constants
import RG_exploder_process as RG_process # Seqrecord Processing
import RG_exploder_builder as RG_builder # For exonplus_lookup recalculation

global mut_types #  Use as a trigger to identify if initialise_global_text_constants() has been called
mut_types=""

global is_append_journalfile # This is a file-open flag for output journal file
is_append_journalfile=False

from RG_exploder_globals import is_htm_journal # Admin set flag for creation of html version of journal file
global is_htm_journal

global is_append_journalfile_htm # This is a file-open flag for output journal html file
is_append_journalfile_htm=False

global is_append_readmefile
is_append_readmefile=False ## This is a file-open flag for output readme file

# =================================================================
# Start of initialising routines for global variables and constants
# =================================================================
            
def general_initiate(is_journal):
    #  Call all the initialising subroutines that define globals in sequence
    initialise_start_time()
    if mut_types=="" : # Double use as a trigger to identify if initialise_global_text_constants() has already been called
        initialise_global_vars() # this needs running once only to set the cross-module (global) shared constants and variables
    initialise_module_outflags()
    set_IO_filenames()
    initialise_module_globals(is_journal)

def initialise_global_vars():
    # Assign global variables - maybe should be command parameters
    '''
    Labels and filename extensions for the mutated (non-reference) sequences generated from variant annotation
    file names are RG_globals.target_locus+RG_globals.mutlabels[n]+".gb"
    '''
    initialise_progver()
    initialise_global_text_constants()
    initialise_global_consts()
    initialise_global_configs()
    initialise_flags_output_configs()
    return
# End of initialise_global_vars()

def initialise_progver():
    global Progver,ProgverDate
    return

def initialise_global_text_constants():
    # Mostly vestigial intent to keep a list of user-unchangeable globals rather than using RG_globals.*
    # Only mut_types remains because it's used as a trigger in general_initiate() - retain mut_types as a global even though it looks unnecessary
    mut_types=RG_globals.mut_types # Don't remove. It's a trigger in general_initiate() 
    return
# end of def initialise_text_constants()

def initialise_module_globals(is_journal):
    # Globals limited to ***this*** module
    global REFSEQ_RECORD
    global MaxRefPos # Maxlength of reference sequence
    global MaxVarPos
    MaxRefPos=0      # read_refseqrecord sets this on sequence read, but initialise here
    MaxVarPos=RG_globals.MaxVarPos # Set to original before resetting in read_refseqrecord    
    initialise_module_counters()
    initialise_module_outfiles()
    if is_journal:
        # Only set up journaling files for exploder function, not others like builder
        initialise_journal()
        initialise_readme()
    return

def initialise_module_counters():
    # Globals limited to ***this*** module
    global var_count, sub_count,  insert_count, delete_count, complex_count #  counted in MutateVarSeq
    zero_varcounts()  #  Above set to zero by zero_varcounts()

    global skip_count
    skip_count=0      # Initialises count of skipped regions (introns) in a varseq by setting skip_count to zero'''

    global mutlistlabels # Variable: builds up from mutlabels in def read_mutrecord.
                         # mutlistlabels is a text string that accumulates as the mutated sequences are read-in
    mutlistlabels=""

    global outfilestring # Accumulates output file annotation for readme file on exit 
    outfilestring=""

    return
# End of initialise_module_counters()


def initialise_flags_output_configs():
    # sequence output configurable flags
    # Hiding global definitions, where they can be amended by GUI module as need to use RG_globals.* instead

    from RG_exploder_globals import bio_parameters # A sledgehammer
    global bio_parameters

    from RG_exploder_globals import is_fastacigar_out,is_onefrag_out,is_muts_only,is_frg_paired_end,is_flip_strand,is_duplex,is_frg_label,is_journal_subs,is_mutate_ref,is_fastq_random,is_make_exome,is_trim_to_gene,is_exome_paired_end,is_pair_monitor
    #global is_fastacigar_out,is_onefrag_out,is_muts_only,is_frg_label,is_journal_subs,is_mutate_ref,is_make_exome

    from RG_exploder_globals import is_show_infilepath  
    global is_show_infilepath # Not intended for GUI module modification

    from RG_exploder_globals import is_fasta_out,is_fastq_out,is_sam_out,is_mut_out
    #global is_fasta_out,is_fastq_out,is_sam_out,is_mut_out

    from RG_exploder_globals import is_write_ref_fasta, is_write_ref_ingb, is_use_absolute
    #global is_write_ref_fasta, is_write_ref_ingb, is_use_absolute

    from RG_exploder_globals import is_CDS
    return
# End of initialise_flags_output_configs()

def initialise_global_configs():
    # No global definitions, where they can be amended by GUI module as need to use RG_globals.* instead

    from RG_exploder_globals import target_locus # defines the target locus name and file-name components.

    from RG_exploder_globals import mutlabels,mutfreqs # names of variants and their 'frequencies' 
 
    from RG_exploder_globals import Fraglen,Fragdepth  # maximum length of a fragment and maximum depth-of-cover

    from RG_exploder_globals import Exome_extend # Extends any splice join boundary by this value where splice_trigger is found
                                                 # eg +3 nt (not negatives) each side of an intron/exon boundary within the "exome" sequence

    # Where should these live?
    from RG_exploder_globals import ensembl_geneid #  variable assigned on reading Refseq
    from RG_exploder_globals import target_transcript_name,target_transcript_id
    return
# End of initialise_global_configs()

def initialise_global_consts():
    # Retain as constants to this module from config file, as these are not changeable by GUI
    from RG_exploder_globals import infilepathroot,outfilepathroot #  CONSTANTS: defines the filepaths for location of input and output files.
    #from RG_exploder_globals import infilepathroot,infilepath,outfilepathroot,outfilepath

    from RG_exploder_globals import Seq_Format,Seq_IO_file_ext # CONSTANTS: defines Genbank or Embl as file format; path to refseq file
    global Seq_Format,Seq_IO_file_ext

    from RG_exploder_globals import Ref_file_name
    global Ref_file_name

    from RG_exploder_globals import journhead,journend,readmehead,readmeend # CONSTANTS: stub names for journal and readme files
    global journhead,journend,readmehead,readmeend

    from RG_exploder_globals import MaxVarPos # The maximum length of Refseq used.
    # Unusual exception: Local global MaxVarPos is reset to RG_globals.MaxVarPos in initialise_module_globals at each new reference-file load

    from RG_exploder_globals import Qualmin,Qualmax
    global Qualmin,Qualmax  #  Define the range of FASTQ quality values that are set at random in get_quality_list

    from RG_exploder_globals import NormaliseUpper
    global NormaliseUpper # Figure for normalising mutfreqs values

    from RG_exploder_globals import Mapq_Min,Mapq_Max # Used in write_frag_samout range mapq=randint(Mapq_Min,Mapq_Max)
    global Mapq_Min,Mapq_Max

    from RG_exploder_globals import DOC_precision # DOC reporting in journal, this is number of decimal places to calculate to 
    global DOC_precision 

    from RG_exploder_globals import is_include_indels,is_vars_to_lower, is_dels_to_dots, is_dels_to_dots_override
    global is_include_indels,is_vars_to_lower, is_dels_to_dots, is_dels_to_dots_override  # Diagnostics

    from RG_exploder_globals import is_print_diagnostics   
    #RG_globals.is_print_diagnostics=True  #  Switch to print more verbose runtime messages for diagnostic purposes. Rarely used now

    from RG_exploder_globals import Paired_insert_Min,Paired_insert_Max
    global Paired_insert_Min, Paired_insert_Max
    
    return
# End of initialise_global_consts()


# ===============================================================
# End of initialising routines for global variables and constants
# ===============================================================


# =================================================
# Start of initialising routines for IO
# =================================================
def set_IO_filenames():
    global Seq_IO_file_ext # Predefined global constant
    global Outfilepath,Infilepath

    # Set the Outfilepath as root, plus gene name.
    Outfilepath=RG_globals.outfilepathroot+RG_globals.target_locus+"/"

    # Set the Infilepath as root, plus gene name.
    Infilepath=RG_globals.infilepathroot+RG_globals.target_locus+"/"

    global locus_transcript # Added May 2021
    locus_transcript=RG_globals.get_locus_transcript()

def initialise_start_time():
    # Initialise global variables for file IO
    global Start_time
    Start_time=time.time()
# End of initialise_global_IO_vars()

def initialise_module_outflags():
    # ''' Initialise all is_append_x (file-handling) flags to False.
    # ''' These are set to True when an output file is opened, ready to be appended-to . Could call these is_open_x I suppose.

    global fastaout,fastqout,fastqout_R1,fastqout_R2,mutout, samout, journout  # ''' These are global file handles '''
    global is_append_fastafile, is_append_fastqfile      # ''' These are file-open flags for fasta & fastq output files '''
    is_append_fastafile=False; is_append_fastqfile=False

    global is_append_samfile, is_append_mutfile, is_append_mutprotfile  # ''' These are file-open flags for SAM & {target}_mut.fasta output files '''
    is_append_samfile=False ; is_append_mutfile=False; is_append_mutprotfile=False

    # COD 6-Mar-2021: These statements also at head of this file for direct execution on first call
    global is_append_journalfile # ''' This is a file-open flag for output journal file '''
    is_append_journalfile=False

    global is_append_journalfile_htm # ''' This is a file-open flag for output journal html file '''
    is_append_journalfile_htm=False

    global is_append_readmefile
    is_append_readmefile=False ## ''' This is a file-open flag for output readme file '''

    return
# End of initialise_module_outflags()

def initialise_module_outfiles():
    #Sequence output filename derivations.
    # Don't know at this stage whether they will be written-to (eg:run may generate zero sequences), so do not open at this point.
    global out_mut_file_head,out_fa_file, out_fa_file_R1,out_fa_file_R2,out_fq_file,out_fq_file_R1,out_fq_file_R2
    global locus_transcript
    global pair_monitor_out,pair_monitor_file

    extralabel,strandlabel=RG_globals.get_strand_extra_label()
    out_mut_file_head="%s%s"%(locus_transcript,extralabel)
    
    out_fa_file="%s%s_%ss.fasta"%(locus_transcript,extralabel,RG_globals.read_annotation)
    out_fa_file_R1="%s%s_%ss_R1.fasta"%(locus_transcript,extralabel,RG_globals.read_annotation)
    out_fa_file_R2="%s%s_%ss_R2.fasta"%(locus_transcript,extralabel,RG_globals.read_annotation)    
    
    out_fq_file="%s%s_%ss.fastq"%(locus_transcript,extralabel,RG_globals.read_annotation)
    out_fq_file_R1="%s%s_%ss_R1.fastq"%(locus_transcript,extralabel,RG_globals.read_annotation)
    out_fq_file_R2="%s%s_%ss_R2.fastq"%(locus_transcript,extralabel,RG_globals.read_annotation)    

    global in_ref_src,in_ref_src_title
    in_ref_src="%s_%s"%(RG_globals.target_locus,Ref_file_name)
    #in_ref_src_title="Primary Source"
    in_ref_src_title="'%s'"%RG_globals.reference_haplotype
     
    global out_sa_file,Out_Ref_Source
    out_sa_file="%s%s_%ss.sam"%(locus_transcript,extralabel,RG_globals.read_annotation)
    plabel="%s"%bio_parameters["is_frg_paired_end"]["label"].split("-")[0].lower()
    #Out_sam_ref="%s_%s"%(RG_globals.target_locus,Ref_file_name)

    #Out_Ref_Source=in_ref_src # Possible alternative
    #Out_Ref_Source="%s-%s"%(RG_globals.target_locus,RG_globals.empty_transcript_name.lower())
    Out_Ref_Source="%s%s_REF"%(locus_transcript,strandlabel)
    if RG_globals.is_frg_paired_end: # Over-riding for SAM format reference
        Out_Ref_Source=in_ref_src
        if RG_globals.is_pair_monitor:
        # Adding pair_monitor
            pair_monitor_file="%s%s_%ss_pair_monitor.txt"%(locus_transcript,extralabel,RG_globals.read_annotation)
            pair_monitor_out = RG_io.open_write("%s%s"%(Outfilepath,pair_monitor_file), 0)
            pair_monitor_out.write("This pair_monitor file %s\n"%pair_monitor_file)
            pair_monitor_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%("Fragcount","start","insert","R1_start","R2_start","R1_end","R2_end","R1_frag","R2_frag"))

    global out_mut,out_mutprot
    #mutxt=RG_globals.bio_parameters["is_mut_out"]["label"].split(" ")[1].lower()
    mutxt="hapvar"
    out_mut="%s%s_%ss"%(locus_transcript,strandlabel,mutxt)

    #outfile="%s_%s%s"%(locus_transcript,id_label,strandlabel)  - see below consistency?
    out_mutprot="%s%s_%ss_prot"%(locus_transcript,strandlabel,mutxt)
    
    return
# End of initialise_module_outfiles()

def initialise_journal():
    #Prior declared globals
    global Progver,ProgverDate
    global journhead,journend,readmehead,readmeend,journal_name,readme_name
    global Outfilepath
    global locus_transcript
    global is_htm_journal

    #Opens the, metadata, journal file
    global journout,is_append_journalfile,is_joindata_in_json
    if not is_append_journalfile:
        extralabel,strandlabel=RG_globals.get_strand_extra_label()
        journal_name="%s%s%s%s"%(locus_transcript,extralabel,journhead,journend)
        readme_name="%s%s%s%s"%(locus_transcript,extralabel,readmehead,readmeend)
        is_append_journalfile=True
        journout = RG_io.open_write("%s%s"%(Outfilepath,journal_name), 0)
        journout.write("This journal file %s"%journal_name)
        if is_htm_journal:
            initialise_journal_htm()
        journals_update(" is created by %s %s starting on %s"%(Progver,ProgverDate,RG_globals.getime()))
        journals_update("Read in conjunction with %s"%readme_name)
        journals_update("User ID:%s"%(RG_globals.CustomerIDText))
        journals_update("Data set:%s"%(RG_globals.DatasetIDText))
        journals_update("Selected %s: %s; Selected %s: %s; Selected %s: %s"%(bio_parameters["target_locus"]["label"],
                                                           RG_globals.target_locus,
                                                           bio_parameters["target_transcript_name"]["label"],
                                                           RG_globals.target_transcript_name,
                                                           bio_parameters["is_CDS"]["label"],
                                                           RG_globals.is_CDS))
        journals_update("If this is the last line, then something has gone wrong reading the source files")
    return
# End of initialise_journal()

def initialise_journal_htm():
    global Outfilepath,journout_htm,is_append_journalfile_htm,journal_name
    journout_htm = RG_io.open_write("%s%s.htm"%(Outfilepath,journal_name), 0)
    is_append_journalfile_htm=True
    outstring="<!DOCTYPE html>\n<html>\n<head><title>%s.htm</title></head>\n<body>\n<pre>"%journal_name
    journout_htm.write(outstring)
    journout_htm.write("This journal file %s.htm"%journal_name)
# End of initialise_journal_htm()

def close_journal_htm():
    global journout_htm,is_append_journalfile_htm
    if is_append_journalfile_htm:
        outstring="</pre></body></html>\n"
        journout_htm.write(outstring)
        journout_htm.close()
        is_append_journalfile_htm=False

def journals_update(instring):
    update_journal_htm(instring)
    update_journal(instring)

def update_journal_htm(instring):
    #Prior declared globals
    global journout_htm,is_append_journalfile_htm
    #Updates the metadata htm file with progress report
    if is_append_journalfile_htm:
        journout_htm.write("%s\n"%instring)
    return
# end of update_journal_htm(instring)

def update_journal(instring):
    #Prior declared globals
    global journout,is_append_journalfile,is_htm_journal
    #Updates the metadata file with progress report
    if is_append_journalfile:
        journout.write("%s\n"%instring)
    return
# end of update_journal(instring)

def initialise_readme():
    # Prior declared globals
    global Progver,ProgverDate
    global journal_name,readme_name
    global Outfilepath
    global locus_transcript
    #Writes a metadata file explaining the files present in output directory
    #New declared globals
    global readmeout,is_append_readmefile
    if not is_append_readmefile:
        extralabel,strandlabel=RG_globals.get_strand_extra_label()
        readmeout = RG_io.open_write("%s%s"%(Outfilepath,readme_name), 0)
        is_append_readmefile=True
        update_readme("This readme file %s is written by Program %s %s starting on %s"%(readme_name,Progver,ProgverDate,RG_globals.getime()))
        update_readme("Read in conjunction with %s"%(journal_name))
        update_readme("Program input files:")
    return
# end of def initialise_readme()

# =================================================
# End of initialising routines for IO
# =================================================

# =================================================
# Sequence-input routines
# =================================================

# Use this instead of the one below to allow for reading of gzip files when using Python
# This is **not** suitable for webapp repository because there is no support for gzip files
#def read_single_record_input(infile,seqform):
#    seq_record,success=RG_io.read_single_record_input_gz_support(infile,seqform)
#    return seq_record,success

# This is safe to use in Python or for webapp, but assumes no input files are gz  
def read_single_record_input(infile,seqform):
    with RG_io.open_read(infile) as f:
        seq_record = SeqIO.read(f,seqform)
        seq_record,success = RG_process.add_extra_objects(seq_record)
        #print("infile %s, seq_record.offset %s"%(infile,seq_record.offset))
    return seq_record,success
# end of read_single_record_input(infile,seqform)

def set_maxrefs(SEQ_record):
    global MaxRefPos,MaxVarPos,in_ref_src_title
    MaxRefPos=len(SEQ_record.seq)
    if MaxVarPos > 0:  # NB: Unresolved bug when setting MaxVarPos > 0
        msgtxt="programmed constant of first"
        SEQ_record.howmany="first"
        SEQ_record=RG_process.modify_seq_in_record(SEQ_record.seq[0:MaxVarPos],SEQ_record)
    else:
        SEQ_record.howmany="all"
        MaxVarPos=MaxRefPos
        msgtxt="full sequence length of"
    journals_update(" MaxVarPos set to %s %s bases from %s %s"%(msgtxt,MaxVarPos,in_ref_src_title,in_ref_src))
    return SEQ_record

def read_refseqrecord(embl_or_genbank):
    global in_ref_src_title,in_ref_src,Ref_file_name,Seq_IO_file_ext # Inherited globals
    # Deliberately made to look like read_mutrecord before calling read_seqrecord - allows hiding of full filepath
    ingb_file=in_ref_src+Seq_IO_file_ext
    SEQ_record,exists = read_seqrecord(ingb_file,embl_or_genbank,Ref_file_name,in_ref_src_title)

    if not exists:
        program_exit("%s file %s not found"%(in_ref_src_title,ingb_file))
    elif len(SEQ_record.seq) < 1:
        exists=False
        program_exit("%s file %s has no sequence"%(in_ref_src_title,ingb_file))
    else:
        journals_update(" Locus range for %s defined as: %s" %(in_ref_src_title,SEQ_record.firstid.replace("chromosome:","")))
        journals_update(" %s  range for %s defined as: %s" %(RG_globals.empty_transcript_name,in_ref_src_title,SEQ_record.locus_range))

        journals_update(" Data taken from %s and retrieved %s"%(RG_globals.Reference_sequences[RG_globals.target_locus]["Release"],RG_globals.Reference_sequences[RG_globals.target_locus]["Retrieval_date"]))

        if RG_globals.is_flip_strand:
            journals_update(" NB: %ss will be reverse-complemented prior to fragmentation because '%s' selected"%(RG_globals.variants_label,bio_parameters["is_flip_strand"]["label"]))
            
        numvarfeats=len(RG_process.get_varfeature_index(SEQ_record))
        if numvarfeats > 0:
            exists=False
            program_exit(" *** Variant feature definitions found in %s ***"%in_ref_src_title)
        else:
            journals_update(" %s %s correctly includes %s variant features"%(in_ref_src,in_ref_src_title,numvarfeats))
            SEQ_record=set_maxrefs(SEQ_record)
            ''' By definition Refseq CIGAR is all Ms'''
            SEQ_record.mutlabel=in_ref_src
            SEQ_record.cigar=str(MaxVarPos)+"M"
            SEQ_record.cigarbox=[0, int(MaxVarPos), 'M', int(MaxVarPos)]
            SEQ_record.mutbox=RG_process.haveamutbox(SEQ_record.cigarbox)
            # Copying in a modified version of RG_globals.bio_parameters["target_build_variant"]["exonplus_lookup"] which will be inherited by each variant entry
            # This is later modified in each variant haplotype to create a lookup table for creating paired-end exomic reads 
            SEQ_record.exonplus_lookup=copy.copy(RG_globals.exonplus_lookup)
            if len(SEQ_record.exonplus_lookup)>1:
                SEQ_record.exonplus_lookup.sort()
                SEQ_record.exonplus_lookup.pop(0)
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
            
            exists,SEQ_record=RG_process.set_seqrec_absolutes(SEQ_record) # must do this even if RG_globals.is_use_absolute== False
            if RG_globals.is_write_ref_ingb:
                write_refseq(SEQ_record,"long")
            else:
                journals_update(" NOTE: Feature-definition files for: %s and '%ss', NOT saved because option '%s' is set to %s"
                               %(in_ref_src_title,RG_globals.variants_label,bio_parameters["is_write_ref_ingb"]["label"],RG_globals.is_write_ref_ingb))
    return SEQ_record,exists
# end of read_refseqrecord(embl_or_genbank)

def read_seqrecord(infile,embl_or_genbank,label,msg):
    global is_show_infilepath
    global Infilepath
    full_filename="%s%s"%(Infilepath,infile)
    if is_show_infilepath:
        msgfilename=full_filename
    else:
        msgfilename=infile

    SEQ_record=""
    exists = RG_io.is_file(full_filename)
    if exists:
        journals_update("\nReading %s file %s"%(msg,msgfilename))
        update_readme("%s - %s"%(msgfilename,msg))
        # second returned parameter "exists" is really a "success" for the file-read, but want to overwrite "exists" to return an overall success/fail condition here
        SEQ_record,exists = read_single_record_input(full_filename,embl_or_genbank)
        if exists:
            SEQ_record.mutlabel=label
        else:
            journals_update("Failure reading ACCESSION line")
    else:
        journals_update(" \n*** %s file not found: %s ***"%(msg,msgfilename))
        exists=False
    return SEQ_record,exists
# end of read_seqrecord(infilename,embl_or_genbank,label,msg)

def grab_mut_annotation(SeqRec,label):
    global in_ref_src
    annotation="Feature definitions for %s '%s'"%(RG_globals.variants_label,label)
    join_ann=""
    introntxt=" "
    if SeqRec.is_ref_spliced:
        if (SeqRec.splicecount-SeqRec.endclipcount) >0:
            introntxt=" and %s introns "%(SeqRec.splicecount-SeqRec.endclipcount)
        ref_splice_ann=", with %s end-trims%sderived from source %s"%(SeqRec.endclipcount,introntxt,in_ref_src)
    else:
        ref_splice_ann=""

    if RG_globals.is_use_absolute:
        abs_ann=" absolute positions"
    else:
        abs_ann=""
    if SeqRec.is_ref_spliced and RG_globals.is_use_absolute:
        join_ann=" and"
    annotation=annotation+ref_splice_ann+join_ann+abs_ann+" added"
    return annotation

def read_mutrecord(label,embl_or_genbank):
    global Seq_IO_file_ext
    global locus_transcript
    inlabel=RG_globals.target_locus+"_"+label
    ingb_file=inlabel+Seq_IO_file_ext
    
    #Mut_record,exists = read_seqrecord(ingb_file,embl_or_genbank,label,"'%s'"%RG_globals.variants_header)
    Mut_record,exists = read_seqrecord(ingb_file,embl_or_genbank,label,"'%s'"%RG_globals.variants_label)
    if exists:
        #Write out input data as its non-processed interpretation
        if RG_globals.is_write_ref_ingb:
            out_file="%s_%s%sin"%(RG_globals.target_locus,label,Seq_IO_file_ext)
            readmetxt="\t- Initial feature list for %s '%s'"%(RG_globals.variants_label,label)
            write_gb_features(Mut_record,out_file,readmetxt,"full")
    return Mut_record,exists
# end of read_mutrecord(label,embl_or_genbank)

def merge_ref_with_mut(REF_record,Mut_record,label):
    global out_mut_file_head
    mergefeatures,accept,messages=RG_process.purl_features(REF_record,Mut_record) # Returns a feature table with only the variant features from Mut_Record
                                                                                  # added to all features (excluding variant features) from REF_record
    journals_update(messages) 
    if accept: # Overwrite the features in Mut_record
        Mut_record.features=mergefeatures    
        Mut_record,messages=RG_process.merge_seqvar_records(REF_record,Mut_record,label)
        if messages !="":
            journals_update(messages) # The return message might not be the same for each merge
        # Set absolutes
        exists=RG_process.set_seqrec_absolutes(Mut_record) # must do this even if RG_globals.is_use_absolute== False
        # Now write out input variant file as its post-merged interpretation, if modified
        # It will only differ if there is a 'gap' feature folded in (which is always since end-trims added), or absolute-positions added to feature table
        if (Mut_record.is_ref_spliced or RG_globals.is_use_absolute) and exists:
            if RG_globals.is_write_ref_ingb:
                annotation=grab_mut_annotation(Mut_record,label)
                out_file="%s_%s%sout"%(out_mut_file_head,label,Seq_IO_file_ext)
                #out_file="%s_%s%sout"%(locus_transcript,label,Seq_IO_file_ext)
                write_gb_features(Mut_record,out_file,"\t- %s"%annotation,"long")
    return Mut_record,accept
# ========================================
# End of sequence-input routines
# ========================================


# =====================================================
# Start of sequence-output writing routines definitions
# =====================================================
def write_gb_features(SeqRec,out_file,readmetxt,style):
    #Prior declared globals
    global Seq_Format
    global Outfilepath

    '''
    When style is "long" or "short":
    a) Clips the sequence to 0 nucleotide and retains the filtered variants list to write out a truncated variant-features-only Genbank file
    b) Save a copy of the variants for each source in the run. Keep these with output.

    otherwise:
    c) saves original file to output directory unchanged (apart from python reformating of location numbering for inserts)

    Included to meet requirement to log full state of program prior to run. '''

    # Several fields in SeqRec.annotations seems to be ignored by SeqIO.write, so 'date' always comes out as '01-JAN-1980'
    # I have tried!
    # Any COMMENT fields, likewise.
    
    if style == "full":
        ''' Point to the full sequence record - the unchanged mutseq'''
        CopySeqRec=SeqRec
    else:
        # Style is "short" or "long", so ...
        #Take clipped sequence and only the base annotations into a new sequence record
        CopySeqRec=RG_process.modify_seq_in_record(SeqRec.seq[0:0],SeqRec)
        if style =="short":
            #Add the filtered-variant features only -> the mutseq
            CopySeqRec.features=RG_process.get_varfeatures(SeqRec)

    out_gbfile="%s%s"%(Outfilepath,out_file)
    journals_update(" Writing %s"%out_file)
    gbout = RG_io.open_write(out_gbfile)
    out_handle=StringIO()
    SeqIO.write(CopySeqRec,out_handle,Seq_Format)
    out_data = out_handle.getvalue()
    gbout.write(out_data)
    out_handle.close()
    gbout.close()
    if readmetxt !="no":
        update_outfilestring("%s %s"%(out_file,readmetxt))
    return
# end of write_gb_features(SeqRec,out_file,readmetxt,style)

def seqwrite(seqrec,outfile,seqtype):
    out_handle=StringIO()
    SeqIO.write(seqrec,out_handle,seqtype)
    out_data = out_handle.getvalue()
    out_handle.close()    
    outfile.write(out_data)

def write_frag_fastout(seqrec,is_R1):
    #Opens or appends a sequence entry to two new files/overwrites in Fasta, Fastq        
    if RG_globals.is_fasta_out:
        write_fastaout(seqrec,is_R1)
    if RG_globals.is_fastq_out:
        write_fastqout(seqrec,is_R1)
    return
# end of def write_frag_fastout # was def write_fragseqout

def write_fastaout(ThisSeqr,is_R1):
    #Writes a sequence entry to Fasta format. First call sets the is_append_fastafile flag
    global fastaout,fastaout_R1,fastaout_R2,is_append_fastafile,Outfilepath,out_fafile,out_fa_file_R1,out_fa_file_R2
    if not is_append_fastafile:
        is_append_fastafile=True
        if RG_globals.is_frg_paired_end:
            fastaout_R1 =  RG_io.open_write("%s%s"%(Outfilepath,out_fa_file_R1))
            fastaout_R2 =  RG_io.open_write("%s%s"%(Outfilepath,out_fa_file_R2))
            journals_update("Writing FASTA paired end %ss to %s and %s"%(RG_globals.read_annotation,out_fa_file_R1,out_fa_file_R2))
        else:
            fastaout = RG_io.open_write("%s%s"%(Outfilepath,out_fa_file))
            journals_update("Writing FASTA %ss to %s"%(RG_globals.read_annotation,out_fa_file))
            
    if not RG_globals.is_frg_paired_end:
        seqwrite(ThisSeqr,fastaout,"fasta")
    elif is_R1:
        seqwrite(ThisSeqr,fastaout_R1,"fasta")
    else:
        seqwrite(ThisSeqr,fastaout_R2,"fasta")
    return
# end of def write_fastaout


def write_fastqout(ThisSeqr,is_R1):
    #Writes a sequence entry to Fastq format. First call sets the is_append_fastqfile flag
    global fastqout,fastqout_R1,fastqout_R2,is_append_fastqfile,Outfilepath,out_fqfile,out_fq_file_R1,out_fq_file_R2
    if not is_append_fastqfile:
        is_append_fastqfile=True
        if RG_globals.is_frg_paired_end:
            fastqout_R1 =  RG_io.open_write("%s%s"%(Outfilepath,out_fq_file_R1))
            fastqout_R2 =  RG_io.open_write("%s%s"%(Outfilepath,out_fq_file_R2))
            journals_update("Writing FASTQ paired end %ss to %s and %s"%(RG_globals.read_annotation,out_fq_file_R1,out_fq_file_R2))
        else:
            fastqout = RG_io.open_write("%s%s"%(Outfilepath,out_fq_file))
            journals_update("Writing FASTQ %ss to %s"%(RG_globals.read_annotation,out_fq_file))

        if RG_globals.Qualmin!=RG_globals.Qualmax:
            if RG_globals.is_fastq_random:
                msgtxt=""
            else: 
                msgtxt="sort-of-"
            msgtxt="%srandomly-assigned a value between %s and %s in each %s"%(msgtxt,RG_globals.Qualmin,RG_globals.Qualmax,RG_globals.read_annotation)
        else:
            msgtxt="identical throughout at a value of %s"%RG_globals.Qualmax

        journals_update("        FASTQ quality scores are %s"%msgtxt)

    if not RG_globals.is_frg_paired_end:
        seqwrite(ThisSeqr,fastqout,"fastq")
    elif is_R1:
        seqwrite(ThisSeqr,fastqout_R1,"fastq")
    else:
        seqwrite(ThisSeqr,fastqout_R2,"fastq")
    return
# end of def write_fastqout


def writemut(seq_record,label):
    global out_mut,Outfilepath,is_append_mutfile, mutout
    #Write mutated reference sequence
    if RG_globals.is_fastacigar_out:
        ''' # Just leaving this here in case it becomes a requirement
        if RG_globals.is_trim_to_gene:
            front_trim=str(Headclip)+"N"
            end_trim=str(Endclip)+"N"
            high=seq_record.cigar.rfind(end_trim)
            #cigar_txt=seq_record.cigar[:high]
            #cigar_txt=cigar_txt[len(front_trim):]
            cigar_txt=seq_record.cigar[len(front_trim):high]
        else:
            cigar_txt=seq_record.cigar
        print("seq_record..endclipcount:%s; front_trim:%s; end_trim:%s; high:%s; cigar_txt:%s; seq_record.cigar:%s"%(seq_record.endclipcount,front_trim,end_trim,high,cigar_txt,seq_record.cigar))

        seq_record.description=seq_record.description+" "+cigar_txt
        '''
        seq_record.description=seq_record.description+" "+seq_record.cigar
    journals_update(" Writing %s_%s FASTA to %s.fasta"%(locus_transcript,label,out_mut))
    if not is_append_mutfile:
        out_mutfile="%s%s.fasta"%(Outfilepath,out_mut)
        mutout= RG_io.open_write(out_mutfile)
        is_append_mutfile=True
    out_handle=StringIO()
    SeqIO.write(seq_record,out_handle,"fasta")
    out_data = out_handle.getvalue()
    mutout.write(out_data)
    out_handle.close()
    if RG_globals.is_CDS and RG_globals.Exome_extend == 0 and (RG_globals.target_transcript_name != RG_globals.empty_transcript_name) and not RG_globals.is_frg_paired_end:
        writeprotmut(seq_record,label)
    return
# end of def writemut

def writeprotmut(seq_record,label):
    global Outfilepath,out_mutprot,is_append_mutprotfile, mutprotout
    is_do_complement = (RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"] and not RG_globals.is_flip_strand) \
                       or (not RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"] and RG_globals.is_flip_strand)
    #Write mutated reference sequence as protein in defined circumstances
    journals_update(" Writing %s FASTA protein to %s"%(label,out_mutprot))
    if not is_append_mutprotfile:
        out_mutprotfile="%s%s.fasta"%(Outfilepath,out_mutprot)
        mutprotout= RG_io.open_write(out_mutprotfile)
        is_append_mutprotfile=True

    mod_seqlen= len(seq_record.seq) % 3
    #VSeq=MutableSeq(str(seq_record.seq),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #VSeq=MutableSeq(str(seq_record.seq))
    VSeq=Biopython_fix.fix_MutableSeq(seq_record.seq)
    
    # If CDS is from a complement, or it isn't and we have switched polarity, then get sense strand
    if is_do_complement:
        #VSeq.reverse_complement(inplace=True)
        Biopython_fix.fix_reverse_complement(VSeq)
    if mod_seqlen >0:
        VSeq=VSeq+"N"*(3-mod_seqlen)
    #ThisSeqr=SeqRecord(Seq(str(VSeq),generic_dna))#generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #ThisSeqr=SeqRecord(Seq(str(VSeq)))
    ThisSeqr=Biopython_fix.fix_SeqRecord(VSeq)
    prot_record=SeqRecord(seq = ThisSeqr.seq.translate(to_stop=True), \
             id = seq_record.id +"_1p", \
             description = "translated sequence")
    
    out_handle=StringIO()
    SeqIO.write(prot_record,out_handle,"fasta")
    out_data = out_handle.getvalue()
    mutprotout.write(out_data)
    out_handle.close()
    return
# end of def writeprotmut

def write_refseq(RefRecord,which):
    #Prior-declared globals
    global Seq_IO_file_ext
    global MaxVarPos
    global Outfilepath
    global Ref_file_name
    global locus_transcript,in_ref_src,in_ref_src_title

    #global VSeqRecREC # testing only
    def writeref_fasta(record,out_file):
        journals_update(" Writing %s.fasta"%(out_file))
        out_filename="%s%s.fasta"%(Outfilepath,out_file)
        refout= RG_io.open_write(out_filename)
        out_handle=StringIO()
        SeqIO.write(record,out_handle,"fasta")
        out_data = out_handle.getvalue()
        refout.write(out_data)
        out_handle.close()
        refout.close()
        return 
        
    #if which=="fasta": # Saving for SAM locseq
    if which==Ref_file_name: # Saving Locseq
        LSeqRec=RG_process.annotate_seq_to_record(RefRecord,RefRecord.seq,in_ref_src,RefRecord.name,Ref_file_name)
        if RG_globals.is_flip_strand:
            journals_update(" Reverse-complementing %s because '%s' selected"%(in_ref_src,bio_parameters["is_flip_strand"]["label"]))
            LSeqRec=RG_process.switch_Rec_polarity(LSeqRec)
            journals_update(" %s %s Range re-defined as: %s" %(in_ref_src,in_ref_src_title,LSeqRec.firstid.replace("chromosome:","")))
            
        LSeqRec.description="%s %i nucleotides from %s"%(LSeqRec.howmany,len(LSeqRec.seq),LSeqRec.firstid)
        writeref_fasta(LSeqRec,Out_Ref_Source)
        journals_update(" %s Location: %s; length: %s %i of %i bases"
                       %(in_ref_src,LSeqRec.firstid.replace("chromosome:",""),LSeqRec.howmany,MaxVarPos,MaxRefPos))
    elif which=="long":
        out_file="%s%sin"%(in_ref_src,Seq_IO_file_ext)
        readmetxt="\t- Feature definitions from the %s %s"%(in_ref_src_title,in_ref_src)
        journals_update("  Feature definitions for %s %s saved as %s"%(in_ref_src_title,in_ref_src,out_file))
        write_gb_features(RefRecord,out_file,readmetxt,"long")

    elif which=="spliced_fasta":
        extralabel,strandlabel=RG_globals.get_strand_extra_label()
        #Create a spliced-sequence version of refseq as the haplotype 'template' when a splice is defined and is_frg_paired_end is False
        if RefRecord.is_ref_spliced and not RG_globals.is_frg_paired_end:
            reftxt="-based %s sequence"%RG_globals.reference_haplotype
            #reftxt2=RG_globals.variants_header
            reftxt2="%ss"%RG_globals.read_annotation
            reftxt3=" in %s, %s and %s"%(out_fa_file,out_fq_file,out_sa_file)
            id_label="refhap"
            if  RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
                reftxt="** Reference Sequence **"
                reftxt2="%ss"%RG_globals.read_annotation
                reftxt3=" in %s, %s and %s"%(out_fa_file,out_fq_file,out_sa_file)
                id_label="REF"
            elif RG_globals.is_CDS:
                reftxt="%s%s"%(RG_globals.CDS_trigger,reftxt)
            else:
                reftxt="%s%s"%(RG_globals.mRNA_trigger,reftxt)
            VSeqRec=make_allvars_in_one_seq(RefRecord,id_label)
            VSeqRec.description=VSeqRec.description+" "+VSeqRec.cigar
            outfile="%s%s_%s"%(locus_transcript,strandlabel,id_label)
            
            #outfile="%s_%s"%(outref,id_label)
            writeref_fasta(VSeqRec,outfile)
            journals_update("  %s, Length: %s bases, is the %s for the %s%s"
                               %(outfile,len(VSeqRec.seq),reftxt,reftxt2,reftxt3))
        else:
            journals_update(" Un-trimmed, and un-spliced, %s %s. Using sequence outside %s boundary for paired ends"%(in_ref_src_title,in_ref_src,bio_parameters["target_transcript_name"]["label"]))

    elif which=="force_Out_Ref_Source": # Forcing save of the correct Reference for mRNA or CDS-based fragments, which is the trimmed Primary Source, never the mRNA or CDS Template
                                        # Should only happen when choosing 'Reads Type' other than paired-ends
        
        # Template and true Reference are the same only when RG_globals.target_transcript_name == RG_globals.empty_transcript_name, so switch locus_transcript to make it happen!
        if RG_globals.target_transcript_name != RG_globals.empty_transcript_name:
            prior_transcript_name=RG_globals.target_transcript_name
            RG_globals.target_transcript_name=RG_globals.empty_transcript_name # swap to force a trimmed-only when calling splice_refseq(RefRecord)
            # NB: Retaining current value for locus_transcript - important part of the force
            journals_update("\n Writing sequence files for both '%s' and '%s' for %s because option '%s' is set to %s and %s '%s' is not '%s'"
                           %(RG_globals.reference_gene,
                             RG_globals.reference_haplotype,
                             locus_transcript,
                             bio_parameters["is_write_ref_fasta"]["label"],
                             RG_globals.is_write_ref_fasta,
                             bio_parameters["target_transcript_name"]["label"],
                             locus_transcript,
                             RG_globals.empty_transcript_name))
            tmp=splice_refseq(RefRecord) # throwaway; adds to journal, but don't need to change it
            RG_globals.target_transcript_name=prior_transcript_name # Now swap-back to restore correct template settings
        
    else:
        journals_update(" Invalid parameter %s to write_refseq"%which)
    return
# end of write_refseq(RefRecord,which)


def write_samheader(ThisSeqr):
    # Writes a header to sam format,unless is_append_samfile already True
    global samout,is_append_samfile,out_sa_file,in_ref_src,Ref_file_name
    global REFSEQ_RECORD,Rname
    Rname="chrX"
    if not is_append_samfile:
        out_samfile="%s%s"%(Outfilepath,out_sa_file)
        samout = RG_io.open_write(out_samfile)
        is_append_samfile=True
        label,genass,chromnum,start,stop,pole= ThisSeqr.firstid.split(":")
        if genass==RG_globals.GRCh37_txt:
            genass="hg19"
        else:
            genass="hg38"
        Rname=in_ref_src # Sets global for reuse by write_frag_samout. Is this correct for RNASeq? Alternative is RG_process.get_outrefname:
        #Rname=RG_process.get_outrefname(REFSEQ_RECORD,Ref_file_name)
        #Rname="chr%s"%chromnum # Sets global for reuse by write_frag_samout
        #Rname="chr%s:%s-%s"%(chromnum,start,stop) # Sets global for reuse by write_frag_samout
        head1="@HD\tVN:1.6\tSO:coordinate\n"
        head2="@SQ\tSN:%s\tLN:%i\tAS:%s\tUR:%s.fasta\n"%(Rname,REFSEQ_RECORD.spliced_length,genass,in_ref_src)
        head3="@CO\tACCESSION:\t%s\n"%(ThisSeqr.firstid)
        head4="@PG\tID:%s\tVN:%s\n"%(Progver,ProgverDate)

        samout.write(head1)
        samout.write(head2)
        samout.write(head3)
        samout.write(head4)
        journals_update("\nWriting %ss in SAM format to %s"%(RG_globals.read_annotation,out_sa_file))
    return
# end of def write_samheader

def write_frag_samout(qname,flag,pos,cigar,rnext,pnext,tlen,forseq):
    # Writes a fragment sequence entry to sam format.
    # pos and pnext parameters sent here are already 1-based; SAM is 1-based
    # Does Rname need to be in format chr:start-end eg: chr2:172936693-172938111 ??
    global samout,is_append_samfile,REFSEQ_RECORD,Ref_file_name
    global Mapq_Min,Mapq_Max,Rname
    #print("ThisSeqr.id %s"%ThisSeqr.id)
    if not is_append_samfile:
        write_samheader(REFSEQ_RECORD)
    if qname =="": qname="barf"
    mapq=randint(Mapq_Min,Mapq_Max)
    qual="*"
    head1="%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\n"%(qname,flag,Rname,pos,mapq,cigar,rnext,pnext,tlen,str(forseq),qual)
    samout.write(head1)
    return
# end of def write_frag_samout

def close_seqout(SeqRec):
    global in_ref_src_title,in_ref_src,out_fa_file,out_fa_file_R1,out_fa_file_R2,out_fq_file,out_fq_file_R1,out_fq_file_R2,out_sa_file,out_mut,out_mutprot,Ref_file_name
    global fastaout,fastqout,fastqout_R1,fastqout_R2,is_append_fastafile,is_append_fastqfile
    global samout,is_append_samfile, mutout, is_append_mutfile, mutprotout, is_append_mutprotfile
    global mutlistlabels,outfilestring
    global locus_transcript,readme_name,journal_name,user_config_file,pair_monitor_file

    extralabel,strandlabel=RG_globals.get_strand_extra_label()
    update_readme("\nProgram output metadata files:\n%s\t- This file"%readme_name)

    update_readme("%s\t- Journal file documenting runtime messages & metadata including %s headers, program parameters"%(journal_name,in_ref_src_title))
    
    if is_htm_journal:
        update_readme("%s.htm\t- html version of Journal file"%journal_name)
    update_readme("%s\t- contains configuration data for this run"%user_config_file)

    if RG_globals.is_frg_paired_end and RG_globals.is_pair_monitor:
        update_readme("%s\t-  contains start-points for each paired read"%(pair_monitor_file))
        pair_monitor_out.close()
        
    sam_fastfile=""
    update_readme(outfilestring)

    update_readme("\nProgram output sequence files:")
    if RG_globals.is_flip_strand:
        flip_extra=", with polarity reversed from '%s'"%RG_globals.reference_gene
        update_readme("NB: polarity of all output sequences reversed from '%s'"%RG_globals.reference_gene)
    else:
        flip_extra=""

    add_extra_ref=" with %s end-trims applied to %s %s sequence"%(SeqRec.endclipcount,in_ref_src_title,in_ref_src)
    add_extra_tem=" with %s end-trims and %s introns applied to %s %s sequence"%(SeqRec.endclipcount,SeqRec.splicecount-SeqRec.endclipcount,in_ref_src_title,in_ref_src)

    if RG_globals.is_write_ref_fasta:
        if RG_globals.is_frg_paired_end:
            update_readme("%s.fasta\t- FASTA file of the un-modified %s %s sequence%s, used to select paired-ends which may lie outside the %s range" %(Out_Ref_Source,in_ref_src_title,in_ref_src,flip_extra,RG_globals.empty_transcript_name))
        else:
            update_readme("%s.fasta\t- FASTA file%s; the Reference Sequence for all %ss%s"%(Out_Ref_Source,add_extra_ref,RG_globals.read_annotation,flip_extra))
        if SeqRec.is_ref_spliced:
            if not RG_globals.empty_transcript_name.lower() in Out_Ref_Source:
                out_spliced=locus_transcript+extralabel+"_refhap.fasta"
                update_readme("%s\t- FASTA file%s; the template sequence for all %s"%(out_spliced,add_extra_tem,RG_globals.variants_label))

    add_txt0="'%s' with frequency >0:"%RG_globals.variants_label

    if is_append_mutfile:
        add_txt1=""
        if SeqRec.is_ref_spliced:
            add_txt1=add_extra_tem
        update_readme("%s.fasta\t- FASTA sequence of each %s%s;%s"%(out_mut,add_txt0,mutlistlabels,add_txt1))
        mutout.close()
        is_append_mutfile = False

    if is_append_mutprotfile:
        update_readme("%s\t- protein translations (frame 1) of sequences in %s"%(out_mutprot,out_mut))
        mutprotout.close()
        is_append_mutprotfile = False

    if SeqRec.is_ref_spliced and ((SeqRec.splicecount-SeqRec.endclipcount) !=0):
        addtxt="spliced "
    else:
        addtxt=""

    sub_select="from each selected %s%s%s"%(addtxt,add_txt0,mutlistlabels)

    if is_append_fastafile:
        is_append_fastafile = False
        
        if RG_globals.is_frg_paired_end:
            sam_fastfile=out_fa_file_R1
            fastaout_R1.close()
            fastaout_R2.close()
            update_readme("%s\t- FASTA sequence R1 %ss %s"%(out_fa_file_R1,RG_globals.read_annotation,sub_select))
            update_readme("%s\t- FASTA sequence R2 %ss %s"%(out_fa_file_R2,RG_globals.read_annotation,sub_select))
        else:
            sam_fastfile=out_fa_file
            fastaout.close()
            update_readme("%s\t- FASTA sequence %ss %s"%(out_fa_file,RG_globals.read_annotation,sub_select))

    if is_append_fastqfile:
        is_append_fastqfile = False

        if RG_globals.is_fastq_random:
            msgtxt=""
        else: 
            msgtxt="sort-of-"
            
        if RG_globals.is_frg_paired_end:
            if sam_fastfile == "":
                sam_fastfile=out_fq_file_R1
                sub_fastq=sub_select
            else:
                sub_fastq= "from %s"%out_fa_file_R1
            fastqout_R1.close()
            fastqout_R2.close()

            update_readme("%s\t- FASTQ sequence R1 %ss %s with a %srandom quality score between %i-%i at each base"
                      %(out_fq_file_R1,RG_globals.read_annotation,sub_fastq,msgtxt,RG_globals.Qualmin,RG_globals.Qualmax))
            update_readme("%s\t- FASTQ sequence R2 %ss %s with a %srandom quality score between %i-%i at each base"
                      %(out_fq_file_R2,RG_globals.read_annotation,sub_fastq,msgtxt,RG_globals.Qualmin,RG_globals.Qualmax))
        else:
            if sam_fastfile == "":
                sam_fastfile=out_fq_file
                sub_fastq=sub_select
            else:
                sub_fastq= "from %s"%out_fa_file
            fastqout.close()
            update_readme("%s\t- FASTQ sequence %ss %s with a %srandom quality score between %i-%i at each base"
                      %(out_fq_file, RG_globals.read_annotation,sub_fastq,msgtxt,RG_globals.Qualmin,RG_globals.Qualmax))   
    
    if is_append_samfile:
        samout.close()
        is_append_samfile = False
        if sam_fastfile != "":
            sam_fastfile=" from %s"%sam_fastfile
        update_readme("%s\t- SAM file of all the sequence %ss%s"%(out_sa_file,RG_globals.read_annotation,sam_fastfile))
        
    if not(RG_globals.is_fasta_out or RG_globals.is_fastq_out or RG_globals.is_sam_out):
        txt=" NB: All 'Output Options: %ss in FASTA/FASTQ/SAM format' are de-selected, so not a lot happened"%RG_globals.read_annotation.capitalize()
        journals_update(txt)
        update_readme(txt)

    elif not REFSEQ_RECORD.Saved_Paired_Fragcount and not REFSEQ_RECORD.Saved_Unpaired_Fragcount:
        txt=" *** No %ss saved to files on this run ***"%RG_globals.read_annotation
        journals_update(txt)
        update_readme(txt)
    return
# end of def close_seqout


# =====================================================
# End of sequence-output writing routines definitions
# =====================================================

# =====================================================
# Start of journal and readme routines definitions
# =====================================================
def update_readme(instring):
    #Updates the metadata file explaining the files present in output directory
    global readmeout,is_append_readmefile
    if is_append_readmefile:
        readmeout.write("%s\n"%instring)
    return
# end of def update_readme(instring)

def update_outfilestring(instring):
    # Accumulates gbin and gbouts
    global outfilestring
    outfilestring+="%s\n"%instring
    return
# end of def update_outfilestrin(instring)

def journalise_dels():
    #Prior declared globals
    global is_dels_to_dots,is_include_indels
    if not is_include_indels:
        journals_update(" WARNING: deletions and insertions will be ** ignored **. Typically test-purposes only")
    if is_dels_to_dots:
        journals_update(" WARNING: deleted bases will be shown as '.' in sequence. Typically test-purposes only, as this destroys position-calculations")
    if RG_globals.is_vars_to_lower:
        journals_update(" WARNING: substituted bases from SNVs, inserts and delins will be shown as lower case in sequence %ss"%RG_globals.read_annotation)
    return
# end of def journalise_dels()


def close_nonseq_files(is_error):
    global is_append_journalfile, is_append_readmefile
    if is_append_journalfile:
        close_journal(is_error)
    if is_append_readmefile:
        close_readme(is_error)
    if is_append_journalfile_htm:
        close_journal_htm()
    return

def journal_final_summary():
    #Prior declared globals
    global Progver,journout,REFSEQ_RECORD
    if not(RG_globals.is_fasta_out or RG_globals.is_fastq_out or RG_globals.is_sam_out):
        journals_update(" ... Rolling tumbleweed ...")
    else:
        readtxt="%s"%RG_globals.read_annotation
        plextxt1="%ss=%s"%(readtxt,REFSEQ_RECORD.Saved_Fragcount)
        endplex1=""
        plextxt2=""
            
        spliceouts=REFSEQ_RECORD.splicecount
        if spliceouts > 0:
            addtxt=", shown in %s CIGAR as 'N'"%readtxt
        else:
            addtxt=""
        journals_update(" %s sections from the Reference Sequence are removed%s"%(spliceouts,addtxt))

        if RG_globals.is_frg_paired_end:
            plextxt=bio_parameters["is_frg_paired_end"]["label"]
            endtxt=""
            plextxt1="%s pairs=%s"%(readtxt,REFSEQ_RECORD.Saved_Paired_Fragcount)
            endplex1="; saved-unpaired %ss = %s, rejected-unpaired %ss = %s"%(readtxt,REFSEQ_RECORD.Saved_Unpaired_Fragcount,readtxt,REFSEQ_RECORD.Unsaved_Unpaired_Fragcount) # Hiding for C Roos version
            plextxt2="; %s=%s"%(bio_parameters["gauss_mean"]["label"],RG_globals.gauss_mean)

        elif RG_globals.is_duplex:
            plextxt=bio_parameters["is_duplex"]["label"]
            endtxt="forward and reverse-complement" 

        else:
            plextxt=bio_parameters["is_simplex"]["label"]
            endtxt="forward only"

        journals_update(" Saved %ss created as %s %s"%(RG_globals.read_annotation,plextxt,endtxt))
        journals_update(" %s=%i%s\n Total number of %s%s"%(bio_parameters["Fraglen"]["label"],RG_globals.Fraglen,plextxt2,
                                                            plextxt1,endplex1))
# end of def journal_final_summary

def close_journal(no_error):
    #Prior declared globals
    global Progver,journout
    if no_error:
        journal_final_summary()
    else:
        journals_update("\n******  Warning or Error report ******")
    journals_update("\nEnding %s at %s "%(Progver,RG_globals.getime()))
    elapsed= time.time()- Start_time
    journals_update("Total time taken:%s"%elapsed)
    journals_update("%s"%RG_globals.CopyrightText)
    journout.close()
    return
# end of def close_journal    

def close_readme(no_error):
    global readmeout,readme_name #  Prior declared filehandle and name for readme file
    global Progver,Start_time
    if not no_error:
        update_readme(" *** Warning or Error report in %s *** "%readme_name)
    update_readme("Ending %s at %s "%(Progver,RG_globals.getime()))
    elapsed= time.time()- Start_time
    update_readme("Total time taken:%s"%elapsed)
    update_readme("%s"%RG_globals.CopyrightText)
    readmeout.close()
    return
# end of def close_readme

# =====================================================
# End of journal and readme routines
# =====================================================

# =====================================================
# Start of basic exit handling
# =====================================================

def program_exit(exitstring):
    # Trapped IO failures should redirect to here
    global is_append_journalfile
    if is_append_journalfile:
        journals_update("\nProgram halted with error: %s "%exitstring)
    else:
        print("\nProgram halted with error: %s "%exitstring)
    # Sending False to the closure routine indicates error
    close_nonseq_files(False)
    if __name__ == "__main__":
        sys.exit()
    else: # Just in case the GUI has a way of handling 
        #sys.exit()
        pass
# end of def program_exit(exitstring)

# =====================================================
# End of basic exit handling
# =====================================================

# =====================================================
# Start of sequence manipulation routines
# =====================================================

def get_mutrecords(REF_record,embl_or_genbank):
    global is_dels_to_dots, is_dels_to_dots_override
    # Protection part 1 - temporarily override is_dels_to_dots flag for this routine
    save_idd=is_dels_to_dots; is_dels_to_dots=is_dels_to_dots_override
    global mutlistlabels # Now used only in close_seqout and journalfrags 

    mutlistlabels=""; mutrecs=[]; mutcount=0 ;is_success=True
    
    def local_add_mutrec(Seq_rec,label,mutfreq):
        global mutlistlabels
        nonlocal mutrecs,mutcount,REF_record
 
        Seq_rec,accept=merge_ref_with_mut(REF_record,Seq_rec,label)
        if accept:
            Seq_rec=make_allvars_in_one_seq(Seq_rec,label)
            if RG_globals.is_mut_out:
                writemut(Seq_rec,label)
            Seq_rec.mutfreq=mutfreq
            Seq_rec.clipped_length=len(Seq_rec.seq)
            mutrecs.append(Seq_rec)
            mutcount+=1
            mutlistlabels=mutlistlabels+" "+RG_globals.target_locus+"_"+label+","
        return accept
   #end of def local_add_mutrec(Seq_rec,label)
        
    #journals_update("\nReading %s (feature) files..."%RG_globals.variants_label)
    journals_update("\nReading '%s' files..."%RG_globals.variants_label)
    if not RG_globals.is_mut_out:
        journals_update(" NOTE: Sequence of '%s'(s) not saved to %s.fasta because option '%s' is set to %s"
                       %(RG_globals.variants_label,out_mut,bio_parameters["is_mut_out"]["label"],RG_globals.is_mut_out))
    
    RG_process.mutfreqs_extend(RG_globals.mutlabels) # First Check, then set, mutfreqs at same length as mutlabels
    addmut_labels= RG_process.get_addmut_labels()
    for label in RG_globals.mutlabels:
        accept= False
        mutfreq=RG_process.get_mutfreq_from_label(label)
        if mutfreq > 0:
            if label in addmut_labels: # Look for user-defined
                Seq_record,exists=RG_process.make_addmut(REF_record,label) 
                if exists:
                    journals_update("\n'%s' %s_%s found as user-defined %s"%(RG_globals.variants_label,RG_globals.target_locus,label,RG_globals.variants_label))
                    accept=local_add_mutrec(Seq_record,label,mutfreq)
                    if not accept:
                        #print("barf")
                        journals_update(" %s definition not accepted"%label)
            else: # Look for data in input directory
                Seq_record,exists = read_mutrecord(label,embl_or_genbank) # Read from input directory
                if exists:
                    accept=local_add_mutrec(Seq_record,label,mutfreq) # Add to mutrecs
                    if not accept:
                        #print("barf")
                        journals_update(" %s definition not accepted"%label)
            if  not exists:
                journals_update(" %s failed to load "%label)
    # end of loop for label in RG_globals.mutlabels:
    #if accept:
    #    mutlistlabels=mutlistlabels[:-1]
    if mutcount < 1:
        #print("special case mutcount %s "%mutcount)
        # Check special case when there is only 1 mutrec and freq is set to zero. Set freq to 1 to avoid nul run. It keeps tripping me up!
        RG_globals.mutfreqs[0]=1
        label=RG_globals.mutlabels[0]
        Seq_record,exists = read_mutrecord(label,embl_or_genbank)
        if exists:
                journals_update(" No haplotypes with frequency >0, so setting the first:%s to %s to get something!"%(RG_globals.mutlabels[0],RG_globals.mutfreqs[0]))
                mutfreq=RG_globals.mutfreqs[0]
                local_add_mutrec(Seq_record,label,mutfreq)
                is_success=exists  
        else:
            program_exit("\t- Unable to read a minimum of 1 %s"%RG_globals.variants_label)
            is_success=False
    mutlistlabels=mutlistlabels[:-1]
    if mutcount > 1:
        addtxt="s"
    else:
        addtxt=""

    journals_update("\nProcessing %s %s%s with frequency >0 : %s"%(mutcount,RG_globals.variants_label,addtxt,mutlistlabels))
    journalise_dels()
    # Protection part 2 - switch back this flag, on leaving this routine
    is_dels_to_dots=save_idd
    return mutrecs,is_success

def MutateVarSeq(VSeq,seq_polarity,this_feature,cigarbox):
    # Take an incoming sequence, and a feature definition: add the feature to the sequence. Update the 'CIGAR box' appropriately
    #This first line of globals would be tidier returned as a list from this function, rather than be globals?
    global var_count, sub_count, insert_count, delete_count, complex_count
    global is_include_indels, skip_count
    global locus_transcript

    #if RG_globals.is_print_diagnostics: print(this_feature)
    #print("Type:%s\nLocation:%s\nQualifiers:%s"%(this_feature.type,this_feature.location,this_feature.qualifiers))
    matchvardef=True  # Default to matching variant definitions are true because it should be set to False when found.
    msgtxt_mvs=""
    is_indel=False; is_insert=False; is_deletion=False; is_complex=False
    # To report Source position for local coord
    feat_trim= -RG_globals.bio_parameters["target_build_variant"]["headclip"]
 
    #xref_ids,reference,replace,replace_string,feat_start,feat_end,feat_polarity,replace_type=RG_process.get_varfeats(this_feature)
    # Longer, but easier-to-follow:
    featurelist=RG_process.get_varfeats2(this_feature)      
    xref_ids=featurelist["xref_ids"];  replace=featurelist["replace"]; replace_string=featurelist["replace_string"]
    feat_start=featurelist["feat_start"]; feat_end=featurelist["feat_end"];feat_polarity=featurelist["feat_polarity"];replace_type=featurelist["replace_type"]
    '''  Try using this instead? (although you'll get the lot) :  for k, v in d.iteritems(): locals()[k]=v'''
    mod_abs_start,mod_abs_end=RG_process.get_absolute_pairs(feat_start,feat_end,REFSEQ_RECORD.offset,REFSEQ_RECORD.polarity)
    
    def do_delete():
        ''' Splice or delete, same thing, apart from 'N' vs 'D'
            First check that feat_end is not beyond end of sequence range for CIGAR build
            Otherwise(Bio)Python copes with subsetting out-of-range. '''
        global skip_count,var_count
        nonlocal VSeq,cigarbox,feat_start,feat_end,replace_string,msgtxt_mvs
        
        # This block is a verification check for delete string match/mis-match
        matchtxt1=""
        matchtxt2=""
        delete_string=VSeq[feat_start:feat_end]
        intended_delete=replace_string.split("/")[0]
        
        if feat_polarity!=seq_polarity:
            #delete_string.reverse_complement(inplace=True)
            Biopython_fix.fix_reverse_complement(delete_string)
        if intended_delete !="N":
            if str(delete_string) != intended_delete:
                matchtxt1="*mis"
                matchtxt2="*"
                matchvardef=False
            #msgtxt_mvs=" Deletion: %s polarity %s, %smatches definition %s polarity +1, loc: %s to %s, abs: %s to %s"%(VSeq[feat_start:feat_end],seq_polarity,matchtxt,delete_string,feat_start+1,feat_end,mod_abs_start,mod_abs_end)
            msgtxt_mvs=" Deletion of: %s polarity %s, %smatches%s sequence %s polarity +1, gen: %s to %s, loc: %s to %s, var: %s to %s"%(intended_delete,seq_polarity,matchtxt1,matchtxt2,delete_string,
                                                                                                                       mod_abs_start,mod_abs_end,feat_start+1+feat_trim,feat_end+feat_trim,feat_start+1,feat_end)
 
        else: # These are the N-gaps / splices
            # was  msgtxt=" Gap: loc: ...            
            #msgtxt_mvs=" %s: loc: %s to %s, abs: %s to %s"%(this_feature.qualifiers['db_xref'],feat_start+1,feat_end,mod_abs_start,mod_abs_end)
            msgtxt_mvs=" %s: length: %s , gen: %s to %s, loc: %s to %s,  var: %s to %s"%(this_feature.qualifiers['db_xref'],feat_end-feat_start,mod_abs_start,mod_abs_end,feat_start+1+feat_trim,feat_end+feat_trim,feat_start+1,feat_end)
                                                                  
        # End of verification check for delete string match/mis-match
        
        if feat_end > cigarbox[len(cigarbox)-1]:
            feat_end = cigarbox[len(cigarbox)-1]
        if is_dels_to_dots and xref_ids[0]!=RG_globals.skip_trigger:
            dots = "."*len(VSeq[feat_start:feat_end] )
        else: dots=""
        VSeq=VSeq[0:feat_start]+dots+VSeq[feat_end:len(VSeq)]

        if xref_ids[0]==RG_globals.skip_trigger:
            # We don't want to call a splice-deletion a variant, but var_count was incremented, so reduce it here
            # N means "reference skipped" & differentiates splice-deletions from variant deletions.
            skip_count+=1; var_count+=-1
            cigarbox=RG_process.fill_cigar(cigarbox,feat_start,feat_end,"N",feat_end-feat_start)
        else:
            cigarbox=RG_process.fill_cigar(cigarbox,feat_start,feat_end,"D",feat_end-feat_start)
        #end of local function do_delete()

    def do_insert():
        '''
        Insertion Looks similar to single substitution, except
        feat_end == feat_start, so force it, rather than use the more obvious string concat
        VSeq=VSeq[0:feat_start]+replace+VSeq[feat_start:len(VSeq)]
        '''
        nonlocal VSeq,cigarbox,feat_start,replace,replace_string,msgtxt_mvs
        # This block is a verification check for insert string match/mis-match
        matchtxt=""
        insert_string=replace

        if feat_polarity!=seq_polarity:
            #insert_string=MutableSeq(insert_string,generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
            #insert_string=MutableSeq(insert_string)
            insert_string=Biopython_fix.fix_MutableSeq(insert_string)
            #insert_string.reverse_complement()# Warning from Python 3.11.0 that this needs replacing for near-future releases 
            #insert_string.reverse_complement(inplace=True)
            Biopython_fix.fix_reverse_complement(insert_string)
            
        if RG_globals.is_vars_to_lower: # Have to convert the replacement back to upper for this comparison
            insert_string=str(insert_string)
            insert_string=insert_string.upper()
 
        if str(insert_string) != replace_string.split("/")[1]:
                matchtxt="mis"
                matchvardef=False
        msgtxt_mvs=" Insert: %s polarity %s, %smatches definition %s polarity +1, gen: %s, loc: %s, var: %s"%(replace,seq_polarity,matchtxt,insert_string,mod_abs_end,feat_end+1+feat_trim,feat_end+1)
        
        # End of verification check for insert string match/mis-match 

        VSeq[feat_start:feat_start]=replace
        cigarbox=RG_process.fill_cigar(cigarbox,feat_start,feat_start,"I",len(replace))
        #end of local function do_insert()

    if RG_globals.is_print_diagnostics: print("start %i end %i length %i"%(feat_start,feat_end,len(VSeq)))
    if (feat_start <= cigarbox[len(cigarbox)-1]) and replace_type!="err":
        # Increase the variant counter
        var_count+=1
        if replace_type=="ins_":
            is_indel=True ; is_insert=True
        elif replace_type=="del_":
            is_indel=True ; is_deletion=True; compseq="-"
        elif replace_type=="delins_":
            is_complex=True
        else:
            # left with replace_type=="sub_"
            # sub_count +=1 is done later, not here
            pass
        if not is_deletion:
            '''
            For non-deletions we need to set the replacement string to reverse complement if the polarities are opposite
            Do reverse-complement because reverse matters when >1 nucleotide
            Where one or both are unset or have the same polarity, we leave them alone.
            '''
            ''' So first, test for identity: '''
            if feat_polarity==seq_polarity:
                compseq=replace
                '''
                The possibilities we are left with are: 1) they are opposite, 2) one set, one unset.
                In this case "unset" is == RG_globals.seq_polarity_none
                So test if either is unset
                '''
            elif seq_polarity==RG_globals.seq_polarity_none or feat_polarity==RG_globals.seq_polarity_none :
                compseq=replace
            else:
                ''' Leaving opposite polarity '''
                #compseq=MutableSeq(replace,generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
                #compseq=MutableSeq(replace)
                compseq=Biopython_fix.fix_MutableSeq(replace)
                #compseq.reverse_complement()# Warning from Python 3.11.0 that this needs replacing for near-future releases 
                #compseq.reverse_complement(inplace=True)
                Biopython_fix.fix_reverse_complement(compseq)
                
            replace=str(compseq)
            if RG_globals.is_vars_to_lower:
                replace=replace.lower()
        # end of conditional: if not is_deletion

        # Now decide which sequence modifications will be made
        if is_complex and is_include_indels:
            # term 'delins' was encountered after this code first written
            # complex variant is expected to be defined like this:
            #   /replace="TGCCAAAACAGT/AG"
            #   requiring a delete followed by an insert at the same position
            # It's just possible there are other 'complex' types that creep in here that might break this
    
            # The insert is performed first, then deletion, to get CIGAR order as nDnI rather than nInD eg: 17D5I instead of 5I17D
            # but to do so, must put the insert 1 base *after* where the deletion will occur
            
            # This line is a pre-verification check for delete string match/mis-match 
            #msgtxt_mvs1=" Delins: %s, loc %s -> %s, abs %s-> %s\n"%(replace_string,feat_start+1,feat_end,mod_abs_start,mod_abs_end)
            msgtxt_mvs1=" Delins: %s: abs %s-> %s, loc %s -> %s\n"%(replace_string,mod_abs_start,mod_abs_end,feat_start+1,feat_end)
            # End of pre-verification check 
            temp_start=feat_start
            feat_start=feat_end
            do_insert()
            msgtxt_mvs1+="         ..."+msgtxt_mvs+"\n"
            feat_start=temp_start
            do_delete()
            msgtxt_mvs1+="         ..."+msgtxt_mvs
            msgtxt_mvs=msgtxt_mvs1
            complex_count +=1
            #if RG_globals.is_print_diagnostics: print("Delins variant processed: %s"%this_feature)
        elif not is_indel:
            #Single substitution. This is where the variant substitution is made in VSeq
            #Expect feat_end == feat_start+1
            
            # Next block is a verification check for SNV string match/mis-match 
            #sub_n=MutableSeq(replace,generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
            #sub_n=MutableSeq(replace)
            sub_n=Biopython_fix.fix_MutableSeq(replace)
            
            current_n=VSeq[feat_start:feat_end]
            if feat_polarity!=seq_polarity: # Have to reverse both because we have done so already and are now turning back to compare to replace_string
                #sub_n.reverse_complement(inplace=True)
                Biopython_fix.fix_reverse_complement(sub_n)
                #current_n.reverse_complement(inplace=True)
                Biopython_fix.fix_reverse_complement(current_n)
                
            if RG_globals.is_vars_to_lower: # Have to convert the replacement back to upper for this comparison
                sub_n=str(sub_n)
                sub_n=sub_n.upper()

            intended_delete=replace_string.split("/")[0]
            matchtxt=""
            if intended_delete !="N":
                #print("replace_string %s, current_n %s, sub_n %s"%(replace_string,current_n,sub_n,))
                if str(current_n) != intended_delete or str(sub_n) != replace_string.split("/")[1] :
                    matchtxt="mis"
                    matchvardef=False
            msgtxt_mvs=" SNV: %s -> %s polarity %s, %smatches definition %s -> %s polarity +1, gen: %s, loc: %s, var:%s"%(VSeq[feat_start:feat_end],replace,seq_polarity,matchtxt,intended_delete,sub_n,mod_abs_end,feat_end+feat_trim,feat_end)
            # End of verification check for SNV string match/mis-match 
            VSeq[feat_start:feat_end]=replace
            #Do CIGAR build
            cigarbox=RG_process.fill_cigar(cigarbox,feat_start,feat_end,"X",1)
            sub_count +=1
        elif is_include_indels:
            # It's an indel: only create the labels & process if is_include_indels==True.
            # When False the labels return empty, which is only important when called from separate_vars_in_multi_seq(now deprecated)
            if is_insert:
                do_insert()
                insert_count+=1 # Has to go here, rather than do_insert(), to avoid clash with complex_count
            else:
                do_delete()
                delete_count+=1 # Has to go here, rather than do_delete(), to avoid clash with complex_count
        else:
            # Looks like we skipped the defined modification - because of is_include_indels=False
            # So dial-down the already-incremented var_count
            var_count+=-1
        # end of series of conditional tests for combination of mutation type and is_include_indels flag
    # end of conditional: if feat_start <= cigarbox[len(cigarbox)-1]
    if RG_globals.is_journal_subs or not matchvardef: # If the variant replacement-definitions don't match the sequence itself, raise note
                                        # is_journal_subs means "document the replacements"
        if not matchvardef:
            msgtxt_mvs="*** WARNING: %s ***"%msgtxt_mvs
        journals_update(msgtxt_mvs)
        
    return VSeq,cigarbox
# end of def MutateVarSeq(VSeq,seq_polarity,this_feature,cigarbox)
#==================================================================
#==================================================================
#==================================================================

def no_splice_refseq(SeqRec):
    # Only called when RG_globals.is_frg_paired_end is True
    # Only modifications to Refseq are setting the clip-lengths and endlocus values, not the feature table and
    # Could actually merge with splice_refseq and add conditionals for RG_globals.is_frg_paired_end
    global out_sa_file,in_ref_src
    
    if RG_globals.is_write_ref_fasta:
        addtxt="is"
        write_refseq(SeqRec,Ref_file_name)
    else:
        addtxt="would be"
    journals_update("  %s.fasta, Length: %s bases, %s the ** Reference Sequence ** for the %ss in %s, %s and %s"\
                   %(Out_Ref_Source,SeqRec.clipped_length,addtxt,RG_globals.read_annotation,out_fa_file,out_fq_file,out_sa_file))

def splice_refseq(SeqRec):
    global in_ref_src_title,in_ref_src,locus_transcript
    
    CopyRec=RG_process.duplicate_record(SeqRec)    
    CopyRec.is_ref_spliced,CopyRec.splicecount,splicedseq,mod_var_features,CopyRec.Headclip,CopyRec.Tailclip=splice_a_sequence(CopyRec,RG_globals.target_locus)

    # splice_a_sequence seems overkill for deriving Headclip and Tailclip values - can get these from the gene & source features
    # Importantly though, splice_a_sequence returns all dbxref=gap: features in  "mod_var_features" that implement the splicing events,
    # of which the two clips are just special cases of splicing
    CopyRec.features=RG_process.get_source_gene_features(SeqRec)+mod_var_features
    
    CopyRec.endclipcount=0                                                                                            
    if CopyRec.is_ref_spliced:
        if CopyRec.Headclip > 0:
            CopyRec.endclipcount+=1
        if CopyRec.Tailclip > 0:
            CopyRec.endclipcount+=1
            
    CopyRec.clipped_length=len(CopyRec.seq)-CopyRec.Tailclip-CopyRec.Headclip
    CopyRec.spliced_length=len(splicedseq)
    CopyRec.endlocus=len(CopyRec.seq)-CopyRec.Tailclip
    CopyRec.rev_endlocus=len(CopyRec.seq)-CopyRec.Headclip
        
    if RG_globals.is_write_ref_fasta:
        write_refseq(CopyRec,"spliced_fasta")
    elif not RG_globals.empty_transcript_name.lower() in locus_transcript:
        if CopyRec.endclipcount !=0:
            addtxt=""
        else:
            addtxt="un"
        journals_update(" NOTE: A reference file for the %ss: '%s_hapvars.fasta', derived from %strimmed %s %s sequence, is not saved because option '%s' is set to %s"
                           %(RG_globals.variants_label,locus_transcript,addtxt,in_ref_src_title,in_ref_src,
                             bio_parameters["is_write_ref_fasta"]["label"],RG_globals.is_write_ref_fasta))
    return CopyRec
# end of def splice_refseq()

def zero_varcounts():
    global var_count, sub_count, insert_count, delete_count, complex_count # counted in MutateVarSeq and probably should not be globals, but returned from that routine
    global skip_count
    skip_count=0
    var_count=0;sub_count=0; insert_count=0; delete_count=0; complex_count=0
    return

def make_allvars_in_one_seq(SeqRec,label):
    # This applies all the variant features in the feature table into the incoming sequence, returning a new seqeunce record with the modified sequence.
    # It does this by applying the variants in reverse order from the tail of the sequence, allowing the original position-definitions to retain relevance.
    # There used to be a companion to this in which each variant was put into its own sequence record. But this was deprecated, there being no requirement for this.
    global REFSEQ_RECORD,MaxVarPos
    global var_count, sub_count, insert_count, delete_count, complex_count
    global skip_count
    global locus_transcript,in_ref_src
    zero_varcounts()

    if SeqRec.firstid != REFSEQ_RECORD.firstid:
        err=True
        msgtxt="WARNING: mis"
    else:
        err=False
        msgtxt="Exact "

    in_label="%s_%s"%(RG_globals.target_locus,label)
    out_label="%s_%s"%(locus_transcript,label)

    #journals_update(" Merging all variants from %s with Source"%(out_label))
    journals_update(" %smatch between Source Ranges for %s and %s"%(msgtxt,in_label,in_ref_src))
    if err:
        journals_update(" %s mismatched Source Range: %s" %(in_label,SeqRec.firstid.replace("chromosome:","")))
    
    #VarSeq=MutableSeq(str(SeqRec.seq),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #VarSeq=MutableSeq(str(SeqRec.seq))
    VarSeq=Biopython_fix.fix_MutableSeq(SeqRec.seq)
    
    cigarbox=[len(VarSeq),len(VarSeq)]
    for (index) in RG_process.get_varfeature_index(SeqRec): # Progressively adds each feature to VarSeq and cigarbox
        VarSeq,cigarbox=MutateVarSeq(VarSeq,SeqRec.polarity,SeqRec.features[index],cigarbox)

    RG_process.pad_cigarbox(cigarbox)# Important - modifies cigarbox!!
    
    subst_txt="End-trim of %s sections; splice-removal of %s sections, from %s. %s variants: %s substitutions; %s inserts; %s deletions; %s delins"\
               %(SeqRec.endclipcount,skip_count-SeqRec.endclipcount,in_ref_src,var_count,sub_count,insert_count,delete_count-skip_count,complex_count)
    VarSeqRec=RG_process.annotate_seq_to_record(REFSEQ_RECORD,VarSeq,out_label,SeqRec.name,label)
    VarSeqRec.description="%s nucleotides from %s: %s"%(len(VarSeq),MaxVarPos,subst_txt)

    #print("IN switch: %s"%VarSeqRec.seq[0:5])
    # Addition of is_flip_strand July 2022.
    if RG_globals.is_flip_strand:
        journals_update(" Reverse-complementing %s because '%s' selected"%(label,bio_parameters["is_flip_strand"]["label"]))
        
        VarSeqRec=RG_process.switch_Rec_polarity(VarSeqRec)
        if label == "REF": # For the case where this is called from write_refseq
            journals_update("  *** %s %s Range re-defined as: %s ***" %(in_ref_src,in_ref_src_title,VarSeqRec.firstid.replace("chromosome:","")))
        cigarbox=RG_process.reverse_cigarbox(cigarbox)

    cigar_label=RG_process.get_untrimmed_cigar(cigarbox) 
    if (cigar_label != ""):
        VarSeqRec.cigar=cigar_label
        VarSeqRec.cigarbox=cigarbox
        VarSeqRec.mutbox=RG_process.haveamutbox(cigarbox) 
    else:
        journals_update(" *** CIGAR is blank - error ***")
        journals_update(" cigar %s, cigarbox %s\n mutbox %s"%(VarSeqRec.cigar,VarSeqRec.cigarbox,VarSeqRec.mutbox))
    journals_update(" %s length: %i bases. %s\n %s CIGAR(wrt %s): %s"%(out_label,len(VarSeq),subst_txt,out_label,in_ref_src,VarSeqRec.cigar))
     #print("OUT switch: %s"%VarSeqRec.seq[0:5])
    return VarSeqRec
# end of def make_allvars_in_one_seq(SeqRec,label)

def refseq_to_frags4(mutrecs):
    global REFSEQ_RECORD
    RG_process.setup_fastq_quality_list()
    is_success= True
    if RG_globals.is_frg_paired_end:
        if RG_globals.gauss_mean < RG_globals.Fraglen:
                journals_update(" *** Processing halted *** Requested %s length %s is higher than the %s of %s"%
                               (RG_globals.read_annotation,RG_globals.Fraglen,RG_globals.bio_parameters["gauss_mean"]["label"],RG_globals.gauss_mean))
                is_success= False
        else:
            generate_multisource_paired_frags(REFSEQ_RECORD,mutrecs,RG_globals.Fraglen,RG_globals.Fragdepth)
    elif RG_globals.is_onefrag_out:
        generate_multisource_all_frags(REFSEQ_RECORD,mutrecs,RG_globals.Fraglen)
    else:
        generate_multisource_random_frags(REFSEQ_RECORD,mutrecs,RG_globals.Fraglen,RG_globals.Fragdepth)
    return is_success
# end of def refseq_to_frags4()
# =====================================================
# End of sequence manipulation routines
# =====================================================

#==========================Start of section dealing with generate_multisource_random_frags and generate_multisource_all_frags===========================

def get_rev_fragseq(inseq):
    # Do a complement on arrivals here
    #outseq=MutableSeq(str(inseq),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #outseq=MutableSeq(str(inseq))
    outseq=Biopython_fix.fix_MutableSeq(inseq)
    #outseq.reverse_complement(inplace=True)
    Biopython_fix.fix_reverse_complement(outseq)
    #print("inseq %s ; outseq %s"%(inseq,outseq))
    outseq=str(outseq)
    return outseq,"r",RG_globals.rev_label

def get_fwd_fragseq(inseq):
    # Just trying symmetry
    return inseq,"f",RG_globals.fwd_label

def get_fragseq(inseq):
    # Do a complement on 50% of arrivals here if RG_globals.is_duplex
    #print("inseq  a %s"%inseq)
    dupstr="f"
    dir_label=RG_globals.fwd_label
    outseq=inseq
    if RG_globals.is_duplex:
        if randint(0,1):#  Use if random() < 0.5 instead?
            outseq,dupstr,dir_label=get_rev_fragseq(inseq)
            #print("outseq b %s"%outseq)
    return outseq,dupstr,dir_label

def generate_multisource_all_frags(RefRec,mutrecs,fraglen):
    # NB: RG_globals.is_onefrag_out must be True to get here AND fraglen must have been tested as OK with validate_fraglength
    # NB: fragtotal is irrelevant for "all_frags"
    
    Generated_Fragcount,Saved_Fragcount,mutfragcount=zero_mutfragcount(mutrecs)

    mutrec_index=0
    mutlabel_count=0
    journals_update(" Relative frequency values ignored because option '%s' selected"%bio_parameters["is_onefrag_out"]["label"])

    for item in mutrecs: # Pick each mutrecs entry in turn to be VarSeq
        for start in range(0,len(item.seq)-fraglen+1):
            #if Generated_Fragcount > 2:
            #    input("enter to continue")
            end=start+fraglen
            Generated_Fragcount+=1
            #if RG_globals.is_duplex: # Increment here because of RG_process.is_mut_cigar filter below
            #   Generated_Fragcount+=1 # Yes: that's +2 for if RG_globals.is_duplex == True
            ref_offset,fwd_cigar_label,rev_cigar_label=RG_process.get_trimmed_cigars(item.cigarbox,item.mutbox,
                                                                             start,fraglen) # Calculate the CIGAR before trimming the sequence
            #if not RG_globals.is_muts_only or (RG_globals.is_muts_only and RG_process.is_mut_cigar(fwd_cigar_label)) : shortens to ...
            if not RG_globals.is_muts_only or RG_process.is_mut_cigar(fwd_cigar_label) : # Only fragment if RG_globals.is_muts_only is False
                                                                                        # or when the fwd_cigar_label includes a variant
                trimseq=item.seq[start:end]  # item.seq[start:end] trims the mutseq to get the fragment
                mutlabel_count+=1
                # Ideally would assign item.Saved_Fragcount+=1 - but this cannot be done in generate_multisource_random_frags, so not doing here
                Saved_Fragcount+=1
                label_and_saveg("f",trimseq,mutlabel_count,start,ref_offset,fwd_cigar_label,rev_cigar_label,item,0,0,True,False,False)
                mutfragcount[mutrec_index]+=1
                if RG_globals.is_duplex: # Just repeat with reverse_complement
                    Saved_Fragcount+=1
                    label_and_saveg("r",trimseq,mutlabel_count,start,ref_offset,fwd_cigar_label,rev_cigar_label,item,0,0,True,False,False)
                    mutfragcount[mutrec_index]+=1                        
        mutrec_index+=1
    journal_frags(mutrecs,Generated_Fragcount,Saved_Fragcount,mutfragcount,0,0)
    return
# end of def generate_multisource_all_frags(RefRec,mutrecs,fraglen)

def zero_mutfragcount(mutrecs):
    # mutfragcount accumulates in def generate_multisource_random_frags or generate_multisource_once_frags as the fragments are generated
    # zeroing it here
    mutfragcount=[]
    for i in range(0, len(mutrecs)):
        if i > len(mutfragcount)-1:
            mutfragcount.append(0)
    Generated_Fragcount=0
    Saved_Fragcount=0
    return Generated_Fragcount,Saved_Fragcount,mutfragcount

def validate_fraglength(mutrecs,fraglen):
    # The variant must be long enough to accommodate fragments
    # In the case of paired-ends, need to account for ends as well
    if RG_globals.is_frg_paired_end:
        mod_fragtest=mutrecs[0].Tailclip+mutrecs[0].Headclip
    else:
        mod_fragtest=0
    longest=0; is_fraglength_OK=True
    shortest=mutrecs[0].clipped_length-mod_fragtest
    short_label=""
    for item in mutrecs:
        thislen=item.clipped_length-mod_fragtest
        #print("thislen %s %s"%(item.mutlabel,thislen))
        if thislen > longest:
            longest = thislen
        elif thislen < shortest:
            shortest=thislen
            short_label=item.mutlabel
    if fraglen*2 > shortest-1:
        is_fraglength_OK=False
    return is_fraglength_OK,shortest,short_label
       
def normalise_mutfreqs(mutrecs):
    global NormaliseUpper,DOC_precision
    doc_mult=10**DOC_precision
    half_adj=0.5/doc_mult

    ''' My first use of Python "list comprehensions" to create a normalised list '''
    ''' mutfreqs defines the ratios of the contributing sequences hap1,hap2, som1 etc when generating fragments.
        The original expectation was two variants hap1 and hap2, representing two parentals, would always be at equal ratio to each other eg: 1:1.
        The other sequences can contribute a ratio of greater or less than one.
        The actual numbers won't matter, as long as ratios are correct.'''
    # Normalise the 'frequencies' to a total of NormaliseUpper
    sumfreq=0; normutfreq=[]
    for item in mutrecs:
        sumfreq=sumfreq+item.mutfreq

    if  sumfreq > 0:
        normutfreq=[int(round(NormaliseUpper*float(item.mutfreq)/sumfreq)) for item in mutrecs]

    min_norm=min(normutfreq)
    mutfreqstring=""; normfreqstring=""  # These strings only for documenting, so can skip zero frequency, as elsewhere
    normratiostring=""
    index=0
    for item in mutrecs:
        mutfreqstring+=str(item.mutfreq)+","
        normfreqstring+=str(normutfreq[index])+","
        ratiofreq=normutfreq[index]/min_norm
        normratiostring+=str(int((ratiofreq+half_adj)*doc_mult)/doc_mult)+","
        index+=1
    journals_update(" Relative frequency values   : [%s]"%(mutfreqstring[:-1]))
    journals_update(" Normalised proportion values: [%s]"%(normfreqstring[:-1]))
    journals_update(" Normalised ratios: [%s]"%(normratiostring[:-1]))
    return normutfreq

def label_and_saveg(in_dupstr,in_fragseq,label_count,start,ref_offset,fwd_cigar_label,rev_cigar_label,RefRec,pnext,tlen,is_read1,is_tworeads,is_paired_end):
    # A common label_and_saveg() for generate_multisource_random_frags, generate_multisource_all_frags, generate_multisource_paired_frags

    id_label="%s%s"%(RG_globals.frag_annotation,label_count)
    flag=3 # +2 aligned (for SAM output); +1 to verify other flag additions (not that I understand this)
    if is_paired_end:
        pnext+=1 # paired-ends
        if is_read1:
            flag+=64 # Read 1 first segment in the template (for SAM output)
        else:
            flag+=128 # Read 2 last segment in the template (for SAM output)

    if not is_tworeads:
        pnext=0 # 
        flag+=8 # next read unmapped, so bitwise for this non=paired read reflects this

    if in_dupstr == "pf":
        loc_fragseq,dupstr,dir_label=get_fwd_fragseq(in_fragseq) # dir_label equivalent to R1 - fwd
        loc_cigar_label=fwd_cigar_label
        if is_tworeads:
            flag+=32  #+32 next segment being reverse complemented (for SAM output)

    elif in_dupstr == "pr":
        loc_fragseq,dupstr,dir_label=get_rev_fragseq(in_fragseq) # dir_label equivalent to R2 - fwd
        loc_cigar_label=rev_cigar_label
        flag+=16  #+16 reverse complement (for SAM output)
         
    elif in_dupstr=="fr":
        loc_fragseq,dupstr,dir_label=get_fragseq(in_fragseq) # 50% chance of returning a "f" or "r" sequence
        if dupstr == "f":
            loc_cigar_label=fwd_cigar_label
        elif dupstr == "r":
            loc_cigar_label=rev_cigar_label
            flag+=16 # reverse complement (for SAM output)
        
    elif in_dupstr == "f":
        loc_fragseq,dupstr,dir_label=get_fwd_fragseq(in_fragseq)
        loc_cigar_label=fwd_cigar_label
        
    elif in_dupstr == "r":
        loc_fragseq,dupstr,dir_label=get_rev_fragseq(in_fragseq)
        loc_cigar_label=rev_cigar_label
        flag+=16 # reverse complement (for SAM output)
    else:
        program_exit("Incorrect call to label_and_saveg")

    if RG_globals.is_frg_label:
        id_label="%s_%s"%(id_label,RefRec.mutlabel)
        FragSeqRec=RG_process.annotate_frag_to_record(loc_fragseq,id_label)
        #ref_offset has been 0-based. Writing out to FASTA & SAM requires 1-based, so +1
        var_begin=start+1
        ref_begin=ref_offset+1
            
        if dupstr == "r" or in_dupstr == "pr":
        # When reverse-complemented: the first base of the fragment is now the last, and the last is the first, so adjust var and reference positions accordingly
            adjust=len(loc_fragseq)-1
            var_begin+=adjust
            ref_begin+=adjust

        var_location=var_begin
        ref_location=ref_begin-RefRec.Headclip
        
        FragSeqRec.description="h:%s%sr:%s"%(var_location,RG_globals.spaceman,ref_location)
        
        if RG_globals.is_use_absolute:
            FragSeqRec.description+="%sg:%s"%(RG_globals.spaceman, RG_process.get_absolute_position(RefRec,ref_begin))

    else:
        FragSeqRec=RG_process.annotate_frag_to_record(loc_fragseq,id_label)
        FragSeqRec.description=""

    if RG_globals.is_fastacigar_out:
        FragSeqRec.description+="%s%s"%(RG_globals.spaceman,loc_cigar_label)

    # Making dir_label constitutive.
    if not RG_globals.is_frg_label:
        FragSeqRec.description+="%s"%(dir_label) # Don't add a space when only the id is present
    else:
        FragSeqRec.description+="%s%s"%(RG_globals.spaceman,dir_label)

    write_frag_fastout(FragSeqRec,is_read1)

    if RG_globals.is_sam_out and in_dupstr != "r":
    # Avoid putting an actual reverse-complement in SAM file!
    # Only the original in_fragseq (forward orientation), is sent to SAM, even if it became reverse-complemented in the get_fragseq call
    # SAM emulates the completed mapping, so only the forward information is retained
        # Do NOT use absolutes in SAM output
        # Make sure to use in_fragseq, not loc_fragseq (*** FragSeqRec.seq holds loc_fragseq***)
        # Use ref_offset+1, for SAM 1-based
        # At this point ref_offset is 0-based but PNEXT is 1-based; 
        if "p" in in_dupstr and is_tworeads:
            write_frag_samout(FragSeqRec.id,flag,ref_offset+1,fwd_cigar_label,"=",pnext,tlen,in_fragseq)
        else:
            write_frag_samout(FragSeqRec.id,flag,ref_offset+1,fwd_cigar_label,"*",pnext,0,in_fragseq)
 
    ###### end of  function label_and_saveg() ########################


def generate_multisource_random_frags(RefRec,mutrecs,fraglen,fragdepth):
    # NB: RG_globals.is_onefrag_out must be False to get here AND fraglen must have been tested as OK with validate_fraglength
    global NormaliseUpper

    Generated_Fragcount,Saved_Fragcount,mutfragcount=zero_mutfragcount(mutrecs)

    # Calculate an expected fragtotal using un-modified spliced-reference
    fragtotal= int(len(mutrecs)*fragdepth*(RefRec.spliced_length/fraglen))
    

    # Get two parallel lists to match mutrecs
    # First a list of normalised frequencies for each mutseq
    normutfreq=normalise_mutfreqs(mutrecs)
    # Next the sequence end of each mutseq, less space for fraglen
    mutseqend=[]
    for item in mutrecs:
        mutseqend.append(len(item.seq)-fraglen)
        #print("mrf: item.id:%s; item.mutlabel:%s; item.cigar:%s; item.mutbox:%s; item.cigarbox:%s"%(item.id,item.mutlabel,item.cigar,item.mutbox,item.cigarbox))
            
    if len(normutfreq)>0:
        # Sum the normalised frequencies progressively to create number-ranges in lottoline
        # lottoline converts the normalised frequencies into a list of cutoff positions for matching a later random number call
        lottoline=[normutfreq[0]]
        for i in range(1, len(normutfreq)):
            lottoline.append(lottoline[len(lottoline)-1]+normutfreq[i])
            
        #  Now do the fragmentation
        while Generated_Fragcount < fragtotal:
            # Pick a random mutrecs entry to be VarSeq, weighted by its normalised values.
            lottowin=randint(1,NormaliseUpper)
            for i in range(0, len(lottoline)):
                if lottowin <= lottoline[i]:
                    mutrec_index=i
                    break
            #print("normalised %i %s %i %i "%(lottowin,RG_globals.mutlabels[mutrec_index],normutfreq[mutrec_index],lottoline[mutrec_index]))
            start=randint(0,mutseqend[mutrec_index])# Lookup of pre-calculated mutseqend -fraglen to save processing
            end=start+fraglen
            Generated_Fragcount+=1
            ref_offset,fwd_cigar_label,rev_cigar_label=RG_process.get_trimmed_cigars(mutrecs[mutrec_index].cigarbox,mutrecs[mutrec_index].mutbox,
                                                                            start,fraglen)  # Calculate the CIGAR before trimming the sequence
            # Only fragment if RG_globals.is_muts_only == False
            #    or, when it's true, when the fwd_cigar_label includes a variant
            # This was originally split into two conditional tests to avoid calling RG_process.is_mut_cigar unnecessarily,
            # but code changes to the first were repeatedly omitted in the second conditional, so reduced to this simpler alternative
            if not RG_globals.is_muts_only or RG_process.is_mut_cigar(fwd_cigar_label) :
                trimseq=mutrecs[mutrec_index].seq[start:end] # Trim the sequence to get a fragment
                mutfragcount[mutrec_index]+=1
                Saved_Fragcount+=1
                label_and_saveg("fr",trimseq,Saved_Fragcount,start,ref_offset,fwd_cigar_label,rev_cigar_label,mutrecs[mutrec_index],0,0,True,False,False)
    journal_frags(mutrecs,Generated_Fragcount,Saved_Fragcount,mutfragcount,0,0)
    return
# end of def generate_multisource_random_frags(RefRec,mutrecs,fraglen,fragdepth)

def generate_multisource_paired_frags(RefRec,mutrecs,fraglen,fragdepth):
    # fraglen must have been tested as OK with validate_fraglength to get here
    global NormaliseUpper,Paired_insert_Max,Paired_insert_Min
    global pair_monitor_out
    
    Generated_Fragcount,Saved_Fragcount,mutfragcount=zero_mutfragcount(mutrecs)
    last_Saved_Fragcount=0
    noread_count=0
    unpaired_count=0
    saved_unpaired=0
    mutseqend=[] # The sequence end of each mutseq
    mutseqendpop=[] # The read boundary
    locus_begin,locus_end=RG_globals.Reference_sequences[RG_globals.target_locus]["Locus_range"].split(":")
    headpop=int(locus_begin) -1 # Adding -1 for zero-based position
    is_one_mut=False
    fraglen_adjust=fraglen-1 # To avoid multiple +1 and -1 recalculations in loop below

    # Generate lookups for each mutseq, such as end-position once rather than recalculate repeatedly
    for item in mutrecs:
        # Length of each mutseq differs, so end-point varies, unlike head
        # Constraining the end position, for a random position, to align with original (locus range - fraglen), not full source sequence
        mutseqend.append(len(item.seq))
        #mutseqendpop.append(len(item.seq)-headpop-fraglen) # Make a lookup to avoid repeat calculations later # oops - wrong for so long
        mutseqendpop.append(len(item.seq)-fraglen) # Make a lookup to avoid repeat calculations later, and zero-based
        #print("mpf: item.mutlabel:%s; item.cigar:%s; item.mutbox:%s; item.cigarbox:%s"%(item.mutlabel,item.cigar,item.mutbox,item.cigarbox))
        #print("exonplus_lookup[0]:%s"%item.exonplus_lookup[0])
        if RG_globals.is_pair_monitor:
            #pair_monitor_out.write("mutseqend:%s\n"%(mutseqend))
            pass
    # End of Generate lookup for end of each mutseq once

    # Generate Booleans once from Globals so do not calculate repeatedly within 'while Generated_Fragcount < fragtotal:' loop
    if (RG_globals.target_transcript_name == RG_globals.empty_transcript_name) or not RG_globals.is_exome_paired_end:
        is_do_exome=False
        #print("is_do_exome %s; head %s; end %s"%(is_do_exome,headpop,mutseqendpop[0]))
    else:
        is_do_exome=True        
        if len(mutrecs) > 1:
            extra="s"
        else:
            extra =""
        journals_update(" *** Adjusting exome parameters on %s variant%s *** "%(len(mutrecs),extra))
        RG_process.make_exonplus_lookups(RefRec.exonplus_lookup,mutrecs)
        #print("is_do_exome %s; head %s; end %s"%(is_do_exome,mutrecs[0].exonplus_lookup[0],mutrecs[0].exonplus_lookup[-1]))
    # End of Generate Booleans once ...

    def label_and_save(pquote,mutrec_index,seq_start,pnext,tlen,is_R1,is_tworeads,is_one_mut):
        nonlocal Generated_Fragcount,Saved_Fragcount,saved_unpaired
        Generated_Fragcount+=1
        ref_offset,fwd_cigar_label,rev_cigar_label=RG_process.get_trimmed_cigars(mutrecs[mutrec_index].cigarbox,mutrecs[mutrec_index].mutbox,seq_start,fraglen)
        if not RG_globals.is_muts_only or RG_process.is_mut_cigar(fwd_cigar_label):
            if not is_tworeads: # We have an unpaired read
                saved_unpaired+=1
        #if not RG_globals.is_muts_only or RG_process.is_mut_cigar(fwd_cigar_label):
            mutfragcount[mutrec_index]+=1
            if last_Saved_Fragcount==Saved_Fragcount: # Not incremented from previous pair or singleton
                Saved_Fragcount+=1             
            inseq=mutrecs[mutrec_index].seq[seq_start:seq_start+fraglen]
            label_and_saveg(pquote,inseq,Saved_Fragcount,seq_start,ref_offset,fwd_cigar_label,rev_cigar_label,mutrecs[mutrec_index],pnext,tlen,is_R1,is_tworeads,True)
        return

    def do_r1_save(R1pnext):# Save the forward sequencing fragment
        if R1_frag:
            label_and_save("pf",mutrec_index,R1_start,R1pnext,insert_len,True,is_tworeads,is_one_mut)
 
    def do_r2_save(R2pnext):# Save the reverse sequencing fragment
        if R2_frag:
            label_and_save("pr",mutrec_index,R2_start,R2pnext,-insert_len,False,is_tworeads,is_one_mut)

    def do_two_reads():
            nonlocal unpaired_count,is_tworeads
            ###########  This section to be skipped if excluding unpaired reads ###########

            # We need to get the reference-based start position, as a pnext for the other read in the pair, which we get by traversing the CIGAR                     
            if R1_frag: # Get pnext and is_cigar_mut for forward sequencing fragment
                R2pnext,R1fwdcigar,revcigar=RG_process.get_trimmed_cigars(mutrecs[mutrec_index].cigarbox,mutrecs[mutrec_index].mutbox,R1_start,fraglen)# pnext and cigar
            else:
                R2pnext=0
                R1fwdcigar=""

            if R2_frag: # Get pnext and is_cigar_mut for reverse sequencing fragment
                R1pnext,R2fwdcigar,revcigar=RG_process.get_trimmed_cigars(mutrecs[mutrec_index].cigarbox,mutrecs[mutrec_index].mutbox,R2_start,fraglen)# pnext and cigar
            else:
                R1pnext=0
                R2fwdcigar=""

            # With the CIGARs, if we are filtering for  variants-only, we need to determine if just one of the reads in the pair contains a variant. 
            # If so, we need to retain both reads, otherwise the paired-end bit is pointless.
            if RG_globals.is_muts_only:
                R1_is_mut=RG_process.is_mut_cigar(R1fwdcigar)
                R2_is_mut=RG_process.is_mut_cigar(R2fwdcigar)
                is_one_mut=R1_is_mut or R2_is_mut
                '''
                if not is_tworeads: # We have an unpaired read
                    # Check values for the unpaired read: R1 expected to be at high end of reefrence, R2 at near end
                    # Use BRCA1_hap2 at 850 insert length, 200 read length, variant-only, >20 DOC, to find cases of unpaired reads that include variants
                    if R1_frag:
                        print("R1 is unpaired with CIGAR %s, starting %s"%(R1fwdcigar,R2pnext))
                    else:
                        print("R2 is unpaired with CIGAR %s, starting %s"%(R2fwdcigar,R1pnext))
                '''
            do_r1_save(R1pnext); do_r2_save(R2pnext) # Always save R1 first, R2 second

            # OR: save forward and reverse in random order, not always R1-first!
            # Also have a suspicion that we need to do one more 50% random: reverse complement both sequences / keep
            # in order to emulate insert-orientation. The question is ... does it make any difference?
            #if random() < 0.5:
            #    # Do a forward, then a reverse, 50% of the time
            #    do_r1_save(); do_r2_save()
            #else:
            #    # Do a reverse, then a forward, 50% of the time
            #    do_r2_save(); do_r1_save()

            ###########  END OF: This section to be skipped if excluding unpaired reads ###########

    ###########  END OF def do_two_reads():

    # Calculate an expected fragtotal using un-modified spliced-reference
    fragtotal= int(len(mutrecs)*fragdepth*(RefRec.spliced_length/fraglen))

    # Get two parallel lists to match mutrecs
    # First a list of normalised frequencies for each mutseq
    normutfreq=normalise_mutfreqs(mutrecs)
           
    if len(normutfreq)>0:
        # Sum the normalised frequencies progressively to create number-ranges in lottoline
        # lottoline converts the normalised frequencies into a list of cutoff positions for matching a later random number call
        lottoline=[normutfreq[0]]
        for i in range(1, len(normutfreq)):
            lottoline.append(lottoline[len(lottoline)-1]+normutfreq[i])
            
        #  Now do the fragmentation
        while Generated_Fragcount < fragtotal:
            last_Saved_Fragcount=Saved_Fragcount
            
            # Pick a random mutrecs entry to be VarSeq, weighted by its normalised values.
            lottowin=randint(1,NormaliseUpper)
            for i in range(0, len(lottoline)):
                if lottowin <= lottoline[i]:
                    mutrec_index=i
                    break               
                
            while True:
                insert_len=int(np.random.normal(loc = RG_globals.gauss_mean, scale=RG_globals.gauss_SD, size=None)) # Gaussian distribution frequencies
                if (Paired_insert_Min <= insert_len <=Paired_insert_Max) and (fraglen <= insert_len):
                    break
            #insert_len=randint(Paired_insert_Min,Paired_insert_Max) # Flat distribution frequencies

            # Constrain the range position to a defined locus range for genomic: headpop & mutseqendpop[mutrec_index] or select a position in exonplus_lookup
            if is_do_exome: # Find a random position in exonplus_lookup. This is zero-based!
                start=choice(mutrecs[mutrec_index].exonplus_lookup)
            else:
                start=randint(headpop,mutseqendpop[mutrec_index]) # For genomic; forced when *not* RG_globals.is_exome_paired_end, even if mRNA or CDS selected
            
            # 50% chance of the position start being the forward read's (R1) start-point or the reverse-read's (R2) end-point
            if random() < 0.5: # Position is the forward-read's start point
                R1_start=start;
                R2_end=R1_start+insert_len-1
            else: # Position is reverse-read's end-point
                R2_end=start;
                R1_start=R2_end-insert_len+1 # Will be < headpop if outside the locus
            #both    
            R1_end=R1_start+fraglen_adjust
            R2_start=R2_end-fraglen_adjust

             # Allow one end to go into pop areas, not both
            R1_frag=(0 < R1_start <= mutseqendpop[mutrec_index]) and (headpop < R1_end <= mutseqend[mutrec_index]) 
            R2_frag=(0 < R2_start <= mutseqendpop[mutrec_index]) and (headpop < R2_end <= mutseqend[mutrec_index])

            if RG_globals.is_pair_monitor:
                if R1_frag or R2_frag:
                    # pair_monitor is done in 1-based
                    pair_monitor_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(Saved_Fragcount+1,start+1,insert_len,R1_start+1,R2_start+1,R1_end+1,R2_end+1,R1_frag,R2_frag))
            # Detect whether two reads are being produced
            is_tworeads=R1_frag and R2_frag
            if is_tworeads: # C Roos version to ignore unpaired reads
                do_two_reads()
            else:
                unpaired_count+=1
    journal_frags(mutrecs,Generated_Fragcount,Saved_Fragcount,mutfragcount,unpaired_count,saved_unpaired)
    #print("noread_count: %s ; unpaired_count %s "%(noread_count,unpaired_count))
    return
# end of def generate_multisource_paired_frags4(RefRec,mutrecs,fraglen,fragdepth)


#==========================End of section dealing with generate_multisource_random_frags and generate_multisource_all_frags===========================

def journal_frags(mutrecs,Generated_Fragcount,Saved_Fragcount,mutfragcount,unpaired,saved_unpaired):
    # This to catch a division-by-zero error. Needs a tidy-up
    min_norm=0
    if len(mutfragcount)>0:
        #min_norm=min(mutfragcount)
        try: # If all are 0, this will fail # Try is not good for pyodide
            min_norm=min([i for i in mutfragcount if i > 0])
        except: # Failed, so make sure it goes nowhere
            min_norm=0
    if min_norm > 0:
        journal_frags2(mutrecs,Generated_Fragcount,Saved_Fragcount,mutfragcount,unpaired,saved_unpaired)
    else:
        journals_update(" No %ss meet the filter criteria eg: option '%s'"%(RG_globals.read_annotation,RG_globals.bio_parameters["is_muts_only"]["label"]))
        
def journal_frags2(mutrecs,Generated_Fragcount,Saved_Fragcount,mutfragcount,unpaired,saved_unpaired):
    global REFSEQ_RECORD
    global DOC_precision
    # DOC: Depth of sequencing should be = (total number of reads * average read length) / total length of reference sequence (which is the genomic region covered, or the exome).
    # The target number of fragments is calculated in generate_multisource_random_frags
    # Generated_Fragcount is the actual number of fragments generated before any filtering when RG_globals.is_muts_only is True 
    # Saved_Fragcount is the number of fragments saved to file; lower then Generated_Fragcount when RG_globals.is_muts_only is True

    fragtallystring=""; fragdocstring="";fragsourcelen=""; fragratio=""
    fragtallycount=0
    doc_mult=10**DOC_precision
    half_adj=0.5/doc_mult
    ref_doc_len=REFSEQ_RECORD.spliced_length
    if len(mutfragcount)>0:
        #min_norm=min(mutfragcount)
        min_norm=min([i for i in mutfragcount if i > 0])
    for item in mutrecs:
        if item.mutfreq >0:
            srclab=RG_globals.target_locus+"_"+item.mutlabel+"("
            fragtallystring=fragtallystring+srclab+str(mutfragcount[fragtallycount])+"), "
            doc=mutfragcount[fragtallycount]*RG_globals.Fraglen/ref_doc_len
            fragdocstring=fragdocstring+srclab+str(int((doc+half_adj)*doc_mult)/doc_mult)+"), "
            fragsourcelen=fragsourcelen+srclab+str(len(item.seq))+"), "
            ratiofreq=mutfragcount[fragtallycount]/min_norm    
            fragratio=fragratio+srclab+str(int((ratiofreq+half_adj)*doc_mult)/doc_mult)+"), "
            fragtallycount+=1
    
    str0="For a '%s' target value of %s, based on template length %s:"%(bio_parameters["Fragdepth"]["label"],RG_globals.Fragdepth,ref_doc_len)
    str1=str(Generated_Fragcount); str2="s";str3=str(RG_globals.Fraglen) ; str4="random";str4a="within" ; str5=str(fragtallycount); str6=RG_globals.variants_label;
    str8="sequence"
    if fragtallycount <2 :
        str7=""
    else:
        str7="s"
    if RG_globals.is_onefrag_out:
        str0=""; str1+=" single"; str2=""; str4="each possible"; str4a="for"

    if RG_globals.is_frg_paired_end:
        str4="paired"
        #str8="%s %ss or un%s"%(str4,RG_globals.read_annotation,str4)
        str8=str4

    journals_update("\n%s\n Generated %s %ss of length %s bases at %s starting position%s %s %s %s%s: %s"
                   %(str0,str1,RG_globals.read_annotation,str3,str4,str2,str4a,str5,str6,str7,mutlistlabels))

    if RG_globals.is_duplex and RG_globals.is_onefrag_out:
        journals_update("\tSelecting option '%s' with '%s' doubles the %s total, by creating a reverse-complement for each %s"
                       %(bio_parameters["is_duplex"]["label"],bio_parameters["is_onefrag_out"]["label"],RG_globals.read_annotation,RG_globals.read_annotation))
    
    if RG_globals.is_muts_only:
        str9="excluded"
        str10=", leaving ..."
    else:
        str9="included"
        str10=""
    journals_update(" Reference-only %s %ss are %s because option '%s' is set to %s%s"%(str8,RG_globals.read_annotation,str9,
                                                                                      bio_parameters["is_muts_only"]["label"],RG_globals.is_muts_only,str10))
    
    journals_update(" source(count):\n\t%s"%(fragtallystring[:-2]))
    journals_update(" source(count-ratio):\n\t%s"%(fragratio[:-2]))
    journals_update(" source(length):\n\t%s"%(fragsourcelen[:-2]))
    journals_update(" source('%s'=count*%s/template_length):\n\t%s"%(bio_parameters["Fragdepth"]["label"],RG_globals.Fraglen,fragdocstring[:-2]))
    REFSEQ_RECORD.Generated_Fragcount=Generated_Fragcount
    REFSEQ_RECORD.Saved_Fragcount=Saved_Fragcount
    REFSEQ_RECORD.Saved_Unpaired_Fragcount=saved_unpaired
    REFSEQ_RECORD.Unsaved_Unpaired_Fragcount=unpaired
    REFSEQ_RECORD.Saved_Paired_Fragcount=Saved_Fragcount-saved_unpaired
# end of journal_frags2()

def split_with_char(instring,splitchar):
    # This routine only created to compensate for Pyodide/Python inconsistency
    # Attempting to split when the character is absent hangs in pyodide/webassembly conversion, Python just carries on
    # Created primarily for transcript ids in format ENST00000311936.8
    if splitchar in instring:
        splitters=instring.split(splitchar)
        outstring=splitters[0]
    else:
        outstring=instring
    return outstring

def splice_a_sequence2(seq_record,target_loc,splice_trigger):
    # This splices any sequence you send it, if join definitions are within the "splice_trigger" feature
    ''' Where "splice_trigger" matches the specific feature field (eg: mRNA, CDS) this creates a feature field of "skip_trigger" (eg:"gap")  '''
    ''' "skip_trigger" is used in MutateVarSeq to introduce a deletion in addition to the usual deletion triggers. '''
    ''' Prior to the actual mutation the object get_mutref_olap is used to identify clashes between these skips and variant definition.'''
    ''' If a splice-extension value is set via RG_globals.Exome_extend, any splice is extended in each direction eg: -3/+3 added to the splice position  '''
    ''' Any feature changes are stored in extra_record_features, not seq_record.features because of the Python passing rules which will change the original'''
    ''' Features are to be merged into REFSEQ on return and then into mut sequences feature tables when calling MutateVarSeq'''
    ''' The spliced sequence is created here, and also passed back as a separate value, not by a change to seq_record.seq  '''
    ''' Development would be to use return_message as in merge_seqvar_records and move this to RG_process'''
    
    ''' from Bio.SeqFeature import SeqFeature, FeatureLocation required'''
    # Set null return parameters to avoid non-assignment crashes
    is_spliced =False
    splicetotal=0; splice_joinlist = [] ; extra_record_features=[];splice_joinlist_txt = []
    #mirror=[] #Hiding mirror: A plain-text representation of SeqFeature that's easier to follow for diagnostics
    spliceseq=""
    target_gene="target_gene"
    message="feature_"
    this_transcript_id="this_transcript_id"
    this_ref_offset=0;this_ref_endclip=0
        
    def skipit(minloc,maxloc,begin,end,second):
        nonlocal splicetotal, splice_joinlist,splice_joinlist_txt,extra_record_features,this_ref_offset,this_ref_endclip
        # Hiding building of the mirror components
        #nonlocal mirror
        #mirrstring="["+str(begin)+":"+str(end)+"](+)"
        #print("skipit: begin %i end %i, minloc %s, maxloc %s"%(begin,end,minloc,maxloc))
        splicetotal+=1
        f=SeqFeature(FeatureLocation(begin,end,strand=1),type="variation")
        if begin ==0:
            dbxref="%s:5-prime upstream trim"%RG_globals.skip_trigger
        elif end == len(seq_record.seq):
            dbxref="%s:3-prime downstream trim"%RG_globals.skip_trigger
        else:
            dbxref="%s:Intron %s-%s"%(RG_globals.skip_trigger,splicetotal-1,splicetotal)
        #dbxref="%s:%snum%s"%(RG_globals.skip_trigger,RG_globals.skip_trigger,str(splicetotal))
        f.qualifiers['replace']=["N/-"]
        f.qualifiers['db_xref']=[dbxref]
        #if RG_globals.is_print_diagnostics: print(begin,first,second)
        if second >= begin:
            # Want highest-position first
            extra_record_features.insert(0,f)
            #mirror.insert(0,mirrstring)
        else:
            extra_record_features.append(f)
            #mirror.append(mirrstring)
        if maxloc != begin:
            splice_joinlist.append(FeatureLocation(minloc,maxloc,strand=+1))
            splice_joinlist_txt.append("%s:%s"%(minloc+1,maxloc))

        #print("mirror %s"%mirror)
        #print("splice_joinlist %s\n"%splice_joinlist)
    # ===== end of def skipit() which is local to splice_a_sequence(seq_record,target_loc,splice_trigger)

    def call_skipit(feature,gap_extend):
                        begin=0
                        final=len(seq_record.seq)
                        parts=feature.location.parts # feature.location is expected to be a BioPython CompoundLocation.
                        #print("Parts: %s"%parts)
                        #print("Parts[0]: %s, Parts[0].strand %s"%(parts[0],parts[0].strand))
        
                        # Potential clash/inconsistency when using this to derive splice_joinlist in get_transcript_data_process
                        if parts[0].strand == -1: # This means that the original Genbank join list is in the format join(complement(a..b,c..d)). rather than join(p..q,r..s)
                                                  # and expect to find the highest numbered locations positioned first.
                                                  # We need the lowest first in this operation, so must reverse the order in which processed
                            rev_parts=[]
                            for item in parts:
                                rev_parts.insert(0,item)
                            parts=rev_parts
                        for item in parts:
                            #print("Part item %s"%item)
                            first=item.start
                            second=item.end
                            
                            # Adjust IE boundaries and correct if exceed limits
                            if gap_extend > 0: # passed-in value gap_extend replaces use of RG_globals.Exome_extend
                                first = first - gap_extend
                                second = second + gap_extend
                                if first < 0:
                                    first =0
                                if second > final:
                                    second = final
                                ''' Check adjusted boundary does not impinge '''
                                if first < begin: first=begin
                                
                            minloc=first
                            maxloc=second
                            end=first
                            if end > begin :
                                skipit(minloc,maxloc,begin,end,second)
                            begin=second
                            
                        if second < final:
                            # There is an intron at the end of the sequence
                            end=final
                        if (final > begin):
                            skipit(minloc,maxloc,begin,end,second)
                        return(item.start,item.end,first,second)
                            
   # ===== end of def call_skipit() which is local to splice_a_sequence2(seq_record,target_loc,splice_trigger)

    '''
    This section locates "splice" feature fields using the "splice-trigger", typically "mRNA" or "CDS".
     mRNA            join(6..186,5542..5663,23525..23703,25164..25323,
                     41026..46148)
                     /gene="ENSG00000133703.11"
                     /standard_name="ENST00000311936.7"
    To be consistent with other definitions, expect to have a minimum of two locations in the definition.
    Using 0..0 is a way to cheat when requiring a single-region slice eg: join(0..0,6..186)
    '''
    
    #settxt=bio_parameters["target_transcript_name"]["label"] # "Template" default
    settxt=RG_globals.reference_haplotype
    if splice_trigger=="gene":
        settxt=RG_globals.reference_gene
        
    journals_update("\n Searching feature table to set '%s' Sequence"%settxt)
    # First match the locus id (format eg: 'KRAS') to determine what the local name is, format eg:ENSG00000133703.11
    
    for (index, features) in enumerate (seq_record.features):
        this_feature=seq_record.features[index]
        if this_feature.type == 'gene': #
            # if str(this_feature.location.strand) == RG_globals.seq_polarity_minus # Come back to this - replace need for is_join_complement in config.json 20/09/2022
            #print("gene feature.location.strand %s"%this_feature.location.strand)
            if this_feature.qualifiers.get("locus_tag")[0] == split_with_char(target_loc,"_"): # eg: trims KRAS_plus to KRAS so do not have to modify source data
                target_gene=str(this_feature.qualifiers.get("gene")[0])
                #print("ensembl_geneid for %s is %s"%(target_loc,RG_globals.ensembl_geneid))
                last_start,last_end,last_first,last_second=call_skipit(this_feature,RG_globals.Exome_extend)
                # call_skipit seems to return 0-based position for last_start. Correction needed for last_start, but I am wary of implication
                journals_update("  Locus %s matches gene id %s from %s Range: %s - %s ; global : %s - %s; length: %s bases"
                               %(target_loc,target_gene,in_ref_src,last_start+1,last_end,RG_process.get_absolute_position(seq_record,last_start+1),
                                 RG_process.get_absolute_position(seq_record,last_end),last_end-last_start))
                this_ref_offset=last_start
                this_ref_endclip=len(seq_record.seq)-last_end
                
                break
    # end of : for (index, features) in enumerate (seq_record.features)

    if RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
        this_transcript_id=target_gene
        if RG_globals.is_trim_to_gene:
            message="trimmed"
            is_spliced=True
        else:
            splicetotal=0; this_ref_offset=0; this_ref_endclip=0 # Re-setting because these become objects
            extra_record_features=[] # Re-setting is necessary to cancel the trim
    else:
        # When left with RG_globals.target_transcript_name != RG_globals.empty_transcript_name, find a splice definition for the same target gene 
        journals_update("  Looking for %s feature in %s to match %s:%s "%(splice_trigger,target_gene,RG_globals.target_transcript_name,RG_globals.target_transcript_id))
        splicetotal=0; splice_joinlist = [] ; splice_joinlist_txt = []; extra_record_features=[] # Re-setting because call_skipit was used for gene feature, above
        #mirror=[] #Hiding mirror. Mirror needs clearing, just like extra_record_features
        for (index, features) in enumerate (seq_record.features):
            this_feature=seq_record.features[index]
            if this_feature.type == splice_trigger:
                if target_gene == str(this_feature.qualifiers.get("gene")[0]):
                    # Found it, so process
                    try: # This to catch Python error for missing objects when feature definitions are missing # Try is not good for pyodide
                        if RG_globals.is_CDS:
                            feature_qualifier="note"
                            message=feature_qualifier
                            get_note=this_feature.qualifiers.get(feature_qualifier)[0]

                            note_tid_qualifier="transcript_id="
                            this_transcript_id=get_note.replace(note_tid_qualifier,"")
                                                        
                            message= "'"+feature_qualifier + "' retrieving "+ this_transcript_id+ " via '" +note_tid_qualifier+"'"

                            feature_qualifier="protein_id"
                            this_protein_id=this_feature.qualifiers.get(feature_qualifier)[0]
                            # If it survives ... fiddle for journal report
                            message= message +" with " + feature_qualifier+" " + this_protein_id
                            #feature_qualifier=" to match note and " + this_protein_id + " with " + this_transcript_id+ " via " +note_tid_qualifier
                        else:
                            feature_qualifier="standard_name"
                            message=feature_qualifier
                            this_transcript_id=this_feature.qualifiers.get(feature_qualifier)[0]
                    except:
                        journals_update(" *** ERROR missing %s qualifier in feature definition ***"%feature_qualifier)
                        
                    # Check if the found transcript id matches the targeted one, ignoring version 
                    if split_with_char(this_transcript_id,'.') == split_with_char(RG_globals.target_transcript_id,'.'):
                        is_spliced=True
                        last_start,last_end,last_first,last_second=call_skipit(this_feature,RG_globals.Exome_extend)
                    if is_spliced:
                        # only do it once: first successful splice-join feature that matches transcript id for this gene, ignore any others, so break
                        break # Exit the loop:  for (index, this_feature) in enumerate (seq_record.features):
                    # end of: if split_with_char(this_transcript_id,'.') == ...
                # end of : if target_gene == str(this_feature.qualifiers.get("gene")[0]):
            # end of: if this_feature.type == splice_trigger
           # end of ...# if RG_globals.target_transcript_name != RG_globals.empty_transcript_name
        # end of: for (index, features) in enumerate (seq_record.features):
    #if RG_globals.is_print_diagnostics:print("mirror:",mirror)

    #print("is_spliced%s, mirror %s"%(is_spliced,mirror))
    if is_spliced:
        #print("splice_joinlist %s \n"%splice_joinlist)
        joinlist=splice_joinlist[0]
        if len(splice_joinlist) > 1:
            for pos in range(1,len(splice_joinlist)):
                joinlist+=splice_joinlist[pos]
                #print("\njoinlist %s\n"%joinlist)
        #print("joinlist %s"%joinlist)
        #spliceseq=MutableSeq(str(joinlist.extract(seq_record.seq)),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
        spliceseq=Biopython_fix.fix_MutableSeq(str(joinlist.extract(seq_record.seq)))
        #print("splice_a_sequence2: len spliceseq %s"%len(spliceseq))
        # Modify the start and end trim-lengths to the addition of Exome_extend
        this_ref_offset-=RG_globals.Exome_extend
        this_ref_endclip-=RG_globals.Exome_extend
    return (is_spliced,splicetotal,spliceseq,extra_record_features,target_gene,message,this_transcript_id,last_first,last_second,this_ref_offset,this_ref_endclip,splice_joinlist_txt)
# end of def splice_a_sequence2(seq_record,target_loc)

def splice_a_sequence(seq_record,target_loc):
    # splice_a_sequence fronts splice_a_sequence2
    global in_ref_src

    def sac_update(template_length,clip_length):
        nonlocal last_first,last_second,splicetotal

        # This documents the addition of nucleotides to each end of the locus range when call_skipit(this_feature,RG_globals.Exome_extend) for this_feature.type == 'gene'
        if RG_globals.Exome_extend==0:
            msgtxt="also"
        else:
            msgtxt="is %s bases longer"%(2*RG_globals.Exome_extend)
            #msgtxt="is now"

        #template_length=last_second-last_first+1
        journals_update("  With option '%s'=%s, %s Range from %s %s: %s - %s ; global : %s - %s; length: %s bases"%(
                    bio_parameters["Exome_extend"]["label"],
                    RG_globals.Exome_extend,
                    bio_parameters["target_transcript_name"]["label"],
                    in_ref_src,
                    msgtxt,
                    last_first+1,last_second,
                    RG_process.get_absolute_position(seq_record,last_first+1),
                    RG_process.get_absolute_position(seq_record,last_second),
                    template_length))
        journals_update("  %s %s has %s bases; %s spliced-out regions compared to %s %s with %s bases"%(bio_parameters["target_transcript_name"]["label"],
                                                                                                       locus_transcript,template_length,splicetotal,in_ref_src_title,
                                                                                                       in_ref_src,clip_length))

    splicetotal=0; is_spliced =False; spliceseq=""; extra_record_features=[]; this_ref_offset=0; this_ref_endclip=0 # Initialise returned parameters although moot with current settings

    if RG_globals.is_make_exome:
        if  RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
            is_spliced,splicetotal,spliceseq,extra_record_features,target_gene,message,gene_id,last_first,last_second,this_ref_offset,this_ref_endclip,splice_joinlist_txt\
                        =splice_a_sequence2(seq_record,target_loc,"gene")
            '''
            if not RG_globals.is_trim_to_gene:
                msgtxt=" NOT"
            else:
                msgtxt=""
            journals_update(" Reference Source sequence is%s trimmed to locus definition for a Reference sequence, because system-config option 'is_trim_to_gene' is set to %s"
                               %(msgtxt,RG_globals.is_trim_to_gene))
            '''
            if RG_globals.is_trim_to_gene:
                sac_update(len(spliceseq),seq_record.unclipped_length)
                
            if RG_globals.is_CDS:
                msgtxt=", even with '%s' selected, "%(bio_parameters["is_CDS"]["label"])
            else:
                msgtxt=" "
            journals_update(" No splicing of exon boundaries%sbecause '%s' is set to '%s'"
                           %(msgtxt,RG_globals.reference_haplotype,RG_globals.target_transcript_name))
        else:
            if RG_globals.is_CDS:
                join_trigger=RG_globals.CDS_trigger # changed from splice_trigger
            else:
                join_trigger=RG_globals.mRNA_trigger # changed from splice_trigger

            is_spliced,splicetotal,spliceseq,extra_record_features,target_gene,message,transcript_id,last_first,last_second,this_ref_offset,this_ref_endclip,splice_joinlist_txt\
                        =splice_a_sequence2(seq_record,target_loc,join_trigger)
 
            if not is_spliced :
                success_text="but"
                within_text=" *not* "
            else:
                success_text="and"
                within_text=" "

            journals_update(" Searched gene %s features for %s, %s found %s which appears%sto match %s:%s"
                           %(target_gene,message,success_text,transcript_id,within_text,RG_globals.target_transcript_name,RG_globals.target_transcript_id))
            sac_update(len(spliceseq),seq_record.clipped_length)
 
    else:
        journals_update(" No trimming nor splicing of exon boundaries because system-config option 'is_make_exome' is set to %s"%(RG_globals.is_make_exome))

    return (is_spliced,splicetotal,spliceseq,extra_record_features,this_ref_offset,this_ref_endclip)
# end of def splice_a_sequence(seq_record,target_loc)

# ===============================================================
# End of main mutation functions definition
# ===============================================================
def make_reference_files1():
    if RG_globals.is_write_ref_fasta:
        if not RG_globals.is_frg_paired_end:
            write_refseq(REFSEQ_RECORD,"force_Out_Ref_Source") # Save the trimmed Source sequence as reference
    else:
        if RG_globals.is_frg_paired_end:
            templatetxt=""
            addtxt=""
        else:
            templatetxt=", and '%s' %s.fasta, "%(RG_globals.reference_haplotype,Out_Ref_Source)
            addtxt="s"
        journals_update(" NOTE: Sequence file%s for: %s %s.fasta %sNOT saved because option '%s' is set to %s"%(addtxt,
                                                                                                               in_ref_src_title,
                                                                                                               Out_Ref_Source,
                                                                                                               templatetxt,
                                                                                                               bio_parameters["is_write_ref_fasta"]["label"],
                                                                                                               RG_globals.is_write_ref_fasta))
# run_processing was formerly the way of starting the batch run.
# It has been adapted into a function to support repeated-calling from a GUI module
def run_processing(is_get_new_refseq):
    global REFSEQ_RECORD,Seq_Format,in_ref_src_title,in_ref_src,Out_Ref_Source,user_config_file

    # Ensure that exonplus_lookup is calculated, which could differ from exonplus_lookup; intially because return from App.vue does not send exonplus_lookup
    #if RG_globals.exonplus_lookup==[0] and (RG_globals.target_transcript_name != RG_globals.empty_transcript_name):
    # Forcing this to recalculate each time. Really need to detect if the value of RG_globals.Exome_extend has changed, or new refseq
    if is_get_new_refseq:
        RG_globals.exonplus_lookup=[0] # Resetting value
    
    if RG_globals.target_transcript_name != RG_globals.empty_transcript_name:
        success,RG_globals.exonplus_lookup=RG_builder.get_muttranscripts2(True) # True is to *extend exome* for exonplus_lookup

    # Save the latest configs first
    user_config_file=RG_globals.save_user_configs()
    is_success=True
    if is_get_new_refseq:
        REFSEQ_RECORD,exists=read_refseqrecord(Seq_Format)
    else:
        exists=True
        REFSEQ_RECORD=set_maxrefs(REFSEQ_RECORD)
        
    if exists:
        make_reference_files1()
        if RG_globals.is_frg_paired_end:
            no_splice_refseq(REFSEQ_RECORD)
        else:
            REFSEQ_RECORD=splice_refseq(REFSEQ_RECORD)
        (Mutrecs,is_success)=get_mutrecords(REFSEQ_RECORD,Seq_Format)
        
        if is_success:
            is_fraglength_OK,shortest,short_label=validate_fraglength(Mutrecs,RG_globals.Fraglen)
            if not is_fraglength_OK:
                journals_update(" *** Processing halted *** Requested %s length %s is more than half the shortest '%s' %s with sequence length %s"%
                               (RG_globals.read_annotation,RG_globals.Fraglen,RG_globals.variants_label,short_label,shortest))
                is_success= False
            elif (RG_globals.is_fasta_out or RG_globals.is_fastq_out or RG_globals.is_sam_out) :
                is_success=refseq_to_frags4(Mutrecs)
                close_seqout(REFSEQ_RECORD)
            else:
                close_seqout(REFSEQ_RECORD)
            close_nonseq_files(is_success)
        run_success=is_success
    else:
        run_success=False    
    return run_success

#########################################
# MAIN
#########################################

# Calls the main exploder function
# Has the potential to behave differently by changing the two booleans.
# There's nothing that currently uses this option. prior versions of GUI did: now all transferred to RG_exploder_builder.py
def exploder_initiate(is_get_new_refseq,is_journal):
    #  This basically runs the job
    #print("arrival: RG_globals.target_transcript_name %s"%RG_globals.target_transcript_name)
    general_initiate(is_journal)
    run_success=run_processing(is_get_new_refseq)
    #print("leaving: RG_globals.target_transcript_name %s"%RG_globals.target_transcript_name)
    return run_success

#########################################
# When called from the python GUI module 
def call_exploder_main():
    run_success=exploder_initiate(True,True)
    return run_success
#########################################
# When run as a batch, as in web version ...
# Current App.js call
if __name__ == "__main__":
    if call_exploder_main():
        print("Done")
    else:
        print("Fail")
    
