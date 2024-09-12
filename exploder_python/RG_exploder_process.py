#!/usr/local/bin/python3
#Progver="RG_exploder_process2"
#ProgverDate="2-Sep-2024"
'''
© author: Cary O'Donnell for Replicon Genetics 2018, 2019, 2020, 2021, 2022, 2023, 2024
'''

# =================================================
# constant unique to this module
# =================================================

global olap_and_base3
# The keys are base3 as strings used to assign a 'type' value for the relationship between the variant feature & reference feature ranges.
# The relationship is the overlap, or lack thereof, between those ranges.
# The type, eg: type 17 for base3 string "20210", is documented outside of this code at
# PDF: https://github.com/Replicon-genetics/public/blob/main/Overlap_types_decription.pdf
# ODP source: https://github.com/Replicon-genetics/public/
olap_and_base3={"00001":1,"00000":2,"00002":3,"00100":4,"10002":5,"00200":6,"00210":7,
                "00220":8,"10101":9,"10200":10,"20102":11,"10210":12,"10220":13,
                "20201":14,"20200":15,"20202":16,"20210":17,"21202":18,"20220":19,
                "21211":20,"21220":21,"22212":22,"22221":23,"22220":24,"22222":25
                }

# =================================================
# python_imports
# =================================================
'''
Python imports. These define global constants, functions, data structures.
As such they must be defined first, being impractical to be called by subroutine and declaring all the globals
'''
from random import randint # COD: Used in get_quality_list # re-instated
from collections import OrderedDict
import copy
import time
'''
BioPython imports
'''
from Bio.SeqFeature import SeqFeature, FeatureLocation
import Biopython_fix  # New Feb 2023 to fix MutableSeq, Seq, SeqRecord, Bio.Alphabet deprecated from Biopython 1.78 (September 2020)

# =================================================
# End of python_imports
# =================================================
'''
RG module imports
'''
import RG_exploder_globals as RG_globals

# =====================================================
# Seqrecord manipulation
# =====================================================

def make_exonplus_lookups(ref_lookup,mutrecs):
    # NB: Incoming ref_lookup and cigarbox are both 1-based sequence counts
    
    def print_counts():
        print("seqlen:%s,seqtype:%s,new_seqcount:%s,old_scount:%s,carry:%s"%(seqlen,seqtype,new_seqcount,old_scount,carry))

    def get_keys(value,first,last):
        #nonlocal elapsed_time
        #print("get_keys: value %s,first %s,last %s"%(value,first,last))
        nonlocal this_lookup
        start_time=time.time()
        begin=0; end=0
        if last > last_lookup:
            last=last_lookup
            #print("last mod: value %s,first %s,last %s"%(value,first,last))
        success1=False; success2=False
        begin_idx=0
        end_idx=0
        for begin in range(first,last+1):
            if begin in this_lookup:
                success1=True
                break
        #print("first %s,last %s, begin %s"%(first,last,begin))    
            
        for end in range(last,begin,-1):
            if end in this_lookup:
                success1=True
                break
        #print(" begin %s,end %s"%(begin,end))
        success=success1 and success2
        if success:
            begin_idx=this_lookup.index(begin)
            end_idx=this_lookup.index(end)
            #print("begin: %s; end: %s"%(this_lookup[begin_idx],this_lookup[end_idx]))
        finish_time=time.time()
        #print("time: %s"%(finish_time-start_time))
        #print("begin:%s; end:%s"%(begin,end))
        return begin_idx,end_idx,success              


    def modify_lists(value):
        #value is "H" for hide; or carry, an integer
        first=last_old_scount+1;last=old_scount
        #print("first %s, last %s, value %s"%(first,last,value))# Hiding unnecessary diagnostics code
        if str(value) in "DX":
            ''' # Hiding unnecessary diagnostics code
            hidden_key="%s,%s"%(first,last)
            hidden_list.setdefault(hidden_key,value)
            '''

            begin_idx,end_idx,success=get_keys(value,first,last)
            if success:
                listpos_key="%s,%s"%(begin_idx,end_idx)
                hidden_listpos.setdefault(listpos_key,value)
                '''
                # Hiding unnecessary diagnostics code
                valpos_key="%s,%s"%(this_lookup[begin_idx],this_lookup[end_idx])
                hidden_valpos.setdefault(valpos_key,value)
                '''
        else:
            ''' # Hiding unnecessary diagnostics code
            carry_key="%s,%s"%(first,last)
            carry_list.setdefault(carry_key,value)
            '''
            
            if value !=0:
                begin_idx,end_idx,success=get_keys(value,first,last)
                if success:
                    listpos_key="%s,%s"%(begin_idx,end_idx)
                    carry_listpos.setdefault(listpos_key,value)
                    '''
                    # Hiding unnecessary diagnostics code
                    valpos_key="%s,%s"%(this_lookup[begin_idx],this_lookup[end_idx])
                    carry_valpos.setdefault(valpos_key,value)
                    '''
                    
    def modify_lookup():
        modify_lookup_carry()
        modify_lookup_hide()
        #print("this_lookup:%s"%this_lookup)
        #print("mutrecs[rec_index].id:%s"%mutrecs[rec_index].id)
        #print("length this_lookup: %s"%(len(this_lookup))) # # Hiding unnecessary diagnostics code
        #print("rec_index:%s,len(mutrecs[rec_index].exonplus_lookup: %s"%(rec_index,len(mutrecs[rec_index].exonplus_lookup))) # # Hiding unnecessary diagnostics code
    
    def modify_lookup_carry():
        nonlocal this_lookup
        for key in carry_listpos:
            #print(" key:%s; value:%s"%(key,carry_listpos[key])) # # Hiding unnecessary diagnostics code
            value=carry_listpos[key]
            start,end=key.split(',')
            for num in range(int(start),int(end)+1):
                this_lookup[num]+=value
                
    def modify_lookup_hide():
        nonlocal elapsed_time
        nonlocal this_lookup
        #print("length this_lookup: %s"%len(this_lookup))
        #begin_time=time.time()
        expanded_hidden_listpos=[]
        for key in hidden_listpos:
            #print(" key:%s; value:%s"%(key,hidden_listpos[key])) # # Hiding unnecessary diagnostics code
            start,end=key.split(',')
            expanded_hidden_listpos = list(range(int(start),int(end)+1))
            expanded_hidden_listpos.reverse()# Reverses.
        for num in expanded_hidden_listpos:# Do it in reverse, so array positions are maintained after deletion
            del(this_lookup[num]) # this_lookup[num]=0 # is faster with remove .reverse() above, but downstream cost?
        #end_time=time.time()
        #elapsed_time=elapsed_time+end_time-begin_time

    def modify_to_zero_based():# Add -1 to each position for zero-based sequence positions
        nonlocal this_lookup
        this_lookup = [item - 1 for item in this_lookup] # list comprehension
        
    #print(" CLEAR\n\n")
    #print(" Reflookup:%s\n" %ref_lookup)
    elapsed_time=0
    begin_time=time.time()
    if len(ref_lookup)>1:
        for mutrec_index in range(len(mutrecs)):
            #print("mutrecs[mutrec_index].id:%s"%mutrecs[mutrec_index].id)
            #print("mutrecs[mutrec_index].exonplus_lookup[2] %s"%mutrecs[mutrec_index].exonplus_lookup[2])
            #print("mutrecs[mutrec_index].id:%s; mutrecs[mutrec_index].mutlabel:%s;\n mutrecs[mutrec_index].cigar:%s;\n mutrecs[mutrec_index].mutbox:%s;\n
            #print("mutrecs[mutrec_index].cigarbox:%s"%mutrecs[mutrec_index].cigarbox)
            this_lookup=copy.copy(ref_lookup)
            first_lookup=this_lookup[0]
            last_lookup=this_lookup[-1]
            old_scount=0; new_seqcount=0; carry=0
            carry_listpos=dict()
            hidden_listpos=dict()
            
            ''' # Hiding unnecessary diagnostics code
            carry_list=dict()
            hidden_list=dict()
            hidden_valpos=dict()
            carry_valpos=dict()
            '''
            for recpos in range(1, len(mutrecs[mutrec_index].cigarbox)-1, 2):
                seqlen=mutrecs[mutrec_index].cigarbox[recpos]
                seqtype=mutrecs[mutrec_index].cigarbox[recpos+1]
                last_old_scount=old_scount
                if seqtype in "D":
                    for item in range(int(seqlen)):
                        carry-=1
                        old_scount+=1
                        #print_counts()
                    modify_lists("D")      
                elif seqtype in "I":
                    for item in range(int(seqlen)):
                        carry+=1
                        new_seqcount+=1
                        #print_counts()
                elif seqtype in "X":
                    for item in range(int(seqlen)):
                        new_seqcount+=1
                        old_scount+=1
                        #print_counts()
                    modify_lists("X")
                else:
                    new_seqcount+=seqlen
                    old_scount+=seqlen
                    modify_lists(carry)
                    #print_counts()
            '''# Hiding unnecessary diagnostics code
            print("hidden_list%s"%hidden_list)
            print("hidden_listpos%s"%hidden_listpos)
            print("hidden_valpos%s"%hidden_valpos)
            print("carry_list%s"%carry_list)
            print("carry_listpos%s"%carry_listpos)
            print("carry_valpos%s"%carry_valpos)
            '''
            modify_lookup()
            #print("this_lookup[0] %s, this_lookup[:-1] %s"%(this_lookup[0],this_lookup[-1]))
            modify_to_zero_based()
            #print("this_lookup[0] %s, this_lookup[:-1] %s"%(this_lookup[0],this_lookup[-1]))
            # copy the modified lookup into relevant mutrecs
            mutrecs[mutrec_index].exonplus_lookup=copy.copy(this_lookup)
            
            '''
            print("\n")
            '''
        end_time=time.time()
        elapsed_time=end_time-begin_time
        #print("elapsed_time:%s"%elapsed_time)

def duplicate_record(SeqRec):
    #CopySeqRec=modify_seq_in_record(SeqRec.seq,SeqRec,SeqRec.howmany)
    CopySeqRec=copy.deepcopy(SeqRec)
    return CopySeqRec

def modify_seq_in_record(inseq,SeqRec):
    # Adaptation of annotate_seq_to_record
    # Just copies across a modified sequence to a new record, including the features
    # annotate_seq_to_record ignores the features
    #print("at modify_seq_in_record: SeqRec.id %s"%(SeqRec.id))

    #CopySeqRec=SeqRecord(Seq(str(inseq),generic_dna)) # Creates an empty record with a given sequence #generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #CopySeqRec=SeqRecord(Seq(str(inseq))) # Creates an empty record with a given sequence
    CopySeqRec=Biopython_fix.fix_SeqRecord(inseq)
    
    # Copies each of the other records from SeqRec into this new record
    # Cannot do this another way because the sequence part of a sequence record is immutable.
    # NB: Have not recently tried CopySeqRec=copy.deepcopy(SeqRec) then CopySeqRec.Seq=Seq(str(inseq),generic_dna)
    CopySeqRec.id=copy.copy(SeqRec.id)
    CopySeqRec.firstid=copy.copy(SeqRec.firstid)
    CopySeqRec.locus_range=copy.copy(SeqRec.locus_range)
    CopySeqRec.name=copy.copy(SeqRec.name)
    CopySeqRec.mutlabel=copy.copy(SeqRec.mutlabel)
    CopySeqRec.description=copy.copy(SeqRec.description)
    CopySeqRec.polarity=copy.copy(SeqRec.polarity)
    CopySeqRec.offset=copy.copy(SeqRec.offset)

    CopySeqRec.GRChver=copy.copy(SeqRec.GRChver)
    CopySeqRec.chrom=copy.copy(SeqRec.chrom)
    CopySeqRec.abs_start=copy.copy(SeqRec.abs_start)
    CopySeqRec.abs_end=copy.copy(SeqRec.abs_end)
    CopySeqRec.strand_mod=copy.copy(SeqRec.strand_mod)
    CopySeqRec.howmany=copy.copy(SeqRec.howmany)
    CopySeqRec.annotations=copy.copy(SeqRec.annotations)
    
    CopySeqRec.features=copy.copy(SeqRec.features)
    CopySeqRec.cigar=copy.copy(SeqRec.cigar)
    CopySeqRec.cigarbox=copy.copy(SeqRec.cigarbox)
    CopySeqRec.mutbox=copy.copy(SeqRec.mutbox)
    #CopySeqRec.exonplus_lookup=copy.copy(SeqRec.exonplus_lookup) # Not needed once make_exonplus_lookups is fully implemented
    CopySeqRec.Headclip=copy.copy(SeqRec.Headclip)
    CopySeqRec.Tailclip=copy.copy(SeqRec.Tailclip)
    CopySeqRec.splicecount=copy.copy(SeqRec.splicecount)
    CopySeqRec.endclipcount=copy.copy(SeqRec.endclipcount)
    CopySeqRec.is_ref_spliced=copy.copy(SeqRec.is_ref_spliced)
    CopySeqRec.unclipped_length=copy.copy(SeqRec.unclipped_length)
    CopySeqRec.clipped_length=copy.copy(SeqRec.clipped_length)
    CopySeqRec.spliced_length=copy.copy(SeqRec.spliced_length)
    CopySeqRec.endlocus=copy.copy(SeqRec.endlocus)
    CopySeqRec.rev_endlocus=copy.copy(SeqRec.rev_endlocus)
    CopySeqRec.Generated_Fragcount=copy.copy(SeqRec.Generated_Fragcount)
    CopySeqRec.Saved_Fragcount=copy.copy(SeqRec.Saved_Fragcount)
    CopySeqRec.Saved_Unpaired_Fragcount=copy.copy(SeqRec.Saved_Unpaired_Fragcount)
    success,CopySeqRec=set_seqrec_absolutes(CopySeqRec) #Must do this even if RG_globals.is_use_absolute == False
    return CopySeqRec

def annotate_seq_to_record(SeqRec,inseq,out_label,seq_name,in_label):
    ThisSeqr=modify_seq_in_record(inseq,SeqRec)
    ThisSeqr.id=out_label
    ThisSeqr.mutlabel=in_label
    ThisSeqr.name=seq_name
    return ThisSeqr

def annotate_seq_to_record1(annotations,inseq,seq_id,seq_name,seq_description,seq_polarity,howmany):
    # Formerly annotate_seqrecord - purpose of this is to omit copying the full feature-table
    #ThisSeqr=SeqRecord(Seq(str(inseq),generic_dna)) #generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #ThisSeqr=SeqRecord(Seq(str(inseq)))
    ThisSeqr=Biopython_fix.fix_SeqRecord(inseq)
    ThisSeqr.id=seq_id
    ThisSeqr.name=seq_name
    ThisSeqr.description=seq_description
    ThisSeqr.polarity=seq_polarity
    ThisSeqr.howmany=howmany
    ThisSeqr.annotations=annotations
     #ThisSeqr=add_extra_objects(ThisSeqr)
    return ThisSeqr
# end of annotate_seq_to_record1(seqrec,seq_id,seq_name,seq_description,seq_polarity)

def annotate_frag_to_record(inseq,seq_id):
    # Only need minimal annotation for a fragment record
    #ThisSeqr=SeqRecord(Seq(str(inseq),generic_dna)) #generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #ThisSeqr=SeqRecord(Seq(str(inseq)))
    ThisSeqr=Biopython_fix.fix_SeqRecord(inseq)

    ThisSeqr.id=seq_id
    #Create a fake quality string for Fastq - only do it if fastq output is selected
    if RG_globals.is_fastq_out:
        ThisSeqr.letter_annotations=get_quality_list(ThisSeqr)
    return ThisSeqr
# end of annotate_frag_to_record(inseq,seq_id)

def setup_fastq_quality_list():
    # Prepare appropriate fastq lookups as necessary. Fastest solutions first
    global shared_phred_dict
    if RG_globals.Qualmax == RG_globals.Qualmin:
        shared_phred_dict=get_quality_list_qualmax()
    elif not RG_globals.is_fastq_random:
        make_fastq_slice()

def get_quality_list(Seqrec):
    global shared_phred_dict
    # Select quality list. Fastest solutions first
    if RG_globals.Qualmax == RG_globals.Qualmin: # Surprisingly not faster than get_quality_list_slice
        phred_dict=shared_phred_dict
    elif not RG_globals.is_fastq_random:
        phred_dict=get_quality_list_slice(Seqrec)
    else:
        phred_dict=get_quality_list_random(Seqrec) 
    return phred_dict
# end of get_quality_list(Seqrec)

def get_quality_list_random(Seqrec):
    #Returns a dictionary suitable for adding to SeqRecord as SeqRecord.letter_annotations
    #Sanger Phred quality scores range from 0 to 93
    #Old Illumina formats include -5 to 62 or 0 to 62
    #Note that quality scores tend to be lower at the ends & higher in the middle - potential refinement
    # NB: The use of randint slows execution, but not noticeably in Python
    #     the conversion of this to browser-enabled equivalent Dec 2020 showed a marked performance hit
    #     Later performance tests: very low sequence lengths or high coverage values are the main problem
    #                            solve the former by imposing a higher loweer-limit; latter by a lower upper-limit
    # Performance tests: the list-comprehension version appears to be slightly-faster, but not significantly
    #                    faster than the loop in the web version.
    # list-comprehenion version:
    return { "phred_quality": [randint(RG_globals.Qualmin,RG_globals.Qualmax) for x in range(len(Seqrec))] }
    '''
    # Explicit-loop version:
    quality_list= []
    for x in range(len(Seqrec)):
        quality_list.append(randint(RG_globals.Qualmin,RG_globals.Qualmax))
    phred_dict = {
        "phred_quality":quality_list
    }
    '''
# end of get_quality_list_random0(Seqrec)

def make_fastq_slice():
    global fastq_slice, slice_range
    slice_factor=9
    slice_range=(slice_factor+1)*RG_globals.Fraglen
    if slice_range > 3*RG_globals.FraglenMax:# Don't let it run away
        slice_range=3*RG_globals.FraglenMax
    fastq_slice= []
    for x in range(slice_range):
        fastq_slice.append(randint(RG_globals.Qualmin,RG_globals.Qualmax))
    #print("slice_range %s"%(slice_range))
    #print("fastq_slice %s"%(fastq_slice))

def get_quality_list_slice(Seqrec):
    global fastq_slice,slice_range
    # Return a slice from a previously-created randomly-generated list.
    # It is measurably quicker than get_quality_list_random: at all fragment-length and DOC values;
    # being comparable to the FASTA-only output speeds until high DOC values
    start=randint(0,slice_range-RG_globals.Fraglen-1)# 0-based
    #print("slice_start=%s"%start)
    return { "phred_quality":fastq_slice[start:start+RG_globals.Fraglen]}
# end of get_quality_list_slice(Seqrec)

def get_quality_list_qualmax():
    # Just return a list with all equal to highest value.
    return { "phred_quality": [RG_globals.Qualmax] * RG_globals.Fraglen }
# end of get_quality_qualmax(Seqrec)

def merge_feats(featlist1,featlist2):
    ''' Merge two sequence-feature lists in descending order, returning one list'''
    ''' Classic sorting-merge of two lists '''
    outfeats=[]
    while len(featlist1) >0 and len(featlist2)>0:
        feat1=featlist1[0]
        f1parts=feat1.location.parts
        for item in f1parts:
            ''' Only one item, so only does this once '''
            f1start=item.start
        feat2=featlist2[0]
        f2parts=feat2.location.parts
        for item in f2parts:
            ''' Only one item, so only does this once '''
            f2start=item.start
        if f1start > f2start:
            outfeats.append(feat1)
            x=featlist1.pop(0)
        else:
            outfeats.append(feat2)
            x=featlist2.pop(0)
    ''' One of the lists now empty '''
    while len(featlist1) >0:
        outfeats.append(featlist1[0])
        x=featlist1.pop(0)

    while len(featlist2) >0:
        outfeats.append(featlist2[0])
        x=featlist2.pop(0)
    return outfeats
# end of def merge_feats(featlist1,featlist2)

def update_journal(instring):
    instring=instring+"\n"
    return instring
    

def purl_features(seqdonor,vardonor):
    #seqdonor is most likely REFSEQ and vardonor is mutseq
    return_message=""
    
    '''
    purl_seqrecords takes the vardonor variant (relative) positions and recalculates that position into the
    seqdonor (relative) positions by first checking the absolute positions are compatible.

    These modified vardonor variant-definitions are then merged with the seqdonor sequence to return a new sequence record
    '''
     
    seqdonor_begin=int(seqdonor.abs_start)
    vardonor_begin=int(vardonor.abs_start)

    local_offset_diff=vardonor_begin-seqdonor_begin

    seqdonor_end=int(seqdonor.abs_end)
    vardonor_index=get_varfeature_index(vardonor)
     
   # print("vardonor_begin %s; seqdonor_begin %s; seqdonor_end %s"%(vardonor_begin,seqdonor_begin,seqdonor_end))

    GRCh_version_match=False
    polarity_match=False
    seqdonor_version=seqdonor.GRChver+":"+seqdonor.chrom
    vardonor_version=vardonor.GRChver+":"+vardonor.chrom

    seqdonor_range=seqdonor.abs_start+":"+seqdonor.abs_end
    vardonor_range=vardonor.abs_start+":"+vardonor.abs_end
   
    if (seqdonor_version == vardonor_version):
        GRCh_version_match=True
    if (seqdonor.polarity == vardonor.polarity):
        polarity_match=True

    seqdonor_version=seqdonor_version+":"+seqdonor.polarity
    vardonor_version=vardonor_version+":"+vardonor.polarity

    if GRCh_version_match and polarity_match:
        return_message+=update_journal(" Compatible GRCh build, chromosome and polarity: "+seqdonor_version)
        part_message="ranges: %s_range: %s; %s_range: %s"%(RG_globals.reference_gene,seqdonor_range,RG_globals.variants_label,vardonor_range)
        if local_offset_diff ==0:
            return_message+=update_journal(" Matching %s"%part_message)
        else:
            return_message+=update_journal(" Differing %s; offset_change: %s"%(part_message,local_offset_diff)) 

        for index in vardonor_index[:-1]:
            xref_vardonor=str(vardonor.features[index].qualifiers.get("db_xref"))
            seq_start=int(vardonor.features[index].location.start)+local_offset_diff 
            seq_end=int(vardonor.features[index].location.end)+local_offset_diff

            var_abs_start=get_absolute_position(seqdonor,seq_start)
            var_abs_end=get_absolute_position(seqdonor,seq_end)
        
            #print(" %s var_abs_start %s var_abs_end %s seqdonor_begin %s seqdonor_end %s"%(xref_vardonor,var_abs_start,var_abs_end,seqdonor_begin,seqdonor_end))
            if (seqdonor_begin<=var_abs_start <=seqdonor_end) and (seqdonor_begin<= var_abs_end <= seqdonor_end):
                f=SeqFeature(FeatureLocation(seq_start,seq_end,strand=1),type="variation")
                #f=SeqFeature(FeatureLocation(seq_start-1,seq_end,strand=1),type="variation")
                f.qualifiers=vardonor.features[index].qualifiers
                #print("purl; f.qualifiers %s"%f.qualifiers)
                vardonor.features[index]=f
                #print(" Feature RETAINED %s DONE"%SeqFeature(FeatureLocation(seq_start,seq_end,strand=1),type="variation"))
            else:
                ''' remove it '''
                del vardonor.features[index]
                #print(" Feature OMITTED %s DONE"%SeqFeature(FeatureLocation(seq_start,seq_end,strand=1),type="variation"))
                return_message+=update_journal(" OMITTED feature %s start: %s, end: %s because it is out-of-range for sequence start: %s, end: %s "
                                               %(xref_vardonor,var_abs_start,var_abs_end,seqdonor_begin,seqdonor_end))                    
    else:
        return_message+=update_journal(" Incompatible GRCh build, chromosome or polarity for sequence: "+seqdonor_version+" vs variants: "+vardonor_version)

        
    newfeatures=knit_features(seqdonor,vardonor)
    return newfeatures,(GRCh_version_match and polarity_match),return_message[:-1]
# end of purl_features(seqdonor,vardonor)


def knit_features(seqdonor,vardonor):
    non_varfeatures=get_filtered_features(seqdonor,False)
    varfeatures=get_varfeatures(vardonor)
    outfeatures=non_varfeatures+varfeatures
    return outfeatures

# =====================================================
# End of Seqrecord manipulation
# =====================================================
def get_outrefname(RefRecord,Ref_file_name): # This could replace similar code in RG_main.write_ref_fasta and RG_main.close_seqout, but I got lost
        extralabel,strandlabel=RG_globals.get_strand_extra_label()
        if RefRecord.is_ref_spliced and not RG_globals.is_frg_paired_end:            
            id_label="refhap"
            if  RG_globals.target_transcript_name == RG_globals.empty_transcript_name:
                id_label="REF"
            out_refname="%s%s_%s"%(RG_globals.get_locus_transcript(),strandlabel,id_label)
        else:
            out_refname="%s_%s"%(RG_globals.target_locus,Ref_file_name)# Defined as in_ref_src= in main...
        return out_refname

def make_addmut(refseqdonor,label):
    success=False; item_count=0
    for item in RG_globals.bio_parameters["target_build_variant"]["AddVars"]:
        #print("AddVars item: %s,target_locus %s "%(item,RG_globals.target_locus))
        if item["locus"]==RG_globals.target_locus and item["hapname"]==label:
            item_count+=1
            #print("Found item %s for label %s ; number %s"%(item,label,item_count))
            #f=SeqFeature(FeatureLocation(item["local_begin"]-1,item["local_end"],strand=1),type="variation")
            # Fixing the above for inconsistent order coming from App.js and RG_exploder_gui Feb 2024; it's RG_exploder_gui that's wrong
            mutstart=item["local_begin"]; mutend=item["local_end"]
            if mutstart > mutend:
                mutend,mutstart=mutstart,mutend
            f=SeqFeature(FeatureLocation(mutstart-1,mutend,strand=1),type="variation") # Turns mutstart from 1-based to 0-based
            qualifiers=OrderedDict()
            qualifiers["replace"]=["%s/%s"%(item["ref_seq"],item["var_seq"])]
            #qualifiers["comment"]=["%s_%s:%s"%(label,item_count,item["varname"])] # comes out bounded by unwanted braces
            qualifiers["comment"]=label+"_%s:"%item_count+"%s"%item["varname"]
            f.qualifiers= qualifiers
            #print("addmut: f.qualifiers %s %s"%(label,f))
            if not success:  # Only do this for the first variant in case there are > 1
                CopySeqRec=duplicate_record(refseqdonor)
            CopySeqRec.features.append(f)
            success=True
    if not success:
        CopySeqRec=""
    return CopySeqRec,success

def merge_seqvar_records(refseqdonor,mutseqrecipient,mutlabel):
    # The main point of writing this function was to merge two variants lists, but parsing out those within 'recipient' that overlap any exon regions in the 'donor'.
    # Needs tidying up:  some unnecessary calls and variable assignments which could be eliminated by using mutidx rather than its derivatives
    # Would be a good idea to separate out get_mutref_olap3 to be stand-alone and re-use for parsing any two lists of variants, not just exons as here

    ''' from Bio.SeqFeature import SeqFeature, FeatureLocation '''
    '''
    This merges the sequence & annotation from refseqdonor and both feature-tables
    This enables the merging of data from a file that has a nominal sequence in it (as resulting from using the function write_gb_features(x,y,z,"short"))
    which also holds a variant subset, originating initially from the refseqdonor. Here they are re-united.

    The second phase is to merge the *variants* list of the two sources.
    The only objective at this point of development is to take the definitions of exome region(s) from a reference sequence feature table
    (eg: "mRNA          join(0..0,156206..156304)" previously processed into the refseq variants table by splice_a_sequence as xref databases "skipnumx"
    and copy it into the mutseq table, so any introns are marked as "skipped" in the CIGAR string

    The only ref variants expected at this point of development are the <splice_trigger> ones
    so refend will always be > refstart.
    To avoid future complications all other variant types are currently excluded and therefore, when present in Refseq, are ignored.
    There may be a case to include this (mutations common to each mutseq) later on, but it's complicated enough as it is now!

    The third phase, merge_feats, attempts a proper sort of these variants into position-order
    '''
    return_message="" # accumulates journal messaging until return of merge_seqvar_records

    def retain_ref():
        save_ref_features=[]
        for refidx in ref_var_index:
            ref_feature=refseqdonor.features[refidx]
            ref_xref0=get_varfeats2(ref_feature)["xref_ids"][0]
            if ref_xref0 != RG_globals.skip_trigger:
                ''' Ignore any ref variant that is NOT a RG_globals.skip_trigger by skipping to next ref '''
                break
            else:
                save_ref_features.append(ref_feature)
        return save_ref_features
    # end of def retain_ref() which is local to merge_seqvar_records(refseqdonor,mutseqrecipient)

    def retain_mut(index):
        nonlocal save_mut_features
        f=mutseqrecipient.features[index]
        save_mut_features.append(f)
    # end of retain_mut(index) which is local to merge_seqvar_records(refseqdonor,mutseqrecipient)

    #def addmut_del(begin,end,style,mrf1,otype):
    def addmut_del(begin,end,style,otype):
        #print("addmut_del: begin %s,end %s,style %s,otype %s"%(begin,end,style,otype))
        nonlocal add_mut_features,mutidx
        featlist=get_varfeats2(mutseqrecipient.features[mutidx])

        #old_replace_string=featlist["replace_string"]
        new_ref_frag_string=str(refseqdonor.seq[begin-1:end])
        #old_ref_frag_string = ''.join(featlist["reference"].split())
        old_mut_frag_string=featlist["replace"]

        old_ref_frag_length=len(''.join(featlist["reference"].split())) # length of old_ref_frag_string
        new_ref_frag_length=len(new_ref_frag_string)
        diff_old_ref=abs(old_ref_frag_length-new_ref_frag_length)
        #old_mut_frag_length=len(old_mut_frag_string)
 
        truncval=str(diff_old_ref)
   
        if old_mut_frag_string == "-" :
            new_mut_frag_string="-"
            
        else:
            if style == "begintrunc":
                #begintrunc: Left-trim, right-retain
                new_mut_frag_string=old_mut_frag_string[diff_old_ref:]
                
            elif style == "endtrunc":
                #endtrunc:  Right-trim, left-retain
                new_mut_frag_string=old_mut_frag_string[:new_ref_frag_length]
                truncval=str(new_ref_frag_length)
            else:
                # Error!!
                #print("Error in style definition '%s' for addmut_del"%style)
                new_mut_frag_string="error"
                truncval="oops"
                
        if new_mut_frag_string=="":
            new_mut_frag_string="-"

        if refseqdonor.polarity==RG_globals.seq_polarity_minus:
            #compseq=MutableSeq(str(new_ref_frag_string),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
            #compseq=MutableSeq(str(new_ref_frag_string))
            compseq=Biopython_fix.fix_MutableSeq(new_ref_frag_string)

            #new_ref_frag_string=compseq.reverse_complement(inplace=True)
            new_ref_frag_string=Biopython_fix.fix_reverse_complement(compseq)
            
            if new_mut_frag_string !="-":
                #compseq=MutableSeq(str(new_mut_frag_string),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
                #compseq=MutableSeq(str(new_mut_frag_string))
                compseq=Biopython_fix.fix_MutableSeq(new_mut_frag_string)
                #new_mut_frag_string=compseq.reverse_complement(inplace=True)
                new_mut_frag_string=Biopython_fix.fix_reverse_complement(compseq)
                
        '''
        print("B: otype %s \n, old_ref_frag_string %s\n,old_mut_frag_string %s\n, new_ref_frag_string %s\n, new_mut_frag_string %s\n"
              %(otype, old_ref_frag_string, old_mut_frag_string, new_ref_frag_string, new_mut_frag_string))
        '''
            
        new_replace_string=str(new_ref_frag_string)+"/"+str(new_mut_frag_string)
        
        #print("new_replace_string %s\n"%new_replace_string)
        #print("Original features %s"%mutseqrecipient.features[mutidx])
        #print("Original replace %s"%mutseqrecipient.features[mutidx].qualifiers['replace'])

        f=SeqFeature(FeatureLocation(begin-1,end,strand=1),type="variation")

        f.qualifiers['replace']=[new_replace_string]

        old_dbxref=mutseqrecipient.features[mutidx].qualifiers.get("db_xref")
        #print("%s"%mutseqrecipient.features[mutidx].qualifiers)
        #print("old_dbxref %s"%old_dbxref)
        extra_dbxref=style+":"+truncval
        #print("extra_dbxref %s"%extra_dbxref)
        if old_dbxref!=None:
            old_dbxref.insert(0,extra_dbxref)
        else:
            old_dbxref=extra_dbxref
        #print("old_dbxref %s"%old_dbxref)
        f.qualifiers["db_xref"]=old_dbxref
        #print("f after: ",f)
        add_mut_features.append(f)
        #print("add_mut_features %s\n"%add_mut_features)
    # end of addmut_del(begin,end,style) which is local to merge_seqvar_records(refseqdonor,mutseqrecipient)


    def get_olaptype(mutstart,mutend,refstart,refend):
        global olap_and_base3
        ''' This assigns a digit to base3string in a specific order of relationships between each pair of the four input parameters
            The digit is based on the comparison operator between each pair of values
            value1  ? value2  (where ? is either <, = , > and assigned 0,1,2 accordingly)
            -----------------
            base3add is called in order with each pair as param1, param2 thus:
            mutstart ? refstart
            mutstart ? refend
            mutend ? refstart
            mutend ? refend
            mutstart ? mutend - least sig position identifies Deletion/complex, SNV, Insert
            the resulting base3string is used to assign a 'type' of relationship between the mutation & reference range using olap_and_base3.
            The type is documented outside of this code - see comment in olap_and_base3
        '''
        def base3add(param1,param2,instring):
            if param1 < param2:
                outstring=instring+"0"
            elif param1 == param2:
                outstring=instring+"1"
            else:
                outstring=instring+"2"
            return outstring

        base3string=base3add(mutstart,refstart,"")
        base3string=base3add(mutstart,refend,base3string)
        base3string=base3add(mutend,refstart,base3string)
        base3string=base3add(mutend,refend,base3string)
        base3string=base3add(mutstart,mutend,base3string)
        olaptype=olap_and_base3[base3string]
        #olapval=int(base3string)- # deprecated
        return olaptype
    # end of get_olaptype(mutstart,mutend,refstart,refend):

    def get_mutref_olap3(showit):
            nonlocal return_message
            ''' skip-trigger ranges (Reference start Rs to Reference end, Re) are introns and so are effectively deletions
                to the reference sequence that are introduced before the resulting "exonic" sequence is fragmented.
                This routine looks for such "skips" in the feature table and determines whether any variant feature ranges
                (Mut-start Ms to Mut-end Me) overlap this same "skip" range.
                The type of overlap is assigned a value 1-25.
                A document exists that defines Rs,Re,Ms,Me relative positions.
                Types 1 to 3 and 23 to 25 are classified as non-overlaps.
                The others are either entirely within Rs-Re or partially overlapping the boundaries between the two ranges.
                Those *within* the Re-Rs range are 'masked', and removed from the feature table.
                Those features overlapping the boundary are trimmed to retain just those parts of the features that are outside Re-Rs
                '''
            ''' Third version, Dec-2020 eliminating the need for complex nested if..then...elif heirarchy is to determine a base3 number '''

            mut_feature=mutseqrecipient.features[mutidx]
            mutparts=mut_feature.location.parts
            mut_replace_type=get_varfeats2(mut_feature)["replace_type"]
            olapcount=0; olaptype=0; mutstart=0; mutend=0; refstart=0; refend=0

            olap_params={"olapcount":olapcount,"olaptype":olaptype,"mutstart":mutstart,"mutend":mutend,"refstart":refstart,
                        "refend":refend,"mut_replace_type":mut_replace_type,"ref_xref0":'',"ref_xref1":''}
            if showit:
                print("mutlabel: %s"%mutlabel)
            for item in mutparts:
                '''Only does this once '''
                # adding +1 to mutstart because things are looking wrong
                mutstart=item.start+1; mutend=item.end
                if showit:
                    print("   mutref: %s, mut_replace_type: %s"%(get_varfeats2(mut_feature)["xref_ids"],mut_replace_type))
                    print("   mutstart: %s, mutend: %s"%(mutstart,mutend))

            for refidx in ref_var_index:
                ref_xref0=get_varfeats2(refseqdonor.features[refidx])["xref_ids"][0]

                if ref_xref0 != RG_globals.skip_trigger:
                    ''' Ignore any ref variant that is NOT a RG_globals.skip_trigger by skipping to next ref '''

                    if RG_globals.is_print_diagnostics:
                        print("ref_xref0 %s is NOT %s"%(ref_xref0,RG_globals.skip_trigger))
                    break
                else:
                    refparts=refseqdonor.features[refidx].location.parts
                    ref_xref1=get_varfeats2(refseqdonor.features[refidx])["xref_ids"][1]

                for item in refparts:
                    ''' Only does this once '''
                    refstart=item.start+1; refend=item.end
                    # NB: this range is the INTRON as it's a skip-region
                    #print("     refstart: %s , refend: %s"%(refstart,refend))

                if refend <= refstart:
                    ''' As with the RG_globals.skip_trigger test above, - only interested in Re > Rs such as introns or subset '''
                    ''' Not expected to occur !!'''
                    return_message+=update_journal("    *** refend <= refstart not expected ***")
                    #print("    *** refend <= refstart not expected ***")
                    break

                # Now have the four relevant parameters, call the overlap type
                olaptype=get_olaptype(mutstart,mutend,refstart,refend)

                # This is a print line to check overlap parameters
                if showit:
                    print("      refstart: %s, refend: %s, olaptype: %s, ref: %s:%s"%(refstart,refend,olaptype,ref_xref0,ref_xref1))

                if (olaptype > 3) and (olaptype < 23) :
                    ''' Retain first overlap type discovered as olap_found. 1,2,3, 23, 24  & 25 are non-overlaps. '''
                    ''' Function returns olaptype:0 if all overlaps outside this range'''
                    '''
                        print("ref_xref[0]",ref_xref0,"refstart ",refstart," refend ",refend,
                          " mutstart ",mutstart," mutend ",mutend," olaptype ", olaptype)
                    '''
                    '''nonlocal olapcount,mutstart,mutend,refstart,refend,mut_replace_type'''
                    olapcount+=1
                    olap_params={"olapcount":olapcount,"olaptype":olaptype,"mutstart":mutstart,"mutend":mutend,"refstart":refstart,
                                 "refend":refend,"mut_replace_type":mut_replace_type,"ref_xref0":ref_xref0,"ref_xref1":ref_xref1}
                    break
            return olap_params
        # end of get_mutref_olap3() which is local to merge_seqvar_records(refseqdonor,mutseqrecipient)

    # This is where merge_seqvar_records(refseqdonor,mutseqrecipient) really starts after the local functions above are declared

    #CopySeqRec=annotate_seq_to_record(refseqdonor.seq,refseqdonor.id,refseqdonor.name,refseqdonor.description,refseqdonor.polarity,refseqdonor.howmany)
    CopySeqRec=duplicate_record(refseqdonor)

    add_mut_features=[]
    save_mut_features=[]
    ref_var_index=get_varfeature_index(refseqdonor)

    # This section needs correcting for it to be used. Never evaluates to < 1
    # if len(ref_var_index) < 1:
    #    ''' No features in ref to process'''
    #    CopySeqRec.features=get_filtered_features(mutseqrecipient,True)
    #    return_message=" No variant features in reference to merge with %s\n"%mutlabel
    # End This section

    if not CopySeqRec.is_ref_spliced:
        CopySeqRec.features=get_filtered_features(mutseqrecipient,True)
        if RG_globals.is_frg_paired_end or (RG_globals.target_transcript_name == RG_globals.empty_transcript_name):
            return_message+=update_journal(" Not splicing %s"%RG_globals.empty_transcript_name)
        else:
            return_message+=update_journal(" No gap features in reference to merge with %s"%mutlabel)
    else:
        return_message+=update_journal(" Merging gap features from %s with variant features from %s"%(RG_globals.bio_parameters["target_transcript_name"]["label"],mutlabel))
        mut_var_index=get_varfeature_index(mutseqrecipient)
        for mutidx in mut_var_index:
            mask=False
            # This became such a large block of code that it was moved to get_mutref_olap() to improve the clarity of this loop
            printit=False
            olap_return_params=get_mutref_olap3(printit) # Third version, uses base3 to determine relationship
            olapcount=olap_return_params["olapcount"]

            #print("Return from get_mutref_olap3 with olaptype at %s"%olap_return_params["olaptype"])
            if olapcount > 1:
                return_message+=update_journal("WARNING: overlap count not expected to exceed 1 but is olapcount %s . for position %s"%(olapcount,mutseqrecipient.features[mutidx]))
                if printit:
                    print(return_message)
                pass

            elif olapcount == 0:
                ''' Should only be true where overlap type has been <4  or > 21 (no overlap) in each case '''
                ''' mask already set to False '''
                '''print("   n  `    No overlap: refstart: %s, refend: %s, olaptype: %s, xref %s:%s\n"
                      %(olap_return_params["refstart"],olap_return_params["refend"],
                        olap_return_params["olaptype"],
                        olap_return_params["ref_xref0"],olap_return_params["ref_xref1"]))
                '''
                if printit:
                    print("   No overlap\n")
                pass

            else:
                ''' leaving us with real overlaps to sort'''
                olaptype=olap_return_params["olaptype"]
                mutstart=olap_return_params["mutstart"]
                mutend=olap_return_params["mutend"]
                refstart=olap_return_params["refstart"]
                refend=olap_return_params["refend"]
                if printit:
                    print("   Overlap type %s\n"%olaptype)
                #print("   Overlap type %s\n"%olaptype)

                ###########  Overlap numbering consistent with get_mutref_olap3() #############
                ''' NB: 1,2 & 3 are non-overlaps and don't get here'''
                if olaptype==4:
                    ''' mutend coincides with refstart, so retain all in the deletion range except the last base '''
                    addmut_del(mutstart,mutend-1,"endtrunc",olaptype) # Adds a deletion feature to replace this one ...
                    mask=True # ... and hide the original one

                elif olaptype==5:
                    ''' mutstart coincides with refstart, but mutend is beforehand, meaning it's an insert just at the boundary '''
                    ''' current judgement is to let this remain mask = False '''
                    pass

                elif olaptype==6 or olaptype==7:
                    ''' mutstart (deletion start point) is prior to start of ref and ends within ref (6) or at end of ref (7).
                        so put in a truncated end-deletion '''
                    #addmut_del(mutstart,refstart,"endtrunc",olaptype) # Adds a deletion feature to replace this one ...
                    addmut_del(mutstart,refstart-1,"endtrunc",olaptype) # Adds a deletion feature to replace this one ...
                    
                    mask=True # ... and hide the original one

                elif olaptype==8:
                    ''' deletion includes entire length of ref and beyond boundary both sides.
                        so create two deletions '''
                    
                    addmut_del(mutstart,refstart,"begintrunc",olaptype)
                    addmut_del(refend+1,mutend,"endtrunc",olaptype)
                    mask=True # ... and hide the original one

                elif olaptype==13 or olaptype==19:
                    ''' mutstart (deletion start point) is at start of ref (13) or after (19), but mutend is beyond refend .
                        so truncate deletion start-point to refend+1 '''
                    addmut_del(refend+1,mutend,"begintrunc",olaptype)
                    mask=True

                elif olaptype==21:
                    ''' mutstart coincides with refend .
                        so increase deletion start-point by one position '''
                    addmut_del(mutstart+1,mutend,"begintrunc",olaptype)
                    mask=True

                else:
                    ''' Leaving 9, 10, 11, 12, 14, 15, 16, 17, 18, 20, 22  which are masked completely'''
                    ''' NB: 23, 24 and 25 are non-overlaps and don't get here'''
                    mask=True

                    ###########  Overlap numbering consistent with get_mutref_olap3() #############

            if mask:
                ''' simply avoid appending the features to the holding place for mut features '''
                #print("Masking olaptype%s feature %s"%(olaptype,mutseqrecipient.features[mutidx]))
                pass
            else:
                ''' hold onto these mut features '''
                #print("Retaining feature %s"%(mutseqrecipient.features[mutidx]))
                retain_mut(mutidx)
                '''replace_mut_features.append(mutseqrecipient.features[mutidx])'''
        # ''' end of loop: for mutidx in mut_var_index '''
        # ''' now build a new CopySeqRec.features from the new features list '''
        #print("save_mut_features: %s \n"%save_mut_features)
        #print("add_mut_features: %s \n"%add_mut_features)
        # Order the features first by end-feature
        sorted_add_mut_features=order_featend(add_mut_features,True)
        # Now order  by start-feature.
        # In that way, we better order those that start in the same place but differ in second
        # Let's just ignore the fact that two features starting in the same place is a bad idea - but works well on test-set
        sorted_add_mut_features=order_featstart(sorted_add_mut_features,True)
        merged_mutseq_features=merge_feats(save_mut_features,sorted_add_mut_features)
        merged_mutseq_ref_features=merge_feats(merged_mutseq_features,retain_ref())
        non_varfeatures=get_filtered_features(mutseqrecipient,False)
        #print("Nonvars: %s"%non_varfeatures)

        CopySeqRec.features=non_varfeatures+merged_mutseq_ref_features
    '''return CopySeqRec'''
    return CopySeqRec,return_message[:-1]
# end of def merge_seqvar_records(refseqdonor,mutseqrecipient,mutlabel)


# ===============================================================
# Start of data-structure query functions
# ===============================================================

def get_varfeature_index(seq_record):
    #print("At get_varfeature_index")
    #print("seq_record.id %s"%seq_record.id)
    var_feature_index,feature_index,source_gene_index=filter_varfeature_index(seq_record)
    #print("At get_varfeature_index %s"%var_feature_index)
    return var_feature_index
# end of get_varfeature_index(seq_record)

def get_varfeatures(SeqRec):
    '''
    Returns a list of all the filtered variants
    '''
    var_feature_index=get_varfeature_index(SeqRec)
    ''' Next line is unncessary as get_varfeature_index has already called order_featindex '''
    '''var_feature_index_sorted=order_featindex(var_feature_index,SeqRec)'''
    features=[]
    for (index) in var_feature_index:
       features.append(SeqRec.features[index])
    return features
# end of get_varfeatures(SeqRec)

def get_filtered_features(seq_record,total):
    '''
    Returns a list of features plus the filtered variants
    Total is Boolean: if True, all the features are returned
                      if False only the non-variant features are returned
    '''
    var_feature_index,other_feature_index,source_gene_index=filter_varfeature_index(seq_record)
    features=[]
    for (index) in source_gene_index:
        features.append(seq_record.features[index])
    for (index) in other_feature_index:
        features.append(seq_record.features[index])
    if total:
        for (index) in var_feature_index:
            features.append(seq_record.features[index])
    return features
# end of get_filtered_features(seq_record):

def get_source_gene_features(seq_record):
    var_feature_index,feature_index,source_gene_index=filter_varfeature_index(seq_record)
    features=[]
    for (index) in source_gene_index:
        features.append(seq_record.features[index])
    return features

def filter_varfeature_index(seq_record):
    '''
    Grab each item in the variants table and list positions of wanted ones
    '''

    #print("At filter_varfeature_index")
    #print("seq_record.id %s"%seq_record.id)
    #print("seq_record.features %s"%seq_record.features)

    feature_index=[]
    varfeat_index=[]
    source_gene_index=[]
    for (index, feature) in enumerate(seq_record.features):
        feature=seq_record.features[index]
        '''
        print(feature.type)
        '''
        if feature.type == 'variation':
            replace=str(feature.qualifiers.get("replace"))
            '''
            print(replace)
            '''
            '''
            We want to filter out where replace is not what we expect eg:HGMD_MUTATION
            '''
            for (item) in RG_globals.unwanted_replace:
                is_keep_feature=True
                if (item in replace):
                    is_keep_feature=False
                if (is_keep_feature):
                    varfeat_index.append(index)
        elif feature.type == 'gene' or feature.type == 'source':
                source_gene_index.append(index)
                # print("seq_record.features[index]%s "%seq_record.features[index])
        else:
            feature_index.append(index)
    varfeat_index=order_featindex(varfeat_index,seq_record,True)
    return varfeat_index,feature_index,source_gene_index
# end of filter_varfeature_index(seq_record)

def get_varfeats(this_feature):
    '''
    Extracts feature information, returns multiple values
    Initially to shorten code in main (became MutateVarSeq), but also used by order_featindex
    '''

    '''
    print(this_feature)
    '''
    feat_start=this_feature.location.start
    feat_end=this_feature.location.end

    if (this_feature.location.strand == 1):feat_polarity=RG_globals.seq_polarity_plus
    elif (this_feature.location.strand == -1):feat_polarity=RG_globals.seq_polarity_minus
    else: feat_polarity=RG_globals.seq_polarity_none

    xref=str(this_feature.qualifiers.get("db_xref"))
    '''
    xref_id=xref.strip("[]'")
    xref_ids=xref_id.split(":")
    '''
    xref2=xref.split(",")
    xref_ids=xref2[0].strip("[]'").split(":")
    replace_string=str(this_feature.qualifiers.get("replace")).strip("[]'")
    # Defend against unintended whitespace of all types in replace string
    replace_string = ''.join(replace_string.split())

    if len(replace_string) < 3:
        if RG_globals.is_print_diagnostics:
            print("problem this_feature %s"%this_feature)
        reference="err"
        replace="err"
    else:
        varts=replace_string.split("/")
        reference=str(varts[0])
        replace=str(varts[1])

    #if RG_globals.is_print_diagnostics:
    #print("xref %s, reference %s, replace %s | varts %s"%(xref,reference,replace, varts))

    replace_type=""
    if (reference == "-"):
        replace_type="ins_"
    elif (replace == "-"):
        replace_type="del_"
    elif (len(reference) !=1 ) or (len(replace) !=1 ):
        replace_type="delins_"
    elif (len(reference) ==1 ) and (len(replace) ==1 ):
        replace_type="sub_"
    else:
        replace_type=replace
    #print("Replace_type %s"%replace_type)
    #return xref_ids,reference,replace,replace_string,feat_start,feat_end,feat_polarity,replace_type
    return xref_ids,reference,replace,replace_string,feat_start,feat_end,feat_polarity,replace_type
# end of get_varfeats(this_feature)

def get_varfeats2(this_feature):
    '''
    Reformats the Genbank feature information into a dictionary list
    '''
    '''
    print(this_feature)
    '''

    xref_ids,reference,replace,replace_string,feat_start,feat_end,feat_polarity,replace_type=get_varfeats(this_feature)
    #print("Replace_type %s"%replace_type)
    featurelist={"xref_ids":xref_ids,"reference":reference,"replace":replace,"replace_string":replace_string,
                 "feat_start":feat_start,"feat_end":feat_end,"feat_polarity":feat_polarity,"replace_type":replace_type}
    #print("Replace_type %s"%featurelist["replace_type"])
    return featurelist
# end of get_varfeats2(this_feature)

def order_featindex(feature_index,seq_record,direction):
    def get_key(in_idx):
        nonlocal seq_record
        return get_varfeats2(seq_record.features[in_idx])["feat_start"]

    feature_index.sort(key=get_key,reverse=direction)
    #sorted_feature_index=sorted(feature_index,key=get_key)
    #print("SORTING feature_index")
    return feature_index

# end of order_featindex(feature_index,seq_record)


def order_featstart(feature_list,direction):
    def get_featkey(feat1):
        f1parts=feat1.location.parts
        for item in f1parts:
            ''' Only one item, so only does this once '''
            f1start=item.start; f1end=item.end
            if f1end < f1start:
                f1start=f1end
        return f1start

    feature_list.sort(key=get_featkey,reverse=direction)
    #sorted_feature_list=sorted(feature_list,key=get_featkey)
    #print("SORTING featstart")
    return feature_list
# end of order_featstart(feature_list)

def order_featend(feature_list,direction):
    def get_featkey(feat1):
        f1parts=feat1.location.parts
        for item in f1parts:
            ''' Only one item, so only does this once '''
            f1start=item.start; f1end=item.end
            if f1start > f1end:
                f1end=f1start
        return f1end

    feature_list.sort(key=get_featkey,reverse=direction)
    #sorted_feature_list=sorted(feature_list,key=get_featkey)
    #print("SORTING featend")
    return feature_list
# end of order_featend(feature_list)

def get_featlimits(this_feature): # Not apparently used
    featurelist=get_varfeats2(this_feature)
    return featurelist["feat_start"],featurelist["feat_end"]
# end of get_featlimits(this_feature)

def get_revcomp(inseq):
    instr=str(inseq)
    if instr !="" and instr !="-":
        #compseq=MutableSeq(instr,generic_dna) #generic_dna # Deprecated from Biopython 1.78 (September 2020)
        #compseq=MutableSeq(instr) #generic_dna # Deprecated from Biopython 1.78 (September 2020)
        compseq=Biopython_fix.fix_MutableSeq(instr)
        #compseq.reverse_complement(inplace=True)
        Biopython_fix.fix_reverse_complement(compseq)
    else:
        compseq=instr
    return str(compseq)

def switch_Rec_polarity(SeqRec):
    #print("IN switch: %s"%SeqRec.seq[0:5])
    #rev_inseq=MutableSeq(get_revcomp(SeqRec.seq),generic_dna)#generic_dna # Deprecated from Biopython 1.78 (September 2020)
    #rev_inseq=MutableSeq(get_revcomp(SeqRec.seq))
    rev_inseq=Biopython_fix.fix_MutableSeq(get_revcomp(SeqRec.seq))
    
    ModSeqRec=modify_seq_in_record(rev_inseq,SeqRec)
    
    ModSeqRec.Headclip,ModSeqRec.Tailclip=ModSeqRec.Tailclip,ModSeqRec.Headclip # Switch head & tail clips (Reference origin)
    ModSeqRec.endlocus,ModSeqRec.rev_endlocus=ModSeqRec.rev_endlocus,ModSeqRec.endlocus # Switch end locus (Reference origin)
    
    if ModSeqRec.polarity == RG_globals.seq_polarity_minus:
        ModSeqRec.polarity = RG_globals.seq_polarity_plus
        ModSeqRec.offset=int(SeqRec.abs_start)
    else:
        ModSeqRec.polarity = RG_globals.seq_polarity_minus
        ModSeqRec.offset=int(SeqRec.abs_end)+1

    ModSeqRec.strand_mod=int(ModSeqRec.polarity)
    
    splits=SeqRec.firstid.split(":")
    ModSeqRec.firstid="%s:%s:%s:%s:%s:%s"%(splits[0],splits[1],splits[2],splits[3],splits[4],str(ModSeqRec.polarity))
    #print("OUT switch: %s"%ModSeqRec.seq[0:5])
    return ModSeqRec

def add_seqrec_objects(SeqRec,splits):
    SeqRec.GRChver=splits[1]
    SeqRec.chrom=splits[2]
    SeqRec.abs_start=splits[3]
    SeqRec.abs_end=splits[4]
    SeqRec.polarity=splits[5]
    SeqRec.strand_mod=int(SeqRec.polarity)
    if SeqRec.polarity == RG_globals.seq_polarity_minus:
        offset=int(SeqRec.abs_end)+1
        #offset=int(SeqRec.abs_end)
    # Empirical discovery that need to add +1 to offset when polarity is -1 to get same position-values as +1 polarity
    else:
        offset=int(SeqRec.abs_start)
    SeqRec.offset=offset # WHOAH!! Different to main.get_muttranscripts(). May explain all the +1/-1 malarkey
    SeqRec.firstid="%s:%s:%s:%s:%s:%s"%(splits[0],splits[1],splits[2],splits[3],splits[4],splits[5])
    locus_begin,locus_end=RG_globals.Reference_sequences[RG_globals.target_locus]["Locus_range"].split(":")
    begin=int(locus_begin)+int(offset)-1; end=int(locus_end)+int(offset)-1
    if RG_globals.Reference_sequences[RG_globals.target_locus]["is_join_complement"]:
        pol="reverse"
    else:
        pol="forward"
    SeqRec.locus_range="%s:%s:%s:%s:%s strand"%(splits[1],splits[2],begin,end,pol)
    return SeqRec
# end of add_seqrec_objects(SeqRec,splits)

def add_extra_objects(SeqRec):
    success=False
    splits=SeqRec.id.split(":")
    if len(splits) == 6:
        success=True
        SeqRec=add_seqrec_objects(SeqRec,splits)
    return success,SeqRec
# end of add_extra_objects(SeqRec)

def get_absolute_position(seq_record,seqpos):
    abs_pos=seq_record.offset+seq_record.strand_mod*seqpos
    # Yet another ... if polarity is +1, need to adjust by -1
    if seq_record.strand_mod == 1:
        abs_pos-=1
    #print("get_absolute_position; strand: %s; pos:%s| abs %s"%(seq_record.strand_mod,seqpos,abs_pos))
    return abs_pos

def get_absolute_pairs(start,end,ref_offset,polarity):
    # called by MutateVarSeq to get the absolute positions from a variation-definition.
    abs_start=ref_offset+int(polarity)*start
    abs_end=ref_offset+int(polarity)*end

    if polarity == RG_globals.seq_polarity_minus:
            # Swap start and end when polarity is minus to derive absolute (forward strand) positions
            abs_start,abs_end=abs_end,abs_start       
    abs_end -=1
    #print("get_absolute_pairs; strand: %s; start: %s; end: %s| abs_start:%s; abs_end: %s"%(polarity,start,end,abs_start,abs_end))
    return abs_start,abs_end

def set_seqrec_absolutes(seq_record):
    # called by read_mutrecord if is_use_absolutes True and only seen if is_write_ref_ingb True, but will follow get_absolute_pairs call
    # based on mashup_gb10.modify_seq_record5(seq_record) - called immediately after read_refseqrecord
    #  Adds two extra feature qualifiers - stating absolute position of feature start and feature end
    #print("At set_seqrec_absolutes")
    #print("seq_record.features %s"%seq_record.features)

    success,seq_record=add_extra_objects(seq_record)
    if success:
        seq_record=set_seqrec_absolutes2(seq_record)
    return success,seq_record

def set_seqrec_absolutes2(seq_record):
    # based on mashup_gb10.modify_seq_record5(seq_record) - called during read_refmutrecord, after muts have been processed
    # therefore after get_absolute_pairs has been prev called
    #  Adds extra feature qualifiers - stating absolute position of feature start and feature end
    for index in get_varfeature_index(seq_record):
        varfeats= get_varfeats2(seq_record.features[index])

        feat_start=int(varfeats['feat_start'])
        feat_end=int(varfeats['feat_end'])

        #print("xref_ids %s, feat_start %s ,feat_end %s"%(varfeats['xref_ids'],feat_start,feat_end))

        feat_abs_start=seq_record.offset+seq_record.strand_mod*feat_start
        feat_abs_end=seq_record.offset+seq_record.strand_mod*feat_end

        if seq_record.polarity == RG_globals.seq_polarity_minus:
            feat_abs_start,feat_abs_end=feat_abs_end,feat_abs_start
        feat_abs_end -=1 # Another normalising-for-polarity +1 step
       
        global_range="%s:%s:%s:%s:1"%(seq_record.GRChver,seq_record.chrom,feat_abs_start,feat_abs_end)
        seq_record.features[index].qualifiers['global_range']=global_range

        if RG_globals.is_print_diagnostics:
            print("index %s , feature %s"%(index,seq_record.features[index]))
    return seq_record
# end of set_seqrec_absolutes2

# ===============================================================
# End of data-structure query functions
# ===============================================================


# ===============================================================
# Start of cigar-manipulation functions definition
# ===============================================================
def reverse_cigarbox(cigarbox):
    # Added 04-July 2022
    #print("input cigarbox: %s"%cigarbox)
    #print("fixed cigarbox: %s, length %s"%(cigarbox,len(cigarbox)))
    #print("cigar_label: %s, rev_cigar_label %s"%(cigar_label,rev_cigar_label))
    if len(cigarbox) > 4:
        out_cigarbox=[0]
        backcount=len(cigarbox)-3
        while backcount >0:
            #print("%s %s"%(cigarbox[backcount+1],cigarbox[backcount-1]))
            if cigarbox[backcount+1]=="I" and (cigarbox[backcount-1]=="D" or cigarbox[backcount-1]=="I"):
                out_cigarbox.append(cigarbox[backcount-2])
                out_cigarbox.append(cigarbox[backcount-1])
                out_cigarbox.append(cigarbox[backcount])
                out_cigarbox.append(cigarbox[backcount+1])
                backcount-=2
            else:
                out_cigarbox.append(cigarbox[backcount])
                out_cigarbox.append(cigarbox[backcount+1])
            backcount-=2
        out_cigarbox.append(cigarbox[-1])
        #print("in_cigarbox %s\nout_cigarbox %s\n"%(cigarbox,out_cigarbox))
    else:
        out_cigarbox=cigarbox    
    return out_cigarbox
 
def fill_cigar(cigarbox,var_start,var_end,mut_type,length):
    '''
    This gets the CIGAR correctly recorded. It is only called from MutateVarSeq, but in several places
    If the sequence has not been modified yet, the first and last items will be the sequence length and so identical.
    (don't compare [0] and [1] as it's possible they will be equal by chance)
    If the sequence has been modified, need to test if the previous modification was the same as the one about to be done.
    If so, just increase the count of nucleotides modified on the reference.
    Otherwise treat as new modification
    '''

    '''
    First Adds Ms from previous mutated position to next-to-be-mutated. Was padM(cigarbox,var_end)
    '''
    ref_gap=cigarbox[0]-var_end
    if ref_gap > 0:
        cigarbox.insert(1,"M")
        cigarbox.insert(1,ref_gap)
    ''' end of old padM '''
    ''' If don't cast to int - this becomes ExactPosition(var_start) for unknown reason.
        It doesn't seem to blow the algorithms, but unclear why this happens'''
    cigarbox[0]=int(var_start)
    is_packed = False
    if cigarbox[1] != cigarbox[len(cigarbox)-1]:
        if cigarbox[2]==mut_type:
            cigarbox[1]=cigarbox[1]+length
            is_packed=True
    if not is_packed:
        cigarbox.insert(1,mut_type)
        cigarbox.insert(1,length)
        if RG_globals.is_print_diagnostics: print("pack %i %s"%(length,mut_type))
    return cigarbox
# end of def fill_cigar(cigarbox,var_start,var_end,mut_type,length)

def pad_cigarbox(cigarbox):
    # Pads the front end of a CIGAR box with Ms
    # Called by RG_main.make_allvars_in_one_seq
    if cigarbox[0] > 0:
        cigarbox.insert(1,"M")
        cigarbox.insert(1,cigarbox[0])
        #Protect it from further manipulation by this function by ending cigarbox with 0
        cigarbox[0]=0 
    return

def get_untrimmed_cigar(cigarbox):
    # Builds a CIGAR string based on contents of the 'cigarbox' list previously built in RG_main.make_allvars_in_one_seq & RG_main.MutateVarSeq
    cigar=""
    index=1; 
    while index<len(cigarbox)-2:
        cigar+="%s%s"%(cigarbox[index],cigarbox[index+1])
        index+=2
    return cigar
# end of get_untrimmed_cigar(cigarbox)

def get_trimmed_cigars(cigarbox,mut_box,start,length):
# Get a new cigar for forward and reverse between the incoming sequence limits   
# NB: cigarbox is returned changed if popped directly
    def popleads(trimbox):
        clear=False
        while not clear and len(trimbox)>2:# Pop leading Ns from cigarbox
            if "N" in trimbox[2]:
                trimbox.pop(1)
                trimbox.pop(1)
            else:
                clear=True

    def popends(trimbox):
        clear=False
        while not clear and len(trimbox)>2:# Pop trailing Ns from cigarbox
            if "N" in trimbox[-2]:
                trimbox.pop(-2)
                trimbox.pop(-2)
            else:
                clear=True
                
    fwd_offset,fwdcigarbox=trimcigarbox(cigarbox,mut_box,start,length)
    revcigarbox=reverse_cigarbox(fwdcigarbox)
    popleads(fwdcigarbox)
    fwd_cigar=get_untrimmed_cigar(fwdcigarbox)
    popleads(revcigarbox)
    popends(revcigarbox)
    rev_cigar=get_untrimmed_cigar(revcigarbox)
    
    #rev_offset,rev_cigar=trimcigar(rev_cigarbox,rev_mut_box,start,length)
    
    monitor=0
    if monitor: #and (is_mut_cigar(fwd_cigar) or is_mut_cigar(rev_cigar)):
        print("cigarbox %s"%cigarbox)
        print("fwdcigarbox %s"%fwdcigarbox)
        print("revcigarbox %s"%revcigarbox)
        print("mutbox %s"%mut_box)
        print("fwd_cigar %s"%fwd_cigar)
        print("fwd_offset %s\n"%fwd_offset)
        
       # print("rev_cigarbox %s"%rev_cigarbox)
       # print("rev_mutbox %s"%rev_mut_box)
        print("rev_cigar %s"%rev_cigar)
       # print("rev_offset %s\n"%rev_offset)
        print("*****")
    #print("fwd_offset %s, rev_offset %s, fwd_cigar %s, "%(fwd_offset,rev_offset,fwd_cigar,rev_cigar))
    return fwd_offset,fwd_cigar,rev_cigar


def haveamutbox(cigarbox):

    '''
    cigarbox states "how-many-bases-to-next-feature-change" within a mutated sequence,stated as a comparison to the Reference
    5M 1D 3I 2X means: 5 Matches, 1 Deletion, 3 Inserts, 2 mismatches(X)
    cigarbox is used to build a CIGAR string: an annotation of how the "sequence read" aligns with the Reference

    mutbox is a list that states the cumulative-offset from the start by successively adding the "how_many-bases..",
    thereby creating a second-order version of the CIGAR string: "how many bases to next-feature change" in the "sequence read"
    5M 5D 8I 10X  - here Deletions & Ns don't move along the sequence.

    The mutbox is used in trimcigar to assist in the correct definition of the CIGAR string when characterising a fragment within the "sequence read".
    Within the purposes of this application: the Reference is mutated based on a list of variants and an entire sequence treated as a
    "sequence read", potentially way in excess of the length of a typical read.

    This "sequence read" is then fragmented to create a realistic read-length.
    The CIGAR needs to be "sliced" correctly to match that fragment, which is done in trimcigar, using both cigarbox and mutbox.

    It might be obvious only to call this routine from trimcigar, but it would mean multiple calls to here during a fragmentation.
    The alternative is to make this routine go away by incorporating it within trimcigar, and it's probably possible,
    but trimcigar is complicated enough as it is!
    Doing it this way (in advance, once) has the advantage of providing a lookup table for trimcigar in the form of mutbox

    '''
    if RG_globals.is_print_diagnostics: print("cigarbox: %s"%(cigarbox))
    mutbox=cigarbox
    if cigarbox!="":
        mut_box=[0]
        mut_pos=0
        '''
        if (cigarbox[2]=="D") or (cigarbox[2]=="N"):
            mut_pos=int(cigarbox[1])
        '''
        index=1
        while index<len(cigarbox)-2:
            length=cigarbox[index]
            mut_type=cigarbox[index+1]
            if RG_globals.is_print_diagnostics: print("mutbox: length %s, type %s"%(length,mut_type))
            if (mut_type != "D") and (mut_type != "N"):
                '''Ds and Ns are Reference-feature offsets of absent nucleotides, so do not contribute to length'''
                mut_pos+=length
            cig=mut_type
            mut_box.append(mut_pos)
            mut_box.append(cig)
            index+=2
        mut_box.append(cigarbox[len(cigarbox)-1])
        mut_box.append("E")
        mutbox=mut_box
    return mutbox
# end of haveamutbox(cigarbox)


def trimcigarbox(cigarbox,mut_box,start,length):
    #Take an existing cigar box and slice it accurately.
    #This is tricksy, begging for an easier-to-follow algorithm
    seqpos=start
    end=start+length
    added=False
    ref_offset=start
    mutboxindex=1
    newcigarbox=[0,0]
    while mutboxindex<len(mut_box)-2:
        length=0
        mut_pos=mut_box[mutboxindex]
        mut_type=mut_box[mutboxindex+1]
        cig_len=cigarbox[mutboxindex]
        if (mut_pos < seqpos):
            ''' Haven't reached target, accumulate reference offsets from omissions ...'''
            if (mut_type == "D") or (mut_type == "N"):
                ref_offset+=cig_len
                ''' ... and also negative reference offsets from inserts'''
            elif mut_type=="I":
                ref_offset-=cig_len
            '''skip to next '''
        else:
            if mut_pos >= end:
                length = end-seqpos
                newcigarbox.insert(len(newcigarbox)-1,length)
                newcigarbox.insert(len(newcigarbox)-1,mut_type)
                ''' Reached end: terminate the while loop '''
                mutboxindex=len(mut_box)
            elif (mut_type == "D") or (mut_type == "N"):
                '''Where CIGAR includes omissions: omit when at start and accumulate in offset '''
                if not added:
                    ref_offset+=cig_len
                    newcigarbox.insert(len(newcigarbox)-1,cig_len)
                    newcigarbox.insert(len(newcigarbox)-1,mut_type)
                else:
                    ''' Append omissions at full length. '''
                    newcigarbox.insert(len(newcigarbox)-1,cig_len)
                    newcigarbox.insert(len(newcigarbox)-1,mut_type)
            else:
                length=mut_pos-seqpos
                if length > 0:
                    ''' Append CIGAR to the length of the (non-omission) feature'''
                    newcigarbox.insert(len(newcigarbox)-1,length)
                    newcigarbox.insert(len(newcigarbox)-1,mut_type)
                    if (not added):
                        added = True
                        if (mut_type == "I"):
                            ''' accumulate negative offset'''
                            ioffset=cig_len-length
                            ref_offset-=ioffset
                elif (mut_type == "I"):
                    ''' Another special case where length=0 at end of I-run, must kill the offset as if feature skipped
                        at mut_pos < seqpos. Omit this and the Is get handled incorrectly'''
                    ref_offset-=cig_len
        if mut_pos > seqpos:
            seqpos=mut_pos
        mutboxindex+=2
        # end of while mutboxindex<len(mut_box)-2 loop
    return ref_offset,newcigarbox
# end of trimcigarbox(cigarbox,mut_box,start,length)


def is_mut_cigar(cigar):
    ''' Mutation types declared in CIGAR: Deletion, Insertion, X - variant '''
    '''          N is a CIGAR item, but ignored as a mutation type'''
    ''' mut_types = "DIX" '''
    return any(i in cigar for i in RG_globals.mut_types)
# end of is_mut_cigar(cigar)
# ===============================================================
# End of cigar-manipulation functions definition
# ===============================================================

# ===============================================================
# Mutfreq and mutlabels functions 
# ===============================================================

def mutfreqs_extend(mutlabels):
    # App.vue should deliver a correctly-sized list. This is to catch any other situations
    #array_extend=0
    for i in range(0, len(mutlabels)):
        # To prevent an out-of-range, append last frequency to end
        if i > len(RG_globals.mutfreqs)-1:
            RG_globals.mutfreqs.append(RG_globals.mutfreqs[-1])
            #array_extend+=1
    #if array_extend > 0:  # Was previously in RG_main; update journal will not work here
    #    update_journal("mutation frequencies array automatically extended by %s to: %s"%(array_extend,RG_globals.mutfreqs))
    return

def get_addmut_labels():
    addmut_labels=[]
    if not RG_globals.bio_parameters["target_build_variant"]["AddVars"]==[]:
        for item in RG_globals.bio_parameters["target_build_variant"]["AddVars"]:
            addmut_labels.append(item["hapname"])
    return addmut_labels

def get_mutfreq_from_label(inlabel):
    # Was used by several, now only by get_mutrecords to avoid reading unselected input files 
    labelcount=0
    for label in RG_globals.mutlabels:
        if label == inlabel:
            break
        else:
            labelcount+=1
    return RG_globals.mutfreqs[labelcount]

