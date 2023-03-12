#!/usr/local/bin/python3
Progver="RG_exploder_Biopython_fix_metro.py"
ProgverDate="21-Feb-2023"
'''
© author: Cary O'Donnell for Replicon Genetics 2023
Sole purpose of this module is to support switching between versions of Biopython to support:

a) Deprecation of 'generic_dna' in MutableSeq and SeqRecord
b) sequence.reverse_complement() change to sequence.reverse_complement(inplace=True)

© author: Cary O'Donnell for Replicon Genetics 2023
'''

'''
BioPython imports
'''

from Bio.Seq import MutableSeq
#from Bio.Alphabet import generic_dna #generic_dna # Deprecated from Biopython 1.78 (September 2020)
from Bio.SeqRecord import SeqRecord # Used in writing protein sequence
from Bio.Seq import Seq

def fix_SeqRecord(inSeq):
    #OutSeqr=SeqRecord(Seq(str(inSeq),generic_dna))#generic_dna # Deprecated from Biopython 1.78 (September 2020)
    OutSeqr=SeqRecord(Seq(str(inSeq)))
    return OutSeqr

def fix_MutableSeq(inSeq):
    #OutSeq=MutableSeq(str(inSeq),generic_dna) #generic_dna # Deprecated from Biopython 1.78 (September 2020)
    OutSeq=MutableSeq(str(inSeq))
    return OutSeq

def fix_reverse_complement(inSeq):
    #inSeq.reverse_complement()# Warning from Python 3.11.0 that this needs replacing for near-future releases 
    inSeq.reverse_complement(inplace=True)
    #return inSeq
