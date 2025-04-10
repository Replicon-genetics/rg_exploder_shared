#!/usr/local/bin/python3
Progver="RG_exploder_Biopython_fix_metro.py"
ProgverDate="21-Feb-2023"
'''

Sole purpose of this module is to support switching between versions of Biopython to support:
a) Deprecation of 'generic_dna' in MutableSeq and SeqRecord
b) sequence.reverse_complement() change to sequence.reverse_complement(inplace=True)

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
