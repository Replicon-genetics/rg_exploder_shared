# Progver="RG_exploder_io.py"
# ProgverDate="18-Feb-2022"
# 

'''
This io module is replaced in the pyodide version with something completely different.
An original version provided, under contract to Replicon Genetics, by Raven Bioseciences

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
'''

from pathlib import Path

def is_file(filename):
    return Path(filename).is_file()

def open_read(filename):
    return open(filename, 'r')

def open_write(filename, buffer_size = 4194304):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    return open(filename, 'w')

# This is currently used solely in the Python GUI module - not in any
# exploder modules used in Pyodide
def list_files(directory,files):
    return(list(Path(directory).glob(files)))
    
    
###!!!### Does not work in webapp repository as no support for gz in web context yet###!!!###
### .............. from here .................... ###!!!###
def read_single_record_input_gz_support(infile,seqform): # Renamed from read_single_record_input to allow web app to use alternative
    # Conditionally read gzip file'''           # call this if want to re-instate in python version 
    # Read in a single record from a genbank or Embl format file
    # Supposed to throw an error if > 1 record found
    
    # Add these imports 
    import binascii  # Used in is_gz_file                            
    import gzip  # Used in read_single_record_input                  
    from functools import partial # Used in read_single_record_input ###!!!### No support for gz in web context ###!!!###
    from Bio import SeqIO #BioPython
    import RG_exploder_process as RG_process

    success=False
    
    
    def is_gz_file(filepath):
    # binary magic detects if the target is a gzip file or not
        with open(filepath, 'rb') as test_f:
            val=binascii.hexlify(test_f.read(2)) == b'1f8b'
        return val
    # end of is_gz_file(filepath)
    is_gz = is_gz_file(infile)
    
    if not is_gz:
        _open = open
    elif is_gz:
        _open = partial(gzip.open, mode='rt')
    else:
        program_exit('Unknown file encoding at %s'%infile)
    with _open(infile) as f:
        seq_record = SeqIO.read(f,seqform)
        seq_record,success = RG_process.add_extra_objects(seq_record)
        #print("infile %s, seq_record.offset %s"%(infile,seq_record.offset))
    return seq_record,success
 #end of read_single_record_input_gz_support(infile,seqform)
### .............. to here .................... ###!!!###
