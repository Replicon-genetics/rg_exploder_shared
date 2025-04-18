# rg_exploder_shared

Python code for generating synthetic sequence data; synthetic DNASEQ or RNASEQ reads

This repository at https://github.com/snowlizardz/rg_exploder_shared/ holds Python code, data and metadata for the fragmentation of DNA sequences intended to emulate NGS style sequencing reads. If you want to know why this is useful, take a look at https://repliconevaluation.com/about/ ; repliconevaluation.com redirects to replicongenomics.com, at least until September 2025.

This code here is in the Public Domain as free software, as defined by https://www.fsf.org/ , specifically the  <a href="https://www.gnu.org/licenses/agpl-3.0.en.html">AGPL-3.0 license</a>. This repository is a tidied-up version of the development repository at https://github.com/snowlizardz/rg_exploder (currently private).

<b>Getting started</b><br>
In <a href="https://github.com/Replicon-genetics/rg_exploder_shared/tree/main/exploder_python">/exploder_python</a>:<br>
a) Execute the Python module <a href="https://github.com/Replicon-genetics/rg_exploder_shared/blob/main/exploder_python/RG_exploder_globals_make.py">RG_exploder_globals_make.py</a> to set the correct values in config.json <br>
b) Execute the Python <a href="https://docs.python.org/3/library/tkinter.html">Tkinter</a> GUI <a href="https://github.com/Replicon-genetics/rg_exploder_shared/blob/main/exploder_python/RG_exploder_gui.py">RG_exploder_gui.py</a>

<b>Dependencies</b><br>
The Python modules require <a href="https://biopython.org/">Biopython</a>, with <a href="https://pillow.readthedocs.io/en/stable/">Pillow(PIL)</a> to support the <a href="https://docs.python.org/3/library/tkinter.html">Tkinter</a> GUI. There may be an X11 dependency on some platforms.

<b>Maintaining and updating genomic data</b><br>
<a href="https://github.com/Replicon-genetics/rg_exploder_shared/tree/main/data_sources">/data_sources</a> holds genomic data downloaded from Ensembl, which has been processed as described:<br>
<a href="https://github.com/Replicon-genetics/rg_exploder_shared/tree/main/documents/">/documents</a> includes a maintenance guide for pre-processing new genomic data using the maintenance scripts in folders <a href="https://github.com/Replicon-genetics/rg_exploder_shared/tree/main/helper_python">/helper_python</a> & <a href="https://github.com/Replicon-genetics/rg_exploder_shared/tree/main/helper_scripts">/helper_scripts</a>. <br>
In <a href="https://github.com/Replicon-genetics/rg_exploder_shared/tree/main/documents/presentations">presentations<a>, there's a document explaining concepts and implementation.

<b>Alternative GUI</b><br>
The same Python code, data etc are behind a Vue.js implementation at https://repliconevaluation.wordpress.com/replicon-genetics (second version released May 2021). The full source for building the Vue.js implementation is held at https://github.com/Replicon-genetics/rg_exploder , currently private. A subset of critical definition files is present in this current repository: folders /pyodide /vue_source and /webdist 

Code and documentation was developed between September 2018 to March 2025 by Cary O'Donnell, originally for Replicon Genetics, a company set up by Dr Gillian Ellison and Jane Theaker in 2018, but de-registered in 2023; IP is claimed by those named. No AI-generation tools were used at any point.

Advice on improving access, offers on collaboration, or other feedback is welcome; please email syrgenreads@gmail.com

Cary O'Donnell 18th April 2025
