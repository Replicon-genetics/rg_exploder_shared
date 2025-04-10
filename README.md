# rg_exploder_shared

Python code for generating synthetic sequence data; synthetic DNASEQ or RNASEQ reads

This repository at https://github.com/snowlizardz/rg_exploder_shared/ holds Python code, data and metadata for the fragmentation of DNA sequences intended to emulate NGS style sequencing reads. It is a tidied-up version of the development repository at https://github.com/snowlizardz/rg_exploder

In /exploder_python there's a Tkinter GUI to generate the reads; /data_sources holds the pre-processed genomic data; /documents include a maintenance guide for pre-processing new genomic data using the maintenance scripts in folders /helper_python & /helper_scripts. There's also a presentation document explaining concept and implementation: simulation5odp.odp

The main objective of the repository is for placement in the Public Domain as free software, as defined by https://www.fsf.org/, the specific licence granted here is the <a href="https://www.gnu.org/licenses/agpl-3.0.en.html">AGPL-3.0 license</a> as indicated in this repository.

The same Python code, data etc are behind a Vue.js implementation at https://repliconevaluation.wordpress.com/replicon-genetics (second version released May 2021). The full source for building the Vue.js implementation is held at https://github.com/Replicon-genetics/rg_exploder. A subset of critical definition files is present in this current repository: folders /pyodide /vue_source and /webdist 

Code and documentation was developed between September 2018 to March 2025 by Cary O'Donnell, originally for Replicon Genetics, a company set up by Dr Gillian Ellison and Jane Theaker in 2018, but de-registered in 2023. IP is claimed by them and Cary O'Donnell. No AI-generation tools were used at any point.

Advice on improving access or other feedback is welcome; please email syrgenreads@gmail.com

Cary O'Donnell 10th April 2025
