self.languagePluginUrl = '/'
importScripts('./pyodide.js')

self.postLines = function(file, lines) {
  self.postMessage({file: file, eof: false, lines: lines});
}

self.postLine = function(file, line) {
  self.postMessage({file: file, eof: false, line: line});
}

self.postEndOfFile = function(file) {
  self.postMessage({file: file, eof: true});
}

var onmessage = function(e) { // eslint-disable-line no-unused-vars
  languagePluginLoader.then(() => {
    console.log(e.data)

    self.pyodide.loadPackage(['biopython','rg_exploder']).then(() => {
      // # Order of these parameters must be coordinated with RG_exploder_webmain.py
      self.pyodide.runPython(
        'import RG_exploder_webmain\n' +
        'import js\n' +
        'RG_exploder_webmain.run('
          + e.data.targetLocus + ', '
          + e.data.transcriptName + ', '
          + (e.data.transcriptId ? e.data.transcriptId : 'None') + ', '
          + (e.data.CDSonly ? 'True' : 'False') + ', '
          + e.data.fragmentLength + ', '

          + e.data.coverageDepth + ', '
          + e.data.exomeExtension + ', '
          + e.data.gaussMean + ', '
          + e.data.gaussSD + ', '
          + e.data.qualityMin + ', '

          + e.data.qualityMax + ', '
          + e.data.mutFreqs + ', '
          + e.data.mutLabels + ', '
          + (e.data.outputSam ? 'True' : 'False') + ', '
          + (e.data.outputFastq ? 'True' : 'False') + ', '

          + (e.data.outputFasta ? 'True' : 'False') + ', '
          + (e.data.refSequenceInFasta ? 'True' : 'False') + ', '
          + (e.data.variantSeqInFasta ? 'True' : 'False') + ', '
          + (e.data.sourceFeatureTables ? 'True' : 'False') + ', '
          + (e.data.eachPossibleRead ? 'True' : 'False') + ', '

          + (e.data.variantReadsOnly ? 'True' : 'False') + ', '
          + (e.data.pairedReads ? 'True' : 'False') + ', '
          + (e.data.duplexReads ? 'True' : 'False') + ', '
          + (e.data.annotateSourcePositions ? 'True' : 'False') + ', '
          + (e.data.plusAbsolutePositions ? 'True' : 'False') + ', '

          + (e.data.cigarAnnotated ? 'True' : 'False') + ', '
          + (e.data.substitutionLowerCase ? 'True' : 'False') + ', '
          + (e.data.journalSubs ? 'True' : 'False') + ', '
          + (e.data.FlipStrand ? 'True' : 'False') + ', '
          + e.data.tbv_AddVars + ', '

          + 'js)\n'
      );
    });
  });
}

//+ e.data.tbv_hap_name + ', '
//+ e.data.tbv_var_name + ', '
//+ e.data.tbv_local_Begin + ', '

//+ e.data.tbv_local_End + ', '
//+ e.data.tbv_ref_subseq + ', '
//+ e.data.tbv_var_subseq + ', '
//+ e.data.tbv_abs_Begin + ', '
//+ e.data.tbv_abs_End + ', '
// tbv_is_save_var:this.tbv_is_save_var,
// tbv_hapname:this.tbv_hap_name.value,
// tbv_varname:this.tbv_var_name.value,
// tbv_local_begin:this.tbv_local_Begin,
// tbv_local_end:this.tbv_local_End,
// tbv_ref_seq:this.tbv_ref_subseq.value,
// tbv_var_seq:this.tbv_var_subseq.value,
//tbv_abs_Begin:this.tbv_abs_Begin.value,
//tbv_abs_End:this.tbv_abs_End.value
