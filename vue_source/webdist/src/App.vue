// This is a GUI for Synthetic Reads Generator developed for Replicon Genetics 2021-24  © Replicon Genetics 
<template>
<div>
  <Header :title="jsonConfig.stringconstants.title_label" />
  <RunButton :running="running" @click="run_method" />

  <div class="details">
    <h1 style="font-size: 1.25rem">
      {{ jsonConfig.stringconstants.CustomerIDText }}
    </h1>
  </div>

  <div class="details">
    <h1 style="font-size: 1.25rem">
      {{ jsonConfig.stringconstants.DatasetIDText }}
    </h1>
  </div>

  <div class="details">
    <a
      v-if="jsonConfig.stringconstants.help_url"
      :href="jsonConfig.stringconstants.help_url"
      target="_blank"
    >
      {{ jsonConfig.stringconstants.help_label }}
    </a>
    <a
      v-if="jsonConfig.stringconstants.about_url"
      :href="jsonConfig.stringconstants.about_url"
      target="_blank"
    >
      {{ jsonConfig.stringconstants.about_label }}
    </a>
    <a
      v-if="jsonConfig.stringconstants.more_url"
      :href="jsonConfig.stringconstants.more_url"
      target="_blank"
    >
      {{ jsonConfig.stringconstants.more_label }}
    </a>
  </div>

  <div class="container">
    <div class="container__child">
      <h3>{{ jsonConfig.stringconstants.reference_gene }}</h3>
      <BaseSelect
        v-if="selectedGene"
        v-model="selectedGene"
        :link="ensembl_locusURL" 
        :label="ensembl_locusURLtxt"
        :options="Object.keys(geneList)"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="FlipStrand"
        v-model="FlipStrand.value"
        :label="FlipStrand.label"
        :disabled="running"
      />
      <h3>{{ jsonConfig.stringconstants.reads_list_label }}</h3>
      <BaseNumberInput
        v-if="fragmentLength"
        v-model="fragmentLength.value"
        :label="fragmentLength.label"
        :min="fragmentLength.min.toString()"
        :max="fragmentLength.max.toString()"
        :disabled="running"
      />
      <BaseNumberInput
        v-if="coverageDepth"
        v-model="coverageDepth.value"
        :label="coverageDepth.label"
        :min="coverageDepth.min.toString()"
        :max="coverageDepth.max.toString()"
        :disabled="running"
      />
       <BaseSelect
        v-if="selectedReadType"
        v-model="selectedReadType"
        :label="jsonConfig.stringconstants.reads_type_label"
        :options="Object(readList)"
        :disabled="running"
      />
      <BaseNumberInput
        v-if="pairedReads.value"
        v-model="gaussMean.value"
        :label="gaussMean.label"
        :min="gaussMean.min.toString()"
        :max="gaussMean.max.toString()"
        :disabled="running"
      />
      <BaseNumberInput
        v-if="pairedReads.value"
        v-model="gaussSD.value"
        :label="gaussSD.label"
        :min="gaussSD.min.toString()"
        :max="gaussSD.max.toString()"
        :disabled="running"
      />
      <img
        :src="require('./assets/replicon_logo.svg')"
        alt="The Replicon Genetics logo"
        style="margin: 2em;"
      />
    </div>
    <div v-if="mutfreqs" class="container__child">
      <h3>{{ jsonConfig.stringconstants.reference_haplotype }}</h3>
      <BaseSelect
        v-if="selectedTranscript"
        v-model="selectedTranscript"
        :link="ensembl_transcriptURL"
        :label="ensembl_transcriptURLtxt"
        :options="transcriptOptions"
        :disabled="running"
      />
      <BaseNumberInput
        v-if="exomeExtension"
        v-model="exomeExtension.value"
        :label="exomeExtension.label"
        :min="exomeExtension.min.toString()"
        :max="exomeExtension.max.toString()"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="CDSonly"
        v-model="CDSonly.value"
        :label="CDSonly.label"
        :disabled="running"
      />
      <div class="details">
         <h4>{{ jsonConfig.stringconstants.variants_label+"s" }}</h4>
      </div>
      <div class="container__child-subtitle">
        <h4>{{ jsonConfig.bio_parameters.target_build_variant.hap_name.label }}</h4>
        <h4>{{ jsonConfig.stringconstants.frequency_label }}</h4>
      </div>
      <BaseSlider
        v-for="(mut, index) of mutfreqs"
        :key="index"
        v-model="mutfreqs[index]"
        :label="selectedGene + ' ' + mutfreqLabels[index]"
        :selected-gene="selectedGene"
        :disabled="running"
      />
    </div>
    <div class="container__child">
      <h3>{{ jsonConfig.stringconstants.options_label }}</h3>
      <BaseCheckbox
        v-if="outputFasta"
        v-model="outputFasta.value"
        :label="outputFasta.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="eachPossibleRead"
        v-model="eachPossibleRead.value"
        :label="eachPossibleRead.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="variantReadsOnly"
        v-model="variantReadsOnly.value"
        :label="variantReadsOnly.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="annotateSourcePositions"
        v-model="annotateSourcePositions.value"
        :label="annotateSourcePositions.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="plusAbsolutePositions"
        v-model="plusAbsolutePositions.value"
        :label="plusAbsolutePositions.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="cigarAnnotated"
        v-model="cigarAnnotated.value"
        :label="cigarAnnotated.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="substitutionLowerCase"
        v-model="substitutionLowerCase.value"
        :label="substitutionLowerCase.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="journalSubs"
        v-model="journalSubs.value"
        :label="journalSubs.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="outputFastq"
        v-model="outputFastq.value"
        :label="outputFastq.label"
        :disabled="running"
      />
      <BaseNumberInput
        v-if="qualityMin"
        v-model="qualityMin.value"
        :label="qualityMin.label"
        :min="qualityMin.min.toString()"
        :max="qualityMin.max.toString()"
        :disabled="running"
      />
      <BaseNumberInput
        v-if="qualityMax"
        v-model="qualityMax.value"
        :label="qualityMax.label"
        :min="qualityMax.min.toString()"
        :max="qualityMax.max.toString()"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="outputSam"
        v-model="outputSam.value"
        :label="outputSam.label"
        :disabled="running"
      />
       <BaseCheckbox
        v-if="refSequenceInFasta"
        v-model="refSequenceInFasta.value"
        :label="refSequenceInFasta.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="variantSeqInFasta"
        v-model="variantSeqInFasta.value"
        :label="variantSeqInFasta.label"
        :disabled="running"
      />
      <BaseCheckbox
        v-if="sourceFeatureTables"
        v-model="sourceFeatureTables.value"
        :label="sourceFeatureTables.label"
        :disabled="running"
      />
      <BaseCheckbox v-model="showMessagePanel" label="Show message panel" />
    </div>
  </div>
  <div class="container">
    <div class="container__child3">
      <div class="details">
    <h1 style="font-size: 1.25rem">
      {{ "Create a New " +this.selectedGene +" "+jsonConfig.stringconstants.variants_label }}
    </h1>
    </div>
      <div class="container__child-subtitle">
        <h3>{{ " Define Location"}}</h3>
      </div>
      <div class="row">
        <div class="column">
          <h4>{{ "Template Coordinates"}}</h4>
        </div>
        <div class="column">
          <h4>{{ "Extension"}}</h4>
        </div>
        <div class="column">
          <h4>{{ "Genome Coordinates"}}</h4>
        </div>
      </div>
      <div class="row">
        <div class="column">
          <BaseNumberInputWide
            v-if="tbv_trans_Begin"
            v-model="tbv_trans_Begin.value"
            :label="build_label+'_Begin'"
            :min="tbv_trans_Begin.min.toString()"
            :max="tbv_trans_End.value.toString()"
            :disabled="running"
          />
        </div>
        <div class="column">
          <BaseNumberInputExt
            v-if="tbv_trans_Begin_ext"
            v-model="tbv_trans_Begin_ext.value"
            :label="'Begin'"
            :min="tbv_trans_Begin_ext.min.toString()"
            :max="tbv_trans_Begin_ext.max.toString()"
            :disabled="running"
          />
        </div>
        <div class="column">
          <BaseStringShow
            v-if="tbv_abs_Begin"
            v-model="tbv_abs_Begin.txt"
            :label="tbv_abs_Begin.label"
            :disabled="true"
          />
        </div>
      </div>
      <div class="row">
        <div class="column">
          <BaseNumberInputWide
            v-if="tbv_trans_End"
            v-model="tbv_trans_End.value"
            :label="build_label+'_End'"
            :min="tbv_trans_Begin.value.toString()"
            :max="tbv_trans_End.max.toString()"
            :disabled="running"
          />
        </div>
        <div class="column">
          <BaseNumberInputExt
            v-if="tbv_trans_End_ext"
            v-model="tbv_trans_End_ext.value"
            :label="'End'"
            :min="tbv_trans_End_ext.min.toString()"
            :max="tbv_trans_End_ext.max.toString()"
            :disabled="running"
          />
        </div>
        <div class="column">
          <BaseStringShow
            v-if="tbv_abs_End"
            v-model="tbv_abs_End.txt"
            :label="tbv_abs_End.label"
            :disabled="true"
          />
        </div>
      </div>  
      <div class="container__child-subtitle">
        <h3>{{ "Define "+tbv_var_name.label.split(' ')[0] }}</h3>
      </div>
        <BaseStringShowWide
        v-if="tbv_ref_subseq"
        v-model="tbv_ref_subseq.viewstring"
        :label="tbv_ref_subseq.label"
        :disabled="true"
      />
        <BaseStringInputDNA
        v-if="tbv_var_subseq"
        v-model="tbv_var_subseq.value"
        :label="tbv_var_subseq.label"
        :disabled="tbv_ref_subseq_notretrieved"
      />
        <BaseStringInput
        v-if="tbv_var_name"
        v-model="tbv_var_name.value"
        :label="tbv_var_name.label"
        :disabled="tbv_ref_subseq_notretrieved"
      />
      <BaseStringInput
        v-if="tbv_hap_name"
        v-model="tbv_hap_name.value"
        :label="tbv_hap_name.label"
        :disabled="tbv_ref_subseq_notretrieved"
      />
      <div class="row">
        <div class="column">
          <EditVarButton :running="running"
                          :label="'Edit'" @click="edit_variant_method" />
        </div>
        <div class="column">
          <SaveButton :running="running"
                      :label="'Save'" @click="save_method" />
        </div>
      </div>
    </div>
  </div>

  <Copyright :copytitle="jsonConfig.stringconstants.CopyrightText" /> 
  <MessagePanel v-if="showMessagePanel" :journal="journal" :readme="readme" />
  </div>
</template>

<script>
// This is a GUI for Synthetic Reads Generator developed for Replicon Genetics 2021-24  © Replicon Genetics 
import Header from '@/components/Header.vue'
import Copyright from '@/components/Copyright.vue'
import MessagePanel from '@/components/MessagePanel.vue'
import BaseNumberInput from '@/components/BaseNumberInput.vue'
import BaseNumberInputExt from '@/components/BaseNumberInputExt.vue'
import BaseNumberInputWide from '@/components/BaseNumberInputWide.vue'
import BaseStringInput from '@/components/BaseStringInput.vue'
import BaseStringInputDNA from '@/components/BaseStringInputDNA.vue'
import BaseStringShow from '@/components/BaseStringShow.vue'
import BaseStringShowWide from '@/components/BaseStringShowWide.vue'
import BaseCheckbox from '@/components/BaseCheckbox.vue'
import RunButton from '@/components/RunButton.vue'
import EditVarButton from '@/components/BaseButton.vue'
import SaveButton from '@/components/BaseButton.vue'
import BaseSelect from '@/components/BaseSelect.vue'
import BaseSlider from '@/components/BaseSlider.vue'

import streamSaver from 'streamsaver'
import config from '../public/input/config.json' // config file

var pyodideWorker = new Worker('./webworker.js')
var encode = TextEncoder.prototype.encode.bind(new TextEncoder())

streamSaver.mitm = './mitm.html'

// The filenames in the /input/ folder, useful for the [mutfreqLabels] & [mutfreq] arrays
// require.context must have literals, so cannot limit this full listing to eg:subdirectory this.selectedGene
const inputGeneKeys = require.context('../public/input/', true, /\.gb$/).keys()

export default {
  name: 'App',
  components: {
    Header,
    Copyright,
    MessagePanel,
    BaseNumberInput,
    BaseNumberInputExt,
    BaseNumberInputWide,
    BaseStringInput,
    BaseStringInputDNA,
    BaseStringShow,
    BaseStringShowWide,
    BaseCheckbox,
    RunButton,
    EditVarButton,
    SaveButton,
    BaseSelect,
    BaseSlider
  },
  data() {
    return {
      running: false,
      inMemFiles: new Map(),
      streamWriters: new Map(),
      journal: '',
      readme: '',

      showMessagePanel: true,

      jsonConfig: config,
      /// Parameters that are used in the pyodide instance
      // Gene list + currently selected gene
      geneList: config?.Reference_sequences ?? {},
      selectedGene: config?.bio_parameters?.target_locus.value ?? 'None',
      
      readList: config.ReadsList,
      selectedReadType:config?.bio_parameters?.is_frg_paired_end.label ?? 'None',

      selectedTranscript:
        config?.bio_parameters?.target_transcript_name.value ?? config.stringconstants.empty_transcript_name,

      // Parameter booleans
      outputFasta: config.bio_parameters.is_fasta_out,
      outputFastq: config.bio_parameters.is_fastq_out,
      outputSam: config.bio_parameters.is_sam_out,
      journalSubs: config.bio_parameters.is_journal_subs,
      refSequenceInFasta: config.bio_parameters.is_write_ref_fasta,
      variantSeqInFasta: config.bio_parameters.is_mut_out,
      sourceFeatureTables: config.bio_parameters.is_write_ref_ingb,
      eachPossibleRead: config.bio_parameters.is_onefrag_out,
      variantReadsOnly: config.bio_parameters.is_muts_only,
      pairedReads: config.bio_parameters.is_frg_paired_end,
      duplexReads: config.bio_parameters.is_duplex,
      annotateSourcePositions: config.bio_parameters.is_frg_label,
      plusAbsolutePositions: config.bio_parameters.is_use_absolute,
      cigarAnnotated: config.bio_parameters.is_fastacigar_out,
      substitutionLowerCase: config.bio_parameters.is_vars_to_lower,
      CDSonly:config.bio_parameters.is_CDS,
      FlipStrand:config.bio_parameters.is_flip_strand,

      // Parameter numbers
      exomeExtension: config?.bio_parameters?.Exome_extend,
      fragmentLength: config?.bio_parameters?.Fraglen,
      coverageDepth: config?.bio_parameters?.Fragdepth,
      qualityMin: config?.bio_parameters?.Qualmin,
      qualityMax: config?.bio_parameters?.Qualmax,
      gaussMean: config?.bio_parameters?.gauss_mean,
      gaussSD: config?.bio_parameters?.gauss_SD,

      // Initialising relative frequency labels and number arrays
      mutfreqLabels: [],
      mutfreqs: [],
      // Can be changed here by Reference_sequences.selectedGene... data - initialised from json
      tbv_GRChver_txt: config?.bio_parameters?.target_build_variant?.GRChver_txt,
      tbv_is_save_var: config?.bio_parameters?.target_build_variant?.is_save_var,
      tbv_headclip:config?.bio_parameters?.target_build_variant?.headclip, 
      //tbv_headclip:0,
      //tbv_tailclip:config?.bio_parameters?.target_build_variant?.tailclip,
      mrnapos_lookup:config?.bio_parameters?.target_build_variant?.mrnapos_lookup,
      tbv_transcript_view:config?.bio_parameters?.target_build_variant?.transcript_view,
      tbv_abs_offset:config?.bio_parameters?.target_build_variant?.abs_offset,
      tbv_ref_strand:config?.bio_parameters?.target_build_variant?.ref_strand,
      tbv_local_Begin: config.bio_parameters.target_build_variant.local_begin,
      tbv_local_End: config.bio_parameters.target_build_variant.local_end,

      // Initialising a variable extracted from Reference_sequences.selectedGene...
      joinlist:config?.bio_parameters.target_build_variant?.joinlist,

      // Changed here dynamically 
      tbv_abs_Begin:config?.bio_parameters.target_build_variant?.abs_Begin,
      tbv_abs_End:config?.bio_parameters.target_build_variant?.abs_End,

      // Changed here dynamically - via user changeable values, save to main via bio_parameters.target_build_variant
      tbv_trans_Begin: config?.bio_parameters.target_build_variant?.trans_Begin,
      tbv_trans_Begin_ext: config?.bio_parameters.target_build_variant?.trans_Begin_ext,
      tbv_trans_End: config?.bio_parameters.target_build_variant?.trans_End,
      tbv_trans_End_ext: config?.bio_parameters.target_build_variant?.trans_End_ext,

      // Calculated value for reference sequence
      tbv_ref_subseq:config?.bio_parameters.target_build_variant?.ref_subseq,
      tbv_ref_subseq_notretrieved: true, // // Implement later as tbv_ref_subseq.noseq ?
      // Initialising as a test value
      REFSEQ:config?.bio_parameters.target_build_variant?.refseq,
      REFSEQ_len:0,
      DNA_COMPLEMENT: {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'},
      tbv_Chrom:"X",

      // User-entered
      tbv_var_subseq:config?.bio_parameters.target_build_variant?.var_subseq,
      tbv_hap_name:config?.bio_parameters.target_build_variant?.hap_name,
      tbv_hap_name_def:config?.bio_parameters.target_build_variant.hap_name?.value, // Saving initial value as a default
      tbv_var_name:config?.bio_parameters.target_build_variant?.var_name,

      //And the object of the build exercise
      tbv_AddVars:config?.bio_parameters.target_build_variant?.AddVars
    }
  },
  watch: {
    selectedGene: {
      /// Finds all the files for the currently selected gene.
      /// Generates two arrays: [mutfreqs] and [mutfreqLabels], and offers these to $data.
      handler() {
        const keys = inputGeneKeys.filter(item => {
          ///if (`${item}`.includes('ref')) {
          if (`${item}`.includes(this.jsonConfig.IOconstants.Ref_file_name)) { 
            //console.log(`selectedGene false item ${item}`)
            return false
          }
          return `${item}`.startsWith(`./${this.selectedGene}/`)
        })
        ///this.mutfreqs = keys.map((e, idx) => (idx <= 1 ? 50 : 30))
        this.mutfreqs = keys.map((e, idx) => (idx < 1 ? 0 : 50))
        ///this.mutfreqs = keys.map((e, idx) => (this.jsonConfig.admin_constants.mutfreqs[idx]))
        this.mutfreqLabels = keys.map(e =>
          e.substring(e.lastIndexOf('_') + 1, e.lastIndexOf('.'))
        )
        this.add_to_mutfreqs()
      },
      immediate: true,
      deep: true,
    },
    selectedReadType: {
      handler(){
        if (this.selectedReadType.includes(config.bio_parameters.is_frg_paired_end.label))
        {this.pairedReads.value = true
         this.duplexReads.value = false
        }
        else if (this.selectedReadType.includes(config.bio_parameters.is_duplex.label))
        {this.pairedReads.value = null
         this.duplexReads.value = true}
        else if (this.selectedReadType.includes(config.bio_parameters.is_simplex.label))
        {this.pairedReads.value = null
         this.duplexReads.value = false}
      },
      immediate: true,
      deep: true,
    },
    transcriptOptions: {
      handler() {
        if (!this.transcriptOptions.includes(this.selectedTranscript)) {
          this.selectedTranscript = this.jsonConfig.stringconstants.empty_transcript_name
        }
        this.join_method()
        //console.log(`TO: joinlist is ${this.joinlist}`)
      },
      immediate: true,
      deep: true,
    },
    selectedTranscript: {
      handler() {
        this.join_method()
        //console.log(` ST: joinlist is ${this.joinlist}`)
      },
      immediate: true,
      deep: true,
    },
    CDSonly:{
      handler() {
        this.join_method()
        //console.log(` CDO: joinlist is ${this.joinlist}`)
      },
      immediate: true,
      deep: true,
    },
    tbv_trans_Begin:{
      handler() {
        this.get_abs_begin()
        //console.log(` tbv_trans_Begin: tbv_abs_Begin.value is ${this.tbv_abs_Begin.value}`)
      },
      immediate: true,
      deep: true,
    },
    tbv_trans_Begin_ext:{
      handler() {
        this.get_abs_begin()
        //console.log(` tbv_trans_Begin_ext: tbv_abs_Begin_ext.value is ${this.tbv_abs_Begin_ext.value}`)
      },
      immediate: true,
      deep: true,
    },
      tbv_trans_End:{
      handler() {
        this.get_abs_end()
        //console.log(` tbv_trans_End: tbv_abs_End.value is ${this.tbv_abs_End.value}`)
      },
      immediate: true,
      deep: true,
    },
    tbv_trans_End_ext:{
      handler() {
        this.get_abs_end()
        //console.log(` tbv_trans_End_ext: tbv_abs_End_ext.value is ${this.tbv_abs_End_ext.value}`)
      },
      immediate: true,
      deep: true,
    },
    joinlist:{
      // Only required values set here are: 
      // this.mrnapos_lookup, this.tbv_ref_strand, this.tbv_abs_offset, this.tbv_trans_Begin.max,this.tbv_trans_End.max

      // Starts, end, exon_length, feature_titles only needed for building transcript_view, to return for journaling 
      handler()
      {
      //console.log(`JL: joinlist change`)
      //console.log(`joinlist: ${this.joinlist}`)

      // Now for the meaty stuff:
      //this.is_join_complement=this.jsonConfig.Reference_sequences[this.selectedGene]["is_join_complement"]
      const Region = this.jsonConfig.Reference_sequences[this.selectedGene]["Region"]
      const abs_start=parseInt(Region.split(":")[2])
      const abs_end=parseInt(Region.split(":")[3])
      this.tbv_ref_strand=parseInt(Region.split(":")[4])
      this.tbv_GRChver_txt=(Region.split(":")[0])
      this.tbv_Chrom=(Region.split(":")[1])
      // Special catch for non-standard chromosome identifiers such as PTEN_a
      if (this.tbv_Chrom.length > 2)
        {this.tbv_Chrom=this.tbv_Chrom.substring(0,3);}

      const maxreflen=Math.abs(abs_end-abs_start)+1

      this.REFSEQ_len=maxreflen

      var offset=abs_end+1

      if (this.tbv_ref_strand == 1)
        {offset=abs_start-1}

      this.tbv_abs_offset=offset*this.tbv_ref_strand
      
      const Locus_Range = this.jsonConfig.Reference_sequences[this.selectedGene]["Locus_range"]
      const locus_begin = parseInt(Locus_Range.split(":")[0])
      const locus_end = parseInt(Locus_Range.split(":")[1])
      const locus_length=locus_end-locus_begin+1

      //console.log(` locus_begin: ${locus_begin}; locus_end: ${locus_end}`) // Validation check
      //console.log(` locus_length: ${locus_length}`) // Validation check


      // Need to store an offset for case of generating variants from *clipped* genomic sequence, so the offset is the Headclip length
      this.tbv_headclip=locus_begin-1
      this.tbv_tailclip=maxreflen-locus_end

      var exon_total=0
      var intron_total=0
      var feature_titles=["Pre-Locus"]
      var starts=[]
      var ends=[]
      var exon_text=""
      var intron_text=""

      //var abs_starts=[]
      //var abs_ends=[]
      var lookup=[0]
      var exon_length=[0]
      var max_seqlength =0 
      var template_length =0 
      var template_cumulative_length=["-"]
      var begin =0
      var end=0
      var this_exon=[]
      var exon_len=0

      // set Pre-Locus range
      if (this.is_join_complement)
      {
        starts.push(maxreflen)
        ends.push(locus_end+1)
        exon_length.push(this.tbv_tailclip)
      }
      else
      {
        starts.push(1)
        ends.push(locus_begin-1)
        exon_length.push(this.tbv_headclip)
      }
      // set Locus ranges
      if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name) // Locus
      {
        feature_titles.push(this.jsonConfig.stringconstants.empty_transcript_name)
        if (this.is_join_complement)
        {
            starts.push(locus_end)
            ends.push(locus_begin)
        }
        else
        {
            starts.push(locus_begin)
            ends.push(locus_end)
        }
        
        exon_length.push(locus_length)
        max_seqlength=locus_length
        template_length+=locus_length
        template_cumulative_length.push(template_length)
      }
      // or set introns & exons ranges
      else{
        feature_titles.push("Upstream")
        exon_length.push(0)
        template_cumulative_length.push("-")
        if (this.is_join_complement){
          starts.push(locus_end)
        }
        else{
          starts.push(locus_begin)
        }
        if (this.CDSonly.value){
          exon_text="CDS "
          intron_text="CDS_Intron "
        }
        else{
          exon_text="Exon "
          intron_text="Intron "
        }
        for (const item of this.joinlist){
          begin=parseInt(item.split(":")[0])
          end=parseInt(item.split(":")[1])

          if (this.is_join_complement){
            if (end== locus_end){
              ends.push(end)
            }
            else {
              ends.push(end+1)
            }
          }
          else if (begin == locus_begin){
            ends.push(begin)
          }
          else {
            ends.push(begin-1)
          }
          exon_len=Math.abs(end-begin+1)
          max_seqlength+=exon_len
          exon_total+=1
          feature_titles.push(exon_text+exon_total)
          this_exon=[]
          for (let n = begin; n < (end+1); n++){
            this_exon.push(n)
          }
          if (this.is_join_complement){
            this_exon.reverse()
          }
          for (const item2 of this_exon){
            lookup.push(item2)
          }
          if (this.is_join_complement){
            starts.push(end)
            ends.push(begin)
            starts.push(begin-1)
          }          
          else{
            starts.push(begin)
            ends.push(end)
            starts.push(end+1)
          }
          exon_length.push(exon_len)
          template_length+=exon_len
          template_cumulative_length.push(template_length)
          //Next intron / downstream start
          intron_total+=1
          feature_titles.push(intron_text+intron_total.toString()+"-"+(intron_total+1).toString())
          template_cumulative_length.push(template_length)
        } // End of: for (const item of this.joinlist)
      feature_titles[-1]="Downstream"
      intron_total-=1
      exon_length.push(exon_len)
      template_cumulative_length[-1]="-"
      if (this.is_join_complement){
        ends.push(locus_begin)
      }
      else{
        ends.push(locus_end)
      }
      } // End of: if RG_globals.target_transcript_name == RG_globals.empty_transcript_name: # Locus
    //Finishing off 
    feature_titles.push("Post-Locus")
    if (this.is_join_complement){
        starts.push(locus_begin-1)
        ends.push(1)
        exon_length.push(this.tbv_headclip)
      }
    else{
        starts.push(locus_end+1)
        ends.push(maxreflen)
        exon_length.push(this.tbv_tailclip)
      }
    template_cumulative_length.push("-")
    // This to catch where there's no downstream because exon boundary coincides with locus boundary
    if (starts[-2]==starts[-1]){
      ends[-2]=starts[-2]
      feature_titles[-2]="No Downstream"
    }
    // This to catch where there's no upstream because exon boundary coincides with locus boundary
    if (starts[1]==ends[1]){
      feature_titles[1]="No Upstream"
    }
      //console.log(` starts: ${starts}`) // Validation check
      //console.log(` ends: ${ends}`) // Validation check
      //console.log(` this.mrnapos_lookup: ${this.mrnapos_lookup}`) // Validation check
      //console.log(` feature_titles: ${feature_titles}`) // Validation check
      
      this.mrnapos_lookup=lookup

      this.tbv_trans_Begin.value=1
      this.tbv_trans_Begin_ext.value=0
      this.tbv_trans_Begin.max=max_seqlength
      this.tbv_trans_End.max=max_seqlength
      this.tbv_trans_End.value=max_seqlength
      this.tbv_trans_End_ext.value=0
      this.exomeExtension.value=0
      },
      immediate: true,
      deep: true,
    }
  },
  created() {
    pyodideWorker.onerror = e => {
      console.log(
        `Error in pyodideWorker at ${e.filename}, Line: ${e.lineno}, ${e.message}`
      )
    }

    pyodideWorker.onmessage = e => {
      if (!e.data) {
        console.error('Empty message received from Pyodide worker')
        return
      }
      if (e.data.error) {
        console.error('Pyodide worker error: %s', e.data.error)
        return
      }
      const message = e.data
      if (message.file && message.eof) {
        const file = message.file
        if (this.streamWriters.has(file)) {
          this.streamWriters.get(file).close()
          this.streamWriters.delete(file)
          console.info('File %s downloaded', file)
          return
        } else if (this.inMemFiles.has(file)) {
          const lines = this.inMemFiles.get(file)
          const writer = streamSaver.createWriteStream(file).getWriter()
          for (const line of lines) {
            writer.write(encode(line))
          }
          writer.close()
          this.inMemFiles.delete(file)
          console.info('File %s downloaded', file)
          if (this.inMemFiles.size === 0) {
            this.running = false
          }
          return
        }
        console.info('File %s is empty -- ignored', file)
        return
      }
      if (message.file && message.line) {
        const { file, line } = message
        if (file.endsWith('_journal')) {
          this.journal += line
        }
        if (file.endsWith('_readme')) {
          this.readme += line
        }
        if (!this.inMemFiles.has(file)) {
          this.inMemFiles.set(file, [])
        }
        this.inMemFiles.get(file).push(line)
        return
      }
      if (message.file && message.lines) {
        const { file, lines } = message
        if (!this.streamWriters.has(file)) {
          this.streamWriters.set(
            file,
            streamSaver.createWriteStream(file).getWriter()
          )
        }
        const writer = this.streamWriters.get(file)
        for (const line of lines) {
          writer.write(encode(line))
        }
        return
      }
      console.error(
        'Pyodide worker message not understood: %s',
        JSON.stringify(e.data)
      )
    }
  },
  beforeUnmount() {
    for (const writer of this.streamWriters.values()) {
      writer.abort()
    }
    this.inMemFiles.clear()
    this.streamWriters.clear()
  },
  methods: {
    edit_variant_method() {
      if (this.running) return
      this.running = true
      this.tbv_var_subseq.value=this.tbv_ref_subseq.value
      if (this.tbv_ref_subseq.length > 0)
      {
      // Stuff
      this.tbv_ref_subseq_notretrieved = false
      }
      // Just do nothing if the sequence is empty!
      this.running = false
    },
    save_method() {
      if (this.running) return
      this.running = true
      if (this.tbv_var_subseq.value!=this.tbv_ref_subseq.value)
      {
      // No repeat hapname acceptable if in mutfreqLabels. 
      var match = false
      for (const hapname of this.mutfreqLabels)
      {
        if (this.tbv_hap_name.value == hapname)
        {match=true}
      }
      if (! match)
      {
        this.add_to_addvars()
        // These next two lines - could be done by this.add_to_mutfreqs() instead, but search there is not needed 
        this.mutfreqLabels.push(this.tbv_hap_name.value)
        this.mutfreqs.push(50)
      }
      else{
        console.log(`Not saving: duplicate ${this.tbv_hap_name.label}`)
        //`Error in pyodideWorker at ${e.filename}, Line: ${e.lineno}, ${e.message}`
      }
      }
      else{
      // Do nothing if the sequence is unchanged!
      console.log(`Not saving: sequence is unmodified`)
      }
      this.running = false
    },
    add_to_addvars(){
      console.log(`Saving modified sequence to Addvars`)
      //console.log(`tbv_local_Begin is: ${this.tbv_local_Begin}`)
      //var varout=this.check_revcomp(this.tbv_var_subseq.value)
      var addstart=this.tbv_local_Begin // Temporary fix
      var addend = this.tbv_local_End   // Temporary fix
      if (addstart>addend ){
        [addstart,addend]=[addend,addstart]
      }
      //console.log(`var_subseq: ${this.tbv_var_subseq.value}`)
      var addvar={"locus":this.selectedGene,
                "hapname":this.tbv_hap_name.value,
                "varname":this.tbv_var_name.value,
                "local_begin":addstart,
                "local_end":addend,
                "ref_seq":this.tbv_ref_subseq.value,
                "var_seq":this.check_revcomp(this.tbv_var_subseq.value),
                "abs_Begin":this.tbv_abs_Begin,
                "abs_End":this.tbv_abs_End
                }
      this.tbv_AddVars.push(addvar)
    },
    add_to_mutfreqs0(){
      // Keep this until certain not required
      var match = false
      if (this.tbv_AddVars.length >0)
      {
        for (const addvar of this.tbv_AddVars)
        {
          if (addvar.locus==this.selectedGene)
          {
            match =false
            for (const hapname of this.mutfreqLabels)
            {
              if (addvar.hapname == hapname)
              {match=true}
            }
            if (! match)
            {
            this.mutfreqLabels.push(addvar.hapname)
            this.mutfreqs.push(50)
            }
          }
        }
      }
    },
    add_to_mutfreqs(){
      var match = false
      if (this.tbv_AddVars.length >0)
      {
        for (const addvar of this.tbv_AddVars)
        {
          if (addvar.locus==this.selectedGene)
          {
          this.mutfreqLabels.push(addvar.hapname)
          this.mutfreqs.push(50)
          }
        }
      }
    },
    run_method() {
      if (this.running) return
      this.journal = ''
      this.readme = ''
      this.running = true
      /// This is received by webworker.js. Parameter order is not significant, names are!
      pyodideWorker.postMessage({
        targetLocus: JSON.stringify(this.selectedGene),

        transcriptName: JSON.stringify(this.selectedTranscript),
        transcriptId: JSON.stringify(
          this.jsonConfig.Reference_sequences[this.selectedGene]?.mRNA[
            this.selectedTranscript
          ]
        ),
        CDSonly:this.CDSonly.value,
        fragmentLength: this.fragmentLength.value,

        coverageDepth: this.coverageDepth.value,
        exomeExtension: this.exomeExtension.value,
        gaussMean: this.gaussMean.value,
        gaussSD: this.gaussSD.value,
        qualityMin: Math.min(this.qualityMin.value, this.qualityMax.value),

        qualityMax: Math.max(this.qualityMin.value, this.qualityMax.value),
        mutFreqs: JSON.stringify(this.mutfreqs),
        mutLabels: JSON.stringify(this.mutfreqLabels),
        outputSam: this.outputSam.value,
        outputFastq: this.outputFastq.value,

        outputFasta: this.outputFasta.value,
        refSequenceInFasta: this.refSequenceInFasta.value,
        variantSeqInFasta: this.variantSeqInFasta.value,
        sourceFeatureTables: this.sourceFeatureTables.value,
        eachPossibleRead: this.eachPossibleRead.value,

        variantReadsOnly: this.variantReadsOnly.value,
        pairedReads: this.pairedReads.value,
        duplexReads: this.duplexReads.value,
        annotateSourcePositions: this.annotateSourcePositions.value,
        plusAbsolutePositions: this.plusAbsolutePositions.value,

        cigarAnnotated: this.cigarAnnotated.value,
        substitutionLowerCase: this.substitutionLowerCase.value,
        journalSubs: this.journalSubs.value,
        FlipStrand:this.FlipStrand.value,
        tbv_AddVars:JSON.stringify(this.tbv_AddVars)
        /// If you add here,must add to webworker.js and RG_exploder_webmain.py
        /// Be sure to JSON.stringify variables that are strings
      })
    },
    join_method() {
        if (!this.transcriptOptions.includes(this.selectedTranscript)) {
          this.selectedTranscript = this.jsonConfig.stringconstants.empty_transcript_name
        }
        // Using this get no change in this.tbv_trans_Begin.label and NaN for genomic calcuations
        if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name)
        {this.joinlist=this.jsonConfig.Reference_sequences[this.selectedGene]["Locus_range"].split(',')
        //this.tbv_trans_Begin.label="Genomic"
        }
        else if (this.CDSonly.value)
        {this.joinlist=this.jsonConfig.Reference_sequences[this.selectedGene]["CDS_join"][this.selectedTranscript].split(',')
        //this.tbv_trans_Begin.label="CDS"
        }
        else
        {this.joinlist=this.jsonConfig.Reference_sequences[this.selectedGene]["mRNA_join"][this.selectedTranscript].split(',')
        //this.tbv_trans_Begin.label="mRNA"
        }
        this.is_join_complement=this.jsonConfig.Reference_sequences[this.selectedGene]["is_join_complement"]
        //console.log(`JM: joinlist is ${this.joinlist}`)
    },
    get_build_data(){
      this.get_tbv_var_name()
      this.get_tbv_ref_subseq()
      //console.log(` tbv_local_Begin ${this.tbv_local_Begin}`)
    },
    get_tbv_var_name(){
      //  Build a default naming string for the variant
      var tbe=""
      var tee=""
      var varfront=""
      var varname=""
      if (this.tbv_trans_Begin_ext.value==0)
        {tbe=""}
      else
        {tbe=this.tbv_trans_Begin_ext.value
        if (this.tbv_trans_Begin_ext.value >0)
          {tbe="+"+tbe}
        }     
      
      if (this.tbv_trans_End_ext.value==0)
        {tee=""}
      else
        {tee=this.tbv_trans_End_ext.value
        if (this.tbv_trans_End_ext.value >0)
        {tee="+"+tee}
        }

      //console.log(`get_tbv_var_name: tee is ${tee} tbe is ${tbe}`)
      if (this.tbv_trans_Begin.value==this.tbv_trans_End.value && (this.tbv_trans_Begin_ext.value==this.tbv_trans_End_ext.value))
         if (this.tbv_trans_Begin_ext.value==0)
            {varname=this.tbv_trans_Begin.value}
          else {varname=this.tbv_trans_Begin.value.toString()+tbe}
      else
        {varname=this.tbv_trans_Begin.value+tbe.toString()+"_"+this.tbv_trans_End.value+tee.toString()}
      //console.log(`get_tbv_var_name: varname is ${varname}`)
      if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name)
      {varfront="g"
      }
      else if (this.CDSonly.value)
      {varfront="c"}
      else
      {varfront="t"}
      //console.log(`get_tbv_var_name: varfront is ${varfront}`)
      this.tbv_var_name.value=varfront+"."+varname
    },
    get_tbv_ref_subseq(){
      var subref=""
      var viewstring=""
      var first=""
      var last =""
      var start = this.tbv_local_Begin
      var end = this.tbv_local_End
      //var is_join_complement=this.jsonConfig.Reference_sequences[this.selectedGene]["is_join_complement"]

      if (start > end){
        [start,end]=[end,start]
      }
      var reflen=Math.abs(end-start+1)

      //console.log(`get_tbv_ref_subseq: tbv_local_End is ${this.tbv_local_End}`)
      //console.log(`get_tbv_ref_subseq: tbv_local_Begin is ${this.tbv_local_Begin}`)
      //console.log(`get_tbv_ref_subseq: reflen is ${reflen}`)
    
     if((reflen <1) || (reflen > this.REFSEQ_len ))
      { // Out of range, should not happen with real sequence. It will with test REFSEQ at top of this code @120nt, or from config.json
        viewstring="OUT OF RANGE AT "+reflen
        //subref="OUT OF RANGE AT "+reflen + "MR " + this.REFSEQ_len
        subref="GATTACA"
        reflen=-1 // Locks out ability to retrieve
      }
      else if (reflen > 15) {
      // Within a range for NNN definition
      //  first = this.REFSEQ.substr(parseInt(start)-1, 3)
      //  last= this.REFSEQ.substr(parseInt(end)-3, 3)
      //  if (this.is_join_complement)
      //  {
      //    first=this.get_revcomp(first)
      //    last=this.get_revcomp(last)
      //  }
      // Cannot implement this properly at present because have no means of reading Reference into App.vue
        //viewstring=first+"..."+last //viewstring="Undefined > 100"
        viewstring="NNN...NNN" //
        //subref="NNN"
        subref="nnn---nnn"

      }
     else { 
      // Within range for subsequence
      // console.log(`**: start is ${start} end is ${end}, length is ${reflen} `)
      //  subref = this.REFSEQ.substr(start-1,reflen) //subref="whatever"
      //  if (this.is_join_complement)
      //  {
      //    subref=this.get_revcomp(subref)
      //  }
      // Cannot implement this at present because have no means of reading Reference into App.vue
        //subref="NNN"
        subref='N'.repeat(reflen)
        //viewstring=subref
     }
      this.tbv_ref_subseq.value=subref
      //viewstring=start+"..."+end // testing
      viewstring=this.tbv_abs_Begin.value+"..."+this.tbv_abs_End.value // testing
      //this.tbv_ref_subseq.viewstring=Math.abs(subref.length)+" bases:"+viewstring+ " b:"+start +" e:"+end
      //this.tbv_ref_subseq.viewstring=reflen+" bases:"+viewstring+ " b:"+start +" e:"+end
      this.tbv_ref_subseq.viewstring=reflen+" bases: "+viewstring
      this.tbv_var_subseq.value=subref
      this.tbv_ref_subseq.length=reflen
      this.tbv_hap_name.value=this.tbv_hap_name_def
      this.tbv_ref_subseq_notretrieved=true
      //console.log(`**: subref is ${subref} reflen is ${reflen}`)
    },
    get_abs_begin(){
      // Sets this.tbv_local_Begin and this.tbv_abs_Begin.value
      var modpos =0
      //console.log(`v1: this.tbv_local_Begin is ${this.tbv_local_Begin}`)
      if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name)
      {this.tbv_local_Begin=this.tbv_trans_Begin.value+this.tbv_trans_Begin_ext.value+this.tbv_headclip
      // console.log(`v2: this.tbv_local_Begin is ${this.tbv_local_Begin}`)
      }
      else
        //{if (this.jsonConfig.Reference_sequences[this.selectedGene]["is_join_complement"])
        {if (this.is_join_complement)
                {modpos= -this.tbv_trans_Begin_ext.value} // alternative to headclip, but does same job here
        else 
                {modpos= this.tbv_trans_Begin_ext.value}
        this.tbv_local_Begin=parseInt(this.mrnapos_lookup[this.tbv_trans_Begin.value])+modpos
        //console.log(`v4: this.tbv_local_Begin is ${this.tbv_local_Begin}`)
        }
        this.tbv_abs_Begin.value=Math.abs(this.tbv_abs_offset+this.tbv_local_Begin)
        this.tbv_abs_Begin.txt= this.tbv_Chrom+":"+this.tbv_abs_Begin.value
      //console.log(`get_abs_begin: modpos is ${modpos}`)
      this.get_build_data()
      //console.log(`get_abs_begin: this.tbv_ref_subseq_notretrieved is ${this.tbv_ref_subseq_notretrieved}`)
    },
    get_abs_end(){
      // Sets this.tbv_local_End and this.tbv_abs_End.value
      var modpos =0
      if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name)
      {this.tbv_local_End=this.tbv_trans_End.value+this.tbv_trans_End_ext.value+this.tbv_headclip}
      else
        //{if (this.jsonConfig.Reference_sequences[this.selectedGene]["is_join_complement"])
        {if (this.is_join_complement)
                {modpos= -this.tbv_trans_End_ext.value} // alternative to headclip, but does same job here
        else 
                {modpos= this.tbv_trans_End_ext.value}
        this.tbv_local_End=this.mrnapos_lookup[this.tbv_trans_End.value]+modpos
        }
        this.tbv_abs_End.value=Math.abs(this.tbv_abs_offset+this.tbv_local_End)
        this.tbv_abs_End.txt= this.tbv_Chrom+":"+this.tbv_abs_End.value
      //console.log(`get_abs_end: modpos is ${modpos}`)
      this.get_build_data()
    },    
    get_complement_seq:function(instring){
      //console.log(`get_complement_seq ${instring}`)
      return this.DNA_COMPLEMENT[instring]
    },
    get_reverse_seq:function(instring){
      //console.log(`get_reverse_seq ${instring}`)
      return instring.split('').reverse().join('')
    },
    get_rev_complement:function(instring){  
      //console.log(`get_rev_complement ${compstring}`)
      return instring.split('').reverse().map(this.get_complement_seq).join('')
    },
    check_revcomp:function(instring){
      // Temporary fix ?? Reverse complement tbv_var_subseq.value only if mRNA or CDS and complement is set 
      var compseq=instring
      if ("build_label" != 'Locus' && this.is_join_complement){
        //console.log(`:function(instring) is triggered: ${this.tbv_var_subseq.value}`)
        compseq= this.get_rev_complement(instring)
      }
      return compseq
    },
  },
  computed: {
    build_label:function() {
        if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name){
          return("Locus")
        }
        else if (this.CDSonly.value){
          return("CDS")
          }
        else {
          return("mRNA")
          }
    },
    ensembl_locusURL: function() { // Supercedes LRG_locusURL as hyperlink from locus-selection
      if  (this.tbv_GRChver_txt=="GRCh38"){
        var GRChver=this.jsonConfig.stringconstants.ensembl38_gene_url
      }
      else{
        var GRChver=this.jsonConfig.stringconstants.ensembl37_gene_url
      }
        return (
          GRChver + 'g=' + this.geneList[this.selectedGene].Ensembl_id.split(".")[0]
          )
    },
    ensembl_transcriptURL: function() {
      var geneid=this.geneList[this.selectedGene].Ensembl_id.split(".")[0]
      if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name) {
        // Same as ensembl_locusURL computed
        if  (this.tbv_GRChver_txt=="GRCh38"){
        var GeneURL=this.jsonConfig.stringconstants.ensembl38_gene_url
      }
      else{
        var GeneURL=this.jsonConfig.stringconstants.ensembl37_gene_url
      }
        return (
          GeneURL + 'g=' + geneid
          )
      }
      else {
        var urltarget="38"
        var transid=this.jsonConfig.Reference_sequences[this.selectedGene]["mRNA"][this.selectedTranscript].split(".")[0]

        if  (this.tbv_GRChver_txt=="GRCh37"){
        // Special case where Ensembl transcript id not present in GRCh37, only in GRCh38, indicated by 'm'
        // Exists for AK2 ENST00000672715m.1 - hand edited into the data file AK2_locseq.gb for this application
        console.log(`transid is ${transid}`)
        
        if (transid.includes("m"))
        {
          transid = transid.slice(0, -1); 
          var urltarget="38"
        }
        else {
          var urltarget="37"
        }
      }
        if  (urltarget=="38")
        {
          var TranscriptURL= this.jsonConfig.stringconstants.ensembl38_transcript_url
        }
        else{
          var TranscriptURL= this.jsonConfig.stringconstants.ensembl37_transcript_url
        }
        {return (
          TranscriptURL +'g=' + geneid + ';t=' + transid
          )
        }
      }
    },
    ensembl_locusURLtxt: function() {
      return (
      this.geneList[this.selectedGene].Ensembl_id.split(".")[0]
      )
    },
    ensembl_transcriptURLtxt: function() {
      if (this.selectedTranscript == this.jsonConfig.stringconstants.empty_transcript_name){
         return (
          this.geneList[this.selectedGene].Ensembl_id.split(".")[0]
         )
        }
      else {
        return(
          this.jsonConfig.Reference_sequences[this.selectedGene]["mRNA"][this.selectedTranscript].split(".")[0]
        )
      }
    },
    transcriptOptions: function() {
      const json = this.jsonConfig.Reference_sequences[this.selectedGene]
      const noTrans = this.jsonConfig.stringconstants.empty_transcript_name
      return [noTrans, ...Object.keys(json?.mRNA ?? {})]
    }
   }
}
</script>

<style>
:root {
  --primary-color: rgb(79, 174, 206);
  --background-color: #ffffff;
  --text-color: #2c3e50;

  --fnt-normal: 16px;
  --fnt-large: 22px;

  --primary-shadow: rgba(0, 0, 0, 0.16) 0px 1px 4px;
}

body {
  margin: 0;
  padding: 0;
  box-sizing: border-box;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
}

#app {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  background-color: var(--background-color);
  color: var(--text-color);
}

.details {
  display: flex;
  align-items: center;
  justify-content: center;
  flex-wrap: wrap;
}

.details > * {
  margin: 0 1em 0;
  padding: 0;
  font-size: 1rem;
  font-weight: bold;
  color: var(--text-color);
}

.container {
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
  justify-content: space-evenly;
  align-items: flex-start;
  width: 80%;
  margin: 2rem auto;
}

.container__child {
  min-width: 18em;
  width: 30%;
  margin: 0.5em;
  padding: 1em;
  padding-right: 0;
  margin-right: 0;
  border-radius: 8px;
  box-shadow: var(--primary-shadow);
}

.container__child2 {
  min-width: 40em;
  width: 30%;
  margin: 0.5em;
  padding: 1em;
  padding-right: 0;
  margin-right: 0;
  border-radius: 8px;
  box-shadow: var(--primary-shadow);
}

.container__child3 {
  min-width: 30em;
  width: 100%;
  margin: 0.5em;
  padding: 1em;
  padding-right: 0;
  margin-right: 0;
  border-radius: 8px;
  box-shadow: var(--primary-shadow);
}

.container__child-subtitle {
  display: flex;
  justify-content: space-between;
  margin-right: 1em;
}

.row {
  display: flex;
}

.column {
  flex: 10%;
}
</style>
