import RG_exploder_globals as RG_globals
import RG_exploder_main as RG_main
import RG_exploder_io as RG_io

# Order of these parameters must be coordinated with webworker.js
def run(target_locus, transcript_name, transcript_id, cdsonly,fragment_length,
        coverage_depth, exome_extension, gauss_mean, gauss_SD, quality_min,
        quality_max, mut_freqs, mut_labels,write_sam, write_fastq,
        write_fasta, refSequenceInFasta, variantSeqInFasta, sourceFeatureTables, eachPossibleRead,
        variantReadsOnly, pairedReads, duplexReads, annotateSourcePositions,plusAbsolutePositions,
        cigarAnnotated, substitutionLowerCase, journal_subs, flip_strand,tbv_AddVars,
        js):
    RG_io.JS = js
    RG_globals.target_locus = target_locus
    if transcript_name == RG_globals.empty_transcript_name:
        RG_globals.target_transcript_name =RG_globals.empty_transcript_name
        RG_globals.target_transcript_id = RG_globals.empty_transcript_id
    else:
        RG_globals.target_transcript_name = transcript_name
        RG_globals.target_transcript_id = transcript_id
    RG_globals.is_CDS = cdsonly
    RG_globals.Fraglen = fragment_length
    RG_globals.Fragdepth = coverage_depth
    RG_globals.Exome_extend = exome_extension
    RG_globals.gauss_mean = gauss_mean
    RG_globals.gauss_SD = gauss_SD
    RG_globals.Qualmin = quality_min
    RG_globals.Qualmax = quality_max
    RG_globals.mutfreqs = mut_freqs
    RG_globals.mutlabels = mut_labels
    RG_globals.is_sam_out = write_sam
    RG_globals.is_fasta_out = write_fasta
    RG_globals.is_fastq_out = write_fastq
    RG_globals.is_write_ref_fasta = refSequenceInFasta
    RG_globals.is_mut_out = variantSeqInFasta
    RG_globals.is_write_ref_ingb = sourceFeatureTables
    RG_globals.is_onefrag_out = eachPossibleRead
    RG_globals.is_muts_only = variantReadsOnly
    RG_globals.is_frg_paired_end = pairedReads
    RG_globals.is_duplex = duplexReads
    RG_globals.is_frg_label = annotateSourcePositions
    RG_globals.is_use_absolute = plusAbsolutePositions
    RG_globals.is_fastacigar_out = cigarAnnotated
    RG_globals.is_vars_to_lower = substitutionLowerCase
    RG_globals.is_journal_subs = journal_subs
    RG_globals.is_flip_strand = flip_strand
    # Extra bits for "build" features
    RG_globals.bio_parameters["target_build_variant"]["AddVars"]=tbv_AddVars
    RG_main.call_exploder_main()


