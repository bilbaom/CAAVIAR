#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================
// MHAAV25 amplicon-seq pipeline
// Usage:
//   nextflow run main.nf --csv SraRunTable.csv --config_file mhaav25_config.sh
//   nextflow run main.nf -profile slurm --csv ... --config_file ...
//   nextflow run main.nf -resume ...
// ============================================================

// ---------- Only these two are passed on the command line ----------
params.csv         = null
params.config_file = null

// ---------- Validate required inputs ----------
if (!params.csv)         error "Please provide --csv <SraRunTable.csv>"
if (!params.config_file) error "Please provide --config_file <mhaav25_config_*.sh>"

// ---------- Parse the shell config file into a Groovy map ----------
def parseCfg(path) {
    def map = [:]
    new File(path).eachLine { line ->
        def m = line =~ /^\s*(\w+)\s*=\s*['"]?([^'"#\n]*)['"]?\s*(?:#.*)?$/
        if (m) map[m[0][1]] = m[0][2].trim()
    }
    return map
}
def cfg = parseCfg(params.config_file)

def require = { map, key ->
    if (!map.containsKey(key) || !map[key])
        error "Required variable '${key}' not found in config file"
    return map[key]
}

// --- Experiment-specific ---
params.fw           = require(cfg, 'fw')
params.rev          = require(cfg, 'rev')
params.amplicon     = require(cfg, 'amplicon')
params.cutSite      = require(cfg, 'cutSite')
params.blat_T_name  = require(cfg, 'blat_T_name')
params.blat_T_start = require(cfg, 'blat_T_start')
params.blat_T_end   = require(cfg, 'blat_T_end')
params.restable     = require(cfg, 'restable')
params.instable     = require(cfg, 'instable')

// --- Paths ---
params.fastq_dir    = require(cfg, 'fastq_dir')
params.genome_fasta = require(cfg, 'genome_fasta')
params.outdir       = require(cfg, 'outdir')

// --- minimap2 scoring params ---
params.minimap_A    = require(cfg, 'minimap_A')
params.minimap_B    = require(cfg, 'minimap_B')
params.minimap_O    = require(cfg, 'minimap_O')
params.minimap_E    = require(cfg, 'minimap_E')

// Convenience label used in publishDir paths across modules
params.param_dir = "A${params.minimap_A}_B${params.minimap_B}_O${params.minimap_O}_E${params.minimap_E}"
                       .replace(',', '_')

// ============================================================
// MODULES
// ============================================================

include { INDEX_GENOME      } from './modules/index_genome.nf'
include { TRIM_READS        } from './modules/trim.nf'
include { MERGE_READS       } from './modules/merge.nf'
include { FILTER_AMPLICONS  } from './modules/filter_amplicons.nf'
include { MAP_READS         } from './modules/map_reads.nf'
include { FILTER_MAPPINGS   } from './modules/filter_mappings.nf'
include { BLAST_AAV         } from './modules/blast_aav.nf'
include { STEP1_CALL_INDELS } from './modules/step1_call_indels.nf'
include { STEP2_INSERTIONS  } from './modules/step2.nf'
include { STEP3_MH          } from './modules/step3.nf'
include { GATHER_RESULTS    } from './modules/gather.nf'

// ============================================================
// WORKFLOW
// ============================================================

workflow {

    // ----------------------------------------------------------
    // 1. Build sample channel from CSV
    //    Emits: tuple(name, group, library)
    // ----------------------------------------------------------
    samples_ch = Channel
        .fromPath(params.csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.Run, row.Group ?: '', row.Library ?: '') }

    // ----------------------------------------------------------
    // 2. Locate FASTQ pairs for each sample
    //    Emits: tuple(name, r1, r2)
    // ----------------------------------------------------------
    fastqs_ch = samples_ch.map { name, group, library ->
        def r1  = file("${params.fastq_dir}/${name}_1.fastq.gz")
        def r2  = file("${params.fastq_dir}/${name}_2.fastq.gz")
        def r1b = file("${params.fastq_dir}/${name}_R1_001.fastq.gz")
        def r2b = file("${params.fastq_dir}/${name}_R2_001.fastq.gz")
        if      (r1.exists())  return tuple(name, r1,  r2)
        else if (r1b.exists()) return tuple(name, r1b, r2b)
        else error "Cannot find FASTQ files for sample: ${name}\n" +
                   "  Looked for: ${r1}\n  Looked for: ${r1b}"
    }

    // ----------------------------------------------------------
    // 3. Trim → Merge → Filter amplicons  (per sample, parallel)
    // ----------------------------------------------------------
    TRIM_READS(fastqs_ch)
    MERGE_READS(TRIM_READS.out.trimmed)
    FILTER_AMPLICONS(MERGE_READS.out.merged)

    // ----------------------------------------------------------
    // 4a. Index genome → Map → Filter BAMs  (per sample, parallel)
    // ----------------------------------------------------------
    INDEX_GENOME( file(params.genome_fasta) )
    MAP_READS(FILTER_AMPLICONS.out.filtered, INDEX_GENOME.out.mmi)
    FILTER_MAPPINGS(MAP_READS.out.bam)

    // ----------------------------------------------------------
    // 4b. BLAST merged reads for AAV counts  (per sample, parallel)
    // ----------------------------------------------------------
    BLAST_AAV(FILTER_AMPLICONS.out.filtered)

    // ----------------------------------------------------------
    // 5. STEP1_CALL_INDELS — collect ALL BAMs, run once
    //    All experiment parameters are passed explicitly as val
    //    so the module has no dependency on params.* internally.
    // ----------------------------------------------------------
    all_bams_ch = FILTER_MAPPINGS.out.bam_clean
        .map { name, bam -> bam }
        .collect()

    all_bais_ch = FILTER_MAPPINGS.out.bam_index
        .map { name, bai -> bai }
        .collect()

    STEP1_CALL_INDELS(
        all_bams_ch,
        all_bais_ch,
        params.csv,
        params.amplicon,
        params.cutSite,
        params.blat_T_name,
        params.blat_T_start,
        params.blat_T_end,
        params.restable,
        params.instable,
        params.param_dir
    )

    // ----------------------------------------------------------
    // 6. Fan out STEP1 results into per-sample (name, restable, instable)
    //    tuples for STEP2.
    //    STEP1 emits results_dir = the <param_dir>/ folder.
    //    We combine it with the sample name channel so each sample
    //    gets its own work dir in STEP2.
    // ----------------------------------------------------------
    sample_names_ch = samples_ch.map { name, group, library -> name }

    per_sample_ch = sample_names_ch
        .combine(STEP1_CALL_INDELS.out.results_dir)
        // emits: tuple(name, results_dir_path)
        .map { name, results_dir ->
            def restable = file("${results_dir}/${name}/${params.restable}")
            def instable = file("${results_dir}/${name}/${params.instable}")
            tuple(name, restable, instable)
        }
        // emits: tuple(name, restable_path, instable_path)

    // ----------------------------------------------------------
    // 7. STEP2 (R) → STEP3 (Python), per sample
    //    Experiment parameters passed explicitly as val.
    // ----------------------------------------------------------
    STEP2_INSERTIONS(per_sample_ch, params.blast_db)
    // emits: tuple(name, all_events_del.tsv, summary_df.tsv)

    STEP3_MH(STEP2_INSERTIONS.out.events, params.amplicon, params.cutSite)
    // emits: tuple(name, <name>_merged_summary.csv)

    // ----------------------------------------------------------
    // 8. Gather — wait for ALL step3 outputs AND all AAV counts
    // ----------------------------------------------------------
    all_summaries_ch = STEP3_MH.out.summary
        .map { name, csv -> csv }
        .collect()

    all_aav_counts_ch = BLAST_AAV.out.aav_count
        .map { name, txt -> txt }
        .collect()

    GATHER_RESULTS(all_summaries_ch, all_aav_counts_ch)

    // Label the terminal output so the DAG renders it correctly
    GATHER_RESULTS.out.gathered_csv
        .view { f -> "Pipeline complete → ${f}" }
}
