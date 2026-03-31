# Configuration file for the MHAAV25 Nextflow pipeline.
# Pass this to the pipeline with:  --config_file

# ── Paths ────────────────────────────────────────────────────────────────────
fastq_dir="/home/bilbaom/mhaav25/mhaav_nf_runs/data_benchmark_2cas9/seq"
genome_fasta="/home/bilbaom/nanopore/ref/GCF_000001635.26_GRCm38.p6_genomic.fna"
outdir="/home/bilbaom/mhaav25/mhaav_nf_runs/data_benchmark_2cas9"

# ── Amplicon primers ─────────────────────────────────────────────────────────
fw='TAAAAGCATCCTAGGAAGGG'
rev='AGACCAATGTTTGTCAGAGG'

# ── Target site (quantification window) ──────────────────────────────────────
# This is the region where indels are quantified
amplicon='TGCTGACTCTCTGTCCTAAAACAGAAGTTGACAGATCGATATCAGCAACGTTGCGAAGCATCCGTGGATAGAGCTTCCATctggaatttaaaataatatt'
cutSite=14
blat_T_name="NC_000068.7"
blat_T_start=134548201
blat_T_end=134548300

# ── CrispRVariants output filenames (no need to change) ──────────────────────
restable="results.csv"
instable="insertions.csv"

# ── minimap2 alignment scoring (default optimised for AAV insertions) ────────
minimap_A=4
minimap_B=27
minimap_O=32
minimap_E=1
