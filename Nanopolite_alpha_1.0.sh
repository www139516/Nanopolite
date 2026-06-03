#!/bin/bash

# ==============================================================================
# Script Name: NanoPolite_alpha_v1.0.sh
# Feature: v5.12 (on top of v5.11):
#   - FIX: Stage D reverted to clean best-per-locus extraction (v5.10 logic).
#          The v5.11 "Surgical INDEL Correction" is now an explicit opt-in
#          only via --indel-fix flag. Default behaviour is unchanged.
#   - NEW: --indel-fix flag enables reference-anchored INDEL correction.
#          Requires -r. Modifies consensus sequences — read the WARNING
#          in the help text before using.
#   - FIX: Bug 1 — Insertion r_pos not incremented correctly.
#   - FIX: Bug 2 — r_pos increments inconsistent across match/del/ins.
#   - FIX: Bug 3 — BLAST called once per (locus,sample) sequentially;
#          now batched into one call per Stage D run (~50x faster).
#   - FIX: Bug 4 — Temp files written to OUTPUT_DIR root without cleanup;
#          now isolated in temp_stageD/ and removed on completion.
#   - FIX: Bug 5 — fixed_indels_log could not distinguish "no minors"
#          from "BLAST failed"; now reports explicit status strings.
#   - FIX: Bug 6 — one-line function defs with no docstrings; restored.
#   - FIX: Bug 7 — v5.11 header missing changelog; fixed.
#   - FIX: Line-concatenation on former line 443; cleaned up.
#
# v5.10 features carried over:
#   - --from-stage: resume from intermediate results (total/oriented/trimmed/best)
#   - --input-total / --input-polished / --input-trimmed-* overrides
#
# v5.9 features carried over:
#   - Stage D: best-per-locus extraction (5_Best_Per_Locus.fasta)
#
# v5.8 features carried over:
#   - Progress bar, per-sample timing, cleanup trap, disk-space check,
#     drop-rate fix, awk portability fix
# ==============================================================================

set -euo pipefail

# ================= Initialization & Parameters =================
export START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

INPUT_LIST=()
INPUT_DIR=""
OUTPUT_DIR="NanoPolite_results_${TIMESTAMP}"
PRIMERS=""
REF_SEQ=""
VARIANT_TABLE=""
MIN_QUAL=10
THREADS=8
HEAD_CROP=0
TAIL_CROP=0
MIN_LEN=20
ID_THRESHOLD=0.9
MIN_CLUSTER_SIZE=3
RACON_ROUNDS=1
FINAL_MIN_LENGTH=500
FINAL_MAX_LENGTH=1200
MISMATCH=5
SEARCH_RANGE="1:-1"
CHIMERA_CHECK=0
POLISHER="racon"
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"

# ── Resume-from-stage parameters ──────────────────────────────────────────────
FROM_STAGE=""               # empty = normal run from FASTQ
INPUT_TOTAL=""              # override path for 1_Total_Clusters.fasta
INPUT_POLISHED=""           # override path for 2_Polished_Clusters.fasta
INPUT_TRIMMED_TOTAL=""      # override path for 3_Total_Clusters_Trimmed.fasta
INPUT_TRIMMED_POLISHED=""   # override path for 4_Polished_Clusters_Trimmed.fasta

# ── Experimental features ──────────────────────────────────────────────────────
INDEL_FIX=0                 # off by default; enable with --indel-fix

# ================= Help / Usage Information =================
usage() {
cat << EOF
NanoPolite alpha v1.0 - Automated Amplicon Processing, Trimming & Variant Calling

Required Input (At least one — unless --from-stage is used):
  -i, --input        Input fastq file(s). Supports .fastq, .fq, .fastq.gz, .fq.gz
                     and wildcards (e.g., -i "*.fastq.gz")
  -d, --dir          Directory containing fastq files (compressed or not).

Resume From Intermediate Results:
  --from-stage       Skip the per-sample loop and re-enter at a later stage.
                     Use this when you already have output files and want to
                     run only downstream steps without reprocessing FASTQs.
                     Accepted values:
                       total    — resume from 1_Total_Clusters.fasta +
                                  2_Polished_Clusters.fasta
                       oriented — resume after Stage A (orientation already done)
                       trimmed  — resume from 3_Total_Clusters_Trimmed.fasta +
                                  4_Polished_Clusters_Trimmed.fasta
                       best     — re-run Stage D + Stage C only, from
                                  4_Polished_Clusters_Trimmed.fasta
  --input-total             Path to existing 1_Total_Clusters.fasta
                            (default: looks inside the -o output directory)
  --input-polished          Path to existing 2_Polished_Clusters.fasta
  --input-trimmed-total     Path to existing 3_Total_Clusters_Trimmed.fasta
  --input-trimmed-polished  Path to existing 4_Polished_Clusters_Trimmed.fasta

  Example — add primer trimming to an existing run:
    NanoPolite_alpha_v1.0.sh --from-stage total -o existing_results/ \\
        -p primers.csv -r ref.fasta

  Example — re-run best-per-locus with a different primer file:
    NanoPolite_alpha_v1.0.sh --from-stage trimmed -o existing_results/ \\
        -p new_primers.csv

Output Options:
  -o, --out          Designate output folder path (Default: NanoPolite_results_TIMESTAMP)

Optional Trimming & Orientation Parameters:
  -p, --primers      Primer CSV file (Format: ID,Forward,Reverse).
                     *If no reference is provided, primers are used to orient sequences.*
  -r, --ref-seq      FASTA file containing reference sequence(s).
                     *If provided, forces all sequences to match the reference forward strand.*
  --amp-min          Min length for retained amplicons (Default: ${FINAL_MIN_LENGTH} bp)
  --amp-max          Max length for retained amplicons (Default: ${FINAL_MAX_LENGTH} bp)
  --mismatch         Allowed mismatches for primer binding (Default: ${MISMATCH})
  --search-range     seqkit amplicon extraction region (Default: ${SEARCH_RANGE})
                     *1:-1 = full amplicon including primer sequences at both ends (default)*
                     *2:-2 = variable region only, primer sequences excluded*
                     *Use 2:-2 for downstream MSA or haplotype analysis*

Optional Variant Calling Parameters:
  -V, --table        Mapping table (Space/Tab separated: RefName Target_Partial_Names)
                     *Requires -r. Generates the HTML SNP/INDEL report.*
                     *NOTE: -t is now reserved as an alias for --table but with strict
                      validation to avoid collision with -T (threads). Prefer -V.*

Optional QC & Clustering Parameters:
  -q, --quality      Min average Phred quality score (Default: ${MIN_QUAL})
  -m, --min-len      Min read length for QC (Default: ${MIN_LEN} bp)
  --id               vsearch clustering identity (Default: ${ID_THRESHOLD})
  --min-cluster      Min read depth required to polish a cluster (Default: ${MIN_CLUSTER_SIZE})
  --chimera-check    Enable vsearch de-novo chimera filtering after QC (Default: off)
  --rounds           Number of polishing rounds (Default: ${RACON_ROUNDS})
  --polisher         Polisher to use: racon | medaka (Default: ${POLISHER})
  --medaka-model     Medaka model name (Default: ${MEDAKA_MODEL})
  -T, --threads      Number of CPU threads (Default: ${THREADS})
  -h, --help         Display this help message and exit

Experimental Options:
  --indel-fix        [EXPERIMENTAL] Enable reference-anchored INDEL correction
                     in Stage D. Requires -r <reference.fasta>.
                     When active, the best-per-locus consensus for each
                     (locus, sample) is compared to the reference via BLAST.
                     Positions where the consensus has an apparent INDEL that
                     is not supported by any minor cluster at the same locus
                     are surgically patched toward the reference base.
                     WARNING: This modifies consensus sequences using the
                     reference genome and may introduce reference bias.
                     It is intended for variety identification workflows only.
                     Do NOT use for variant discovery or SNP calling.
                     Off by default. Produces an INDELs_Fixed column in the
                     5_Best_Per_Locus_summary.tsv reporting every correction.
EOF
    exit 0
}

if [ $# -eq 0 ]; then usage; fi

# ================= Parse Command Line Arguments =================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input)
            shift
            while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do INPUT_LIST+=("$1"); shift; done
            continue ;;
        -d|--dir) INPUT_DIR="$2"; shift ;;
        -o|--out) OUTPUT_DIR="$2"; shift ;;
        -p|--primers) PRIMERS="$2"; shift ;;
        -r|--ref-seq) REF_SEQ="$2"; shift ;;
        -V|--table) VARIANT_TABLE="$2"; shift ;;
        -t)
            # Backward-compat: accept -t but validate to avoid -T (threads) confusion.
            if [[ "$2" =~ ^[0-9]+$ ]]; then
                echo -e "[ERROR] -t was passed a numeric value ('$2'). This used to mean threads in earlier drafts." >&2
                echo -e "        -t/--table now expects a path to a variant mapping table file." >&2
                echo -e "        Use -T or --threads for the thread count, and -V/--table for the table." >&2
                exit 1
            fi
            VARIANT_TABLE="$2"; shift ;;
        -q|--quality) MIN_QUAL="$2"; shift ;;
        -m|--min-len) MIN_LEN="$2"; shift ;;
        --headcrop) HEAD_CROP="$2"; shift ;;
        --tailcrop) TAIL_CROP="$2"; shift ;;
        --id) ID_THRESHOLD="$2"; shift ;;
        --min-cluster) MIN_CLUSTER_SIZE="$2"; shift ;;
        --rounds) RACON_ROUNDS="$2"; shift ;;
        --mismatch) MISMATCH="$2"; shift ;;
        --amp-min) FINAL_MIN_LENGTH="$2"; shift ;;
        --amp-max) FINAL_MAX_LENGTH="$2"; shift ;;
        --search-range) SEARCH_RANGE="$2"; shift ;;
        --chimera-check) CHIMERA_CHECK=1 ;;
        --indel-fix) INDEL_FIX=1 ;;
        --polisher) POLISHER="$2"; shift ;;
        --medaka-model) MEDAKA_MODEL="$2"; shift ;;
        -T|--threads) THREADS="$2"; shift ;;
        --from-stage) FROM_STAGE="$2"; shift ;;
        --input-total) INPUT_TOTAL="$2"; shift ;;
        --input-polished) INPUT_POLISHED="$2"; shift ;;
        --input-trimmed-total) INPUT_TRIMMED_TOTAL="$2"; shift ;;
        --input-trimmed-polished) INPUT_TRIMMED_POLISHED="$2"; shift ;;
        -h|--help) usage ;;
        *) echo -e "[ERROR] Unknown parameter: $1" >&2; exit 1 ;;
    esac
    shift
done

# Validate polisher choice
if [[ "${POLISHER}" != "racon" && "${POLISHER}" != "medaka" ]]; then
    echo -e "[ERROR] --polisher must be 'racon' or 'medaka' (got '${POLISHER}')" >&2
    exit 1
fi

# Validate --indel-fix dependencies
if [[ "${INDEL_FIX}" -eq 1 ]]; then
    if [[ -z "${REF_SEQ}" ]]; then
        echo -e "[ERROR] --indel-fix requires a reference sequence (-r <ref.fasta>)." >&2
        exit 1
    fi
    if [[ -z "${PRIMERS}" ]]; then
        echo -e "[ERROR] --indel-fix requires a primer file (-p <primers.csv>)." >&2
        echo -e "        INDEL correction operates on Stage D output, which requires Stage B." >&2
        exit 1
    fi
fi

# ================= Setup Output Structure =================
mkdir -p "${OUTPUT_DIR}"
export OUTPUT_DIR
export LOG_FILE="${OUTPUT_DIR}/NanoPolite_run_${TIMESTAMP}.log"
export STATS_FILE="$(pwd)/${OUTPUT_DIR}/NanoPolite_stats_${TIMESTAMP}.tsv"
export REPORT_HTML="$(pwd)/${OUTPUT_DIR}/NanoPolite_Report_${TIMESTAMP}.html"
export VARIANT_HTML="$(pwd)/${OUTPUT_DIR}/NanoPolite_Variant_Report_${TIMESTAMP}.html"

CURRENT_TOTAL="${OUTPUT_DIR}/1_Total_Clusters.fasta"
CURRENT_POLISHED="${OUTPUT_DIR}/2_Polished_Clusters.fasta"

# Only truncate the output FASTAs for a fresh run.
# In resume mode the files are populated by the file-injection block below.
if [[ -z "${FROM_STAGE}" ]]; then
    > "${CURRENT_TOTAL}"
    > "${CURRENT_POLISHED}"
fi

# ================= Logging Function =================
log() {
    local TYPE=$1
    shift
    local MSG="$@"
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [${TYPE}] ${MSG}" | tee -a "${LOG_FILE}"
}

# ================= Safe count helper (fixes grep -c || echo bug) =================
# Returns 0 cleanly when file is empty/missing instead of concatenating outputs.
count_fasta_seqs() {
    local F="$1"
    if [ ! -s "$F" ]; then
        echo 0
        return
    fi
    grep -c "^>" "$F" 2>/dev/null || echo 0
}

# ================= Progress bar (terminal-aware) =================
# Renders a single-line progress bar to stderr if stderr is a TTY,
# otherwise prints a one-line milestone.
# Args: current_count total_count optional_label
draw_progress() {
    local cur=$1
    local total=$2
    local label="${3:-}"
    if [ "$total" -le 0 ]; then return; fi
    local pct=$(( cur * 100 / total ))
    local bar_width=30
    local filled=$(( cur * bar_width / total ))
    local empty=$(( bar_width - filled ))
    local bar=""
    local i
    for ((i=0; i<filled; i++)); do bar="${bar}█"; done
    for ((i=0; i<empty; i++)); do bar="${bar}░"; done

    # Compute ETA from PROGRESS_START_TIME (set by caller)
    local eta_str=""
    if [ -n "${PROGRESS_START_TIME:-}" ] && [ "$cur" -gt 0 ]; then
        local now=$(date +%s)
        local elapsed=$(( now - PROGRESS_START_TIME ))
        if [ "$elapsed" -gt 0 ] && [ "$cur" -lt "$total" ]; then
            local total_est=$(( elapsed * total / cur ))
            local remaining=$(( total_est - elapsed ))
            local eta_min=$(( remaining / 60 ))
            local eta_sec=$(( remaining % 60 ))
            eta_str=$(printf "  ETA %dm%02ds" "$eta_min" "$eta_sec")
        fi
    fi

    if [ -t 2 ]; then
        # TTY — overwrite the line
        printf "\r  [%s] %3d%%  %d/%d%s  %s\033[K" \
            "$bar" "$pct" "$cur" "$total" "$eta_str" "$label" >&2
        if [ "$cur" -eq "$total" ]; then echo "" >&2; fi
    else
        # Non-TTY (log file, pipe) — print milestones at every 10%
        local prev_pct
        if [ -n "${LAST_PROGRESS_PCT:-}" ]; then
            prev_pct="${LAST_PROGRESS_PCT}"
        else
            prev_pct=-1
        fi
        if [ "$pct" -ge $(( prev_pct / 10 * 10 + 10 )) ] || [ "$cur" -eq "$total" ]; then
            printf "  Progress: %d%%  %d/%d%s\n" "$pct" "$cur" "$total" "$eta_str" >&2
            export LAST_PROGRESS_PCT=$pct
        fi
    fi
}

# ================= Cleanup trap (handles Ctrl+C and crashes) =================
cleanup_temp() {
    if [ -n "${OUTPUT_DIR:-}" ] && [ -d "${OUTPUT_DIR}" ]; then
        find "${OUTPUT_DIR}" -maxdepth 1 -type d -name "temp_work_*" -exec rm -rf {} + 2>/dev/null || true
    fi
}

# Error trap so failures don't die silently under set -e
trap 'log "ERROR" "Script failed at line $LINENO (exit code $?). See ${LOG_FILE}."; cleanup_temp' ERR
trap 'echo ""; log "WARNING" "Interrupted by user. Cleaning up..."; cleanup_temp; exit 130' INT TERM
trap 'cleanup_temp' EXIT

# ================= Pre-run Validation & Dependency Checks =================
log "INFO" ">>> Initializing NanoPolite alpha v1.0 <<<"
log "INFO" "Output directory: ${OUTPUT_DIR}"

# ── Validate --from-stage value ───────────────────────────────────────────────
VALID_STAGES=("total" "oriented" "trimmed" "best")
if [[ -n "${FROM_STAGE}" ]]; then
    STAGE_VALID=0
    for s in "${VALID_STAGES[@]}"; do
        [[ "${FROM_STAGE}" == "$s" ]] && STAGE_VALID=1 && break
    done
    if [[ "${STAGE_VALID}" -eq 0 ]]; then
        log "ERROR" "--from-stage '${FROM_STAGE}' is not valid. Choose: total | oriented | trimmed | best"
        exit 1
    fi
    log "INFO" "Resume mode: --from-stage ${FROM_STAGE} — per-sample FASTQ loop will be skipped."
fi

# Core dependencies always required
for cmd in NanoFilt vsearch minimap2 seqkit python3 awk; do
    command -v $cmd &> /dev/null || { log "ERROR" "Missing dependency: $cmd"; exit 1; }
done

# Polisher-specific dependencies
if [[ "${POLISHER}" == "racon" ]]; then
    command -v racon &> /dev/null || { log "ERROR" "Missing dependency: racon (required for --polisher racon)"; exit 1; }
elif [[ "${POLISHER}" == "medaka" ]]; then
    command -v medaka_consensus &> /dev/null || { log "ERROR" "Missing dependency: medaka_consensus (required for --polisher medaka)"; exit 1; }
fi

# ================= Dependency Version Logging =================
log "INFO" "--- Dependency versions ---"
log "INFO" "NanoFilt:  $(NanoFilt --version 2>&1 | head -1 || echo 'unknown')"
log "INFO" "vsearch:   $(vsearch --version 2>&1 | head -1 || echo 'unknown')"
log "INFO" "minimap2:  $(minimap2 --version 2>&1 | head -1 || echo 'unknown')"
log "INFO" "seqkit:    $(seqkit version 2>&1 | head -1 || echo 'unknown')"
log "INFO" "python3:   $(python3 --version 2>&1 | head -1 || echo 'unknown')"
if [[ "${POLISHER}" == "racon" ]]; then
    log "INFO" "racon:     $(racon --version 2>&1 | head -1 || echo 'unknown')"
else
    log "INFO" "medaka:    $(medaka --version 2>&1 | head -1 || echo 'unknown')"
fi
if [[ -n "${REF_SEQ}" ]] || [[ -n "${VARIANT_TABLE}" ]]; then
    if command -v blastn &> /dev/null; then
        log "INFO" "blastn:    $(blastn -version 2>&1 | head -1 || echo 'unknown')"
    fi
fi
log "INFO" "---------------------------"

if [[ -n "${REF_SEQ}" ]]; then
    command -v blastn &> /dev/null || { log "ERROR" "Missing dependency for orientation: blastn"; exit 1; }
    if [[ ! -f "${REF_SEQ}" ]]; then
        log "ERROR" "Reference sequence file not found: ${REF_SEQ}"; exit 1
    fi
    log "INFO" "Reference sequence provided. Output sequences will be forced into the forward orientation."
fi

if [[ -n "${VARIANT_TABLE}" ]]; then
    if [[ -z "${REF_SEQ}" ]]; then
        log "ERROR" "Variant calling (-V/--table) requires a reference sequence (-r)."; exit 1
    fi
    if [[ ! -f "${VARIANT_TABLE}" ]]; then
        log "ERROR" "Mapping table file not found: ${VARIANT_TABLE}"; exit 1
    fi
    log "INFO" "Variant calling table provided. BLASTn HTML report will be generated."
fi

# ================= Build File List (with gzip support) =================
# Skipped entirely when --from-stage is provided.
if [[ -z "${FROM_STAGE}" ]]; then
    FILES=()
    for f in "${INPUT_LIST[@]}"; do [ -f "$f" ] && FILES+=("$f"); done
    if [ -n "$INPUT_DIR" ] && [ -d "$INPUT_DIR" ]; then
        shopt -s nullglob
        for f in "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq "$INPUT_DIR"/*.fastq.gz "$INPUT_DIR"/*.fq.gz; do
            [ -f "$f" ] && FILES+=("$f")
        done
        shopt -u nullglob
    fi

    # Safely handle filenames with spaces during deduplication
    if [ ${#FILES[@]} -gt 0 ]; then
        mapfile -t FILES < <(printf "%s\n" "${FILES[@]}" | sort -u)
    else
        log "ERROR" "No valid fastq files found via -i or -d."; exit 1
    fi
fi

echo -e "Sample\tRaw_Reads\tClean_Reads\tClusters\tPolished_Seqs\tUnpolished_Seqs" > "${STATS_FILE}"

# ================= Resume-from-stage: inject existing files =================
# When --from-stage is set, locate the user's existing FASTA files and point
# CURRENT_TOTAL / CURRENT_POLISHED at them so downstream stages run correctly.
# Files are COPIED into the output directory so the originals are not modified.
if [[ -n "${FROM_STAGE}" ]]; then

    # Helper: resolve input file path.
    # Sets global RESOLVED_PATH; exits with helpful message if not found.
    _resolve_input() {
        local override="$1" default_path="$2" label="$3"
        if [[ -n "${override}" ]]; then
            if [[ ! -f "${override}" ]]; then
                echo -e "[ERROR] Resume mode: --input-${label} file not found: ${override}" | tee -a "${LOG_FILE}" >&2
                exit 1
            fi
            RESOLVED_PATH="${override}"
        elif [[ -f "${default_path}" ]] && [[ -s "${default_path}" ]]; then
            RESOLVED_PATH="${default_path}"
        else
            echo -e "[ERROR] Resume mode (--from-stage ${FROM_STAGE}) needs the '${label}' file." | tee -a "${LOG_FILE}" >&2
            echo -e "[ERROR] Expected at: ${default_path}" | tee -a "${LOG_FILE}" >&2
            echo -e "[ERROR] Either copy it there or use: --input-${label} /path/to/file" | tee -a "${LOG_FILE}" >&2
            exit 1
        fi
    }

    case "${FROM_STAGE}" in

        total|oriented)
            # Need 1_Total and 2_Polished
            _resolve_input "${INPUT_TOTAL}"    "${OUTPUT_DIR}/1_Total_Clusters.fasta"   "total"
            SRC_TOTAL="${RESOLVED_PATH}"
            _resolve_input "${INPUT_POLISHED}" "${OUTPUT_DIR}/2_Polished_Clusters.fasta" "polished"
            SRC_POLISHED="${RESOLVED_PATH}"
            # Copy to canonical names in output dir if they are external files
            if [[ "${SRC_TOTAL}"    != "${CURRENT_TOTAL}"    ]]; then cp "${SRC_TOTAL}"    "${CURRENT_TOTAL}"; fi
            if [[ "${SRC_POLISHED}" != "${CURRENT_POLISHED}" ]]; then cp "${SRC_POLISHED}" "${CURRENT_POLISHED}"; fi
            N_TOTAL=$(count_fasta_seqs "${CURRENT_TOTAL}")
            N_POL=$(count_fasta_seqs "${CURRENT_POLISHED}")
            log "INFO" "Resume from stage '${FROM_STAGE}': loaded ${N_TOTAL} total, ${N_POL} polished sequences."
            log "WARNING" "Stats TSV will be empty (per-sample loop was skipped). Re-run from FASTQ to regenerate stats."
            # If oriented, Stage A is already done — skip it by unsetting REF_SEQ guard
            if [[ "${FROM_STAGE}" == "oriented" ]]; then
                log "INFO" "Stage '${FROM_STAGE}': Stage A orientation will be skipped (already applied)."
                _SKIP_STAGE_A=1
            else
                _SKIP_STAGE_A=0
            fi
            ;;

        trimmed)
            # Need 3_Total_Trimmed and 4_Polished_Trimmed
            _resolve_input "${INPUT_TRIMMED_TOTAL}"    "${OUTPUT_DIR}/3_Total_Clusters_Trimmed.fasta"   "trimmed-total"
            SRC_TRIM_TOTAL="${RESOLVED_PATH}"
            _resolve_input "${INPUT_TRIMMED_POLISHED}" "${OUTPUT_DIR}/4_Polished_Clusters_Trimmed.fasta" "trimmed-polished"
            SRC_TRIM_POL="${RESOLVED_PATH}"            OUT_TOTAL_TRIMMED="${OUTPUT_DIR}/3_Total_Clusters_Trimmed.fasta"
            OUT_POLISHED_TRIMMED="${OUTPUT_DIR}/4_Polished_Clusters_Trimmed.fasta"
            if [[ "${SRC_TRIM_TOTAL}" != "${OUT_TOTAL_TRIMMED}" ]]; then cp "${SRC_TRIM_TOTAL}" "${OUT_TOTAL_TRIMMED}"; fi
            if [[ "${SRC_TRIM_POL}"   != "${OUT_POLISHED_TRIMMED}" ]]; then cp "${SRC_TRIM_POL}" "${OUT_POLISHED_TRIMMED}"; fi
            export CURRENT_TOTAL="${OUT_TOTAL_TRIMMED}"
            export CURRENT_POLISHED="${OUT_POLISHED_TRIMMED}"
            N_TRIM=$(count_fasta_seqs "${CURRENT_TOTAL}")
            log "INFO" "Resume from stage 'trimmed': loaded ${N_TRIM} trimmed sequences."
            log "INFO" "Stages A and B will be skipped. Running Stage D and Stage C only."
            log "WARNING" "Stats TSV will be empty (per-sample loop was skipped)."
            _SKIP_STAGE_A=1
            _SKIP_STAGE_B=1
            ;;

        best)
            # Need 4_Polished_Trimmed only — re-run Stage D + C
            _resolve_input "${INPUT_TRIMMED_POLISHED}" "${OUTPUT_DIR}/4_Polished_Clusters_Trimmed.fasta" "trimmed-polished"
            SRC_TRIM_POL="${RESOLVED_PATH}"
            OUT_POLISHED_TRIMMED="${OUTPUT_DIR}/4_Polished_Clusters_Trimmed.fasta"
            if [[ "${SRC_TRIM_POL}" != "${OUT_POLISHED_TRIMMED}" ]]; then cp "${SRC_TRIM_POL}" "${OUT_POLISHED_TRIMMED}"; fi
            export CURRENT_POLISHED="${OUT_POLISHED_TRIMMED}"
            N_BEST=$(count_fasta_seqs "${CURRENT_POLISHED}")
            log "INFO" "Resume from stage 'best': loaded ${N_BEST} trimmed polished sequences."
            log "INFO" "Running Stage D (best-per-locus) and Stage C (QC report) only."
            _SKIP_STAGE_A=1
            _SKIP_STAGE_B=1
            ;;
    esac
else
    _SKIP_STAGE_A=0
    _SKIP_STAGE_B=0
fi

# ================= Disk space pre-check =================
# Each sample creates a temp dir; warn if available space looks tight.
AVAIL_GB=$(df -BG "${OUTPUT_DIR}" 2>/dev/null | awk 'NR==2 {gsub(/G/, "", $4); print $4}')
if [ -n "${AVAIL_GB}" ] && [ "${AVAIL_GB}" -lt 5 ]; then
    log "WARNING" "Available disk space in output directory is only ${AVAIL_GB} GB. Consider freeing space if processing many samples."
fi

# Helper: stream a fastq file (handles .gz transparently)
stream_fastq() {
    local F="$1"
    if [[ "$F" == *.gz ]]; then
        zcat -- "$F"
    else
        cat -- "$F"
    fi
}
export -f stream_fastq

# Helper: count reads in a fastq file (handles .gz)
count_fastq_reads() {
    local F="$1"
    local LINES
    if [[ "$F" == *.gz ]]; then
        LINES=$(zcat -- "$F" | wc -l)
    else
        LINES=$(wc -l < "$F")
    fi
    echo $((LINES / 4))
}

# ================= Resume support =================
# Build set of already-completed sample names from the stats TSV.
# Note: STATS_FILE was just truncated above with a header; resume is enabled
# only if the file has more than just the header (i.e. a previous partial run
# was preserved by the user via a manual --resume flag in the future).
# For now, this scaffolding is in place but the file is fresh on every run.
declare -A COMPLETED_SAMPLES=()

# ================= Core Processing Logic =================
# Skipped when --from-stage is provided (user supplies existing FASTA files).
if [[ -z "${FROM_STAGE}" ]]; then
TOTAL_FILES=${#FILES[@]}
log "INFO" "Processing ${TOTAL_FILES} files..."
export PROGRESS_START_TIME=$(date +%s)
SAMPLE_INDEX=0

for INPUT_FILE in "${FILES[@]}"; do
    SAMPLE_INDEX=$(( SAMPLE_INDEX + 1 ))

    # Strip .gz first if present, then strip the fastq extension
    BASE_NAME=$(basename "$INPUT_FILE")
    BASE_NAME="${BASE_NAME%.gz}"
    FILENAME="${BASE_NAME%.*}"

    # Skip if already completed (resume support)
    if [ -n "${COMPLETED_SAMPLES[$FILENAME]:-}" ]; then
        log "INFO" "Skipping (already completed): ${FILENAME}"
        draw_progress "${SAMPLE_INDEX}" "${TOTAL_FILES}" "skipped: ${FILENAME}"
        continue
    fi

    SAMPLE_DIR="${OUTPUT_DIR}/temp_work_${FILENAME}"
    SAMPLE_START=$(date +%s)
    log "INFO" "[${SAMPLE_INDEX}/${TOTAL_FILES}] Processing: ${FILENAME}"

    RAW_READS=$(count_fastq_reads "$INPUT_FILE")
    rm -rf "${SAMPLE_DIR}"; mkdir -p "${SAMPLE_DIR}/clean_data" "${SAMPLE_DIR}/clusters"

    # QC filtering with gzip-aware streaming
    stream_fastq "$INPUT_FILE" \
        | NanoFilt -q ${MIN_QUAL} -l ${MIN_LEN} --headcrop ${HEAD_CROP} --tailcrop ${TAIL_CROP} \
        > "${SAMPLE_DIR}/clean_data/cleaned.fastq"

    if [ ! -s "${SAMPLE_DIR}/clean_data/cleaned.fastq" ]; then
        echo -e "${FILENAME}\t${RAW_READS}\t0\t0\t0\t0" >> "${STATS_FILE}"
        rm -rf "${SAMPLE_DIR}"; continue
    fi

    # Optional chimera filtering (de-novo, on cleaned reads)
    if [ "${CHIMERA_CHECK}" -eq 1 ]; then
        log "INFO" "Running de-novo chimera detection (vsearch --uchime3_denovo)..."
        # vsearch chimera detection works on FASTA; convert temporarily, then re-pair to FASTQ.
        seqkit fq2fa "${SAMPLE_DIR}/clean_data/cleaned.fastq" > "${SAMPLE_DIR}/clean_data/cleaned.fasta"
        vsearch --uchime3_denovo "${SAMPLE_DIR}/clean_data/cleaned.fasta" \
                --nonchimeras "${SAMPLE_DIR}/clean_data/nonchim.fasta" \
                --threads ${THREADS} --quiet 2>> "${LOG_FILE}" || true
        if [ -s "${SAMPLE_DIR}/clean_data/nonchim.fasta" ]; then
            # Filter the FASTQ to keep only non-chimeric IDs
            seqkit seq -n -i "${SAMPLE_DIR}/clean_data/nonchim.fasta" \
                > "${SAMPLE_DIR}/clean_data/keep.ids"
            seqkit grep -f "${SAMPLE_DIR}/clean_data/keep.ids" \
                "${SAMPLE_DIR}/clean_data/cleaned.fastq" \
                > "${SAMPLE_DIR}/clean_data/cleaned.nonchim.fastq" 2>/dev/null || true
            if [ -s "${SAMPLE_DIR}/clean_data/cleaned.nonchim.fastq" ]; then
                mv "${SAMPLE_DIR}/clean_data/cleaned.nonchim.fastq" "${SAMPLE_DIR}/clean_data/cleaned.fastq"
            fi
        fi
        rm -f "${SAMPLE_DIR}/clean_data/cleaned.fasta" \
              "${SAMPLE_DIR}/clean_data/nonchim.fasta" \
              "${SAMPLE_DIR}/clean_data/keep.ids"
    fi

    CLEAN_READS=$(($(wc -l < "${SAMPLE_DIR}/clean_data/cleaned.fastq") / 4))
    vsearch --cluster_fast "${SAMPLE_DIR}/clean_data/cleaned.fastq" \
            --id ${ID_THRESHOLD} --strand both \
            --clusters "${SAMPLE_DIR}/clusters/cluster_" \
            --minseqlength ${MIN_LEN} --threads ${THREADS} --quiet

    CLUSTER_COUNT=$(find "${SAMPLE_DIR}/clusters" -name "cluster_*" | wc -l)
    POLISHED_COUNT=0; UNPOLISHED_COUNT=0

    while read -r CLUSTER_FILE; do
        mv "${CLUSTER_FILE}" "${CLUSTER_FILE}.fasta"
        FILE_EXT="${CLUSTER_FILE}.fasta"
        CLUSTER_NAME=$(basename "${FILE_EXT}" .fasta)
        READ_COUNT=$(grep -c "^>" "${FILE_EXT}" || true)

        if [ "${READ_COUNT}" -lt "${MIN_CLUSTER_SIZE}" ]; then
            # FIX (v5.4): Write the single representative sequence directly to CURRENT_TOTAL,
            # rather than appending to a growing unpolished.fasta and then cat-ing the whole
            # file each iteration (which caused exponential duplication).
            seqkit sort -l -r "${FILE_EXT}" > "${SAMPLE_DIR}/clusters/sorted_tmp.fasta" 2>/dev/null || true
            seqkit head -n 1 "${SAMPLE_DIR}/clusters/sorted_tmp.fasta" \
                | seqkit replace -p ".+" -r "${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_Unpolished" \
                >> "${CURRENT_TOTAL}"
            rm -f "${SAMPLE_DIR}/clusters/sorted_tmp.fasta"
            UNPOLISHED_COUNT=$((UNPOLISHED_COUNT + 1))
        else
            seqkit sort -l -r "${FILE_EXT}" > "${SAMPLE_DIR}/clusters/sorted_tmp.fasta" 2>/dev/null || true
            seqkit head -n 1 "${SAMPLE_DIR}/clusters/sorted_tmp.fasta" > "${SAMPLE_DIR}/clusters/draft.fasta"
            rm -f "${SAMPLE_DIR}/clusters/sorted_tmp.fasta"
            DRAFT="${SAMPLE_DIR}/clusters/draft.fasta"

            if [[ "${POLISHER}" == "racon" ]]; then
                for ((i=1; i<=${RACON_ROUNDS}; i++)); do
                    minimap2 -x map-ont -t ${THREADS} "${DRAFT}" "${FILE_EXT}" \
                        > "${SAMPLE_DIR}/clusters/temp.paf" 2>/dev/null
                    if [ -s "${SAMPLE_DIR}/clusters/temp.paf" ]; then
                        racon -t ${THREADS} "${FILE_EXT}" "${SAMPLE_DIR}/clusters/temp.paf" "${DRAFT}" \
                            > "${SAMPLE_DIR}/clusters/round${i}.fasta" 2>/dev/null || true
                        if [ -s "${SAMPLE_DIR}/clusters/round${i}.fasta" ]; then
                            DRAFT="${SAMPLE_DIR}/clusters/round${i}.fasta"
                        else
                            log "WARNING" "Racon round ${i} produced empty output for ${CLUSTER_NAME}. Retaining previous draft."
                        fi
                    fi
                done
            else
                # Medaka polishing. medaka_consensus needs reads in FASTQ form.
                # Convert the cluster FASTA reads to a placeholder FASTQ for medaka input.
                MEDAKA_OUT="${SAMPLE_DIR}/clusters/medaka_${CLUSTER_NAME}"
                # medaka_consensus accepts fasta/fastq for -i; the FILE_EXT here is FASTA already.
                medaka_consensus -i "${FILE_EXT}" -d "${DRAFT}" \
                    -o "${MEDAKA_OUT}" -t ${THREADS} -m "${MEDAKA_MODEL}" \
                    >> "${LOG_FILE}" 2>&1 || log "WARNING" "Medaka failed for ${CLUSTER_NAME}; retaining draft."
                if [ -s "${MEDAKA_OUT}/consensus.fasta" ]; then
                    DRAFT="${MEDAKA_OUT}/consensus.fasta"
                fi
            fi

            seqkit replace -p ".+" -r "${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_Polished" "${DRAFT}" \
                > "${SAMPLE_DIR}/clusters/final_draft.fasta"
            cat "${SAMPLE_DIR}/clusters/final_draft.fasta" >> "${CURRENT_TOTAL}"
            cat "${SAMPLE_DIR}/clusters/final_draft.fasta" >> "${CURRENT_POLISHED}"
            POLISHED_COUNT=$((POLISHED_COUNT + 1))
        fi
    done < <(find "${SAMPLE_DIR}/clusters" -name "cluster_*" | sort)

    echo -e "${FILENAME}\t${RAW_READS}\t${CLEAN_READS}\t${CLUSTER_COUNT}\t${POLISHED_COUNT}\t${UNPOLISHED_COUNT}" >> "${STATS_FILE}"
    rm -rf "${SAMPLE_DIR}"

    # Per-sample timing for performance diagnostics
    SAMPLE_ELAPSED=$(( $(date +%s) - SAMPLE_START ))
    log "INFO" "[${SAMPLE_INDEX}/${TOTAL_FILES}] Done: ${FILENAME} (clusters=${CLUSTER_COUNT}, polished=${POLISHED_COUNT}, unpolished=${UNPOLISHED_COUNT}, ${SAMPLE_ELAPSED}s)"
    draw_progress "${SAMPLE_INDEX}" "${TOTAL_FILES}" "${FILENAME}"
done

# Total processing summary
PROCESSING_TOTAL_SEC=$(( $(date +%s) - PROGRESS_START_TIME ))
PROCESSING_MIN=$(( PROCESSING_TOTAL_SEC / 60 ))
PROCESSING_SEC=$(( PROCESSING_TOTAL_SEC % 60 ))
log "INFO" "All ${TOTAL_FILES} samples processed in ${PROCESSING_MIN}m${PROCESSING_SEC}s"

fi  # end of: if [[ -z "${FROM_STAGE}" ]]

# ================= Stage A: Reference-Based Orientation =================
if [[ -n "${REF_SEQ}" ]] && [[ "${_SKIP_STAGE_A:-0}" -eq 0 ]]; then
    export CURRENT_TOTAL CURRENT_POLISHED REF_SEQ
    log "INFO" "Aligning sequences to match the Reference forward strand..."

    blastn -query "${CURRENT_TOTAL}" -subject "${REF_SEQ}" \
        -outfmt "6 qseqid sstrand bitscore" \
        | sort -k1,1 -k3,3nr > "${OUTPUT_DIR}/orientation.txt"

    python3 << 'EOF'
import os
def rc(seq):
    trans = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')
    return seq.translate(trans)[::-1]

orient_file = os.path.join(os.environ['OUTPUT_DIR'], 'orientation.txt')
orient = {}
if os.path.exists(orient_file):
    with open(orient_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                if parts[0] not in orient: orient[parts[0]] = parts[1]

def process_fasta(input_fa, output_fa):
    if not os.path.exists(input_fa): return
    with open(input_fa, 'r') as fin, open(output_fa, 'w') as fout:
        name = ''
        seq = []
        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    final_seq = ''.join(seq)
                    if orient.get(name) == 'minus': final_seq = rc(final_seq)
                    fout.write(f">{name}\n{final_seq}\n")
                name = line[1:].split()[0]
                seq = []
            else: seq.append(line)
        if name:
            final_seq = ''.join(seq)
            if orient.get(name) == 'minus': final_seq = rc(final_seq)
            fout.write(f">{name}\n{final_seq}\n")

process_fasta(os.environ['CURRENT_TOTAL'], os.environ['CURRENT_TOTAL'] + ".tmp")
process_fasta(os.environ['CURRENT_POLISHED'], os.environ['CURRENT_POLISHED'] + ".tmp")
EOF
    mv "${CURRENT_TOTAL}.tmp" "${CURRENT_TOTAL}"
    mv "${CURRENT_POLISHED}.tmp" "${CURRENT_POLISHED}"
    rm -f "${OUTPUT_DIR}/orientation.txt"
    log "INFO" "Reference Orientation Complete."
fi

# ================= Stage B: Primer Trimming (Smart Fallback) =================
if [[ -n "${PRIMERS}" ]] && [[ "${_SKIP_STAGE_B:-0}" -eq 0 ]]; then
    export CURRENT_TOTAL CURRENT_POLISHED REF_SEQ
    log "INFO" "Primer file provided. Trimming target amplicons..."

    # Count input sequences before trimming for drop-rate reporting
    SEQS_BEFORE_TRIM=$(count_fasta_seqs "${CURRENT_TOTAL}")

    OUT_TOTAL_TRIMMED="${OUTPUT_DIR}/3_Total_Clusters_Trimmed.fasta"
    OUT_POLISHED_TRIMMED="${OUTPUT_DIR}/4_Polished_Clusters_Trimmed.fasta"
    > "${OUT_TOTAL_TRIMMED}"
    > "${OUT_POLISHED_TRIMMED}"

    sed 's/\r//g' "${PRIMERS}" | tr -d '"' > "${OUTPUT_DIR}/temp_primers.csv"

    python3 << 'EOF'
import sys, os
def rc(seq):
    trans = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')
    return seq.translate(trans)[::-1]

temp_in = os.path.join(os.environ['OUTPUT_DIR'], 'temp_primers.csv')
clean_out = os.path.join(os.environ['OUTPUT_DIR'], 'clean_primers.csv')

with open(temp_in, 'r') as fin, open(clean_out, 'w') as fout:
    lines = fin.readlines()
    if lines:
        fout.write(lines[0])
        for line in lines[1:]:
            line = line.strip()
            if not line: continue
            fout.write(line + '\n')
            parts = line.split(',')
            if len(parts) >= 3:
                new_id = parts[0] + '_Rev'
                new_fwd = rc(parts[2].strip())
                new_rev = rc(parts[1].strip())
                fout.write(f"{new_id},{new_fwd},{new_rev}\n")
EOF

    if [ -s "${CURRENT_TOTAL}" ]; then
        sed '1d' "${OUTPUT_DIR}/clean_primers.csv" | while IFS=, read -r P_ID P_FWD P_REV; do
            [ -z "$P_ID" ] && continue

            seqkit amplicon -F "${P_FWD}" -R "${P_REV}" -m ${MISMATCH} --region ${SEARCH_RANGE} "${CURRENT_TOTAL}" 2>/dev/null \
            | seqkit replace -p "(.+)" -r "${P_ID}_\$1" \
            | seqkit seq -g -m ${FINAL_MIN_LENGTH} -M ${FINAL_MAX_LENGTH} > "${OUTPUT_DIR}/temp_amp.fasta"

            if [[ "$P_ID" == *"_Rev" ]]; then
                TRUE_ID=${P_ID%_Rev}
                if [[ -z "${REF_SEQ}" ]]; then
                    # Fallback: No reference provided, flip using the Rev primer match
                    seqkit seq -r -p "${OUTPUT_DIR}/temp_amp.fasta" 2>/dev/null \
                        | seqkit replace -p "${P_ID}_" -r "${TRUE_ID}_" >> "${OUT_TOTAL_TRIMMED}"
                else
                    seqkit replace -p "${P_ID}_" -r "${TRUE_ID}_" "${OUTPUT_DIR}/temp_amp.fasta" >> "${OUT_TOTAL_TRIMMED}"
                fi
            else
                cat "${OUTPUT_DIR}/temp_amp.fasta" >> "${OUT_TOTAL_TRIMMED}"
            fi
        done
        seqkit rmdup -n "${OUT_TOTAL_TRIMMED}" -o "${OUTPUT_DIR}/tmp1.fasta" 2>/dev/null \
            && mv "${OUTPUT_DIR}/tmp1.fasta" "${OUT_TOTAL_TRIMMED}"
        rm -f "${OUTPUT_DIR}/temp_amp.fasta"
    fi

    if [ -s "${OUT_TOTAL_TRIMMED}" ]; then
        # Portable awk syntax (works in mawk and BSD awk, not just gawk)
        awk '/^>/ { keep = 0; if (/_Polished/) keep = 1 } keep { print }' \
            "${OUT_TOTAL_TRIMMED}" > "${OUT_POLISHED_TRIMMED}"
    fi

    # Drop-count reporting (v5.4)
    SEQS_AFTER_TRIM=$(count_fasta_seqs "${OUT_TOTAL_TRIMMED}")
    DROPPED=$(( SEQS_BEFORE_TRIM - SEQS_AFTER_TRIM ))
    if [ "${SEQS_BEFORE_TRIM}" -gt 0 ]; then
        DROP_PCT=$(awk -v d="$DROPPED" -v t="$SEQS_BEFORE_TRIM" 'BEGIN{ printf "%.1f", (d/t)*100 }')
    else
        DROP_PCT="0.0"
    fi
    log "INFO" "Primer trimming summary: ${SEQS_BEFORE_TRIM} input -> ${SEQS_AFTER_TRIM} retained, ${DROPPED} dropped (${DROP_PCT}%)."
    if [ "${DROPPED}" -gt 0 ] && [ "${SEQS_BEFORE_TRIM}" -gt 0 ]; then
        log "INFO" "Note: dropped sequences did not match any primer pair within --mismatch=${MISMATCH} and --search-range=${SEARCH_RANGE}."
    fi

    export CURRENT_TOTAL="${OUT_TOTAL_TRIMMED}"
    export CURRENT_POLISHED="${OUT_POLISHED_TRIMMED}"
    log "INFO" "Trimming complete."
    rm -f "${OUTPUT_DIR}/temp_primers.csv" "${OUTPUT_DIR}/clean_primers.csv"
fi

# ================= Stage D: Best-Per-Locus Extraction =================
# For every (primer_locus, sample) pair in the trimmed polished FASTA,
# retain only the sequence with the highest Reads=N value.
# This gives one representative consensus per locus per sample —
# the best-supported sequence for downstream genotyping.
#
# Header format expected (produced by Stage B):
#   >PrimerID_SampleName_cluster_N_Reads=N_Polished
#
# Outputs:
#   5_Best_Per_Locus.fasta       — one sequence per (locus, sample) pair
#   5_Best_Per_Locus_summary.tsv — QC table: Sample, Locus, Cluster,
#                                   Reads, Sequence_Length, [INDELs_Fixed]
#
# Stage D is skipped when Stage B was not run (no -p flag), because
# without primer trimming, sequences carry no locus-prefix label.
#
# Optional: --indel-fix enables experimental reference-anchored INDEL
# correction AFTER best-per-locus selection. See help text for WARNING.
if [[ -n "${PRIMERS}" ]] || [[ "${FROM_STAGE}" == "best" ]] || [[ "${FROM_STAGE}" == "trimmed" ]]; then
    log "INFO" "Running Stage D: best-per-locus sequence extraction..."
    if [[ "${INDEL_FIX}" -eq 1 ]]; then
        log "WARNING" "--indel-fix is ACTIVE. Consensus sequences will be modified using reference-anchored INDEL correction. This may introduce reference bias. Do not use for variant discovery."
    fi
    export BEST_FASTA="$(pwd)/${OUTPUT_DIR}/5_Best_Per_Locus.fasta"
    export BEST_TSV="$(pwd)/${OUTPUT_DIR}/5_Best_Per_Locus_summary.tsv"
    export QUERY_FOR_BEST="${CURRENT_POLISHED}"
    export INDEL_FIX_EXPORT="${INDEL_FIX}"
    export REF_SEQ_EXPORT="${REF_SEQ}"

    python3 << 'EOF'
import os
import re
import subprocess
import tempfile

query_file   = os.environ.get('QUERY_FOR_BEST', '')
out_fasta    = os.environ.get('BEST_FASTA', '')
out_tsv      = os.environ.get('BEST_TSV', '')
indel_fix    = os.environ.get('INDEL_FIX_EXPORT', '0') == '1'
ref_file     = os.environ.get('REF_SEQ_EXPORT', '')
out_dir      = os.environ.get('OUTPUT_DIR', '.')

# ── Parsing helpers ──────────────────────────────────────────────────────────

def parse_reads(header):
    """Extract integer read count from 'Reads=N' tag in header. Returns 0 if absent."""
    m = re.search(r'Reads=(\d+)', header)
    return int(m.group(1)) if m else 0

def parse_locus(header):
    """
    Extract locus ID = first token before the first underscore.
    e.g. 'M3613_1-A-10_TSM..._cluster_9_Reads=108_Polished' -> 'M3613'
    """
    return header.split('_')[0]

def parse_sample(header):
    """
    Extract sample name = everything between locus_id and _cluster_.
    Uses rfind to handle sample names that contain the word 'cluster'.
    e.g. 'M3613_1-A-10_TSM..._cluster_9_Reads=108_Polished'
         -> '1-A-10_TSM...'
    """
    after_locus = header.split('_', 1)[1] if '_' in header else header
    idx = after_locus.rfind('_cluster_')
    if idx >= 0:
        return after_locus[:idx]
    return after_locus

def parse_cluster(header):
    """Extract cluster number from _cluster_N_Reads=N tag."""
    m = re.search(r'_cluster_(\d+)_', header)
    return m.group(1) if m else 'unknown'

# ── Read all entries from the trimmed polished FASTA ────────────────────────
entries = []
if not os.path.exists(query_file) or os.path.getsize(query_file) == 0:
    print(f"[Stage D] No trimmed polished FASTA found or file is empty: {query_file}")
    print("[Stage D] Skipping best-per-locus extraction.")
    open(out_fasta, 'w').close()
    with open(out_tsv, 'w') as f:
        header_cols = "Sample\tLocus\tCluster\tReads\tSequence_Length"
        if indel_fix:
            header_cols += "\tINDELs_Fixed"
        f.write(header_cols + "\tHeader\n")
    raise SystemExit(0)

with open(query_file) as f:
    current_header = None
    seq_lines = []
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if current_header is not None:
                entries.append((current_header, seq_lines))
            current_header = line[1:].split()[0]
            seq_lines = []
        elif line:
            seq_lines.append(line)
    if current_header is not None:
        entries.append((current_header, seq_lines))

# ── Group by (locus, sample) keeping all clusters per group ─────────────────
grouped = {}   # (locus, sample) -> list of (reads, header, seq_str)
for header, seq_lines in entries:
    locus  = parse_locus(header)
    sample = parse_sample(header)
    reads  = parse_reads(header)
    seq    = ''.join(seq_lines)
    key    = (locus, sample)
    grouped.setdefault(key, []).append((reads, header, seq))

# Sort each group descending by reads; the first entry is the major consensus
for key in grouped:
    grouped[key].sort(key=lambda x: x[0], reverse=True)

# ── Optional INDEL correction ─────────────────────────────────────────────────
# Runs ONLY when --indel-fix is active.
#
# Strategy (Bug-fixed version of v5.11):
#   1. Batch ALL sequences (major + minors from every group) into one FASTA
#      and run a SINGLE blastn call against the reference (Bug 3 fix).
#   2. For each group, parse the major and minor alignments.
#   3. Build a per-reference-position map from minor clusters showing which
#      positions have MATCH support.
#   4. Walk the major alignment column by column, maintaining ONE consistent
#      r_pos counter (Bug 1+2 fix: r_pos advances only on reference-consuming
#      columns — matches and deletions — not on insertions).
#   5. Patch only positions where the major has an INDEL that is NOT supported
#      by any minor cluster.
#   6. Report all corrections in the summary TSV (Bug 5 fix).

alignments = {}   # header -> {qseq, sseq, qstart, qend, sstart, send}

if indel_fix and ref_file and os.path.exists(ref_file):
    # ── Bug 4 fix: use a dedicated temp directory inside OUTPUT_DIR ──────────
    tmp_dir = os.path.join(out_dir, 'temp_stageD')
    os.makedirs(tmp_dir, exist_ok=True)
    batch_fasta = os.path.join(tmp_dir, 'all_sequences.fa')

    with open(batch_fasta, 'w') as bf:
        for key, cluster_list in grouped.items():
            for reads, header, seq in cluster_list:
                bf.write(f">{header}\n{seq}\n")

    # ── Bug 3 fix: one batched blastn call for all sequences ─────────────────
    cmd = [
        'blastn',
        '-query',   batch_fasta,
        '-subject', ref_file,
        '-outfmt',  '6 qseqid sseqid qseq sseq qstart qend sstart send',
        '-max_hsps', '1',      # keep only top HSP per query
        '-max_target_seqs', '1'
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)

    for line in p.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 8:
            continue
        qid = parts[0]
        if qid not in alignments:   # first hit wins (best score)
            alignments[qid] = {
                'qseq':   parts[2],
                'sseq':   parts[3],
                'qstart': int(parts[4]),
                'qend':   int(parts[5]),
                'sstart': int(parts[6]),
                'send':   int(parts[7]),
            }

def apply_indel_fix(major_head, major_seq, minors, alignments):
    """
    Surgically patch the major consensus sequence at INDEL positions that
    are not supported by any minor cluster.

    Returns (patched_seq, status_string).

    r_pos tracks the current reference position. It advances by 1 for every
    reference-consuming column (match or deletion in the query), and does NOT
    advance for insertion columns (where the reference has '-').
    This is the standard pairwise-alignment reference-position convention
    and fixes Bugs 1 and 2 from v5.11.
    """
    if major_head not in alignments:
        # Bug 5 fix: distinguish BLAST-miss from no-minors case
        return major_seq, "Skipped:no_BLAST_hit"

    if not minors:
        return major_seq, "Skipped:no_minor_clusters"

    maj_aln = alignments[major_head]

    # Build reference-position support map from all minor clusters
    # Key: int r_pos for match/deletion support; str f"{r_pos}_INS" for ins.
    minor_ref_support = {}   # r_pos -> set of ('MATCH', 'DEL', 'INS')
    for _, m_head, _ in minors:
        if m_head not in alignments:
            continue
        m_aln = alignments[m_head]
        r_pos = m_aln['sstart']
        for mq, ms in zip(m_aln['qseq'], m_aln['sseq']):
            if ms != '-':               # reference-consuming column
                if mq == ms:
                    minor_ref_support.setdefault(r_pos, set()).add('MATCH')
                elif mq == '-':
                    minor_ref_support.setdefault(r_pos, set()).add('DEL')
                r_pos += 1              # advance ref pos only here
            else:                       # insertion in query (ref has gap)
                ins_key = f"{r_pos}_INS"
                minor_ref_support.setdefault(ins_key, set()).add('INS')
                # r_pos does NOT advance for insertions

    # Walk the major alignment and patch
    patched_seq = []
    r_pos = maj_aln['sstart']
    fixed_del = 0
    fixed_ins = 0

    for q_base, s_base in zip(maj_aln['qseq'], maj_aln['sseq']):
        if s_base != '-':               # reference-consuming column
            if q_base == s_base:        # match
                patched_seq.append(q_base)
            elif q_base == '-':         # deletion in major consensus
                # Fix only if at least one minor cluster shows MATCH here
                if 'MATCH' in minor_ref_support.get(r_pos, set()):
                    patched_seq.append(s_base)  # restore reference base
                    fixed_del += 1
                # else: keep the deletion (don't append anything)
            else:                       # SNP — never alter
                patched_seq.append(q_base)
            r_pos += 1                  # advance ref pos for all ref-consuming cols

        else:                           # insertion in major (ref has '-')
            ins_key = f"{r_pos}_INS"
            # Keep insertion only if at least one minor also shows it
            if 'INS' in minor_ref_support.get(ins_key, set()):
                patched_seq.append(q_base)
            else:
                fixed_ins += 1          # drop the inserted base
            # r_pos does NOT advance (reference did not consume a position)

    if fixed_del == 0 and fixed_ins == 0:
        return major_seq, "None"

    # Reconstruct: unaligned left flank + patched aligned region + right flank
    left_clip  = major_seq[:maj_aln['qstart'] - 1]
    right_clip = major_seq[maj_aln['qend']:]
    patched    = left_clip + ''.join(patched_seq) + right_clip
    status     = f"Fixed:{fixed_del}DEL,{fixed_ins}INS"
    return patched, status

# ── Select best and optionally correct ───────────────────────────────────────
out_data = []
sorted_keys = sorted(grouped.keys(), key=lambda k: (k[1], k[0]))

for (locus, sample) in sorted_keys:
    cluster_list = grouped[(locus, sample)]
    major_reads, major_head, major_seq = cluster_list[0]
    minors = cluster_list[1:]

    if indel_fix:
        final_seq, fix_status = apply_indel_fix(
            major_head, major_seq, minors, alignments
        )
    else:
        final_seq  = major_seq
        fix_status = None    # column not written when indel_fix is off

    out_data.append({
        'header':  major_head,
        'seq':     final_seq,
        'sample':  sample,
        'locus':   locus,
        'cluster': parse_cluster(major_head),
        'reads':   major_reads,
        'fixed':   fix_status,
    })

# ── Write outputs ─────────────────────────────────────────────────────────────
n_written = 0
with open(out_fasta, 'w') as ff, open(out_tsv, 'w') as ft:
    tsv_header = "Sample\tLocus\tCluster\tReads\tSequence_Length"
    if indel_fix:
        tsv_header += "\tINDELs_Fixed"
    ft.write(tsv_header + "\tHeader\n")

    for d in out_data:
        seq = d['seq']
        ff.write(f">{d['header']}\n")
        for i in range(0, len(seq), 60):
            ff.write(seq[i:i+60] + '\n')

        row = (f"{d['sample']}\t{d['locus']}\t{d['cluster']}\t"
               f"{d['reads']}\t{len(seq)}")
        if indel_fix:
            row += f"\t{d['fixed']}"
        row += f"\t{d['header']}"
        ft.write(row + '\n')
        n_written += 1

# ── Bug 4 fix: clean up the temp directory ───────────────────────────────────
if indel_fix:
    import shutil
    tmp_dir_path = os.path.join(out_dir, 'temp_stageD')
    if os.path.exists(tmp_dir_path):
        shutil.rmtree(tmp_dir_path, ignore_errors=True)

print(f"[Stage D] {len(entries)} sequences -> {n_written} best-per-locus sequences written.")
print(f"[Stage D] Unique loci:    {len(set(k[0] for k in grouped))}")
print(f"[Stage D] Unique samples: {len(set(k[1] for k in grouped))}")
if indel_fix:
    n_fixed = sum(1 for d in out_data
                  if d['fixed'] and d['fixed'] not in ('None','Skipped:no_BLAST_hit',
                                                         'Skipped:no_minor_clusters'))
    print(f"[Stage D] INDEL-fix corrections applied to: {n_fixed} sequences")
EOF

    BEST_COUNT=$(count_fasta_seqs "${BEST_FASTA}")
    log "INFO" "Stage D complete: ${BEST_COUNT} best-per-locus sequences written to 5_Best_Per_Locus.fasta"
    log "INFO" "Stage D summary table: 5_Best_Per_Locus_summary.tsv"
else
    log "INFO" "Stage D skipped: no primer file provided (-p). Re-run with -p <primers.csv> to enable."
fi

# ================= Stage C: BLAST Variant Calling =================
if [[ -n "${VARIANT_TABLE}" ]]; then
    log "INFO" "Executing Variant Calling HTML Generation..."
    export VARIANT_TABLE
    export QUERY_FASTA="${CURRENT_POLISHED}"

    python3 << 'EOF'
import sys, os, subprocess
import pandas as pd

ref_file = os.environ.get('REF_SEQ')
table_file = os.environ.get('VARIANT_TABLE')
query_file = os.environ.get('QUERY_FASTA')
out_html = os.environ.get('VARIANT_HTML')
out_dir = os.environ.get('OUTPUT_DIR')

tmp_ref = os.path.join(out_dir, 'tmp_ref.fa')
tmp_q = os.path.join(out_dir, 'tmp_q.fa')

def read_fasta(file_path):
    seqs = {}
    if not os.path.exists(file_path): return seqs
    with open(file_path, 'r') as f:
        name = ''
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name: seqs[name] = ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else: seq.append(line)
        if name: seqs[name] = ''.join(seq)
    return seqs

try:
    df_map = pd.read_csv(table_file, sep=r'\s+', header=None, dtype=str)
    ref_seqs = read_fasta(ref_file)
    query_seqs = read_fasta(query_file)
    mapping_pairs = []
    missing_pairs = []  # (ref_name, sample_partial) with no matching cluster

    for _, row in df_map.iterrows():
        ref_name = str(row[0]).strip()
        if ref_name not in ref_seqs: continue
        q_partials = [q.strip() for q in str(row[1]).split(',') if q.strip()]
        # Restrict candidate clusters to those whose header starts with
        # "<ref_name>_" — Stage B trimming prepends the primer ID, so
        # Locus01-trimmed clusters look like "Locus01_<sample>_cluster_...".
        # Without this prefix check, a substring like "sample01" would also
        # match Locus02_sample01_*, Locus03_sample01_*, etc. and the first hit
        # would be returned regardless of which locus row we're processing.
        ref_prefix = ref_name + "_"
        candidates = [q for q in query_seqs.keys() if q.startswith(ref_prefix)]
        for partial in q_partials:
            matches = [q for q in candidates if partial in q]
            if matches:
                # Pick the highest-depth cluster if multiple match (heuristic:
                # 'Reads=N' tag in header). Fall back to first if not parseable.
                def depth(name):
                    try:
                        return int(name.split('Reads=')[1].split('_')[0])
                    except (IndexError, ValueError):
                        return 0
                matches.sort(key=depth, reverse=True)
                mapping_pairs.append((ref_name, matches[0]))
            else:
                missing_pairs.append((ref_name, partial))

    results = []
    # Emit a row for every (locus, sample) pairing that produced no cluster
    # so users can see which loci failed to amplify or were dropped.
    for ref_name, partial in missing_pairs:
        results.append({
            'Reference_Name': ref_name,
            'Target_Polished_Sample': f"(no cluster matching '{partial}')",
            'Identity': 'No Cluster',
            'Variants_Detected': 'Locus not recovered for this sample'
        })

    for ref_name, q_name in mapping_pairs:
        with open(tmp_ref, 'w') as f: f.write(f">{ref_name}\n{ref_seqs[ref_name]}\n")
        with open(tmp_q, 'w') as f: f.write(f">{q_name}\n{query_seqs[q_name]}\n")

        cmd = ['blastn', '-query', tmp_q, '-subject', tmp_ref, '-outfmt', '6 qseqid sseqid qstart qend sstart send pident btop']
        p = subprocess.run(cmd, capture_output=True, text=True)
        blast_lines = p.stdout.strip().split('\n')

        if not p.stdout.strip():
            results.append({'Reference_Name': ref_name, 'Target_Polished_Sample': q_name, 'Identity': 'No Alignment', 'Variants_Detected': 'No Alignment Found'})
            continue

        fields = blast_lines[0].split('\t')
        sstart = int(fields[4])
        send = int(fields[5])
        pident = fields[6] + "%"
        btop = fields[7]

        ref_pos = sstart
        step = 1 if sstart < send else -1

        # Robust character-iteration BTOP parser
        tokens = []
        idx = 0
        while idx < len(btop):
            if btop[idx].isdigit():
                num_str = ""
                while idx < len(btop) and btop[idx].isdigit():
                    num_str += btop[idx]
                    idx += 1
                tokens.append(num_str)
            else:
                tokens.append(btop[idx:idx+2])
                idx += 2

        variants = []
        for token in tokens:
            if token.isdigit():
                ref_pos += (int(token) * step)
            else:
                q_base, r_base = token[0], token[1]
                if q_base == '-':
                    variants.append(f"DEL [Pos {ref_pos}]: -{r_base.upper()}")
                    ref_pos += step
                elif r_base == '-':
                    variants.append(f"INS [Pos {ref_pos}]: +{q_base.upper()}")
                else:
                    variants.append(f"SNP [Pos {ref_pos}]: {r_base.upper()} &rarr; {q_base.upper()}")
                    ref_pos += step

        var_str = "<br>".join(variants) if variants else "<i>Perfect Match (0 Variants)</i>"
        results.append({'Reference_Name': ref_name, 'Target_Polished_Sample': q_name, 'Identity': pident, 'Variants_Detected': var_str})

    if os.path.exists(tmp_ref): os.remove(tmp_ref)
    if os.path.exists(tmp_q): os.remove(tmp_q)

    if results:
        df_res = pd.DataFrame(results)
        html_content = f"""
        <!DOCTYPE html>
        <html><head><style>
            body {{ font-family: sans-serif; margin: 40px; color: #333; }}
            table {{ border-collapse: collapse; width: 100%; margin-top: 20px; text-align: left; }}
            th, td {{ border: 1px solid #BDC3C7; padding: 12px; vertical-align: top; }}
            th {{ background-color: #ECF0F1; }}
        </style></head><body>
            <h1>NanoPolite Variant Discovery Report</h1>
            {df_res.to_html(index=False, classes='table', escape=False)}
        </body></html>
        """
        with open(out_html, 'w') as f: f.write(html_content)
except Exception as e:
    print(f"Variant Calling Error: {e}")
EOF
fi

# ================= Reporting =================
log "INFO" "Generating main HTML report..."
export END_TIME=$(date +"%Y-%m-%d %H:%M:%S")

python3 << 'EOF'
import pandas as pd
import os

try:
    df = pd.read_csv(os.environ['STATS_FILE'], sep='\t')
    total_raw = df['Raw_Reads'].sum()
    total_clean = df['Clean_Reads'].sum()
    pass_rate = (total_clean / total_raw * 100) if total_raw > 0 else 0

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8"><title>NanoPolite Report</title>
        <style>
            body {{ font-family: sans-serif; margin: 40px; color: #333; }}
            h1 {{ color: #2C3E50; border-bottom: 2px solid #3498DB; padding-bottom: 10px; }}
            table {{ border-collapse: collapse; width: 100%; margin-top: 20px; text-align: left; }}
            th, td {{ border: 1px solid #BDC3C7; padding: 12px; }}
            th {{ background-color: #ECF0F1; }}
            .summary-box {{ background-color: #F8F9F9; padding: 20px; border-left: 5px solid #3498DB; margin-bottom: 30px; }}
        </style>
    </head>
    <body>
        <h1>NanoPolite Analysis Report</h1>
        <div class="summary-box">
            <h2>Run Summary</h2>
            <p><b>Period:</b> {os.environ['START_TIME']} to {os.environ['END_TIME']}</p>
            <p><b>Total Samples:</b> {len(df)} | <b>Average Retention Rate:</b> {pass_rate:.2f}%</p>
        </div>
        <h2>Sample Detailed Statistics</h2>
        {df.to_html(index=False, classes='table')}
    </body>
    </html>
    """
    with open(os.environ['REPORT_HTML'], 'w') as f:
        f.write(html_content)

except Exception as e:
    print(f"Reporting Error: {e}")
EOF

log "INFO" "🎉 Analysis Complete. 📂 Output Folder: ${OUTPUT_DIR}"
