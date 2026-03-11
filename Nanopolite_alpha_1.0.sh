#!/bin/bash

# ==============================================================================
# Script Name: Nanopolite_alpha_1.0.sh
# Feature: Exposed --amp-min and --amp-max for strict amplicon length filtering.
# ==============================================================================

# ================= Initialization & Parameters =================
export START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

INPUT_LIST=()
INPUT_DIR=""
OUTPUT_DIR="Nanopolite_results_${TIMESTAMP}"
PRIMERS=""
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
MISMATCH=3
SEARCH_RANGE="1:-1"

# ================= Help / Usage Information =================
usage() {
cat << EOF
Nanopolite v3.2 - Automated Nanopore Amplicon Processing Pipeline

Required Input (At least one):
  -i, --input        Input fastq file(s). Supports wildcards (e.g., -i "*.fastq")
  -d, --dir          Directory containing fastq files.

Output Options:
  -o, --out          Designate output folder path (Default: Nanopolite_results_TIMESTAMP)

Optional Primer & Trimming Parameters:
  -p, --primers      Primer CSV file path (Format: ID,Forward,Reverse)
  --amp-min          Min length for retained trimmed amplicons (Default: ${FINAL_MIN_LENGTH} bp)
  --amp-max          Max length for retained trimmed amplicons (Default: ${FINAL_MAX_LENGTH} bp)
  --mismatch         Allowed mismatches for primer binding (Default: ${MISMATCH})

Optional QC & Clustering Parameters:
  -q, --quality      Min average Phred quality score (Default: ${MIN_QUAL})
  -m, --min-len      Min read length for QC (Default: ${MIN_LEN} bp)
  --id               vsearch clustering identity (Default: ${ID_THRESHOLD})
  --min-cluster      Min reads for Racon polishing (Default: ${MIN_CLUSTER_SIZE})
  -t, --threads      Number of CPU threads (Default: ${THREADS})
  -h, --help         Display this help message and exit
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
        -t|--threads) THREADS="$2"; shift ;;
        -h|--help) usage ;;
        *) echo -e "[ERROR] Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# ================= Setup Output Structure =================
mkdir -p "${OUTPUT_DIR}"
export LOG_FILE="${OUTPUT_DIR}/nanopolite_run_${TIMESTAMP}.log"
export STATS_FILE="$(pwd)/${OUTPUT_DIR}/nanopolite_stats_${TIMESTAMP}.tsv"
export REPORT_HTML="$(pwd)/${OUTPUT_DIR}/Nanopolite_Report_${TIMESTAMP}.html"

OUT_TOTAL_UNTRIMMED="${OUTPUT_DIR}/1_Total_Clusters_Untrimmed.fasta"
OUT_POLISHED_UNTRIMMED="${OUTPUT_DIR}/2_Polished_Clusters_Untrimmed.fasta"

> "${OUT_TOTAL_UNTRIMMED}"
> "${OUT_POLISHED_UNTRIMMED}"

# ================= Logging Function =================
log() {
    local TYPE=$1
    shift
    local MSG="$@"
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [${TYPE}] ${MSG}" | tee -a "${LOG_FILE}"
}

# ================= Pre-run Checks =================
log "INFO" ">>> Initializing Nanopolite v3.2 <<<"
log "INFO" "Output directory: ${OUTPUT_DIR}"

for cmd in NanoFilt vsearch minimap2 racon seqkit python3 awk; do
    command -v $cmd &> /dev/null || { log "ERROR" "Missing dependency: $cmd"; exit 1; }
done

export HAS_WEASYPRINT="true"
if ! python3 -c "import weasyprint" &> /dev/null; then
    export HAS_WEASYPRINT="false"
fi

FILES=()
for f in "${INPUT_LIST[@]}"; do [ -f "$f" ] && FILES+=("$f"); done
if [ -d "$INPUT_DIR" ]; then
    for f in "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq; do [ -f "$f" ] && FILES+=("$f"); done
fi
FILES=($(echo "${FILES[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

if [ ${#FILES[@]} -eq 0 ]; then
    log "ERROR" "No valid fastq files found via -i or -d."; exit 1
fi

echo -e "Sample\tRaw_Reads\tClean_Reads\tClusters\tPolished_Seqs\tUnpolished_Seqs" > "${STATS_FILE}"

# ================= Core Processing Logic =================
log "INFO" "Processing ${#FILES[@]} files..."

for INPUT_FILE in "${FILES[@]}"; do
    FILENAME=$(basename "${INPUT_FILE%.*}")
    SAMPLE_DIR="${OUTPUT_DIR}/temp_work_${FILENAME}"
    log "INFO" "------------------------------------------------"
    log "INFO" "Processing: ${FILENAME}"
    
    RAW_READS=$(($(wc -l < "$INPUT_FILE") / 4))
    rm -rf "${SAMPLE_DIR}"; mkdir -p "${SAMPLE_DIR}/clean_data" "${SAMPLE_DIR}/clusters"
    
    NanoFilt -q ${MIN_QUAL} -l ${MIN_LEN} --headcrop ${HEAD_CROP} --tailcrop ${TAIL_CROP} < "$INPUT_FILE" > "${SAMPLE_DIR}/clean_data/cleaned.fastq"
    
    if [ ! -s "${SAMPLE_DIR}/clean_data/cleaned.fastq" ]; then
        echo -e "${FILENAME}\t${RAW_READS}\t0\t0\t0\t0" >> "${STATS_FILE}"
        rm -rf "${SAMPLE_DIR}"; continue
    fi
    
    CLEAN_READS=$(($(wc -l < "${SAMPLE_DIR}/clean_data/cleaned.fastq") / 4))
    vsearch --cluster_fast "${SAMPLE_DIR}/clean_data/cleaned.fastq" --id ${ID_THRESHOLD} --strand both --clusters "${SAMPLE_DIR}/clusters/cluster_" --minseqlength ${MIN_LEN} --threads ${THREADS} --quiet
    
    CLUSTER_COUNT=$(find "${SAMPLE_DIR}/clusters" -name "cluster_*" | wc -l)
    POLISHED_COUNT=0; UNPOLISHED_COUNT=0

    while read -r CLUSTER_FILE; do
        mv "${CLUSTER_FILE}" "${CLUSTER_FILE}.fasta"
        FILE_EXT="${CLUSTER_FILE}.fasta"
        CLUSTER_NAME=$(basename "${FILE_EXT}" .fasta)
        READ_COUNT=$(grep -c "^>" "${FILE_EXT}" || true)

        if [ "${READ_COUNT}" -lt "${MIN_CLUSTER_SIZE}" ]; then
            seqkit sort -l -r "${FILE_EXT}" | seqkit head -n 1 | seqkit replace -p ".+" -r "${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_Unpolished" >> "${OUT_TOTAL_UNTRIMMED}"
            ((UNPOLISHED_COUNT++))
        else
            seqkit sort -l -r "${FILE_EXT}" | seqkit head -n 1 > "${SAMPLE_DIR}/clusters/draft.fasta"
            DRAFT="${SAMPLE_DIR}/clusters/draft.fasta"
            
            for ((i=1; i<=${RACON_ROUNDS}; i++)); do
                minimap2 -x map-ont -t ${THREADS} "${DRAFT}" "${FILE_EXT}" > "${SAMPLE_DIR}/clusters/temp.paf" 2>/dev/null
                if [ -s "${SAMPLE_DIR}/clusters/temp.paf" ]; then
                    racon -t ${THREADS} "${FILE_EXT}" "${SAMPLE_DIR}/clusters/temp.paf" "${DRAFT}" > "${SAMPLE_DIR}/clusters/round${i}.fasta" 2>/dev/null
                    [ -s "${SAMPLE_DIR}/clusters/round${i}.fasta" ] && DRAFT="${SAMPLE_DIR}/clusters/round${i}.fasta"
                fi
            done
            
            seqkit replace -p ".+" -r "${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_Polished" "${DRAFT}" > "${SAMPLE_DIR}/clusters/final_draft.fasta"
            
            cat "${SAMPLE_DIR}/clusters/final_draft.fasta" >> "${OUT_TOTAL_UNTRIMMED}"
            cat "${SAMPLE_DIR}/clusters/final_draft.fasta" >> "${OUT_POLISHED_UNTRIMMED}"
            ((POLISHED_COUNT++))
        fi
    done < <(find "${SAMPLE_DIR}/clusters" -name "cluster_*" | sort)
    
    echo -e "${FILENAME}\t${RAW_READS}\t${CLEAN_READS}\t${CLUSTER_COUNT}\t${POLISHED_COUNT}\t${UNPOLISHED_COUNT}" >> "${STATS_FILE}"
    rm -rf "${SAMPLE_DIR}"
done

# ================= Primer Extraction =================
if [ -n "${PRIMERS}" ]; then
    log "INFO" "Primer file provided. Extracting target regions..."
    log "INFO" "Length Filtering Applied: Min ${FINAL_MIN_LENGTH} bp | Max ${FINAL_MAX_LENGTH} bp"
    
    OUT_TOTAL_TRIMMED="${OUTPUT_DIR}/3_Total_Clusters_Trimmed.fasta"
    OUT_POLISHED_TRIMMED="${OUTPUT_DIR}/4_Polished_Clusters_Trimmed.fasta"
    > "${OUT_TOTAL_TRIMMED}"
    > "${OUT_POLISHED_TRIMMED}"

    sed 's/\r//g' "${PRIMERS}" | tr -d '"' > "${OUTPUT_DIR}/clean_primers.csv"
    
    # Generate File 3 (Total Trimmed) with strict Length Filtering (-m and -M)
    if [ -s "${OUT_TOTAL_UNTRIMMED}" ]; then
        sed '1d' "${OUTPUT_DIR}/clean_primers.csv" | while IFS=, read -r P_ID P_FWD P_REV; do
            [ -z "$P_ID" ] && continue
            seqkit amplicon -F "${P_FWD}" -R "${P_REV}" -m ${MISMATCH} --region ${SEARCH_RANGE} "${OUT_TOTAL_UNTRIMMED}" 2>/dev/null \
            | seqkit replace -p "(.+)" -r "${P_ID}_\$1" \
            | seqkit seq -g -m ${FINAL_MIN_LENGTH} -M ${FINAL_MAX_LENGTH} >> "${OUT_TOTAL_TRIMMED}"
        done
        seqkit rmdup -n "${OUT_TOTAL_TRIMMED}" -o "${OUTPUT_DIR}/tmp1.fasta" 2>/dev/null && mv "${OUTPUT_DIR}/tmp1.fasta" "${OUT_TOTAL_TRIMMED}"
    fi

    # Generate File 4 (Polished Trimmed) via explicit AWK extraction from File 3
    log "INFO" "Extracting strictly Polished sequences from Trimmed pool..."
    if [ -s "${OUT_TOTAL_TRIMMED}" ]; then
        awk '/^>/ {keep = /_Polished/ ? 1 : 0} keep {print}' "${OUT_TOTAL_TRIMMED}" > "${OUT_POLISHED_TRIMMED}"
    fi
    
    log "INFO" "Trimming complete."
    rm -f "${OUTPUT_DIR}/clean_primers.csv"
else
    log "INFO" "No primer file provided. Skipping trimming steps."
fi

# ================= Reporting =================
log "INFO" "Generating HTML report..."
export END_TIME=$(date +"%Y-%m-%d %H:%M:%S")

python3 << 'EOF'
import pandas as pd
import os, sys

try:
    df = pd.read_csv(os.environ['STATS_FILE'], sep='\t')
    total_raw = df['Raw_Reads'].sum()
    total_clean = df['Clean_Reads'].sum()
    pass_rate = (total_clean / total_raw * 100) if total_raw > 0 else 0
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8"><title>Nanopolite Report</title>
        <style>
            body {{ font-family: sans-serif; margin: 40px; color: #333; }}
            h1 {{ color: #2C3E50; border-bottom: 2px solid #3498DB; padding-bottom: 10px; }}
            table {{ border-collapse: collapse; width: 100%; margin-top: 20px; text-align: left; }}
            th, td {{ border: 1px solid #BDC3C7; padding: 12px; }}
            th {{ background-color: #ECF0F1; }}
            .summary-box {{ background-color: #F8F9F9; padding: 20px; border-left: 5px solid #3498DB; margin-bottom: 30px; }}
            .warning-box {{ background-color: #FFF3CD; padding: 10px; border: 1px solid #FFEEBA; margin-bottom: 20px; display: {'none' if os.environ['HAS_WEASYPRINT'] == 'true' else 'block'}; }}
        </style>
    </head>
    <body>
        <h1>Nanopolite Analysis Report</h1>
        <div class="warning-box">&#9888; Note: PDF generation skipped (weasyprint missing).</div>
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

log "INFO" "🎉 Analysis Complete."
log "INFO" "📂 Output Folder: ${OUTPUT_DIR}"
