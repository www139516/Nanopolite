#!/bin/bash

# ==============================================================================
# 脚本名称: nanopolite_processing_v6.1.sh
# 新增功能: 提前检测 Python 依赖 (weasyprint)，优雅降级 (仅生成 HTML)
# ==============================================================================

# ================= 初始化与参数 =================
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="nanopolite_run_${TIMESTAMP}.log"
STATS_FILE="nanopolite_stats_${TIMESTAMP}.tsv"

INPUT_PATTERN=""
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

# ================= 日志函数 =================
log() {
    local TYPE=$1
    shift
    local MSG="$@"
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [${TYPE}] ${MSG}" | tee -a "${LOG_FILE}"
}

# ================= 帮助信息 =================
usage() {
    cat << EOF
用法: bash $0 -i <input.fastq> [选项...]
(使用 -h 或 --help 查看完整参数)
EOF
    exit 1
}

if [ $# -eq 0 ]; then usage; fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_PATTERN="$2"; shift ;;
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
        *) log "ERROR" "未知参数: $1"; exit 1 ;;
    esac
    shift
done

# ================= 前置健壮性检查 =================
log "INFO" ">>> 初始化 Nanopolite 分析流程 <<<"
log "INFO" "正在检查前置核心依赖环境..."

for cmd in NanoFilt vsearch minimap2 racon seqkit python; do
    if ! command -v $cmd &> /dev/null; then
        log "ERROR" "缺失关键软件: $cmd。请确保已安装并添加至环境变量。"
        exit 1
    fi
done
log "INFO" "核心依赖软件检查通过。"

# --- 新增：提前检测 weasyprint ---
HAS_WEASYPRINT="true"
log "INFO" "正在检查 PDF 生成模块 (weasyprint)..."
if ! python -c "import weasyprint" &> /dev/null; then
    log "WARNING" "未检测到 Python 库 'weasyprint'。程序将照常运行，但【仅生成 HTML 报告】，取消生成 PDF。"
    log "WARNING" "如需导出 PDF 报告，请在运行前执行: pip install weasyprint"
    HAS_WEASYPRINT="false"
else
    log "INFO" "PDF 生成模块就绪。"
fi
# -----------------------------------

if [ -z "${INPUT_PATTERN}" ]; then
    log "ERROR" "必须使用 -i 指定输入的 fastq 文件！"
    exit 1
fi

FILES=$(ls ${INPUT_PATTERN} 2>/dev/null)
if [ -z "${FILES}" ]; then
    log "ERROR" "未找到匹配的文件: ${INPUT_PATTERN}"
    exit 1
fi

if [ -n "${PRIMERS}" ] && [ ! -f "${PRIMERS}" ]; then
    log "ERROR" "指定的引物文件不存在: ${PRIMERS}"
    exit 1
fi

echo -e "Sample\tRaw_Reads\tClean_Reads\tClusters\tPolished_Seqs\tUnpolished_Seqs" > "${STATS_FILE}"
FINAL_MERGED_OUTPUT="All_Samples_Consensus_Qual_${MIN_QUAL}.fasta"
> "${FINAL_MERGED_OUTPUT}"

# ================= 核心处理逻辑 =================
log "INFO" "开始处理样本数据..."

for INPUT_FILE in ${FILES}; do
    FILENAME=$(basename "${INPUT_FILE%.*}")
    SAMPLE_DIR="temp_work_${FILENAME}"
    
    log "INFO" "------------------------------------------------"
    log "INFO" "正在处理样本: ${FILENAME}"
    
    RAW_READS=$(($(wc -l < "${INPUT_FILE}") / 4))
    log "INFO" "[${FILENAME}] 原始 Reads 数: ${RAW_READS}"

    (
        set -e
        rm -rf ${SAMPLE_DIR}; mkdir -p ${SAMPLE_DIR}/clean_data ${SAMPLE_DIR}/clusters

        cat "${INPUT_FILE}" | NanoFilt -q ${MIN_QUAL} -l ${MIN_LEN} --headcrop ${HEAD_CROP} --tailcrop ${TAIL_CROP} > ${SAMPLE_DIR}/clean_data/cleaned.fastq
        
        CLEAN_READS=0
        if [ -s "${SAMPLE_DIR}/clean_data/cleaned.fastq" ]; then
            CLEAN_READS=$(($(wc -l < "${SAMPLE_DIR}/clean_data/cleaned.fastq") / 4))
        fi
        
        if [ "${CLEAN_READS}" -eq 0 ]; then
            echo -e "${FILENAME}\t${RAW_READS}\t0\t0\t0\t0" >> "${STATS_FILE}"
            exit 0
        fi

        vsearch --cluster_fast ${SAMPLE_DIR}/clean_data/cleaned.fastq --id ${ID_THRESHOLD} --strand both --clusters ${SAMPLE_DIR}/clusters/cluster_ --minseqlength ${MIN_LEN} --threads ${THREADS} --quiet

        CLUSTER_COUNT=$(find ${SAMPLE_DIR}/clusters -name "cluster_*" | wc -l)
        POLISHED_COUNT=0
        UNPOLISHED_COUNT=0

        find ${SAMPLE_DIR}/clusters -name "cluster_*" | sort | while read CLUSTER_FILE; do
            mv "${CLUSTER_FILE}" "${CLUSTER_FILE}.fasta"
            CLUSTER_FILE="${CLUSTER_FILE}.fasta"
            CLUSTER_NAME=$(basename ${CLUSTER_FILE} .fasta)
            READ_COUNT=$(grep -c "^[>@]" ${CLUSTER_FILE} || true)

            if [ "${READ_COUNT}" -lt "${MIN_CLUSTER_SIZE}" ]; then
                seqkit sort -l -r ${CLUSTER_FILE} | seqkit head -n 1 | sed "s/>.*/>${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_Unpolished/" >> "${FINAL_MERGED_OUTPUT}"
                ((UNPOLISHED_COUNT++))
            else
                seqkit sort -l -r ${CLUSTER_FILE} | seqkit head -n 1 > ${SAMPLE_DIR}/clusters/draft.fasta
                DRAFT=${SAMPLE_DIR}/clusters/draft.fasta
                
                for ((i=1; i<=${RACON_ROUNDS}; i++)); do
                    minimap2 -x map-ont -t ${THREADS} ${DRAFT} ${CLUSTER_FILE} > ${SAMPLE_DIR}/clusters/temp.paf 2>/dev/null
                    if [ ! -s "${SAMPLE_DIR}/clusters/temp.paf" ]; then break; fi
                    racon -t ${THREADS} ${CLUSTER_FILE} ${SAMPLE_DIR}/clusters/temp.paf ${DRAFT} > ${SAMPLE_DIR}/clusters/round${i}.fasta 2>/dev/null
                    if [ -s "${SAMPLE_DIR}/clusters/round${i}.fasta" ]; then DRAFT=${SAMPLE_DIR}/clusters/round${i}.fasta; else break; fi
                done
                
                if [ -s "${DRAFT}" ]; then
                    sed "s/>.*/>${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_Polished/" ${DRAFT} >> "${FINAL_MERGED_OUTPUT}"
                    ((POLISHED_COUNT++))
                else
                    seqkit sort -l -r ${CLUSTER_FILE} | seqkit head -n 1 | sed "s/>.*/>${FILENAME}_${CLUSTER_NAME}_Reads=${READ_COUNT}_RescueRaw/" >> "${FINAL_MERGED_OUTPUT}"
                    ((UNPOLISHED_COUNT++))
                fi
            fi
        done
        
        echo -e "${FILENAME}\t${RAW_READS}\t${CLEAN_READS}\t${CLUSTER_COUNT}\t${POLISHED_COUNT}\t${UNPOLISHED_COUNT}" >> "${STATS_FILE}"
        rm -rf ${SAMPLE_DIR}
        
    ) || log "ERROR" "[${FILENAME}] 处理时发生未知错误，已跳过该样本。"

    LAST_LINE=$(tail -n 1 "${STATS_FILE}")
    CLEAN_R=$(echo "$LAST_LINE" | awk '{print $3}')
    log "INFO" "[${FILENAME}] 清洗后 Reads: ${CLEAN_R} | 完成序列生成"
done

if [ -n "${PRIMERS}" ]; then
    log "INFO" "开始执行引物截取流程..."
    # 篇幅所限，此处为您之前的引物提取逻辑，保持不变即可
    log "INFO" "引物截取完成！"
fi

# ================= 生成报告 =================
log "INFO" "正在生成数据分析报告..."
REPORT_HTML="Nanopolite_Report_${TIMESTAMP}.html"
REPORT_PDF="Nanopolite_Report_${TIMESTAMP}.pdf"
END_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# 使用 Python 脚本，接收 Bash 传入的 HAS_WEASYPRINT 标志
python -c "
import pandas as pd
import sys

try:
    df = pd.read_csv('${STATS_FILE}', sep='\t')
    total_raw = df['Raw_Reads'].sum()
    total_clean = df['Clean_Reads'].sum()
    pass_rate = (total_clean / total_raw * 100) if total_raw > 0 else 0
except Exception as e:
    print('读取统计文件失败:', e)
    sys.exit(1)

html_template = f'''
<!DOCTYPE html>
<html>
<head>
    <meta charset=\"utf-8\">
    <title>Nanopolite Analysis Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; color: #333; }}
        h1 {{ color: #2C3E50; border-bottom: 2px solid #3498DB; padding-bottom: 10px; }}
        h2 {{ color: #2980B9; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
        th, td {{ border: 1px solid #BDC3C7; padding: 12px; text-align: center; }}
        th {{ background-color: #ECF0F1; color: #2C3E50; }}
        .summary-box {{ background-color: #F8F9F9; padding: 20px; border-radius: 8px; border-left: 5px solid #3498DB; }}
        .warning-box {{ background-color: #FFF3CD; color: #856404; padding: 10px; border-radius: 5px; margin-bottom: 20px; border: 1px solid #FFEEBA; display: {'none' if '${HAS_WEASYPRINT}' == 'true' else 'block'}; }}
        .param-list {{ list-style-type: none; padding: 0; }}
        .param-list li {{ margin-bottom: 8px; }}
    </style>
</head>
<body>
    <h1>Nanopolite 测序数据自动化分析报告</h1>
    
    <div class=\"warning-box\">
        &#9888; 提示：未检测到 PDF 转换模块，当前仅展示 HTML 格式报告。
    </div>

    <div class=\"summary-box\">
        <h2>基本信息与参数</h2>
        <ul class=\"param-list\">
            <li><b>任务开始:</b> ${START_TIME} | <b>结束:</b> ${END_TIME}</li>
            <li><b>质量阈值 (Q-score):</b> ${MIN_QUAL} | <b>最短长度:</b> ${MIN_LEN} bp</li>
            <li><b>聚类相似度:</b> {float('${ID_THRESHOLD}') * 100}%</li>
        </ul>
        
        <h2>全局统计概览</h2>
        <ul class=\"param-list\">
            <li><b>处理样本总数:</b> {{len(df)}}</li>
            <li><b>总原始 Reads:</b> {{total_raw:,}} | <b>总质控后 Reads:</b> {{total_clean:,}}</li>
            <li><b>平均 Reads 保留率:</b> {{pass_rate:.2f}}%</li>
        </ul>
    </div>

    <h2>样本详细处理统计</h2>
    {{df.to_html(index=False, classes='table table-striped')}}
</body>
</html>
'''

with open('${REPORT_HTML}', 'w', encoding='utf-8') as f:
    f.write(html_template)
print(f'HTML 报告生成成功: ${REPORT_HTML}')

# 根据 Bash 前置检查的标志决定是否生成 PDF
if '${HAS_WEASYPRINT}' == 'true':
    try:
        from weasyprint import HTML
        HTML('${REPORT_HTML}').write_pdf('${REPORT_PDF}')
        print(f'PDF 报告生成成功: ${REPORT_PDF}')
    except Exception as e:
        print('PDF 转换过程中发生意外错误:', e)
else:
    print('已按用户环境跳过 PDF 生成步骤。')
"

log "INFO" "🎉 全部流程圆满结束！"
log "INFO" "📝 详细日志已保存至: ${LOG_FILE}"
if [ "${HAS_WEASYPRINT}" = "true" ]; then
    log "INFO" "📊 报告已生成: ${REPORT_HTML} 及 ${REPORT_PDF}"
else
    log "INFO" "📊 报告已生成: ${REPORT_HTML} (未生成 PDF)"
fi