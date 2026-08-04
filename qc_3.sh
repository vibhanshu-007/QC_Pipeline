#!/usr/bin/env bash
#set -euo pipefail

# Pipeline build by ---> Vibhanshu singh

# Load Conda
source "$(conda info --base)/etc/profile.d/conda.sh"

Activate environment
conda activate qc_pipeline || {
    echo "Failed to activate Conda environment."
    exit 1
}
################################################################################
# FAST Pre-QC PIPELINE
################################################################################

# ── Configuration ──────────────────────────────────────────────────────────────
FASTQ_DIR="$(realpath "${1:-$(pwd)}")"
OUTDIR="${FASTQ_DIR}/qc_output"
THREADS=16
FASTQC_JOBS=16
CSV_FILE=""

declare -A SIZE_THRESHOLDS=(
    ["TARGET_Indiegene_Liquid"]=46.9
    ["TARGET_Indiegene"]=4.9
    ["Germline_plus_"]=0.5
    ["TARGET_First_Solid_Lite"]=0.5
    ["TARGET_First_Liquid"]=7
    ["TARGT_FIRST_LIQUID_LITE"]=6
    ["TARGET_First_Solid"]=1
    ["TARGET_Absolute"]=20
    ["TA_Germline"]=6
)

declare -A MIN_SIZE_THRESHOLDS=(
    ["TARGET_Indiegene_Liquid"]=4.7
    ["TARGET_Indiegene"]=0.5
    ["Germline_plus_"]=0.05
    ["TARGET_First_Solid_Lite"]=0.05
    ["TARGET_First_Liquid"]=0.7
    ["TARGT_FIRST_LIQUID_LITE"]=0.6
    ["TARGET_First_Solid"]=0.1
    ["TARGET_Absolute"]=2
    ["TA_Germline"]=0.6
)

declare -A PANEL_SIZE_MB=(
    ["TARGET_Indiegene"]=4.5
    ["TARGET_Indiegene_Liquid"]=4.5
    ["TARGET_First_Solid"]=0.45
    ["TARGET_First_Liquid"]=0.45
    ["TARGET_First_Solid_Lite"]=0.33
    ["Germline_plus_"]=0.56
    ["TARGT_FIRST_LIQUID_LITE"]=0.33
    ["TA_Germline"]=41
    ["TARGET_Absolute"]=41
)

# NEW: Expected depth thresholds
declare -A EXPECTED_DEPTH=(
    ["TARGET_Indiegene"]=1000
    ["TARGET_Indiegene_Liquid"]=10000
    ["TARGET_First_Solid"]=1500
    ["TARGET_First_Liquid"]=10000
    ["TARGET_First_Solid_Lite"]=1500
    ["Germline_plus_"]=300
    ["TARGT_FIRST_LIQUID_LITE"]=15000
    ["TA_Germline"]=1000
    ["TARGET_Absolute"]=1000
)

Q30_THRESHOLD=90
GC_PERCENT_THRESHOLD=40
DUPLICATION_THRESHOLD=40
LOG_FILE="$OUTDIR/pipeline.log"

OUTFILE=""
READCOUNT_CSV=""
DATASIZE_CSV=""
FINAL_REPORT=""
FASTQC_DIR=""

declare -A MD5_STATUS
declare -A READCOUNT_STATUS
declare -A Min_Data_Size_Status
declare -A DATA_SIZE_STATUS
declare -A ACTUAL_DEPTH_STATUS  
declare -A SAMPLE_TO_TEST

# ── Logging ────────────────────────────────────────────────────────────────────
setup_logging() {
    mkdir -p "$OUTDIR"
    LOG_FILE="$OUTDIR/pipeline.log"
    exec 1> >(tee -a "$LOG_FILE") 2>&1
}

log_info()    { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO   ] $*"; }
log_error()   { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR  ] $*"; }
log_success() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUCCESS] $*"; }
log_warn()    { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN   ] $*"; }

# ── Dependency check ───────────────────────────────────────────────────────────
validate_requirements() {
    log_info "Checking required tools..."
    local missing=()
    for tool in parallel md5sum perl python3 fastqc stat multiqc; do
        command -v "$tool" &>/dev/null || missing+=("$tool")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing[*]}"
        exit 1
    fi
    log_success "All required tools present"
}

# ── Input validation ───────────────────────────────────────────────────────────
validate_inputs() {
    log_info "Validating inputs in: $FASTQ_DIR"
    if [[ -z "$CSV_FILE" ]]; then
        mapfile -t csv_candidates < <(find "$FASTQ_DIR" -maxdepth 1 -name "*.csv" 2>/dev/null | sort)
        if [[ ${#csv_candidates[@]} -eq 0 ]]; then
            log_error "No CSV manifest found in $FASTQ_DIR"
            exit 1
        fi
        if [[ ${#csv_candidates[@]} -gt 1 ]]; then
            log_warn "Multiple CSV files found – using first: $(basename "${csv_candidates[0]}")"
        fi
        CSV_FILE="${csv_candidates[0]}"
    fi
    [[ ! -f "$CSV_FILE" ]] && { log_error "CSV not found: $CSV_FILE"; exit 1; }
    local fq_count
    fq_count=$(find "$FASTQ_DIR" -maxdepth 1 -name "*.fastq.gz" 2>/dev/null | wc -l)
    [[ $fq_count -eq 0 ]] && { log_error "No .fastq.gz files found in $FASTQ_DIR"; exit 1; }
    local header
    header=$(head -1 "$CSV_FILE" | tr -d '\r')
    SAMPLE_ID_COL=$(echo "$header" | tr ',' '\n' | grep -n "^Sample_ID$" | cut -d: -f1)
    TEST_NAME_COL=$(echo "$header" | tr ',' '\n' | grep -n "^Test_Name$" | cut -d: -f1)
    if [[ -z "$SAMPLE_ID_COL" ]]; then
        log_error "Column 'Sample_ID' not found in CSV header: $header"
        exit 1
    fi
    if [[ -z "$TEST_NAME_COL" ]]; then
        log_error "Column 'Test_Name' not found in CSV header: $header"
        exit 1
    fi
    log_success "CSV: $(basename "$CSV_FILE") | FASTQs: $fq_count"
    log_info "  Sample_ID col: $SAMPLE_ID_COL | Test_Name col: $TEST_NAME_COL"
}

# ── Populate SAMPLE_TO_TEST map from CSV ───────────────────────────────────────
load_sample_test_map() {
    log_info "Loading sample → test mapping from CSV..."
    local line_num=0
    while IFS=',' read -r -a fields; do
        (( line_num++ ))
        [[ $line_num -eq 1 ]] && continue
        local sample_id test_name
        sample_id=$(echo "${fields[$((SAMPLE_ID_COL - 1))]}" | tr -d '\r')
        test_name=$(echo "${fields[$((TEST_NAME_COL - 1))]}" | tr -d '\r')
        [[ -z "$sample_id" ]] && continue
        SAMPLE_TO_TEST["$sample_id"]="$test_name"
    done < "$CSV_FILE"
    log_success "Loaded ${#SAMPLE_TO_TEST[@]} sample(s) from CSV"
}
validate_sample_fastq_count() {

    log_info "Validating FASTQ pairs..."

    local failed=0

    for sample in "${!SAMPLE_TO_TEST[@]}"; do

        r1="$FASTQ_DIR/${sample}_R1.fastq.gz"
        r2="$FASTQ_DIR/${sample}_R2.fastq.gz"

        echo "Checking:"
        echo "$r1"
        echo "$r2"

        [[ ! -f "$r1" ]] && {
            log_error "Missing R1 for $sample"
            failed=1
        }

        [[ ! -f "$r2" ]] && {
            log_error "Missing R2 for $sample"
            failed=1
        }
    done

    [[ $failed -eq 1 ]] && exit 1

    log_success "All FASTQ pairs are present."
}
# ── Helper functions ─────────────────────────
get_total_sequences() {
    local fastqc_zip="$1"
    unzip -p "$fastqc_zip" "*/fastqc_data.txt" | awk -F'\t' '/^Total Sequences/{print $2}'
}

get_size_threshold() {
    local test_name="$1"
    [[ -v SIZE_THRESHOLDS["$test_name"] ]] || { log_error "Threshold not defined for Test_Name: '$test_name'"; exit 1; }
    echo "${SIZE_THRESHOLDS[$test_name]}"
}

get_min_size_threshold() {
    local test_name="$1"
    [[ -v MIN_SIZE_THRESHOLDS["$test_name"] ]] || { log_error "Minimum threshold not defined for Test_Name: '$test_name'"; exit 1; }
    echo "${MIN_SIZE_THRESHOLDS[$test_name]}"
}

get_panel_size_mb() {
    local test_name="$1"
    [[ -v PANEL_SIZE_MB["$test_name"] ]] || { log_error "Panel size not defined for Test_Name: '$test_name'"; exit 1; }
    echo "${PANEL_SIZE_MB[$test_name]}"
}

get_expected_depth() {
    local test_name="$1"
    [[ -v EXPECTED_DEPTH["$test_name"] ]] || { log_error "Expected depth not defined for Test_Name: '$test_name'"; exit 1; }
    echo "${EXPECTED_DEPTH[$test_name]}"
}

calculate_data_size() {
    local sequences="$1"
    awk -v seq="$sequences" 'BEGIN { printf "%.2f", (seq * 2 * 150) / 1000000000 }'
}

calculate_actual_depth() {
    local data_gb="$1"
    local panel_size_mb="$2"
    awk -v d="$data_gb" -v p="$panel_size_mb" 'BEGIN { printf "%.2f", d / (p / 1000) }'
}


calculate_q_counts() {
    local r1="$1"
    local r2="$2"

    pigz -p "$THREADS" -dc "$r1" "$r2" | perl -ne '
        next unless $. % 4 == 0;

        for (unpack("C*", $_)) {
            $total++;
            $q20++ if ($_ >= 53);
            $q30++ if ($_ >= 63);
            $q40++ if ($_ >= 73);
        }

        END {
            print "$total,$q20,$q30,$q40\n";
        }
    '
}
# ── Step 1: Renaming ───────────────────────────────────────────────────────────
run_renaming() {
    log_info "Running renaming scripts..."
    #python3 /path/to/Renaming.py
    #python3 /path/to/Renaming-S1.py
    log_success "Renaming complete"
}

# ── Step 2: MD5 checksum validation ───────────────────────────────────────────
run_md5_check() {
    log_info "Running MD5 checksum validation..."

    OUTFILE="$OUTDIR/md5_check.csv"
    echo "Filename,MD5,Status" > "$OUTFILE"

    local TMP_MD5="$OUTDIR/md5_parallel.tmp"

    # Calculate MD5 hashes in parallel
    find "$FASTQ_DIR" -maxdepth 1 -name "*.fastq.gz" | \
    parallel -j "$THREADS" --no-notice '
        md5=$(md5sum {} | awk "{print \$1}")
        echo "$(basename {})|$md5"
    ' > "$TMP_MD5"

    local -A seen_md5

    while IFS='|' read -r sample_file md5; do

        local sample
        sample="${sample_file%%_R1*}"
        sample="${sample%%_R2*}"

        if [[ -n "${seen_md5[$md5]:-}" ]]; then
            MD5_STATUS["$sample"]="FAIL"

            local other
            other=$(basename "${seen_md5[$md5]}")

            echo "$sample_file,$md5,DUPLICATE_OF:$other" >> "$OUTFILE"

            log_warn "Duplicate MD5 for sample '$sample' — matches $other"
        else
            MD5_STATUS["$sample"]="PASS"
            seen_md5[$md5]="$sample_file"

            echo "$sample_file,$md5,UNIQUE" >> "$OUTFILE"
        fi

    done < "$TMP_MD5"

    rm -f "$TMP_MD5"

    log_success "MD5 check complete. Results: $OUTFILE"
}

# ── Step 3: FastQC with GNU Parallel ──────────────────────────────────────────
run_fastqc() {
    log_info "Running FastQC on MD5-passed samples using GNU Parallel..."

    FASTQC_DIR="$OUTDIR/fastqc"
    mkdir -p "$FASTQC_DIR"

    local tmp_fastq_list="$OUTDIR/fastqc_files.list"
    > "$tmp_fastq_list"

    # Build FastQC input list
    for sample in "${!MD5_STATUS[@]}"; do
        [[ "${MD5_STATUS[$sample]}" != "PASS" ]] && continue

        local r1 r2

        r1=$(find "$FASTQ_DIR" -maxdepth 1 -name "${sample}*_R1*.fastq.gz" -print -quit)
        r2=$(find "$FASTQ_DIR" -maxdepth 1 -name "${sample}*_R2*.fastq.gz" -print -quit)

        if [[ -z "$r1" || -z "$r2" ]]; then
            log_warn "R1/R2 files not found for sample '$sample' — skipping FastQC"
            continue
        fi

        printf "%s\n%s\n" "$r1" "$r2" >> "$tmp_fastq_list"
    done

    if [[ ! -s "$tmp_fastq_list" ]]; then
        log_warn "No FASTQ files found for FastQC"
        rm -f "$tmp_fastq_list"
        return
    fi

    export FASTQC_DIR

    local file_count
    file_count=$(wc -l < "$tmp_fastq_list")

    log_info "Starting FastQC on $file_count FASTQ files using $FASTQC_JOBS parallel jobs..."

    parallel \
        --jobs "$FASTQC_JOBS" \
        --halt soon,fail=1 \
        --joblog "$OUTDIR/fastqc_parallel.log" \
        'fastqc -t 1 -o "$FASTQC_DIR" {}' \
        :::: "$tmp_fastq_list" >>"$LOG_FILE" 2>&1

    local rc=$?

    rm -f "$tmp_fastq_list"

    if [[ $rc -eq 0 ]]; then
        log_success "FastQC parallel processing complete"
    else
        log_warn "FastQC finished with exit code $rc"
        log_info "See $OUTDIR/fastqc_parallel.log"
    fi

    log_success "FastQC complete. Results: $FASTQC_DIR"
}
# ── Step 4: Read-count check ───────────────────────────────────────────────────
run_readcount_check() {
    log_info "Checking R1/R2 read counts..."
    READCOUNT_CSV="$OUTDIR/readcounts.csv"
    echo "Sample_ID,R1_Count,R2_Count,Status" > "$READCOUNT_CSV"
    for sample in "${!MD5_STATUS[@]}"; do
        [[ "${MD5_STATUS[$sample]}" != "PASS" ]] && continue
        local r1_zip r2_zip
        r1_zip=$(find "$FASTQC_DIR" -name "${sample}*R1*_fastqc.zip" | head -1)
        r2_zip=$(find "$FASTQC_DIR" -name "${sample}*R2*_fastqc.zip" | head -1)
        if [[ -z "$r1_zip" || -z "$r2_zip" ]]; then
            log_warn "FastQC zip not found for sample '$sample' — skipping read count check"
            READCOUNT_STATUS["$sample"]="FAIL"
            echo "$sample,NA,NA,FAIL" >> "$READCOUNT_CSV"
            continue
        fi
        local r1_count r2_count
        r1_count=$(get_total_sequences "$r1_zip")
        r2_count=$(get_total_sequences "$r2_zip")
        if [[ "$r1_count" == "$r2_count" ]]; then
            READCOUNT_STATUS["$sample"]="PASS"
        else
            READCOUNT_STATUS["$sample"]="FAIL"
            log_warn "Read count mismatch for '$sample': R1=$r1_count R2=$r2_count"
        fi
        echo "$sample,$r1_count,$r2_count,${READCOUNT_STATUS[$sample]}" >> "$READCOUNT_CSV"
    done
    log_success "Read count check complete. Results: $READCOUNT_CSV"
}

# ── Step 5: MultiQC ───────────────────────────────────────────────────────────
run_multiqc() {
    log_info "Running MultiQC..."
    multiqc "$FASTQC_DIR" -o "$OUTDIR/multiqc"
    log_success "MultiQC complete. Results: $OUTDIR/multiqc"
}

# ── Step 6: Q20/Q30/Q40 metrics ───────────────────────────────────────────────
run_quality_metrics() {
    log_info "Calculating Q20/Q30/Q40 metrics..."
    echo "Sample_ID,PCTQ20,PCTQ30,PCTQ40" > "$OUTDIR/qmetrics.csv"

    for sample in "${!READCOUNT_STATUS[@]}"; do
        # Skip if ReadCount failed
        [[ "${READCOUNT_STATUS[$sample]}" != "PASS" ]] && continue
        # skip if actual depth failed
        [[ "${ACTUAL_DEPTH_STATUS[$sample]:-PASS}" == "FAIL" ]] && continue
        local r1_fastq r2_fastq

        r1_fastq=$(find "$FASTQ_DIR" -maxdepth 1 -name "${sample}*_R1*.fastq.gz" -print -quit)
        r2_fastq=$(find "$FASTQ_DIR" -maxdepth 1 -name "${sample}*_R2*.fastq.gz" -print -quit)

        if [[ -z "$r1_fastq" || -z "$r2_fastq" ]]; then
            log_warn "FASTQ files not found for '$sample' — skipping Q-metrics"
            continue
        fi

        local total q20 q30 q40
        IFS=',' read -r total q20 q30 q40 <<< "$(calculate_q_counts "$r1_fastq" "$r2_fastq")"
        echo "DEBUG:"
        echo "sample=$sample"
        echo "total=$total"
        echo "q20=$q20"
        echo "q30=$q30"
        echo "q40=$q40"

        local pct20 pct30 pct40
       pct20=$(awk -v q="$q20" -v t="$total" '
BEGIN{
    if(t==0 || q==""){
        print "0.00"
    }else{
        printf "%.2f",100*q/t
    }
}')

pct30=$(awk -v q="$q30" -v t="$total" '
BEGIN{
    if(t==0 || q==""){
        print "0.00"
    }else{
        printf "%.2f",100*q/t
    }
}')

pct40=$(awk -v q="$q40" -v t="$total" '
BEGIN{
    if(t==0 || q==""){
        print "0.00"
    }else{
        printf "%.2f",100*q/t
    }
}')

        echo "$sample,$pct20,$pct30,$pct40" >> "$OUTDIR/qmetrics.csv"
    done

    log_success "Quality metrics complete. Results: $OUTDIR/qmetrics.csv"
}

# ── Step 7: Data size & Actual depth check ───────────────────────────────────────────
run_datasize_check() {
    log_info "Checking data size and depth thresholds..."
    DATASIZE_CSV="$OUTDIR/datasize.csv"
    echo "Sample_ID,Test_Name,Data_GB,Min_Data_Threshold_GB,Expected_Data_Threshold_GB,Min_Data_Size_Status,Expected_Data_Size_Status,Panel_MB,Actual_Depth,Expected_Depth,Actual_Depth_Status" > "$DATASIZE_CSV"
    
    for sample in "${!READCOUNT_STATUS[@]}"; do
        [[ "${READCOUNT_STATUS[$sample]}" != "PASS" ]] && continue
        
        local test_name="${SAMPLE_TO_TEST[$sample]:-}"
        if [[ -z "$test_name" ]]; then
            log_error "Sample '$sample' not found in CSV manifest"
            exit 1
        fi
        
        local r1_zip
        r1_zip=$(find "$FASTQC_DIR" -name "${sample}*R1*_fastqc.zip" | head -1)
        if [[ -z "$r1_zip" ]]; then
            log_warn "FastQC zip not found for '$sample' — skipping data size check"
            continue
        fi
        
        local sequences
        sequences=$(get_total_sequences "$r1_zip")
        local data_gb
        data_gb=$(calculate_data_size "$sequences")
        local min_threshold exp_threshold
        min_threshold=$(get_min_size_threshold "$test_name")
        exp_threshold=$(get_size_threshold "$test_name")
        
        if awk "BEGIN{exit !($data_gb >= $min_threshold)}"; then
            Min_Data_Size_Status["$sample"]="PASS"
        else
            Min_Data_Size_Status["$sample"]="FAIL"
            log_warn "'$sample' data size ${data_gb} GB < min threshold ${min_threshold} GB"
        fi
        
        if awk "BEGIN{exit !($data_gb >= $exp_threshold)}"; then
            DATA_SIZE_STATUS["$sample"]="PASS"
        else
            DATA_SIZE_STATUS["$sample"]="FAIL"
            log_warn "'$sample' data size ${data_gb} GB < expected threshold ${exp_threshold} GB"
        fi
        
        local panel_size_mb actual_depth min_depth expected_depth
        panel_size_mb=$(get_panel_size_mb "$test_name")
        actual_depth=$(calculate_actual_depth "$data_gb" "$panel_size_mb")
        expected_depth=$(get_expected_depth "$test_name")
        
       
        if awk "BEGIN{exit !($actual_depth >= $expected_depth)}"; then
            ACTUAL_DEPTH_STATUS["$sample"]="PASS"
            log_info "'$sample' Actual Depth: ${actual_depth}x >= Expected: ${expected_depth}x - PASS"
        else
            ACTUAL_DEPTH_STATUS["$sample"]="FAIL"
            log_warn "'$sample' Actual Depth: ${actual_depth}x < Expected: ${expected_depth}x - FAIL"
        fi
        
        echo "$sample,$test_name,$data_gb,$min_threshold,$exp_threshold,${Min_Data_Size_Status[$sample]},${DATA_SIZE_STATUS[$sample]},$panel_size_mb,$actual_depth,$expected_depth,${ACTUAL_DEPTH_STATUS[$sample]},$min_depth" >> "$DATASIZE_CSV"
    done
    log_success "Data size check complete. Results: $DATASIZE_CSV"
}

# ── Merge R1 and R2 module statuses ──────────────────────────────────────────
merge_r1_r2_modules() {
    local r1="$1"
    local r2="$2"

    IFS=',' read -ra A <<< "$r1"
    IFS=',' read -ra B <<< "$r2"

    local merged=()

    for ((i=0; i<${#A[@]}; i++)); do
        if [[ "${A[$i]}" == "FAIL" || "${B[$i]}" == "FAIL" ]]; then
            merged+=("FAIL")
        elif [[ "${A[$i]}" == "WARN" || "${B[$i]}" == "WARN" ]]; then
            merged+=("WARN")
        else
            merged+=("PASS")
        fi
    done

    IFS=','
    echo "${merged[*]}"
}

# ── Step 8: Parse ALL FastQC module statuses ──────────────────────────────────
run_parse_fastqc_modules() {
    module_csv="$OUTDIR/fastqc_module_status.csv"
    echo "Sample_ID,Basic_Statistics,Per_Base_Sequence_Quality,Per_Sequence_Quality_Scores,Per_Base_Sequence_Content,Per_Sequence_GC_Content,Per_Base_N_Content,Sequence_Length_Distribution,Sequence_Duplication_Levels,Overrepresented_Sequences,Adapter_Content" > "$module_csv"

    declare -A R1_STATUS
    declare -A R2_STATUS

    # Define the exact order of modules as they appear in FastQC summary.txt
    local module_order=(
        "Basic Statistics"
        "Per base sequence quality"
        "Per sequence quality scores"
        "Per base sequence content"
        "Per sequence GC content"
        "Per base N content"
        "Sequence Length Distribution"
        "Sequence Duplication Levels"
        "Overrepresented sequences"
        "Adapter Content"
    )

    for zip in "$FASTQC_DIR"/*_fastqc.zip; do
        [[ -f "$zip" ]] || continue

        sample=$(basename "$zip" "_fastqc.zip")
        sample=${sample%_R1}
        sample=${sample%_R2}

        # Skip if Actual Depth failed
        #[[ "${ACTUAL_DEPTH_STATUS[$sample]:-PASS}" == "FAIL" ]] && continue

        # Read summary.txt lines into an associative array
        local -A mod_status

        while IFS=$'\t' read -r status module filename || [[ -n "$status" ]]; do
            status=$(echo "$status" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            module=$(echo "$module" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            mod_status["$module"]="$status"
        done < <(unzip -p "$zip" "*/summary.txt")

        # Build status line in the defined order
        local status_line=""
        for mod in "${module_order[@]}"; do
            local st="${mod_status[$mod]:-NA}"
            if [[ -n "$status_line" ]]; then
                status_line="${status_line},${st}"
            else
                status_line="${st}"
            fi
        done

        # Store R1 and R2 separately
        if [[ "$zip" == *"_R1_"* || "$zip" == *"_R1_fastqc.zip" ]]; then
            R1_STATUS["$sample"]="$status_line"
        else
            R2_STATUS["$sample"]="$status_line"
        fi
    done

    # Merge R1 and R2
    for sample in "${!R1_STATUS[@]}"; do
        r1="${R1_STATUS[$sample]}"
        r2="${R2_STATUS[$sample]:-}"
        if [[ -z "$r2" ]]; then
            echo "$sample,$r1" >> "$module_csv"
        else
            merged=$(merge_r1_r2_modules "$r1" "$r2")
            echo "$sample,$merged" >> "$module_csv"
        fi
    done

    log_success "FastQC module status extraction complete"
    log_info "Module status: $module_csv"
}
run_extract_gc_dup_metrics() {
    log_info "Extracting average GC and Duplication metrics from FastQC..."

    local gc_dup_csv="$OUTDIR/gc_dup_metrics.csv"

    echo "Sample_ID,GC_Content_Percent,GC_FLAG,Duplication_Percent,DUPLICATION_FLAG" > "$gc_dup_csv"

    for sample in "${!READCOUNT_STATUS[@]}"; do

        [[ "${READCOUNT_STATUS[$sample]}" != "PASS" ]] && continue

        local r1_zip r2_zip
        r1_zip=$(find "$FASTQC_DIR" -name "${sample}*R1*_fastqc.zip" -print -quit)
        r2_zip=$(find "$FASTQC_DIR" -name "${sample}*R2*_fastqc.zip" -print -quit)

        if [[ -z "$r1_zip" || -z "$r2_zip" ]]; then
            log_warn "FastQC reports missing for sample '$sample'"
            continue
        fi

########################
# GC%
########################

local gc1 gc2 avg_gc

gc1=$(unzip -p "$r1_zip" "*/fastqc_data.txt" | \
      awk -F'\t' '$1=="%GC"{printf "%.2f",$2; exit}')

gc2=$(unzip -p "$r2_zip" "*/fastqc_data.txt" | \
      awk -F'\t' '$1=="%GC"{printf "%.2f",$2; exit}')

avg_gc=$(awk -v a="$gc1" -v b="$gc2" \
          'BEGIN{printf "%.2f",(a+b)/2}')
        ########################
        # Duplication %
        ########################

        local dup1 dup2 avg_dup

        dup1=$(unzip -p "$r1_zip" "*/fastqc_data.txt" | awk -F'\t' '
            /^#Total Deduplicated Percentage/{
                printf "%.2f",100-$2
                exit
            }')

        dup2=$(unzip -p "$r2_zip" "*/fastqc_data.txt" | awk -F'\t' '
            /^#Total Deduplicated Percentage/{
                printf "%.2f",100-$2
                exit
            }')

        avg_dup=$(awk -v a="$dup1" -v b="$dup2" 'BEGIN{printf "%.2f",(a+b)/2}')

        ########################
        # Flags
        ########################

        local gc_flag dup_flag

        if awk -v gc="$avg_gc" -v th="$GC_PERCENT_THRESHOLD" 'BEGIN{exit !(gc > th)}'
        then
            gc_flag="FAIL"
        else
            gc_flag="PASS"
        fi

        if awk -v dup="$avg_dup" -v th="$DUPLICATION_THRESHOLD" 'BEGIN{exit !(dup > th)}'
        then
            dup_flag="FAIL"
        else
            dup_flag="PASS"
        fi

        echo "$sample,$avg_gc,$gc_flag,$avg_dup,$dup_flag" >> "$gc_dup_csv"

    done

    log_success "GC/Duplication extraction complete."
    log_info "Results: $gc_dup_csv"
}

# ── Step 10: Final per-sample conclusion ──────────────────────────────────────
run_final_report() {
    log_info "Generating final QC report..."
    FINAL_REPORT="$OUTDIR/final_report.csv"
    echo "Sample_ID,MD5_Status,ReadCount_Status,MinSize_Status,DataSize_Status,Actual_Depth_Status,PCTQ30,Q30_Status,GC_Content_Percent,GC_FLAG,Duplication_Percent,DUPLICATION_FLAG,Basic_Statistics,Per_Base_Sequence_Quality,Per_Sequence_Quality_Scores,Per_Base_Sequence_Content,Per_Sequence_GC_Content,Per_Base_N_Content,Sequence_Length_Distribution,Sequence_Duplication_Levels,Overrepresented_Sequences,Adapter_Content,Conclusion,Pre_QC_Status_MD5_Readcount_Actual_Depth,Pre_QC_Status" > "$FINAL_REPORT"

    declare -A PCTQ30_MAP
    if [[ -f "$OUTDIR/qmetrics.csv" ]]; then
        while IFS=',' read -r sid pct20 pct30 pct40; do
            [[ "$sid" == "Sample_ID" ]] && continue
            PCTQ30_MAP["$sid"]="$pct30"
        done < "$OUTDIR/qmetrics.csv"
    fi

    declare -A MOD_BASIC MOD_PBQ MOD_PSQ MOD_PBC MOD_GC MOD_PBN MOD_SLD MOD_SDL MOD_OS MOD_ADP
    if [[ -f "$OUTDIR/fastqc_module_status.csv" ]]; then
        while IFS=',' read -r sid bs pbq psq pbc gc nc ld dup or adp; do
            [[ "$sid" == "Sample_ID" ]] && continue
            MOD_BASIC["$sid"]=$(echo "$bs" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_PBQ["$sid"]=$(echo "$pbq" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_PSQ["$sid"]=$(echo "$psq" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_PBC["$sid"]=$(echo "$pbc" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_GC["$sid"]=$(echo "$gc" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_PBN["$sid"]=$(echo "$nc" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_SLD["$sid"]=$(echo "$ld" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_SDL["$sid"]=$(echo "$dup" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_OS["$sid"]=$(echo "$or" | xargs | tr '[:lower:]' '[:upper:]')
            MOD_ADP["$sid"]=$(echo "$adp" | xargs | tr '[:lower:]' '[:upper:]')
        done < "$OUTDIR/fastqc_module_status.csv"
    fi

    declare -A GC_PERCENT GC_FLAG DUP_PERCENT DUP_FLAG
    if [[ -f "$OUTDIR/gc_dup_metrics.csv" ]]; then
        while IFS=',' read -r sid gc_pct gc_flag dup_pct dup_flag; do
            [[ "$sid" == "Sample_ID" ]] && continue
            GC_PERCENT["$sid"]="$gc_pct"
            GC_FLAG["$sid"]=$(echo "$gc_flag" | xargs | tr '[:lower:]' '[:upper:]')
            DUP_PERCENT["$sid"]="$dup_pct"
            DUP_FLAG["$sid"]=$(echo "$dup_flag" | xargs | tr '[:lower:]' '[:upper:]')
        done < "$OUTDIR/gc_dup_metrics.csv"
    fi

    local -A all_samples
    for s in "${!MD5_STATUS[@]}"       ; do all_samples["$s"]=1; done
    for s in "${!READCOUNT_STATUS[@]}" ; do all_samples["$s"]=1; done

    for sample in "${!all_samples[@]}"; do
        local md5="${MD5_STATUS[$sample]:-FAIL}"
        local readcount="${READCOUNT_STATUS[$sample]:-NOT_RUN}"
        local min_size="${Min_Data_Size_Status[$sample]:-NOT_RUN}"
        local data_size="${DATA_SIZE_STATUS[$sample]:-NOT_RUN}"
        local actual_depth_status="${ACTUAL_DEPTH_STATUS[$sample]:-NOT_RUN}"
        local pctq30="${PCTQ30_MAP[$sample]:-0}"

        local gc_pct="${GC_PERCENT[$sample]:-0}"
        local gc_flag="${GC_FLAG[$sample]:-NA}"
        local dup_pct="${DUP_PERCENT[$sample]:-0}"
        local dup_flag="${DUP_FLAG[$sample]:-NA}"

        local bs="${MOD_BASIC[$sample]:-NA}"
        local pbq="${MOD_PBQ[$sample]:-NA}"
        local psq="${MOD_PSQ[$sample]:-NA}"
        local pbc="${MOD_PBC[$sample]:-NA}"
        local gc_mod="${MOD_GC[$sample]:-NA}"
        local pbn="${MOD_PBN[$sample]:-NA}"
        local sld="${MOD_SLD[$sample]:-NA}"
        local sdl="${MOD_SDL[$sample]:-NA}"
        local os="${MOD_OS[$sample]:-NA}"
        local adp="${MOD_ADP[$sample]:-NA}"

        local q30_status="FAIL"
        if awk "BEGIN{exit !(${pctq30} >= ${Q30_THRESHOLD})}"; then
            q30_status="PASS"
        fi

        # Build pre-QC status comments for the three core checks
        local pre_qc_comments=()
        local pre_qc_overall="PASS"
        
        if [[ "$md5" != "PASS" ]]; then
            pre_qc_comments+=("Duplicate MD5 checksum failed")
            pre_qc_overall="FAIL"
        fi
        
        if [[ "$readcount" != "PASS" ]]; then
            pre_qc_comments+=("R1 and R2 read counts are mismatched")
            pre_qc_overall="FAIL"
        fi
        
        if [[ "$actual_depth_status" == "FAIL" ]]; then
            pre_qc_comments+=("Actual depth below threshold")
            pre_qc_overall="FAIL"
        fi
        
        # Join comments with "; " if multiple, otherwise single comment or empty
        local pre_qc_comment_str=""
        if [[ ${#pre_qc_comments[@]} -gt 0 ]]; then
            pre_qc_comment_str=$(printf "; %s" "${pre_qc_comments[@]}")
            pre_qc_comment_str="${pre_qc_comment_str:2}"  # Remove leading "; "
        fi

        # Build detailed conclusion
        local conclusion=""
        local fail_list=()
        local warn_list=()

        # Check if any of the three core pre-QC checks failed
        local pre_qc_failed=false
        if [[ "$md5" != "PASS" ]] || [[ "$readcount" != "PASS" ]] || [[ "$actual_depth_status" == "FAIL" ]]; then
            pre_qc_failed=true
        fi

        if [[ "$pre_qc_failed" == true ]]; then
            # If any pre-QC check failed, conclusion is just the comments
            conclusion="$pre_qc_comment_str"
        else
            # All three pre-QC checks passed - use existing downstream QC logic
            # MinSize (size_flag)
            if [[ "$min_size" != "PASS" ]]; then fail_list+=("size_flag"); fi
            # DataSize (Size_Default_Flag)
            if [[ "$data_size" != "PASS" ]]; then fail_list+=("Size_Default_Flag"); fi
            
            # Q30
            if [[ "$q30_status" != "PASS" ]]; then fail_list+=("Q30"); fi
            # GC flag
            if [[ "$gc_flag" == "FAIL" ]]; then
                fail_list+=("GC_flag")
            elif [[ "$gc_flag" == "WARN" ]]; then
                warn_list+=("GC_flag")
            fi
            # Duplication flag
            if [[ "$dup_flag" == "FAIL" ]]; then
                fail_list+=("Duplication_flag")
            elif [[ "$dup_flag" == "WARN" ]]; then
                warn_list+=("Duplication_flag")
            fi

            # FastQC modules
            if [[ "$pbc" == "FAIL" ]]; then
                fail_list+=("Per_base_sequence_content")
            elif [[ "$pbc" == "WARN" ]]; then
                warn_list+=("Per_base_sequence_content")
            fi

            if [[ "$gc_mod" == "FAIL" ]]; then
                fail_list+=("Per_sequence_GC_content")
            elif [[ "$gc_mod" == "WARN" ]]; then
                warn_list+=("Per_sequence_GC_content")
            fi

            if [[ "$sdl" == "FAIL" ]]; then
                fail_list+=("Sequence_Duplication_Levels")
            elif [[ "$sdl" == "WARN" ]]; then
                warn_list+=("Sequence_Duplication_Levels")
            fi

            if [[ "$adp" == "FAIL" ]]; then
                fail_list+=("Adapter_Content")
            elif [[ "$adp" == "WARN" ]]; then
                warn_list+=("Adapter_Content")
            fi

            # Build conclusion string
            if [[ ${#fail_list[@]} -gt 0 ]]; then
                local fail_str
                fail_str=$(printf ", %s" "${fail_list[@]}")
                conclusion="FAIL (${fail_str:2})"
            fi

            if [[ ${#warn_list[@]} -gt 0 ]]; then
                local warn_str
                warn_str=$(printf ", %s" "${warn_list[@]}")
                if [[ -n "$conclusion" ]]; then
                    conclusion="${conclusion}; WARN (${warn_str:2})"
                else
                    conclusion="WARN (${warn_str:2})"
                fi
            fi

            [[ -z "$conclusion" ]] && conclusion="PASS"
        fi

        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s","%s","%s"\n' \
"$sample" \
"$md5" \
"$readcount" \
"$min_size" \
"$data_size" \
"$actual_depth_status" \
"$pctq30" \
"$q30_status" \
"$gc_pct" \
"$gc_flag" \
"$dup_pct" \
"$dup_flag" \
"$bs" \
"$pbq" \
"$psq" \
"$pbc" \
"$gc_mod" \
"$pbn" \
"$sld" \
"$sdl" \
"$os" \
"$adp" \
"${conclusion//\"/\"\"}" \
"${pre_qc_comment_str//\"/\"\"}" \
"$pre_qc_overall" >> "$FINAL_REPORT"
    done

    log_success "Final report: $FINAL_REPORT"
}

# ── Step 11: Generate master XLSX ─────────────────────────────────────────────
run_generate_xlsx() {
    log_info "Generating master XLSX report..."
    python3 - <<'PYEOF'
import csv, sys, os
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

OUTDIR = os.environ.get("OUTDIR", ".")
FINAL_REPORT   = os.path.join(OUTDIR, "final_report.csv")
DATASIZE_CSV   = os.path.join(OUTDIR, "datasize.csv")
QMETRICS_CSV   = os.path.join(OUTDIR, "qmetrics.csv")
READCOUNT_CSV  = os.path.join(OUTDIR, "readcounts.csv")
MD5_CSV        = os.path.join(OUTDIR, "md5_check.csv")
OUTPUT_XLSX    = os.path.join(OUTDIR, "QC_Master_Report.xlsx")
FASTQC_MODULE_CSV = os.path.join(OUTDIR, "fastqc_module_status.csv")
GC_DUP_CSV     = os.path.join(OUTDIR, "gc_dup_metrics.csv")

PASS_FILL    = PatternFill("solid", fgColor="C6EFCE")
FAIL_FILL    = PatternFill("solid", fgColor="FFC7CE")
WARN_FILL    = PatternFill("solid", fgColor="FFEB9C")
NA_FILL      = PatternFill("solid", fgColor="EFEFEF")
HEADER_FILL  = PatternFill("solid", fgColor="4472C4")
HEADER_FONT  = Font(bold=True, color="FFFFFF", size=11)
TITLE_FONT   = Font(bold=True, size=13)
THIN_BORDER  = Border(left=Side(style="thin"), right=Side(style="thin"), top=Side(style="thin"), bottom=Side(style="thin"))

def read_csv(path):
    if not os.path.exists(path):
        return [], []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return fieldnames, rows

def cell_fill(cell, value):
    v = str(value).strip().upper()
    if v.startswith("FAIL") or "FAIL" in v:
        cell.fill = FAIL_FILL
    elif v.startswith("WARN") or "WARN" in v:
        cell.fill = WARN_FILL
    elif v in ("PASS",):
        cell.fill = PASS_FILL
    elif v in ("NA", "NOT_RUN", ""):
        cell.fill = NA_FILL

def write_header(ws, headers, row=1, col_offset=1):
    for ci, h in enumerate(headers, start=col_offset):
        cell = ws.cell(row=row, column=ci, value=h)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
        cell.border = THIN_BORDER

def auto_col_width(ws):
    for col in ws.columns:
        max_len = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[get_column_letter(col[0].column)].width = min(max_len + 4, 50)

_, final_rows   = read_csv(FINAL_REPORT)
_, ds_rows      = read_csv(DATASIZE_CSV)
_, qm_rows      = read_csv(QMETRICS_CSV)
_, rc_rows      = read_csv(READCOUNT_CSV)
_, gc_dup_rows  = read_csv(GC_DUP_CSV)

# Determine which samples passed pre-QC
pre_qc_pass_samples = set()
for r in final_rows:
    if str(r.get("Pre_QC_Status", "")).upper() == "PASS":
        pre_qc_pass_samples.add(r.get("Sample_ID", ""))

wb = Workbook()
ws_summary = wb.active
ws_summary.title = "QC Summary"

# NEW: Updated metrics list with Pre-QC columns - split into two groups
pre_qc_metrics = [
    "MD5_Status",
    "ReadCount_Status",
    "Pre_QC_Status_MD5_Readcount_Actual_Depth",
    "Pre_QC_Status"
]

downstream_metrics = [
    "Q30_Status",
    "GC_FLAG",
    "DUPLICATION_FLAG",
    "Basic_Statistics",
    "Per_Base_Sequence_Quality",
    "Per_Sequence_Quality_Scores",
    "Per_Base_Sequence_Content",
    "Per_Sequence_GC_Content",
    "Per_Base_N_Content",
    "Sequence_Length_Distribution",
    "Sequence_Duplication_Levels",
    "Overrepresented_Sequences",
    "Adapter_Content"
]

ws_summary["A1"] = "QC Summary — Pass / Warn / Fail Counts"
ws_summary["A1"].font = TITLE_FONT
ws_summary["A1"].fill = PatternFill("solid", fgColor="DCE6F1")
ws_summary.merge_cells("A1:D1")
headers_sum = ["Metric", "PASS", "FAIL / WARN", "TOTAL"]
write_header(ws_summary, headers_sum, row=2)

row_idx = 3
# Pre-QC metrics - count all samples
for metric in pre_qc_metrics:
    pass_cnt = sum(1 for r in final_rows if str(r.get(metric,"")).upper() in ("PASS",))
    fail_cnt = sum(1 for r in final_rows if str(r.get(metric,"")).upper() not in ("PASS","NA","NOT_RUN",""))
    total = len(final_rows)
    ws_summary.cell(row=row_idx, column=1, value=metric).border = THIN_BORDER
    pc = ws_summary.cell(row=row_idx, column=2, value=pass_cnt)
    pc.fill = PASS_FILL; pc.border = THIN_BORDER; pc.alignment = Alignment(horizontal="center")
    fc = ws_summary.cell(row=row_idx, column=3, value=fail_cnt)
    fc.fill = FAIL_FILL if fail_cnt > 0 else PASS_FILL; fc.border = THIN_BORDER; fc.alignment = Alignment(horizontal="center")
    tc = ws_summary.cell(row=row_idx, column=4, value=total)
    tc.border = THIN_BORDER; tc.alignment = Alignment(horizontal="center")
    row_idx += 1

# Downstream metrics - count only samples where Pre_QC_Status == PASS
for metric in downstream_metrics:
    pass_cnt = sum(1 for r in final_rows 
                   if r.get("Sample_ID", "") in pre_qc_pass_samples 
                   and str(r.get(metric,"")).upper() in ("PASS",))
    fail_cnt = sum(1 for r in final_rows 
                   if r.get("Sample_ID", "") in pre_qc_pass_samples 
                   and str(r.get(metric,"")).upper() not in ("PASS","NA","NOT_RUN",""))
    total = len(pre_qc_pass_samples)
    ws_summary.cell(row=row_idx, column=1, value=metric).border = THIN_BORDER
    pc = ws_summary.cell(row=row_idx, column=2, value=pass_cnt)
    pc.fill = PASS_FILL; pc.border = THIN_BORDER; pc.alignment = Alignment(horizontal="center")
    fc = ws_summary.cell(row=row_idx, column=3, value=fail_cnt)
    fc.fill = FAIL_FILL if fail_cnt > 0 else PASS_FILL; fc.border = THIN_BORDER; fc.alignment = Alignment(horizontal="center")
    tc = ws_summary.cell(row=row_idx, column=4, value=total)
    tc.border = THIN_BORDER; tc.alignment = Alignment(horizontal="center")
    row_idx += 1

auto_col_width(ws_summary)

ws_samples = wb.create_sheet("Sample-Wise Report")
_, fastqc_rows = read_csv(FASTQC_MODULE_CSV)
fastqc_map = {r["Sample_ID"]: r for r in fastqc_rows if "Sample_ID" in r}
ds_map = {r["Sample_ID"]: r for r in ds_rows if "Sample_ID" in r}
qm_map = {r["Sample_ID"]: r for r in qm_rows if "Sample_ID" in r}
rc_map = {r["Sample_ID"]: r for r in rc_rows if "Sample_ID" in r}
gc_dup_map = {r["Sample_ID"]: r for r in gc_dup_rows if "Sample_ID" in r}

# UPDATED: New column order with Pre-QC columns after Conclusion
COMBINED_HEADERS = [
    "Sample_ID","Test_Name","MD5_Status","R1_Count","R2_Count","ReadCount_Status",
    "Data_GB","Min_Data_Threshold_GB","Expected_Data_Threshold_GB","Min_Data_Size_Status","Expected_Data_Size_Status",
    "Panel_MB","Actual_Depth","Expected_Depth","Actual_Depth_Status","Min_Depth",
    "PCTQ20","PCTQ30","PCTQ40","Q30_Status",
    "GC_Content_Percent","GC_FLAG","Duplication_Percent","DUPLICATION_FLAG",
    "Basic_Statistics","Per_Base_Sequence_Quality","Per_Sequence_Quality_Scores",
    "Per_Base_Sequence_Content","Per_Sequence_GC_Content","Per_Base_N_Content",
    "Sequence_Length_Distribution","Sequence_Duplication_Levels",
    "Overrepresented_Sequences","Adapter_Content","Conclusion",
    "Pre_QC_Status_MD5_Readcount_Actual_Depth","Pre_QC_Status"
]
write_header(ws_samples, COMBINED_HEADERS, row=1)

for ri, row in enumerate(final_rows, start=2):
    sid = row.get("Sample_ID", "")
    ds = ds_map.get(sid, {})
    qm = qm_map.get(sid, {})
    rc = rc_map.get(sid, {})
    fq = fastqc_map.get(sid, {})
    gc_dup = gc_dup_map.get(sid, {})
    q30_val = qm.get("PCTQ30", "0") or "0"
    try:
        q30_status = "PASS" if float(q30_val) >= 90 else "FAIL"
    except:
        q30_status = "FAIL"
    combined = {
        "Sample_ID": sid,
        "Test_Name": ds.get("Test_Name", ""),
        "MD5_Status": row.get("MD5_Status", ""),
        "R1_Count": rc.get("R1_Count", ""),
        "R2_Count": rc.get("R2_Count", ""),
        "ReadCount_Status": row.get("ReadCount_Status", ""),
        "Data_GB": ds.get("Data_GB", ""),
        "Min_Data_Threshold_GB": ds.get("Min_Data_Threshold_GB", ""),
        "Expected_Data_Threshold_GB": ds.get("Expected_Data_Threshold_GB", ""),
        "Min_Data_Size_Status": row.get("MinSize_Status", ""),
        "Expected_Data_Size_Status": row.get("DataSize_Status", ""),
        "Panel_MB": ds.get("Panel_MB", ""),
        "Actual_Depth": ds.get("Actual_Depth", ""),
        "Expected_Depth": ds.get("Expected_Depth", ""),
        "Actual_Depth_Status": row.get("Actual_Depth_Status", ""),
        "Min_Depth": ds.get("Min_Depth", ""),
        "PCTQ20": qm.get("PCTQ20", ""),
        "PCTQ30": qm.get("PCTQ30", ""),
        "PCTQ40": qm.get("PCTQ40", ""),
        "Q30_Status": q30_status,
        "GC_Content_Percent": gc_dup.get("GC_Content_Percent", ""),
        "GC_FLAG": gc_dup.get("GC_FLAG", ""),
        "Duplication_Percent": gc_dup.get("Duplication_Percent", ""),
        "DUPLICATION_FLAG": gc_dup.get("DUPLICATION_FLAG", ""),
        "Basic_Statistics": fq.get("Basic_Statistics", ""),
        "Per_Base_Sequence_Quality": fq.get("Per_Base_Sequence_Quality", ""),
        "Per_Sequence_Quality_Scores": fq.get("Per_Sequence_Quality_Scores", ""),
        "Per_Base_Sequence_Content": fq.get("Per_Base_Sequence_Content", ""),
        "Per_Sequence_GC_Content": fq.get("Per_Sequence_GC_Content", ""),
        "Per_Base_N_Content": fq.get("Per_Base_N_Content", ""),
        "Sequence_Length_Distribution": fq.get("Sequence_Length_Distribution", ""),
        "Sequence_Duplication_Levels": fq.get("Sequence_Duplication_Levels", ""),
        "Overrepresented_Sequences": fq.get("Overrepresented_Sequences", ""),
        "Adapter_Content": fq.get("Adapter_Content", ""),
        "Conclusion": row.get("Conclusion", ""),
        # FIX: Use the correct column name from final_report.csv (with parentheses)
        "Pre_QC_Status_MD5_Readcount_Actual_Depth": row.get("Pre_QC_Status_MD5_Readcount_Actual_Depth", ""),
        "Pre_QC_Status": row.get("Pre_QC_Status", "")
    }
    for ci, h in enumerate(COMBINED_HEADERS, start=1):
        val = combined.get(h, "")
        cell = ws_samples.cell(row=ri, column=ci, value=val)
        cell_fill(cell, val)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
        cell.border = THIN_BORDER

ws_samples.freeze_panes = "A2"
auto_col_width(ws_samples)

# NEW: Create PreQC worksheet
ws_preqc = wb.create_sheet("PreQC")

# Build data from final_report.csv and datasize.csv
preqc_data = []
for row in final_rows:
    sid = row.get("Sample_ID", "")
    ds = ds_map.get(sid, {})
    preqc_data.append({
        "Sample_ID": sid,
        "Test_Name": ds.get("Test_Name", ""),
        # FIX: Use the correct column name from final_report.csv (with parentheses)
        "Pre_QC_Status_MD5_Readcount_Actual_Depth": row.get("Pre_QC_Status_MD5_Readcount_Actual_Depth", ""),
        "Pre_QC_Status": row.get("Pre_QC_Status", "")
    })

PREQC_HEADERS = ["Sample_ID", "Test_Name", "Pre_QC_Status_MD5_Readcount_Actual_Depth", "Pre_QC_Status"]
write_header(ws_preqc, PREQC_HEADERS, row=1)

for ri, data in enumerate(preqc_data, start=2):
    for ci, h in enumerate(PREQC_HEADERS, start=1):
        val = data.get(h, "")
        cell = ws_preqc.cell(row=ri, column=ci, value=val)
        cell_fill(cell, val)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
        cell.border = THIN_BORDER

ws_preqc.freeze_panes = "A2"
auto_col_width(ws_preqc)

wb.save(OUTPUT_XLSX)
print(f"XLSX saved: {OUTPUT_XLSX}")
PYEOF
    log_success "XLSX report: $OUTDIR/QC_Master_Report.xlsx"
}
# ── Main orchestration ─────────────────────────────────────────────────────────
main() {
    setup_logging
    log_info "========== QC PIPELINE START =========="
    log_info "FASTQ dir : $FASTQ_DIR"
    log_info "Output dir: $OUTDIR"
    validate_requirements
    validate_inputs
    load_sample_test_map
    validate_sample_fastq_count
    run_renaming
    run_md5_check
    run_fastqc
    run_readcount_check
    run_multiqc
    run_datasize_check
    run_quality_metrics
    run_parse_fastqc_modules
    #run_extract_multiqc_metrics
    run_extract_gc_dup_metrics
    run_final_report
    OUTDIR="$OUTDIR" run_generate_xlsx
    log_info "========== QC PIPELINE COMPLETE =========="
    log_success "All outputs in: $OUTDIR"
}

main "$@"
conda deactivate
