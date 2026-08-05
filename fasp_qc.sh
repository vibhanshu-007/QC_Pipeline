#!/usr/bin/env bash
#set -euo pipefail

# Pipeline build by ---> Vibhanshu singh

# Load Conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate environment
# conda activate qc_pipeline || {
#     echo "Failed to activate Conda environment."
#     exit 1
# }
################################################################################
# FAST Pre-QC PIPELINE
################################################################################

# ── Configuration ──────────────────────────────────────────────────────────────
FASTQ_DIR="$(realpath "${1:-$(pwd)}")"
OUTDIR="${FASTQ_DIR}/qc_output"
THREADS=16
FAST_P=8
PAIRSYNC_JOBS=$FAST_P
CSV_FILE=""

declare -A SIZE_THRESHOLDS=(
    ["TARGET_Indiegene_Liquid"]=46.9
    ["TARGET_Indiegene"]=4.9
    ["Germline_plus_"]=0.35
    ["TARGET_First_Solid_Lite"]=0.5
    ["TARGET_First_Liquid"]=7
    ["TARGT_FIRST_LIQUID_LITE"]=6
    ["TARGET_First_Solid"]=0.7
    ["TARGET_Absolute"]=20
    ["TA_Germline"]=6
)

declare -A MIN_SIZE_THRESHOLDS=(
    ["TARGET_Indiegene_Liquid"]=4.7
    ["TARGET_Indiegene"]=0.5
    ["Germline_plus_"]=0.035
    ["TARGET_First_Solid_Lite"]=0.05
    ["TARGET_First_Liquid"]=0.7
    ["TARGT_FIRST_LIQUID_LITE"]=0.6
    ["TARGET_First_Solid"]=0.07
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
FASTP_DIR=""
Read_Mismatch=""

declare -A MD5_STATUS
declare -A READCOUNT_STATUS
declare -A Min_Data_Size_Status
declare -A DATA_SIZE_STATUS
declare -A ACTUAL_DEPTH_STATUS  
declare -A SAMPLE_TO_TEST
declare -A PAIRSYNC_STATUS

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
    for tool in parallel md5sum python3 fastp stat; do
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

# Reads a single field out of a fastp JSON report.

get_fastp_field() {
    local json_file="$1"
    local field_path="$2"

    if [[ ! -f "$json_file" ]]; then
        echo "NA"
        return
    fi

    python3 - "$json_file" "$field_path" <<'PYEOF'
import json
import sys

json_file = sys.argv[1]
field_path = sys.argv[2]

try:
    with open(json_file) as f:
        data = json.load(f)
    val = data
    for key in field_path.split('.'):
        val = val[key]
    print(val)
except Exception:
    print("NA")
PYEOF
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

    TMP_MD5="$OUTDIR/md5_parallel.tmp"

    #############################################
    # Calculate MD5 for every FASTQ
    #############################################

    find "$FASTQ_DIR" -maxdepth 1 -name "*.fastq.gz" | \
    parallel -j "$THREADS" --no-notice '
        if [[ ! -s {} ]]; then
            echo "$(basename {})|NA|EMPTY"
        else
            echo "$(basename {})|$(md5sum {} | cut -d" " -f1)|OK"
        fi
    ' | sort > "$TMP_MD5"

    #############################################
    # Initialize every sample as PASS
    #############################################

    unset MD5_STATUS

    declare -gA MD5_STATUS
    declare -A MD5_COUNT
    declare -A MD5_FIRSTFILE

    while IFS='|' read -r file md5 status
    do
        [[ "$status" == "EMPTY" ]] && continue

        ((MD5_COUNT[$md5]++))

        if [[ -z "${MD5_FIRSTFILE[$md5]}" ]]; then
            MD5_FIRSTFILE[$md5]="$file"
        fi

    done < "$TMP_MD5"

    #############################################
    # Second pass
    #############################################

    while IFS='|' read -r file md5 status
    do

        sample="${file%%_R1*}"
        sample="${sample%%_R2*}"

        [[ -z "${MD5_STATUS[$sample]}" ]] && MD5_STATUS[$sample]="PASS"

        ####################################
        # Empty file
        ####################################

        if [[ "$status" == "EMPTY" ]]; then

            MD5_STATUS[$sample]="FAIL"

            echo "$file,NA,EMPTY_FILE" >> "$OUTFILE"

            continue
        fi

        ####################################
        # Duplicate
        ####################################

        if (( MD5_COUNT[$md5] > 1 ))
        then

            MD5_STATUS[$sample]="FAIL"

            echo "$file,$md5,DUPLICATE_OF:${MD5_FIRSTFILE[$md5]}" >> "$OUTFILE"

        else

            echo "$file,$md5,UNIQUE" >> "$OUTFILE"

        fi

    done < "$TMP_MD5"

    rm -f "$TMP_MD5"

    #############################################
    # Print sample summary
    #############################################

    log_info "Sample MD5 summary"

    for s in "${!MD5_STATUS[@]}"
    do
        log_info "$s : ${MD5_STATUS[$s]}"
    done

    log_success "MD5 check complete."
}

# ── Step 3: fastp with GNU Parallel ───────────────────────────────────────────

run_fastp() {

    log_info "Running fastp on MD5-passed samples using GNU Parallel..."

    FASTP_DIR="$OUTDIR/fastp"
    mkdir -p "$FASTP_DIR"

    local tmp_sample_list="$OUTDIR/fastp_samples.list"
    > "$tmp_sample_list"

    ############################################################
    # Calculate fastp threads per job
    ############################################################
    local fastp_threads

    if [[ -z "$FAST_P" || "$FAST_P" -lt 1 ]]; then
        FAST_P=1
    fi

    fastp_threads=$(( THREADS / FAST_P ))

    [[ "$fastp_threads" -lt 1 ]] && fastp_threads=1

    log_info "Total threads : $THREADS"
    log_info "Parallel jobs : $FAST_P"
    log_info "Threads/job   : $fastp_threads"

    ############################################################
    # Build FASTQ lookup (single directory scan)
    ############################################################
    declare -A R1_FILES
    declare -A R2_FILES

    for fq in "$FASTQ_DIR"/*.fastq.gz; do

        [[ ! -e "$fq" ]] && continue

        fname=$(basename "$fq")

        if [[ "$fname" =~ ^(.+)_R1.*\.fastq\.gz$ ]]; then
            R1_FILES["${BASH_REMATCH[1]}"]="$fq"

        elif [[ "$fname" =~ ^(.+)_R2.*\.fastq\.gz$ ]]; then
            R2_FILES["${BASH_REMATCH[1]}"]="$fq"
        fi

    done

    ############################################################
    # Build sample list
    ############################################################
    for sample in "${!MD5_STATUS[@]}"; do

        [[ "${MD5_STATUS[$sample]}" != "PASS" ]] && continue

        r1="${R1_FILES[$sample]}"
        r2="${R2_FILES[$sample]}"

        if [[ -z "$r1" || -z "$r2" ]]; then
            log_warn "FASTQ pair not found for sample '$sample' - skipping"
            continue
        fi

        printf "%s\t%s\t%s\n" \
            "$sample" \
            "$r1" \
            "$r2" >> "$tmp_sample_list"

    done

    if [[ ! -s "$tmp_sample_list" ]]; then
        log_warn "No FASTQ files found for fastp"
        rm -f "$tmp_sample_list"
        return
    fi

    export FASTP_DIR

    local file_count
    file_count=$(wc -l < "$tmp_sample_list")

    log_info "Starting fastp on $file_count sample(s)..."

    parallel \
        --jobs "$FAST_P" \
        --colsep '\t' \
        --halt soon,fail=1 \
        --joblog "$OUTDIR/fastp_parallel.log" \
        "fastp \
            -i {2} \
            -I {3} \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            --thread $fastp_threads \
            --json \"$FASTP_DIR/{1}.json\" \
            -h /dev/null \
            > \"$FASTP_DIR/{1}.fastp.log\" 2>&1" \
        :::: "$tmp_sample_list"

    local rc=$?

    rm -f "$tmp_sample_list"

    if [[ $rc -eq 0 ]]; then
        log_success "fastp parallel processing complete"
    else
        log_warn "fastp finished with exit code $rc"
        log_info "See $OUTDIR/fastp_parallel.log"
    fi

    log_success "fastp complete. Results: $FASTP_DIR"
}


#-------------------Read pair sync------------
run_pair_sync_check() {

    log_info "Running paired-end read synchronization check..."

    local outfile="$OUTDIR/pair_sync.csv"
    echo "Sample_ID,Reads_Checked,Status" > "$outfile"

    export FASTQ_DIR

    check_pair_sync() {

        local sample="$1"

        local r1="$FASTQ_DIR/${sample}_R1.fastq.gz"
        local r2="$FASTQ_DIR/${sample}_R2.fastq.gz"

        if [[ ! -f "$r1" || ! -f "$r2" ]]; then
            echo "$sample,0,FAIL"
            return
        fi

        perl -MIO::Uncompress::Gunzip=gunzip -e '

use strict;
use warnings;

my ($r1,$r2,$sample)=@ARGV;

my $fh1=IO::Uncompress::Gunzip->new($r1)
    or do { print "$sample,0,FAIL\n"; exit; };

my $fh2=IO::Uncompress::Gunzip->new($r2)
    or do { print "$sample,0,FAIL\n"; exit; };

my $count=0;

while (1){

    my $h1=<$fh1>;
    my $h2=<$fh2>;

    last if !defined($h1) && !defined($h2);

    if(!defined($h1) || !defined($h2)){
        print "$sample,$count,FAIL\n";
        exit;
    }

    chomp($h1);
    chomp($h2);

    my ($id1)=split(/\s+/,$h1);
    my ($id2)=split(/\s+/,$h2);

    $id1 =~ s/\/[12]$//;
    $id2 =~ s/\/[12]$//;

    if($id1 ne $id2){
        print "$sample,$count,FAIL\n";
        exit;
    }

    <$fh1>; <$fh1>; <$fh1>;
    <$fh2>; <$fh2>; <$fh2>;

    $count++;
}

print "$sample,$count,PASS\n";

' "$r1" "$r2" "$sample"
    }

    export -f check_pair_sync

    printf "%s\n" "${!READCOUNT_STATUS[@]}" |
    while read -r sample; do
        [[ "${READCOUNT_STATUS[$sample]}" == "PASS" ]] && echo "$sample"
    done |
    parallel -j "$PAIRSYNC_JOBS" --no-notice check_pair_sync {} >> "$outfile"

    while IFS=',' read -r sample reads status; do

        [[ "$sample" == "Sample_ID" ]] && continue

        PAIRSYNC_STATUS["$sample"]="$status"

        if [[ "$status" == "PASS" ]]; then
            log_success "Pair synchronization PASS for '$sample' ($reads reads)"
        else
            log_warn "Pair synchronization FAILED for '$sample' after $reads reads"
        fi

    done < "$outfile"

    log_success "Pair synchronization check complete. Results: $outfile"
}
# ── Step 4: Read-count check ───────────────────────────────────────────────────
run_readcount_check() {
    log_info "Checking R1/R2 read counts..."
    READCOUNT_CSV="$OUTDIR/readcounts.csv"
    echo "Sample_ID,R1_Count,R2_Count,Status" > "$READCOUNT_CSV"
    for sample in "${!MD5_STATUS[@]}"; do
        [[ "${MD5_STATUS[$sample]}" != "PASS" ]] && continue
        local json_file="$FASTP_DIR/${sample}.json"
        if [[ ! -f "$json_file" ]]; then
            log_warn "fastp JSON not found for sample '$sample' — skipping read count check"
            READCOUNT_STATUS["$sample"]="FAIL"
            echo "$sample,NA,NA,FAIL" >> "$READCOUNT_CSV"
            continue
        fi
        local r1_count r2_count
        r1_count=$(get_fastp_field "$json_file" "read1_before_filtering.total_reads")
        r2_count=$(get_fastp_field "$json_file" "read2_before_filtering.total_reads")
        if [[ "$r1_count" == "NA" || "$r2_count" == "NA" ]]; then
            log_warn "Malformed or incomplete fastp JSON for '$sample' — skipping read count check"
            READCOUNT_STATUS["$sample"]="FAIL"
            echo "$sample,$r1_count,$r2_count,FAIL" >> "$READCOUNT_CSV"
            continue
        fi
if [[ "$r1_count" -eq 0 || "$r2_count" -eq 0 ]]; then
    READCOUNT_STATUS["$sample"]="FAIL"
    log_warn "Zero reads detected for '$sample': R1=$r1_count R2=$r2_count"

elif [[ "$r1_count" -eq "$r2_count" ]]; then
    READCOUNT_STATUS["$sample"]="PASS"

else
    READCOUNT_STATUS["$sample"]="FAIL"
    log_warn "Read count mismatch for '$sample': R1=$r1_count R2=$r2_count"
fi

echo "$sample,$r1_count,$r2_count,${READCOUNT_STATUS[$sample]}" >> "$READCOUNT_CSV"
    done
    log_success "Read count check complete. Results: $READCOUNT_CSV"
}

# ── Step 5: generate_pre_qc_review ─────────────────────────────────────────────
generate_pre_qc_review() {

    log_info "Generating Pre-QC_review.csv..."

    python3 - <<PY
import csv
import os

outdir = "$OUTDIR"

md5_file = os.path.join(outdir, "md5_check.csv")
read_file = os.path.join(outdir, "readcounts.csv")
depth_file = os.path.join(outdir, "datasize.csv")
outfile = os.path.join(outdir, "Pre_QC_review.csv")

# ------------------------------------------------------------------
# Load MD5
# ------------------------------------------------------------------
md5 = {}

with open(md5_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample = row["Filename"].split("_R")[0]

        if row["Status"] == "UNIQUE":
            md5[sample] = "PASS"
        else:
            md5[sample] = "FAIL"

# ------------------------------------------------------------------
# Load Read Counts
# ------------------------------------------------------------------
readcount = {}

with open(read_file) as f:
    reader = csv.DictReader(f)

    for row in reader:
        readcount[row["Sample_ID"]] = row["Status"]

# ------------------------------------------------------------------
# Load Datasize / Depth
# ------------------------------------------------------------------
rows = []

with open(depth_file) as f:

    reader = csv.DictReader(f)

    for row in reader:

        sample = row["Sample_ID"]

        md5_status = md5.get(sample, "FAIL")
        rc_status = readcount.get(sample, "FAIL")
        depth_status = row["Actual_Depth_Status"]

        remarks = []

        if md5_status != "PASS":
            remarks.append("MD5 failed")

        if rc_status != "PASS":
            remarks.append("Read count mismatch")

        if depth_status != "PASS":
            remarks.append("Actual depth below threshold")

        pre_qc = "PASS" if not remarks else "FAIL"

        rows.append({
            "Sample_ID": sample,
            "Test_Name": row["Test_Name"],
            "MD5_Status": md5_status,
            "ReadCount_Status": rc_status,
            "Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth": "; ".join(remarks),
            "Actual_Depth_Status": depth_status,
            "Pre_QC_Status": pre_qc
        })

# ------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------
with open(outfile, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Sample_ID",
            "Test_Name",
            "MD5_Status",
            "ReadCount_Status",
            "Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth",
            "Actual_Depth_Status",
            "Pre_QC_Status"
        ]
    )

    writer.writeheader()
    writer.writerows(rows)

print(f"Generated {outfile}")

PY

    log_success "Pre_QC_review.csv generated."

}
# ── Generate list of PreQC PASS samples ───────────────────────────────────────
generate_preqc_pass_samples() {

    log_info "Generating PreQC PASS sample list..."

    local input_csv="$OUTDIR/Pre_QC_review.csv"
    local output_csv="$OUTDIR/PreQC_PASS_Samples.csv"

    [[ ! -f "$input_csv" ]] && {
        log_error "Pre_QC_review.csv not found."
        return 1
    }

    awk -F',' '
    BEGIN { OFS="," }

    NR==1 {
        print "Sample_ID","Test_Name","Pre_QC_Status"
        next
    }

    {
        gsub(/\r/, "", $7)
    }

    $7 == "PASS" {
        print $1, $2, $7
    }
    ' "$input_csv" > "$output_csv"

    local count
    count=$(( $(wc -l < "$output_csv") - 1 ))

    log_success "Generated $output_csv"
    log_info "PreQC PASS samples: $count"
}
# ── Step 6: Q20/Q30/Q40 metrics (from fastp JSON) ─────────────────────────────

run_quality_metrics() {

    log_info "Extracting Q20/Q30/Q40 metrics from fastp JSON..."

    echo "Sample_ID,PCTQ20,PCTQ30,PCTQ40" > "$OUTDIR/qmetrics.csv"

    for sample in "${!READCOUNT_STATUS[@]}"; do

        # Skip failed samples
        [[ "${READCOUNT_STATUS[$sample]}" != "PASS" ]] && continue
        #[[ "${ACTUAL_DEPTH_STATUS[$sample]:-PASS}" == "FAIL" ]] && continue

        local json_file="$FASTP_DIR/${sample}.json"

        if [[ ! -f "$json_file" ]]; then
            log_warn "fastp JSON not found for '$sample' — skipping Q-metrics"
            continue
        fi

        ########################################################
        # Read values from fastp JSON
        ########################################################

        local q20_rate q30_rate
        local total_bases
        local q40_r1 q40_r2 q40_bases
        local pct20 pct30 pct40

        q20_rate=$(get_fastp_field "$json_file" "summary.before_filtering.q20_rate")
        q30_rate=$(get_fastp_field "$json_file" "summary.before_filtering.q30_rate")
        total_bases=$(get_fastp_field "$json_file" "summary.before_filtering.total_bases")

        # q40_bases are available only in read1/read2 sections
        q40_r1=$(get_fastp_field "$json_file" "read1_before_filtering.q40_bases")
        q40_r2=$(get_fastp_field "$json_file" "read2_before_filtering.q40_bases")

        ########################################################
        # Q20 %
        ########################################################

        if [[ "$q20_rate" == "NA" ]]; then
            pct20="NA"
            log_warn "Missing q20_rate for '$sample'"
        else
            pct20=$(awk -v r="$q20_rate" 'BEGIN{printf "%.2f", r*100}')
        fi

        ########################################################
        # Q30 %
        ########################################################

        if [[ "$q30_rate" == "NA" ]]; then
            pct30="NA"
            log_warn "Missing q30_rate for '$sample'"
        else
            pct30=$(awk -v r="$q30_rate" 'BEGIN{printf "%.2f", r*100}')
        fi

        ########################################################
        # Q40 %
        ########################################################

        if [[ "$q40_r1" == "NA" || "$q40_r2" == "NA" ]]; then
            pct40="NA"
            log_warn "Missing q40_bases for '$sample'"
        else
            q40_bases=$((q40_r1 + q40_r2))

            if [[ "$total_bases" == "NA" || "$total_bases" -eq 0 ]]; then
                pct40="NA"
                log_warn "Missing total_bases for '$sample'"
            else
                pct40=$(awk -v q="$q40_bases" -v t="$total_bases" \
                    'BEGIN{printf "%.2f", (q/t)*100}')
            fi
        fi

        echo "$sample,$pct20,$pct30,$pct40" >> "$OUTDIR/qmetrics.csv"

    done

    log_success "Q20/Q30/Q40 extraction complete. Results: $OUTDIR/qmetrics.csv"
}



# ── Step 7: Data size & Actual depth check ────────────────────────────────────

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
        
        local json_file="$FASTP_DIR/${sample}.json"
        if [[ ! -f "$json_file" ]]; then
            log_warn "fastp JSON not found for '$sample' — skipping data size check"
            continue
        fi
        
        local sequences
        sequences=$(get_fastp_field "$json_file" "read1_before_filtering.total_reads")
        if [[ "$sequences" == "NA" ]]; then
            log_warn "Could not extract read count from fastp JSON for '$sample' — skipping data size check"
            continue
        fi

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

# ── Step 8: GC% and Duplication% from fastp JSON ──────────────────────────────
run_extract_gc_dup_metrics() {
    log_info "Extracting average GC and Duplication metrics from fastp JSON..."

    local gc_dup_csv="$OUTDIR/gc_dup_metrics.csv"

    echo "Sample_ID,GC_Content_Percent,GC_FLAG,Duplication_Percent,DUPLICATION_FLAG" > "$gc_dup_csv"

    for sample in "${!READCOUNT_STATUS[@]}"; do

        [[ "${READCOUNT_STATUS[$sample]}" != "PASS" ]] && continue

        local json_file="$FASTP_DIR/${sample}.json"
        if [[ ! -f "$json_file" ]]; then
            log_warn "fastp JSON missing for sample '$sample'"
            continue
        fi

########################################################
# GC %
########################################################

local gc_content avg_gc

gc_content=$(get_fastp_field "$json_file" "summary.before_filtering.gc_content")

if [[ "$gc_content" == "NA" ]]; then
    log_warn "Missing GC content in fastp JSON for '$sample'"
    avg_gc="NA"
else
    avg_gc=$(awk -v gc="$gc_content" 'BEGIN{printf "%.2f", gc*100}')
fi

        ########################
        # Duplication %
        ########################

        local dup_rate avg_dup

        dup_rate=$(get_fastp_field "$json_file" "duplication.rate")

        if [[ "$dup_rate" == "NA" ]]; then
            log_warn "Missing duplication rate in fastp JSON for '$sample'"
            avg_dup="NA"
        else
            avg_dup=$(awk -v r="$dup_rate" 'BEGIN{printf "%.2f", r*100}')
        fi

        ########################
        # Flags
        ########################

        local gc_flag dup_flag

        if [[ "$avg_gc" == "NA" ]]; then
            gc_flag="NA"
        elif awk -v gc="$avg_gc" -v th="$GC_PERCENT_THRESHOLD" 'BEGIN{exit !(gc > th)}'
        then
            gc_flag="FAIL"
        else
            gc_flag="PASS"
        fi

        if [[ "$avg_dup" == "NA" ]]; then
            dup_flag="NA"
        elif awk -v dup="$avg_dup" -v th="$DUPLICATION_THRESHOLD" 'BEGIN{exit !(dup > th)}'
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

# ── Step 9: Final per-sample conclusion ───────────────────────────────────────
run_final_report() {
    log_info "Generating final QC report..."
    FINAL_REPORT="$OUTDIR/final_report.csv"
    echo "Sample_ID,MD5_Status,ReadCount_Status,Pair_Synchronization_Status,MinSize_Status,DataSize_Status,Actual_Depth_Status,PCTQ30,Q30_Status,GC_Content_Percent,GC_FLAG,Duplication_Percent,DUPLICATION_FLAG,Conclusion,Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth,Pre_QC_Status" > "$FINAL_REPORT" 

    declare -A PCTQ30_MAP
    if [[ -f "$OUTDIR/qmetrics.csv" ]]; then
        while IFS=',' read -r sid pct20 pct30 pct40; do
            [[ "$sid" == "Sample_ID" ]] && continue
            PCTQ30_MAP["$sid"]="$pct30"
        done < "$OUTDIR/qmetrics.csv"
    fi
declare -A PAIRSYNC_MAP

if [[ -f "$OUTDIR/pair_sync.csv" ]]; then
    while IFS=',' read -r sid reads_checked status; do
        [[ "$sid" == "Sample_ID" ]] && continue
        PAIRSYNC_MAP["$sid"]="$status"
    done < "$OUTDIR/pair_sync.csv"
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
        local pair_sync="${PAIRSYNC_MAP[$sample]:-NOT_RUN}"
        local min_size="${Min_Data_Size_Status[$sample]:-NOT_RUN}"
        local data_size="${DATA_SIZE_STATUS[$sample]:-NOT_RUN}"
        local actual_depth_status="${ACTUAL_DEPTH_STATUS[$sample]:-NOT_RUN}"
        local pctq30="${PCTQ30_MAP[$sample]:-0}"

        local gc_pct="${GC_PERCENT[$sample]:-0}"
        local gc_flag="${GC_FLAG[$sample]:-NA}"
        local dup_pct="${DUP_PERCENT[$sample]:-0}"
        local dup_flag="${DUP_FLAG[$sample]:-NA}"

        local q30_status="FAIL"
        if awk "BEGIN{exit !(${pctq30} >= ${Q30_THRESHOLD})}" 2>/dev/null; then
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

if [[ "$pair_sync" != "PASS" ]]; then
    pre_qc_comments+=("Read pair synchronization failed")
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

        # Check if any of the four core pre-QC checks failed
        if [[ "$pre_qc_overall" == "FAIL" ]]; then
            # If any pre-QC check failed, conclusion is just the pre-QC comments
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

printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s","%s","%s"\n' \
"$sample" \
"$md5" \
"$readcount" \
"$pair_sync" \
"$min_size" \
"$data_size" \
"$actual_depth_status" \
"$pctq30" \
"$q30_status" \
"$gc_pct" \
"$gc_flag" \
"$dup_pct" \
"$dup_flag" \
"${conclusion//\"/\"\"}" \
"${pre_qc_comment_str//\"/\"\"}" \
"$pre_qc_overall" >> "$FINAL_REPORT"
    done

    log_success "Final report: $FINAL_REPORT"
}

# ── Step 10: Generate master XLSX ─────────────────────────────────────────────
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

# Metrics list with Pre-QC columns - split into two groups
pre_qc_metrics = [
    "MD5_Status",
    "ReadCount_Status",
    "Pair_Synchronization_Status",
    "Pre_QC_Status"
]

downstream_metrics = [
    "Q30_Status",
    "GC_FLAG",
    "DUPLICATION_FLAG"
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
ds_map = {r["Sample_ID"]: r for r in ds_rows if "Sample_ID" in r}
qm_map = {r["Sample_ID"]: r for r in qm_rows if "Sample_ID" in r}
rc_map = {r["Sample_ID"]: r for r in rc_rows if "Sample_ID" in r}
gc_dup_map = {r["Sample_ID"]: r for r in gc_dup_rows if "Sample_ID" in r}

# Column order (FastQC module columns removed; everything else unchanged)
COMBINED_HEADERS = [
    "Sample_ID","Test_Name","MD5_Status","R1_Count","R2_Count","ReadCount_Status",
    "Data_GB","Min_Data_Threshold_GB","Expected_Data_Threshold_GB","Min_Data_Size_Status","Expected_Data_Size_Status",
    "Panel_MB","Actual_Depth","Expected_Depth","Actual_Depth_Status",
    "PCTQ20","PCTQ30","PCTQ40","Q30_Status",
    "GC_Content_Percent","GC_FLAG","Duplication_Percent","DUPLICATION_FLAG",
    "Conclusion","Pair_Synchronization_Status",
    "Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth","Pre_QC_Status"
]
write_header(ws_samples, COMBINED_HEADERS, row=1)

for ri, row in enumerate(final_rows, start=2):
    sid = row.get("Sample_ID", "")
    ds = ds_map.get(sid, {})
    qm = qm_map.get(sid, {})
    rc = rc_map.get(sid, {})
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
        "PCTQ20": qm.get("PCTQ20", ""),
        "PCTQ30": qm.get("PCTQ30", ""),
        "PCTQ40": qm.get("PCTQ40", ""),
        "Q30_Status": q30_status,
        "GC_Content_Percent": gc_dup.get("GC_Content_Percent", ""),
        "GC_FLAG": gc_dup.get("GC_FLAG", ""),
        "Duplication_Percent": gc_dup.get("Duplication_Percent", ""),
        "DUPLICATION_FLAG": gc_dup.get("DUPLICATION_FLAG", ""),
        "Conclusion": row.get("Conclusion", ""),
        "Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth": row.get("Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth", ""),
        "Pre_QC_Status": row.get("Pre_QC_Status", ""),
        "Pair_Synchronization_Status": row.get("Pair_Synchronization_Status",""),
        
    }
    for ci, h in enumerate(COMBINED_HEADERS, start=1):
        val = combined.get(h, "")
        cell = ws_samples.cell(row=ri, column=ci, value=val)
        cell_fill(cell, val)
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
        cell.border = THIN_BORDER

ws_samples.freeze_panes = "A2"
auto_col_width(ws_samples)

# PreQC worksheet
ws_preqc = wb.create_sheet("PreQC")

preqc_data = []
for row in final_rows:
    sid = row.get("Sample_ID", "")
    ds = ds_map.get(sid, {})
    preqc_data.append({
        "Sample_ID": sid,
        "Test_Name": ds.get("Test_Name", ""),
        "Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth": row.get("Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth", ""),
        "Pre_QC_Status": row.get("Pre_QC_Status", "")
    })

PREQC_HEADERS = ["Sample_ID", "Test_Name", "Pre_QC_Status_MD5_Readcount_PairSync_Actual_Depth", "Pre_QC_Status"]
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
    run_fastp
    run_readcount_check
    run_pair_sync_check
    run_datasize_check
    generate_pre_qc_review
    generate_preqc_pass_samples
    run_quality_metrics
    run_extract_gc_dup_metrics
    run_final_report
    OUTDIR="$OUTDIR" run_generate_xlsx
    log_info "========== QC PIPELINE COMPLETE =========="
    log_success "All outputs in: $OUTDIR"
}

main "$@"
# conda deactivate
