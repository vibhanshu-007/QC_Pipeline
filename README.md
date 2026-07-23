# QC Pipeline Formulas, Tools, and Threading

## Mathematical Formulas

### 1. Sequencing Data Size (GB)

**Formula**

    Data Size (GB) = (Total Reads × 2 × 150) / 1,000,000,000

-   2 = paired-end reads (R1 + R2)
-   150 = read length (bp)

### 2. Actual Sequencing Depth

    Actual Depth = Data Size (GB) / (Panel Size (MB) / 1000)

### 3. Minimum Depth

    Minimum Depth = Actual Depth × 0.70

### 4. Q20 Percentage

    Q20% = (Reads ≥ Q20 / Total Reads) × 100

### 5. Q30 Percentage

    Q30% = (Reads ≥ Q30 / Total Reads) × 100

### 6. Q40 Percentage

    Q40% = (Reads ≥ Q40 / Total Reads) × 100

### 7. GC Content

Weighted average from the FastQC GC distribution.

### 8. Duplication Percentage

    Duplication = 100 − %Deduplicated

### 9. Q30 Status

-   PASS if Q30 ≥ 90
-   FAIL otherwise

### 10. GC Status

-   PASS if GC ≤ 40%
-   FAIL if GC \> 40%

### 11. Duplication Status

-   PASS if Duplication ≤ 40%
-   FAIL otherwise

### 12. Read Count Check

-   PASS if R1 Reads = R2 Reads
-   FAIL otherwise

### 13. Minimum Data Size Check

-   PASS if Data Size ≥ Minimum Threshold

### 14. Expected Data Size Check

-   PASS if Data Size ≥ Expected Threshold

### 15. Expected Depth Check

-   PASS if Actual Depth ≥ Expected Depth

## Tools Used

  Tool           Purpose
  -------------- --------------------------
  bash           Pipeline scripting
  conda          Environment management
  GNU Parallel   Parallel execution
  md5sum         MD5 checksum
  FastQC         FASTQ QC
  MultiQC        Aggregate QC reports
  unzip          Read FastQC ZIPs
  awk            Parsing and calculations
  sed            Text processing
  grep           Search columns
  find           File discovery
  wc             Count files/lines
  head           Read headers
  tee            Logging
  stat           File statistics
  python3        Excel generation
  openpyxl       XLSX report creation
  realpath       Absolute paths

## Threading

Configuration:

``` bash
THREADS=16
FASTQC_JOBS=16
```

### MD5

    parallel -j 16

Runs 16 MD5 jobs concurrently.

### FastQC

    parallel --jobs 16
    fastqc -t 2

-   16 FastQC jobs
-   2 CPU threads per job
-   Up to 32 CPU threads may be used simultaneously.

### Single-threaded Steps

-   Read count extraction
-   MultiQC
-   Data size/depth calculations
-   Q20/Q30/Q40 calculations
-   FastQC module parsing
-   GC/Duplication extraction
-   Excel report generation

## Parallelization Summary

  Step                Parallel   Threads
  ------------------- ---------- ---------------------
  MD5 Check           Yes        16 jobs
  FastQC              Yes        16 jobs × 2 threads
  Read Count          No         1
  MultiQC             No         1
  Quality Metrics     No         1
  Data Size & Depth   No         1
  Module Parsing      No         1
  GC/Duplication      No         1
  Excel Report        No         1

