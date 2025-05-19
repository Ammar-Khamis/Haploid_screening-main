#!/bin/bash
# Enhanced Haploid Screening Analysis Workflow
# For analyzing ONT data from eHAP1 cells with PB-DGTV vectors
# Includes enhanced gene frequency analysis and condition comparisons

set -e  # Exit on any error
set -u  # Treat unset variables as errors

##################################
# USER CONFIGURABLE PARAMETERS
##################################
# Set ONT_BARCODE to the barcode you want to analyze (e.g., "bc07" or "bc08")
ONT_BARCODE="barcode07"
THREADS=4
MIN_MAPQ=10                                     # Minimum mapping quality for filtering
MIN_LENGTH=50                                   # Minimum read length after trimming
HIGH_FREQ_THRESHOLD=5                           # Threshold for high-frequency insertion genes (≥5)

# Path configuration
BASE_DIR="/users/paxak6/Haploid_screening_ONT"
GENOME_REF="/users/paxak6/Haploid_screening_ONT/Human_Ref_GTF/Ensembls/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENE_ANNOTATIONS="/users/paxak6/Haploid_screening_ONT/Human_Ref_GTF/Ensembls/Homo_sapiens.GRCh38.113.chr.gtf"

# Illumina barcode treatments - mapping between barcode and treatment condition
# You can modify these based on your experimental design
declare -A TREATMENT_MAP
TREATMENT_MAP["01"]="hyPB_Control_DMSO"
TREATMENT_MAP["02"]="hyPB_Cordycepin_200uM"
TREATMENT_MAP["03"]="hyPB_Cordycepin_150uM"
TREATMENT_MAP["04"]="suPB_Control_DMSO"
TREATMENT_MAP["05"]="suPB_Cordycepin_200uM"
TREATMENT_MAP["06"]="suPB_Cordycepin_150uM"
TREATMENT_MAP["07"]="hyPB_Library_1"
TREATMENT_MAP["08"]="suPB_Library_2"

# Create a timestamped output directory for this run
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="$BASE_DIR/${ONT_BARCODE}_analysis_${TIMESTAMP}"
SCRIPT_DIR="$OUTPUT_DIR/scripts"  # Store copies of scripts used

##################################
# SETUP AND UTILITY FUNCTIONS
##################################
# Create directory structure
mkdir -p $OUTPUT_DIR/logs
mkdir -p $OUTPUT_DIR/working
mkdir -p $OUTPUT_DIR/demultiplexed
mkdir -p $OUTPUT_DIR/trimmed
mkdir -p $OUTPUT_DIR/mapped
mkdir -p $OUTPUT_DIR/insertion_sites
mkdir -p $OUTPUT_DIR/annotation
mkdir -p $OUTPUT_DIR/final_results
mkdir -p $OUTPUT_DIR/condition_comparison
mkdir -p $OUTPUT_DIR/high_frequency_genes
mkdir -p $SCRIPT_DIR

# Log function
log() {
    local level="INFO"
    if [ "$#" -eq 2 ]; then
        level="$1"
        shift
    fi
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" | tee -a "$OUTPUT_DIR/logs/workflow.log"
}

# Initialize conda
source /users/paxak6/miniconda3/etc/profile.d/conda.sh

# Define a function to run commands that require samtools, bedtools, or minimap2
run_biotools_cmd() {
    (conda activate biotools && "$@")
}

# Define a function to run commands that require LAST or perl scripts
run_lastenv_cmd() {
    (conda activate last_env && "$@")
}

# Function to check if a command succeeded
check_success() {
    local exit_code=$?
    local message="$1"
    if [ $exit_code -ne 0 ]; then
        log "ERROR" "Command failed: $message (exit code: $exit_code)"
        return 1
    else
        log "SUCCESS" "$message"
        return 0
    fi
}

# Function to check file existence and size
check_file() {
    local file="$1"
    local description="$2"
    
    if [ ! -f "$file" ]; then
        log "ERROR" "$description file not found: $file"
        return 1
    elif [ ! -s "$file" ]; then
        log "WARNING" "$description file is empty: $file"
        return 2
    else
        local size=$(du -h "$file" | cut -f1)
        log "File check: $description ($size): $file"
        return 0
    fi
}

# Save a copy of this script
cp "$0" "$SCRIPT_DIR/$(basename "$0")"
log "Starting workflow for $ONT_BARCODE in $OUTPUT_DIR"

##################################
# STEP 1: GENOME PREPARATION
##################################
log "=== STEP 1: GENOME PREPARATION ==="

# Create genome index if it doesn't exist
if [ ! -f "$GENOME_REF.fai" ]; then
    log "Creating genome index..."
    run_biotools_cmd samtools faidx "$GENOME_REF"
    check_success "Genome indexing"
else
    log "Using existing genome index: $GENOME_REF.fai"
fi

##################################
# STEP 2: INITIAL CONCATENATION
##################################
log "=== STEP 2: INITIAL CONCATENATION ==="

if [ ! -f "$OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz" ]; then
    log "Concatenating $ONT_BARCODE files..."
    cat $BASE_DIR/$ONT_BARCODE/*fastq.gz > "$OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz"
    check_success "Concatenating reads"
    check_file "$OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz" "Concatenated reads"
else
    log "Using existing concatenated reads: $OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz"
fi

##################################
# STEP 3: DEMULTIPLEXING
##################################
log "=== STEP 3: DEMULTIPLEXING ==="

# Check if barcodes file exists
if [ ! -f "$BASE_DIR/barcodes_and_adapters.fa" ]; then
    log "ERROR" "Barcodes and adapters file not found: $BASE_DIR/barcodes_and_adapters.fa"
    exit 1
fi

# Generate LAST index for barcodes/adapters if it doesn't exist
if [ ! -f "$BASE_DIR/barcodes_and_adapters.prj" ]; then
    log "Building LAST index for barcodes and adapters..."
    run_lastenv_cmd lastdb -uRY4 -R01 "$BASE_DIR/barcodes_and_adapters" "$BASE_DIR/barcodes_and_adapters.fa"
    check_success "Building LAST index"
else
    log "Using existing LAST index for barcodes and adapters"
fi

# Create training matrix if it doesn't exist
if [ ! -f "$OUTPUT_DIR/${ONT_BARCODE}.mat" ]; then
    log "Training LAST matrix for $ONT_BARCODE..."
    run_lastenv_cmd last-train -Q1 -P$THREADS "$BASE_DIR/barcodes_and_adapters" \
        "$OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz" > "$OUTPUT_DIR/${ONT_BARCODE}.mat"
    check_success "Training LAST matrix"
    check_file "$OUTPUT_DIR/${ONT_BARCODE}.mat" "LAST matrix"
else
    log "Using existing training matrix: $OUTPUT_DIR/${ONT_BARCODE}.mat"
fi

# Demultiplex if not already done
if [ ! -f "$OUTPUT_DIR/${ONT_BARCODE}_assignments.txt.gz" ]; then
    log "Demultiplexing $ONT_BARCODE reads..."
    run_lastenv_cmd lastal --split -p "$OUTPUT_DIR/${ONT_BARCODE}.mat" -P$THREADS "$BASE_DIR/barcodes_and_adapters" \
        <(zcat "$OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz") | \
        run_lastenv_cmd maf-convert -n tab | \
        cut -f 2,7 | \
        sort | \
        uniq | \
        gzip > "$OUTPUT_DIR/${ONT_BARCODE}_assignments.txt.gz"
    check_success "Creating barcode assignments"
    
    log "Running fastx-fetch.pl for demultiplexing..."
    run_lastenv_cmd $BASE_DIR/fastx-fetch.pl -lengths -demultiplex "$OUTPUT_DIR/${ONT_BARCODE}_assignments.txt.gz" \
        -prefix "$OUTPUT_DIR/demultiplexed/reads_" \
        <(zcat "$OUTPUT_DIR/reads_${ONT_BARCODE}.fastq.gz") > "$OUTPUT_DIR/${ONT_BARCODE}_counts.txt"
    check_success "Running fastx-fetch.pl"
    
    log "Demultiplexing complete for $ONT_BARCODE. Summary:"
    cat "$OUTPUT_DIR/${ONT_BARCODE}_counts.txt"
else
    log "Using existing demultiplexed data in $OUTPUT_DIR/demultiplexed/"
fi

##################################
# STEP 4: PROCESS EACH ILLUMINA BARCODE
##################################
log "=== STEP 4: PROCESS EACH ILLUMINA BARCODE ==="

process_barcode() {
    local illum_bc=$1
    local treatment="${TREATMENT_MAP[$illum_bc]}"
    
    log "Processing $ONT_BARCODE illumina barcode $illum_bc ($treatment)"
    
    # Merge left and right reads if not already done
    if [ ! -f "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all.fq.gz" ]; then
        log "Merging left+right reads for ${ONT_BARCODE}_bc${illum_bc}..."
        left_file="$OUTPUT_DIR/demultiplexed/reads__left_bc${illum_bc}-P5-${illum_bc}-PB5.fq.gz"
        right_file="$OUTPUT_DIR/demultiplexed/reads__right_bc${illum_bc}-P5-${illum_bc}-PB3.fq.gz"
        
        # Check if both files exist
        if [ -f "$left_file" ] && [ -f "$right_file" ]; then
            zcat "$left_file" "$right_file" > "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all.fq"
            gzip "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all.fq"
            check_success "Merging left+right reads"
        else
            log "WARNING" "Missing left or right files for bc${illum_bc}, skipping"
            return 1
        fi
    fi
    
    # Align to get adapter positions
    if [ ! -f "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_alignments.tsv" ]; then
        log "Aligning ${ONT_BARCODE}_bc${illum_bc}_all to find adapter positions..."
        run_lastenv_cmd lastal --split -p "$OUTPUT_DIR/${ONT_BARCODE}.mat" -P$THREADS "$BASE_DIR/barcodes_and_adapters" \
            <(zcat "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all.fq.gz") | \
            run_lastenv_cmd maf-convert -n tab > "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_alignments.tsv"
        check_success "Aligning reads to adapters"
    fi
    
    # Extract best alignment for each read
    if [ ! -f "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_for_trimdata.tsv" ]; then
        log "Finding best alignments for ${ONT_BARCODE}_bc${illum_bc}..."
        cat "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_alignments.tsv" | \
            cut -f7,8,9,10,13 | \
            awk 'BEGIN{FS="\t";OFS="\t"}{print $5,$2,$3,$4,$1}' | \
            cut -c8- | \
            sort -g -k1,1 | \
            nl | \
            sort -k6,6 -k1,1n | \
            uniq -f5 > "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_for_trimdata.tsv"
        check_success "Extracting best alignments"
    fi
    
    # Copy trimming scripts to our script directory if not already there
    if [ ! -f "$SCRIPT_DIR/Improved_trim_fastq.py" ]; then
        cp "$BASE_DIR/Improved_trim_fastq.py" "$SCRIPT_DIR/"
    fi
    
    if [ ! -f "$SCRIPT_DIR/Filter_TTAA.py" ]; then
        cp "$BASE_DIR/Filter_TTAA.py" "$SCRIPT_DIR/"
    fi
    
    # Trim reads based on alignment
    if [ ! -f "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed.fq.gz" ]; then
        log "Trimming reads for ${ONT_BARCODE}_bc${illum_bc}..."
        python "$SCRIPT_DIR/Improved_trim_fastq.py" \
            --fastq "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all.fq.gz" \
            --trimdata "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_for_trimdata.tsv" \
            --out "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed.fq" \
            --report "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_trim_report.txt"
        check_success "Trimming reads"
        
        log "Filtering for TTAA and minimum length ($MIN_LENGTH bp)..."
        python "$SCRIPT_DIR/Filter_TTAA.py" \
            --fastq "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed.fq" \
            --min_length $MIN_LENGTH \
            --out "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed_TTAA.fq" \
            --report "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_ttaa_report.txt"
        check_success "TTAA filtering"
        
        gzip "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed.fq"
        gzip "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed_TTAA.fq"
    fi
    
    # Verify FASTQ has reads
    check_file "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed_TTAA.fq.gz" "Trimmed TTAA FASTQ"
    read_count=$(zcat "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed_TTAA.fq.gz" | grep -c "^@" || echo "0")
    if [ "$read_count" -eq 0 ]; then
        log "WARNING" "No reads found in trimmed FASTQ, skipping mapping"
        return 1
    else
        log "Trimmed FASTQ contains $read_count reads"
    fi
    
    # Map to reference genome
    if [ ! -f "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_mapped.bam" ]; then
        log "Mapping ${ONT_BARCODE}_bc${illum_bc} to reference genome..."
        run_biotools_cmd minimap2 -ax map-ont --secondary=no -t $THREADS "$GENOME_REF" \
            <(zcat "$OUTPUT_DIR/trimmed/${ONT_BARCODE}_bc${illum_bc}_all_trimmed_TTAA.fq.gz") 2>"$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_minimap2.log" | \
            run_biotools_cmd samtools view -h -b -q $MIN_MAPQ | \
            run_biotools_cmd samtools sort > "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_mapped.bam"
        check_success "Mapping reads"
        
        run_biotools_cmd samtools index "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_mapped.bam"
        run_biotools_cmd samtools flagstat "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_mapped.bam" > "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_flagstat.txt"
        
        # Check if mapping produced results
        mapped_reads=$(grep "mapped (" "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_flagstat.txt" | head -1 | awk '{print $1}')
        log "${mapped_reads} reads mapped for ${ONT_BARCODE}_bc${illum_bc}"
    fi
    
    # Copy insertion site extraction script if not already there
    if [ ! -f "$SCRIPT_DIR/Extract_Insertion_Sites.py" ]; then
        cp "$BASE_DIR/Extract_Insertion_Sites.py" "$SCRIPT_DIR/"
    fi
    
    # Extract insertion sites
    if [ ! -f "$OUTPUT_DIR/insertion_sites/${ONT_BARCODE}_bc${illum_bc}_sites.bed" ]; then
        log "Extracting insertion sites for ${ONT_BARCODE}_bc${illum_bc}..."
        run_biotools_cmd samtools view "$OUTPUT_DIR/mapped/${ONT_BARCODE}_bc${illum_bc}_mapped.bam" | \
            python "$SCRIPT_DIR/Extract_Insertion_Sites.py" \
                --output "$OUTPUT_DIR/insertion_sites/${ONT_BARCODE}_bc${illum_bc}_sites.bed" \
                --min_mapq $MIN_MAPQ \
                --merge_window 5
        check_success "Extracting insertion sites"
        
        if [ -f "$OUTPUT_DIR/insertion_sites/${ONT_BARCODE}_bc${illum_bc}_sites.bed" ]; then
            site_count=$(wc -l < "$OUTPUT_DIR/insertion_sites/${ONT_BARCODE}_bc${illum_bc}_sites.bed")
            log "Extracted $site_count insertion sites for ${ONT_BARCODE}_bc${illum_bc}"
        else
            log "WARNING" "No insertion sites extracted for ${ONT_BARCODE}_bc${illum_bc}"
        fi
    fi
}

# Process all illumina barcodes
for illum_bc in {01..08}; do
    process_barcode $illum_bc
done

log "Processing of all illumina barcodes complete"

##################################
# STEP 5: PREPARE GTF FILE
##################################
log "=== STEP 5: PREPARE GTF FILE ==="

log "Preparing GTF file for annotation..."

# Extract gene entries from GTF
SIMPLE_GTF="$OUTPUT_DIR/annotation/simplified.gtf"
if [ ! -f "$SIMPLE_GTF" ]; then
    log "Extracting gene entries from GTF file..."
    # Using broader pattern match for gene entries as suggested
    grep -E '[[:space:]]gene[[:space:]]' "$GENE_ANNOTATIONS" > "$SIMPLE_GTF.tmp"
    check_success "Extracting gene entries"
    
    # Check if we found any entries
    if [ -s "$SIMPLE_GTF.tmp" ]; then
        log "Sorting GTF file..."
        sort -k1,1 -k4,4n "$SIMPLE_GTF.tmp" > "$SIMPLE_GTF"
        rm "$SIMPLE_GTF.tmp"
        log "Created simplified GTF with $(wc -l < "$SIMPLE_GTF") gene entries"
    else
        log "ERROR" "No gene entries found in GTF file"
        # Create a minimal dummy GTF as fallback
        {
            echo -e "1\t.\tgene\t1000\t2000\t.\t+\t.\tgene_id \"DUMMY1\""
            echo -e "2\t.\tgene\t1000\t2000\t.\t+\t.\tgene_id \"DUMMY2\""
            echo -e "X\t.\tgene\t1000\t2000\t.\t+\t.\tgene_id \"DUMMYX\""
        } > "$SIMPLE_GTF"
        log "WARNING" "Created minimal dummy GTF as fallback"
    fi
fi

##################################
# STEP 6: PREPARE CHROMOSOME FILES
##################################
log "=== STEP 6: PREPARE CHROMOSOME FILES ==="

# Extract all chromosomes from insertion sites
log "Collecting chromosomes from insertion sites..."
SITES_CHROMS="$OUTPUT_DIR/working/sites_chromosomes.txt"
find "$OUTPUT_DIR/insertion_sites" -name "*sites.bed" -type f -exec cut -f1 {} \; | sort -u > "$SITES_CHROMS"
log "Found $(wc -l < "$SITES_CHROMS") chromosomes in insertion sites"

# Extract all chromosomes from GTF file
log "Extracting chromosomes from GTF file..."
GTF_CHROMS="$OUTPUT_DIR/working/gtf_chromosomes.txt"
cut -f1 "$SIMPLE_GTF" | sort -u > "$GTF_CHROMS"
log "Found $(wc -l < "$GTF_CHROMS") chromosomes in GTF file"

# Combine all chromosomes
log "Creating comprehensive chromosome list..."
ALL_CHROMS="$OUTPUT_DIR/working/all_chromosomes.txt"
cat "$SITES_CHROMS" "$GTF_CHROMS" | sort -u > "$ALL_CHROMS"
log "Combined list has $(wc -l < "$ALL_CHROMS") unique chromosomes"

# Create comprehensive genome file
COMPREHENSIVE_GENOME="$OUTPUT_DIR/annotation/comprehensive_genome.txt"
awk '{print $1 "\t100000000"}' "$ALL_CHROMS" > "$COMPREHENSIVE_GENOME"
log "Created comprehensive genome file: $COMPREHENSIVE_GENOME"

# Filter GTF to only include chromosomes from insertion sites
log "Creating filtered GTF with only relevant chromosomes..."
FILTERED_GTF="$OUTPUT_DIR/annotation/filtered_simplified.gtf"
grep -f "$SITES_CHROMS" "$SIMPLE_GTF" > "$FILTERED_GTF"
log "Created filtered GTF with $(wc -l < "$FILTERED_GTF") entries"

# If filtered GTF is empty, use a fallback approach
if [ ! -s "$FILTERED_GTF" ]; then
    log "WARNING" "Filtered GTF is empty! Creating a minimal GTF file..."
    # Create a minimal valid GTF
    {
        while read chrom; do
            echo -e "${chrom}\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id \"DUMMY_${chrom}\"; gene_name \"DUMMY_${chrom}\""
        done < "$SITES_CHROMS"
    } > "$FILTERED_GTF"
    log "Created minimal GTF with $(wc -l < "$FILTERED_GTF") entries"
fi

##################################
# STEP 7: FINAL ANALYSIS
##################################
log "=== STEP 7: FINAL ANALYSIS ==="

# Create master sites file
log "Creating master sites file..."
MASTER_SITES="$OUTPUT_DIR/final_results/${ONT_BARCODE}_master_sites.bed"
ALL_SITES="$OUTPUT_DIR/working/all_sites_for_merge.bed"
mkdir -p "$OUTPUT_DIR/final_results"

# Check if we have any insertion site files
SITE_FILES=($(ls $OUTPUT_DIR/insertion_sites/${ONT_BARCODE}_bc*_sites.bed 2>/dev/null || echo ""))
if [ ${#SITE_FILES[@]} -gt 0 ]; then
    log "Combining sites from ${#SITE_FILES[@]} files..."
    cat "$OUTPUT_DIR"/insertion_sites/${ONT_BARCODE}_bc*_sites.bed > "$ALL_SITES"
    
    # Sort using the comprehensive genome file
    run_biotools_cmd bedtools sort -i "$ALL_SITES" -g "$COMPREHENSIVE_GENOME" > "$ALL_SITES.sorted"
    
    # Merge overlapping sites
    run_biotools_cmd bedtools merge -i "$ALL_SITES.sorted" -c 4,6 -o sum,distinct > "$MASTER_SITES"
    log "Created master sites file with $(wc -l < "$MASTER_SITES") entries"
else
    log "ERROR" "No insertion site files found!"
    exit 1
fi

# Annotate master sites
log "Annotating master sites..."
MASTER_ANNOTATED="$OUTPUT_DIR/final_results/${ONT_BARCODE}_master_annotated.bed"
run_biotools_cmd bedtools closest -a "$MASTER_SITES" -b "$FILTERED_GTF" -g "$COMPREHENSIVE_GENOME" -D ref > "$MASTER_ANNOTATED"
check_success "Annotating master sites"

# Process individual site files
log "Annotating individual condition files..."
for illum_bc in {01..08}; do
    treatment="${TREATMENT_MAP[$illum_bc]}"
    SITES_BED="$OUTPUT_DIR/insertion_sites/${ONT_BARCODE}_bc${illum_bc}_sites.bed"
    FINAL_SITES="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_sites.bed"
    ANNOTATED_BED="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_annotated.bed"
    COVERAGE_BG="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_coverage.bedgraph"
    
    if [ -f "$SITES_BED" ] && [ -s "$SITES_BED" ]; then
        log "Processing ${ONT_BARCODE}_bc${illum_bc} ($treatment)..."
        # Copy and sort the sites
        run_biotools_cmd bedtools sort -i "$SITES_BED" -g "$COMPREHENSIVE_GENOME" > "$FINAL_SITES"
        
        # Create coverage track (bedgraph)
        run_biotools_cmd bedtools genomecov -bg -i "$FINAL_SITES" -g "$COMPREHENSIVE_GENOME" > "$COVERAGE_BG"
        
        # Annotate with the filtered GTF
        run_biotools_cmd bedtools closest -a "$FINAL_SITES" -b "$FILTERED_GTF" -g "$COMPREHENSIVE_GENOME" -D ref > "$ANNOTATED_BED"
        log "Completed annotation and coverage track for ${ONT_BARCODE}_bc${illum_bc} ($treatment)"
    fi
done

# Create master coverage track
MASTER_COVERAGE="$OUTPUT_DIR/final_results/${ONT_BARCODE}_master_coverage.bedgraph"
run_biotools_cmd bedtools genomecov -bg -i "$MASTER_SITES" -g "$COMPREHENSIVE_GENOME" > "$MASTER_COVERAGE"
log "Created master coverage track: $MASTER_COVERAGE"

##################################
# STEP 8: IDENTIFY HIGH-FREQUENCY GENES
##################################
log "=== STEP 8: IDENTIFY HIGH-FREQUENCY GENES ==="

# Extract genes with high insertion frequencies (≥5)
HIGH_FREQ_GENES="$OUTPUT_DIR/high_frequency_genes/high_frequency_genes.txt"
HIGH_FREQ_BED="$OUTPUT_DIR/high_frequency_genes/high_frequency_genes.bed"

log "Identifying genes with ≥$HIGH_FREQ_THRESHOLD insertions..."

if grep -q "gene_name" "$MASTER_ANNOTATED"; then
    # Extract gene hits and count frequencies using read counts
    awk -F'\t' 'BEGIN {OFS="\t"}
    {
        if ($0 !~ /\t-1$/ && match($0, /gene_name "([^"]+)"/, name) && match($0, /gene_id "([^"]+)"/, id)) {
            # Use the actual read count from column 4
            genes[name[1]] += $4;  # Add the read count
            gene_ids[name[1]] = id[1];
        }
    } 
    END {
        for (gene in genes) {
            if (genes[gene] >= '$HIGH_FREQ_THRESHOLD') {
                print genes[gene], gene, gene_ids[gene];
            }
        }
    }' "$MASTER_ANNOTATED" | sort -nr > "$HIGH_FREQ_GENES"
    
    log "Found $(wc -l < "$HIGH_FREQ_GENES") genes with ≥$HIGH_FREQ_THRESHOLD insertions"
    
    # Create a BED file with locations of high-frequency insertions
    if [ -s "$HIGH_FREQ_GENES" ]; then
        # Extract gene names for grep
        cut -f2 "$HIGH_FREQ_GENES" > "$OUTPUT_DIR/high_frequency_genes/high_freq_gene_names.txt"
        
        # Extract all insertions in these genes
        grep -f "$OUTPUT_DIR/high_frequency_genes/high_freq_gene_names.txt" "$MASTER_ANNOTATED" > "$HIGH_FREQ_BED"
        
        log "Created BED file with high-frequency gene insertions: $HIGH_FREQ_BED"
        
        # Create individual files for each high-frequency gene with accurate counts
        mkdir -p "$OUTPUT_DIR/high_frequency_genes/individual_genes"
        
        while read line; do
            total_count=$(echo "$line" | cut -f1)
            gene=$(echo "$line" | cut -f2)
            gene_id=$(echo "$line" | cut -f3)
            
            # Create a sanitized filename
            safe_gene=$(echo "$gene" | tr ' /:*?"<>|' '_')
            
            # Extract all insertions for this gene
            grep "gene_name \"$gene\"" "$MASTER_ANNOTATED" > "$OUTPUT_DIR/high_frequency_genes/individual_genes/${safe_gene}_insertions.bed"
            
            # Count insertions per condition for this gene using read counts
            {
                echo "Gene: $gene (ID: $gene_id)"
                echo "Total insertions: $total_count"
                echo ""
                echo "Condition | Insertions"
                echo "--------- | ----------"
                
                for illum_bc in {01..08}; do
                    treatment="${TREATMENT_MAP[$illum_bc]}"
                    cond_file="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_annotated.bed"
                    
                    if [ -f "$cond_file" ]; then
                        # Sum the read counts (column 4) for this gene
                        cond_count=$(awk -v gene="$gene" '
                            $0 ~ "gene_name \"" gene "\"" {sum += $4} 
                            END {print sum}
                        ' "$cond_file")
                        echo "${treatment} | $cond_count"
                    else
                        echo "${treatment} | N/A"
                    fi
                done
            } > "$OUTPUT_DIR/high_frequency_genes/individual_genes/${safe_gene}_summary.txt"
        done < "$HIGH_FREQ_GENES"
        
        log "Created individual files for each high-frequency gene with accurate counts"
    fi
else
    log "WARNING" "Cannot identify genes by name (no gene_name attribute in annotations)"
    # Try using gene_id instead with the same read count-based approach
    awk -F'\t' 'BEGIN {OFS="\t"}
    {
        if ($0 !~ /\t-1$/ && match($0, /gene_id "([^"]+)"/, id)) {
            genes[id[1]] += $4;  # Add the read count
        }
    } 
    END {
        for (gene in genes) {
            if (genes[gene] >= '$HIGH_FREQ_THRESHOLD') {
                print genes[gene], gene;
            }
        }
    }' "$MASTER_ANNOTATED" | sort -nr > "$HIGH_FREQ_GENES"
    
    log "Found $(wc -l < "$HIGH_FREQ_GENES") genes with ≥$HIGH_FREQ_THRESHOLD insertions (using gene_id)"
fi

##################################
# STEP 9: CONDITION COMPARISON
##################################
log "=== STEP 9: CONDITION COMPARISON ==="

# Create a matrix for comparing insertion patterns across conditions
COMP_MATRIX="$OUTPUT_DIR/condition_comparison/condition_matrix.txt"
GENE_MATRIX="$OUTPUT_DIR/condition_comparison/gene_matrix.txt"

# Create chromosome-level comparison matrix
{
    printf "Chromosome"
    for illum_bc in {01..08}; do
        printf "\t%s" "${TREATMENT_MAP[$illum_bc]}"
    done
    printf "\n"
    
    # Get all unique chromosomes
    cut -f1 "$MASTER_SITES" | sort -u | while read chrom; do
        printf "%s" "$chrom"
        
        for illum_bc in {01..08}; do
            sites_file="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_sites.bed"
            
            if [ -f "$sites_file" ] && [ -s "$sites_file" ]; then
                count=$(grep -c "^$chrom" "$sites_file" || echo "0")
                printf "\t%s" "$count"
            else
                printf "\t0"
            fi
        done
        printf "\n"
    done
} > "$COMP_MATRIX"

# Fix any potential tab formatting issues in condition matrix
temp_file=$(mktemp)
awk '{gsub(/\\t/, "\t"); print}' "$COMP_MATRIX" > "$temp_file"
mv "$temp_file" "$COMP_MATRIX"

log "Created condition comparison matrix: $COMP_MATRIX"

# Create gene-level comparison matrix for top genes
if [ -s "$HIGH_FREQ_GENES" ]; then
    {
        printf "Gene\tGeneID\tTotal"
        for illum_bc in {01..08}; do
            printf "\t%s" "${TREATMENT_MAP[$illum_bc]}"
        done
        printf "\n"
        
        # Process each high-frequency gene
        cat "$HIGH_FREQ_GENES" | while read line; do
            gene=$(echo "$line" | cut -f2)
            gene_id=$(echo "$line" | cut -f3)
            
            # Store all condition counts to calculate total
            condition_counts=()
            
            # Get counts for each condition using read counts
            for illum_bc in {01..08}; do
                cond_file="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_annotated.bed"
                
                if [ -f "$cond_file" ]; then
                    # Sum the read counts (column 4) for this gene
                    cond_count=$(awk -v gene="$gene" '
                        $0 ~ "gene_name \"" gene "\"" {sum += $4} 
                        END {print sum+0}  # +0 ensures 0 for no matches
                    ' "$cond_file")
                    condition_counts+=("$cond_count")
                else
                    condition_counts+=("0")
                fi
            done
            
            # Calculate total from condition counts
            total_count=0
            for count in "${condition_counts[@]}"; do
                total_count=$((total_count + count))
            done
            
            # Print gene info and total
            printf "%s\t%s\t%d" "$gene" "$gene_id" "$total_count"
            
            # Print individual condition counts
            for count in "${condition_counts[@]}"; do
                printf "\t%d" "$count"
            done
            printf "\n"
        done
    } > "$GENE_MATRIX"
    
    log "Created gene-level comparison matrix: $GENE_MATRIX"
fi

# Create a heatmap-compatible matrix for visualization
HEATMAP_MATRIX="$OUTPUT_DIR/condition_comparison/heatmap_matrix.txt"
if [ -s "$GENE_MATRIX" ]; then
    # Create proper header for heatmap matrix
    {
        # Get the header line from gene matrix but skip the "Total" column
        head -1 "$GENE_MATRIX" | cut -f1-2,4- > "$HEATMAP_MATRIX"
        
        # Add the data rows, skipping the "Total" column
        tail -n +2 "$GENE_MATRIX" | cut -f1-2,4- >> "$HEATMAP_MATRIX"
    }
    log "Created matrix for heatmap visualization: $HEATMAP_MATRIX"
    
    # Fix any potential tab formatting issues
    for matrix_file in "$GENE_MATRIX" "$HEATMAP_MATRIX"; do
        # Create temporary file
        temp_file=$(mktemp)
        # Replace any literal "\t" with actual tabs
        awk '{gsub(/\\t/, "\t"); print}' "$matrix_file" > "$temp_file"
        # Replace original file with fixed version
        mv "$temp_file" "$matrix_file"
    done
    log "Verified tab formatting in matrix files"
fi

##################################
# STEP 10: GENERATE COMPREHENSIVE REPORTS
##################################
log "=== STEP 10: GENERATE COMPREHENSIVE REPORTS ==="

# Create a comprehensive report
REPORT="$OUTPUT_DIR/final_results/${ONT_BARCODE}_analysis_report.txt"
{
    echo "# Haploid Screening Insertion Site Analysis"
    echo "Date: $(date)"
    echo ""
    echo "## Summary Statistics"
    echo "Total unique insertion sites: $(wc -l < "$MASTER_SITES")"
    echo ""
    echo "## Condition-specific Counts"
    echo "Condition | Treatment | Sites"
    echo "--------- | --------- | -----"
    for illum_bc in {01..08}; do
        FINAL_SITES="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_sites.bed"
        if [ -f "$FINAL_SITES" ] && [ -s "$FINAL_SITES" ]; then
            COUNT=$(wc -l < "$FINAL_SITES")
            echo "${ONT_BARCODE}_bc${illum_bc} | ${TREATMENT_MAP[$illum_bc]} | $COUNT"
        else
            echo "${ONT_BARCODE}_bc${illum_bc} | ${TREATMENT_MAP[$illum_bc]} | N/A"
        fi
    done
    echo ""
    echo "## Chromosome Distribution"
    cut -f1 "$MASTER_SITES" | sort | uniq -c | sort -nr | head -20 | awk '{print $2 ": " $1 " sites"}'
    echo ""
    echo "## High-Frequency Genes (≥$HIGH_FREQ_THRESHOLD insertions)"
    if [ -s "$HIGH_FREQ_GENES" ]; then
        echo "The following genes have high insertion frequencies:"
        echo ""
        echo "Gene | Gene ID | Insertions"
        echo "---- | ------- | ----------"
        cat "$HIGH_FREQ_GENES" | awk '{print $2 " | " $3 " | " $1}'
        echo ""
        echo "See $OUTPUT_DIR/high_frequency_genes/ for detailed analysis of these genes."
    else
        echo "No genes with ≥$HIGH_FREQ_THRESHOLD insertions found."
    fi
    echo ""
    echo "## Condition Comparison"
    echo "A comparison matrix of insertions across conditions is available at:"
    echo "$COMP_MATRIX"
    echo ""
    echo "A gene-level comparison matrix is available at:"
    echo "$GENE_MATRIX"
    echo ""
    echo "## Analysis Parameters"
    echo "ONT Barcode: $ONT_BARCODE"
    echo "Reference genome: $GENOME_REF"
    echo "Gene annotations: $GENE_ANNOTATIONS"
    echo "Minimum mapping quality: $MIN_MAPQ"
    echo "Minimum read length: $MIN_LENGTH"
    echo "High-frequency threshold: $HIGH_FREQ_THRESHOLD insertions"
    echo ""
    echo "## Output Files"
    echo "Master sites: $MASTER_SITES"
    echo "Master annotated: $MASTER_ANNOTATED"
    echo "Master coverage: $MASTER_COVERAGE"
    echo ""
    echo "Individual condition files are in: $OUTPUT_DIR/final_results/"
    echo "High-frequency gene analysis: $OUTPUT_DIR/high_frequency_genes/"
    echo "Condition comparison files: $OUTPUT_DIR/condition_comparison/"
} > "$REPORT"

log "Created comprehensive analysis report: $REPORT"

# Create a simple HTML report for easier viewing
HTML_REPORT="$OUTPUT_DIR/final_results/${ONT_BARCODE}_report.html"
{
    echo "<!DOCTYPE html>"
    echo "<html><head><title>Haploid Screening Analysis Report</title>"
    echo "<style>"
    echo "body { font-family: Arial, sans-serif; line-height: 1.6; margin: 40px; }"
    echo "h1 { color: #2c3e50; }"
    echo "h2 { color: #3498db; margin-top: 30px; }"
    echo "h3 { color: #2980b9; }"
    echo "table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }"
    echo "th, td { text-align: left; padding: 12px; border-bottom: 1px solid #ddd; }"
    echo "th { background-color: #f2f2f2; }"
    echo "tr:hover {background-color: #f5f5f5;}"
    echo ".info { background-color: #d1ecf1; border: 1px solid #bee5eb; padding: 10px; border-radius: 5px; }"
    echo ".highlight { background-color: #fff3cd; }"
    echo "</style></head><body>"
    
    echo "<h1>Haploid Screening Insertion Site Analysis</h1>"
    echo "<p>Date: $(date)</p>"
    echo "<p class='info'>ONT Barcode: $ONT_BARCODE | Min MAPQ: $MIN_MAPQ | Min Length: $MIN_LENGTH | High-freq threshold: ≥$HIGH_FREQ_THRESHOLD</p>"
    
    echo "<h2>Summary Statistics</h2>"
    echo "<p>Total unique insertion sites: <strong>$(wc -l < "$MASTER_SITES")</strong></p>"
    
    echo "<h2>Condition-specific Counts</h2>"
    echo "<table><tr><th>Condition</th><th>Treatment</th><th>Sites</th></tr>"
    for illum_bc in {01..08}; do
        FINAL_SITES="$OUTPUT_DIR/final_results/${ONT_BARCODE}_bc${illum_bc}_sites.bed"
        if [ -f "$FINAL_SITES" ] && [ -s "$FINAL_SITES" ]; then
            COUNT=$(wc -l < "$FINAL_SITES")
            echo "<tr><td>${ONT_BARCODE}_bc${illum_bc}</td><td>${TREATMENT_MAP[$illum_bc]}</td><td>$COUNT</td></tr>"
        else
            echo "<tr><td>${ONT_BARCODE}_bc${illum_bc}</td><td>${TREATMENT_MAP[$illum_bc]}</td><td>N/A</td></tr>"
        fi
    done
    echo "</table>"
    
    echo "<h2>Chromosome Distribution</h2>"
    echo "<table><tr><th>Chromosome</th><th>Sites</th></tr>"
    # Don't use >> operator here - keep the data flow inside the {...} block
    cut -f1 "$MASTER_SITES" | sort | uniq -c | sort -nr | head -20 | \
    awk '{print "<tr><td>" $2 "</td><td>" $1 "</td></tr>"}'
    echo "</table>"
    
    if [ -s "$HIGH_FREQ_GENES" ]; then
        echo "<h2>High-Frequency Genes (≥$HIGH_FREQ_THRESHOLD insertions)</h2>"
        echo "<table><tr><th>Gene</th><th>Gene ID</th><th>Insertions</th></tr>"
        # Don't use >> operator here - keep the data flow inside the {...} block
        cat "$HIGH_FREQ_GENES" | \
        awk '{print "<tr><td>" $2 "</td><td>" $3 "</td><td>" $1 "</td></tr>"}'
        echo "</table>"
    fi
    
    if [ -s "$GENE_MATRIX" ]; then
        echo "<h2>Gene-Condition Matrix</h2>"
        echo "<p>Distribution of insertions per condition for high-frequency genes:</p>"
        # Create a table with headers from the matrix file
        echo "<table>"
        
        # Read the header row and convert to HTML
        head -1 "$GENE_MATRIX" | awk -F'\t' '{
            printf "<tr>";
            for(i=1; i<=NF; i++) {
                printf "<th>%s</th>", $i;
            }
            printf "</tr>\n";
        }'
        
        # Process the data rows
        tail -n +2 "$GENE_MATRIX" | awk -F'\t' '{
            printf "<tr>";
            for(i=1; i<=NF; i++) {
                printf "<td>%s</td>", $i;
            }
            printf "</tr>\n";
        }'
        echo "</table>"
    fi
    
    echo "<h2>Analysis Details</h2>"
    echo "<p>Reference genome: $GENOME_REF</p>"
    echo "<p>Gene annotations: $GENE_ANNOTATIONS</p>"
    echo "<p>Analysis directory: $OUTPUT_DIR</p>"
    
    echo "<h2>Output Files</h2>"
    echo "<ul>"
    echo "<li>Master sites: $MASTER_SITES</li>"
    echo "<li>Master annotated: $MASTER_ANNOTATED</li>"
    echo "<li>Master coverage: $MASTER_COVERAGE</li>"
    echo "<li>Individual condition files: $OUTPUT_DIR/final_results/</li>"
    echo "<li>High-frequency gene analysis: $OUTPUT_DIR/high_frequency_genes/</li>"
    echo "<li>Condition comparison matrices: $OUTPUT_DIR/condition_comparison/</li>"
    echo "</ul>"
    
    echo "</body></html>"
} > "$HTML_REPORT"

log "Created HTML report: $HTML_REPORT"
log "Workflow completed successfully."
